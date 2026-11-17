#!/usr/bin/env python3
"""Summarize input and output synapses across brain regions for PAG neurons.

Read frozen per-neuron synapse tables and neuron metadata, assign synapses to
MECE atlas regions, and export regional proportions, counts, figures, and
Neuroglancer annotation layers. This analysis runs offline with the bundled
inputs; it does not query CAVE.

Paths default to this analysis directory and can be overridden on the command
line. Run ``python scripts/all_cells_synaptome.py --help`` from synapse_regions
for options. See ../README.md for methods, data provenance, and usage.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

try:
    from .mece_region_assignment import (
        HINDBRAIN_MECE_LEVEL_1_REGIONS,
        coarsen_hindbrain_to_mece_level_1,
    )
except ImportError:  # Direct execution: python scripts/all_cells_synaptome.py
    from mece_region_assignment import (
        HINDBRAIN_MECE_LEVEL_1_REGIONS,
        coarsen_hindbrain_to_mece_level_1,
    )

# Configuration

VOX_SIZE_UM_EM = (0.016, 0.016, 0.030)
NG_VOX_M = (16e-9, 16e-9, 30e-9)
MIDLINE_POSITION = (27909, 16358, 6034)
SEGMENTATION_SOURCE = "graphene://middleauth+https://pcgv3local.brain-wire-test.org/segmentation/table/fish1_v250915/"
PLOT_DPI = 300

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
H5_MASKS = ANALYSIS_DIR / "data" / "fish1.0_MECE_labels.hdf5"
CELLS_BASE = ANALYSIS_DIR / "data" / "per_neuron"
DEFAULT_CELLS_CSV = ANALYSIS_DIR / "data" / "neurons.csv"
RESULTS_DIR = ANALYSIS_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
})

# Path helpers

def ensure_dirs() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    (RESULTS_DIR / "ng_layers").mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def figure_pdf_path(name: str) -> Path:
    if name.endswith(".pdf"):
        name = name[:-4]
    return FIGURES_DIR / f"{name}.pdf"


def cell_dir(soma_id: int) -> Path:
    return CELLS_BASE / f"cell_{soma_id}"


def cell_csv_inputs(soma_id: int) -> Path:
    return cell_dir(soma_id) / f"cell_{soma_id}_inputs.csv"


def cell_csv_outputs(soma_id: int) -> Path:
    return cell_dir(soma_id) / f"cell_{soma_id}_outputs.csv"

# Coordinate and region helpers

def px_to_um_arrays(x_px, y_px, z_px, vox_um) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sx, sy, sz = vox_um
    return x_px.astype(float) * sx, y_px.astype(float) * sy, z_px.astype(float) * sz


def load_label_volume(h5_path: str) -> Tuple[np.ndarray, List[str], Tuple[float, float, float]]:
    with h5py.File(h5_path, "r") as f:
        if "label_volume" not in f:
            raise RuntimeError("HDF5 does not contain /label_volume")
        label = np.ascontiguousarray(f["label_volume"][()])
        if label.dtype != np.uint16:
            label = label.astype(np.uint16, copy=False)

        names_attr = f.attrs.get("region_names", None)
        if names_attr is None:
            region_names = [f"region_{i}" for i in range(1, int(label.max()) + 1)]
        else:
            raw = names_attr[()] if hasattr(names_attr, "__getitem__") else names_attr
            region_names = [
                x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x)
                for x in raw
            ]

        vox = f.attrs.get("voxel_size_um", (0.512, 0.512, 0.480))
        voxel_size_um = tuple(float(v) for v in vox)
        return label, region_names, voxel_size_um


def um_to_voxel(x_um, y_um, z_um, voxel_size_um) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    vx, vy, vz = voxel_size_um
    return (
        np.round(np.asarray(x_um, float) / vx).astype(np.int64),
        np.round(np.asarray(y_um, float) / vy).astype(np.int64),
        np.round(np.asarray(z_um, float) / vz).astype(np.int64),
    )


def map_df_to_regions(
    df: pd.DataFrame,
    prefix: str,
    label: np.ndarray,
    region_names: List[str],
    voxel_size_um: Tuple[float, float, float],
) -> pd.DataFrame:
    x_col, y_col, z_col = f"{prefix}_x_um", f"{prefix}_y_um", f"{prefix}_z_um"
    if not all(c in df.columns for c in [x_col, y_col, z_col]):
        raise ValueError(f"Missing required columns: {x_col}, {y_col}, {z_col}")

    X, Y, Z = label.shape
    xi, yi, zi = um_to_voxel(df[x_col], df[y_col], df[z_col], voxel_size_um)
    regions = []
    for a, b, c in zip(xi, yi, zi):
        if np.any(np.isnan([a, b, c])):
            regions.append([])
            continue
        a, b, c = int(a), int(b), int(c)
        if 0 <= a < X and 0 <= b < Y and 0 <= c < Z:
            lab = int(label[a, b, c])
            regions.append(
                [coarsen_hindbrain_to_mece_level_1(region_names[lab - 1])]
                if lab > 0
                else []
            )
        else:
            regions.append([])

    out = df.copy()
    out[f"{prefix}_regions"] = regions
    out[f"{prefix}_n_regions"] = out[f"{prefix}_regions"].apply(len)
    return out


def counts_from_region_lists(series_of_lists: pd.Series) -> pd.DataFrame:
    exploded = series_of_lists.explode()
    exploded = exploded[exploded.notna() & (exploded != "")]
    exploded = exploded.map(display_region_name)
    return (
        exploded.value_counts()
        .rename_axis("region")
        .reset_index(name="count")
        .sort_values("region")
        .reset_index(drop=True)
    )


def normalize(counts: Dict[str, int]) -> Dict[str, float]:
    total = float(sum(counts.values()) or 1.0)
    return {k: v / total for k, v in counts.items()}


HYPOTHALAMUS_THALAMUS_DISPLAY_NAME = "Hypothalamus/Thalamus"
PREGLOMERULAR_DISPLAY_NAME = "Preglomerular complex"
POSTERIOR_TUBERCULUM_DISPLAY_NAME = "Posterior tuberculum"
CAUDAL_HINDBRAIN_LEGEND_NAME = "Rhombomere 8"
HINDBRAIN_GROUP_NAME = "Hindbrain"


# Retinal arborization / optic tract merged display name
RETINAL_OPTIC_DISPLAY_NAME = "Retinal arborization/optic tract"

# Related color families make membership clear while preserving regional
# separation. Tegmentum and Other are intentionally fixed by the user.
REGION_GROUP_ORDER = (
    "Forebrain",
    "Diencephalon",
    "Midbrain",
    HINDBRAIN_GROUP_NAME,
)
HINDBRAIN_REGION_COLORS = dict(zip(
    HINDBRAIN_MECE_LEVEL_1_REGIONS,
    (
        "#FEE0D2",
        "#FCBBA1",
        "#FC9272",
        "#FB6A4A",
        "#EF3B2C",
        "#CB181D",
        "#A50F15",
        "#67000D",
    ),
))
REGION_LEGEND_ORDER = {
    "Forebrain": ["Forebrain"],
    "Diencephalon": [
        "Retinal arborization/optic tract",
        HYPOTHALAMUS_THALAMUS_DISPLAY_NAME,
        PREGLOMERULAR_DISPLAY_NAME,
        "Pretectum",
        POSTERIOR_TUBERCULUM_DISPLAY_NAME,
        "Hypocretin_Neuron_Cluster",
        "Preoptic_Area",
    ],
    "Midbrain": [
        "Tectum",
        "Tegmentum",
        "Torus_Semicircularis",
        "Oculomotor_Nucleus_nIII",
    ],
    HINDBRAIN_GROUP_NAME: list(HINDBRAIN_MECE_LEVEL_1_REGIONS),
}
DIENCEPHALON_REGION_COLORS = {
    "Retinal arborization/optic tract": "#FFFFFF",
    HYPOTHALAMUS_THALAMUS_DISPLAY_NAME: "#C8E2C2",
    PREGLOMERULAR_DISPLAY_NAME: "#96C98B",
    "Pretectum": "#5FA65B",
    POSTERIOR_TUBERCULUM_DISPLAY_NAME: "#174F2B",
    "Hypocretin_Neuron_Cluster": "#238B45",
    "Preoptic_Area": "#005824",
}
MIDBRAIN_REGION_COLORS = {
    "Tectum": "#90D5FF",
    "Tegmentum": "#6CA0DF",  # Preserve the established Tegmentum blue.
    "Torus_Semicircularis": "#2C7FB8",
    "Oculomotor_Nucleus_nIII": "#000080",
}
REGION_GROUP_PALETTES = {
    "Forebrain": ("#332288",),
    # Fixed anatomical palettes shared by all plots and Neuroglancer layers.
    "Diencephalon": tuple(DIENCEPHALON_REGION_COLORS.values()),
    "Midbrain": tuple(MIDBRAIN_REGION_COLORS.values()),
    HINDBRAIN_GROUP_NAME: tuple(HINDBRAIN_REGION_COLORS.values()),
}

# Matplotlib builds stacked bars from bottom to top.  These membership sets
# support the requested visual order (top to bottom): Diencephalon, Midbrain,
# Hindbrain, Other.
DIENCEPHALON_STACK_REGIONS = frozenset(REGION_LEGEND_ORDER["Diencephalon"])
MIDBRAIN_STACK_REGIONS = frozenset(REGION_LEGEND_ORDER["Midbrain"])
HINDBRAIN_STACK_BOTTOM_TO_TOP = (
    "Caudal_Hindbrain",
    "Rhombomere_7",
    "Rhombomere_6",
    "Rhombomere_5",
    "Rhombomere_4",
    "Rhombomere_3",
    "Rhombomere_2",
    "Rhombomere_1",
)
HINDBRAIN_STACK_RANK = {
    region: rank
    for rank, region in enumerate(HINDBRAIN_STACK_BOTTOM_TO_TOP)
}

# Desired visual order in the plot (top to bottom)
DIENCEPHALON_STACK_TOP_TO_BOTTOM = tuple(
    REGION_LEGEND_ORDER["Diencephalon"][::-1]
)
DIENCEPHALON_STACK_BOTTOM_TO_TOP = (
    DIENCEPHALON_STACK_TOP_TO_BOTTOM[::-1]
)
DIENCEPHALON_STACK_RANK = {
    region: rank
    for rank, region in enumerate(DIENCEPHALON_STACK_TOP_TO_BOTTOM)
}

MIDBRAIN_STACK_TOP_TO_BOTTOM = tuple(
    REGION_LEGEND_ORDER["Midbrain"][::-1]
)
MIDBRAIN_STACK_BOTTOM_TO_TOP = (
    MIDBRAIN_STACK_TOP_TO_BOTTOM[::-1]
)
MIDBRAIN_STACK_RANK = {
    region: rank
    for rank, region in enumerate(MIDBRAIN_STACK_TOP_TO_BOTTOM)
}

REGION_COLORS = {"Other": "#000000"}
REGION_SPECIFIC_COLORS = {
    **DIENCEPHALON_REGION_COLORS,
    **MIDBRAIN_REGION_COLORS,
    **HINDBRAIN_REGION_COLORS,
}

def display_region_name(region: str) -> str:
    region_str = coarsen_hindbrain_to_mece_level_1(region)
    key = re.sub(r"[^a-z0-9]+", "_", region_str.strip().lower()).strip("_")
    if "retinal_arborization" in key or "optic_tract_and_accessory_optic_system" in key:
        return RETINAL_OPTIC_DISPLAY_NAME
    if key in {"rostral_hypothalamus", "intermediate_hypothalamus", "ventral_thalamus", "thalamus"}:
        return HYPOTHALAMUS_THALAMUS_DISPLAY_NAME
    if key == "preglomerular_complex_approximate_area":
        return PREGLOMERULAR_DISPLAY_NAME
    if key == "posterior_tuberculum":
        return POSTERIOR_TUBERCULUM_DISPLAY_NAME
    return region_str


def build_region_color_map(regions: List[str]) -> Dict[str, str]:
    """Assign related colors to regions and keep a stable legend order."""
    colors: Dict[str, str] = {}
    counters = {group: 0 for group in REGION_GROUP_ORDER}
    for region in regions:
        if region == "Other":
            colors[region] = "#000000"
            continue
        if region in REGION_SPECIFIC_COLORS:
            colors[region] = REGION_SPECIFIC_COLORS[region]
            continue
        key = re.sub(r"[^a-z0-9]+", "_", str(region).strip().lower()).strip("_")
        if key == "forebrain":
            group = "Forebrain"
        elif key in {
            "hypothalamus", "retinal_arborization_optic_tract",
            "preglomerular_complex", "pretectum", "thalamus",
            "posterior_tuberculum",
            "hypocretin_neuron_cluster",
        }:
            group = "Diencephalon"
        elif key in {"tectum", "tegmentum", "torus_semicircularis", "oculomotor_nucleus_niii"}:
            group = "Midbrain"
        elif region in HINDBRAIN_MECE_LEVEL_1_REGIONS:
            group = HINDBRAIN_GROUP_NAME
        else:
            colors[region] = "#000000"
            continue

        if key == "tegmentum":
            colors[region] = "#6ca0df"
            continue
        palette = REGION_GROUP_PALETTES[group]
        colors[region] = palette[counters[group] % len(palette)]
        counters[group] += 1
    return colors


def legend_region_name(region: str) -> str:
    """Return the user-facing name for a plotted region."""
    if region == "Caudal_Hindbrain":
        return CAUDAL_HINDBRAIN_LEGEND_NAME
    if region.startswith("Rhombomere_"):
        return region.replace("_", " ")
    return region


def grouped_region_legend(regions: List[str], region_colors: Dict[str, str]):
    """Build legend handles with anatomical membership section headers."""
    handles = []
    labels = []
    for group in REGION_GROUP_ORDER:
        group_regions = [
            region for region in REGION_LEGEND_ORDER[group]
            if region in regions
        ]
        if not group_regions:
            continue
        handles.append(Line2D([], [], linestyle="none", marker="", color="none"))
        labels.append(group)
        for region in group_regions:
            handles.append(
                plt.Rectangle(
                    (0, 0),
                    1,
                    1,
                    facecolor=region_colors[region],
                    edgecolor="black",
                    linewidth=0.3,
                )
            )
            labels.append(f"  {legend_region_name(region)}")
    if "Other" in regions:
        handles.append(Line2D([], [], linestyle="none", marker="", color="none"))
        labels.append("Other")
        handles.append(plt.Rectangle((0, 0), 1, 1, color=region_colors["Other"]))
        labels.append("  Other")
    return handles, labels

# Neuroglancer helpers

def um_to_ng_vox(x_um, y_um, z_um, ng_vox_m=NG_VOX_M):
    x_vox = (np.asarray(x_um, dtype=float) * 1e-6) / ng_vox_m[0]
    y_vox = (np.asarray(y_um, dtype=float) * 1e-6) / ng_vox_m[1]
    z_vox = (np.asarray(z_um, dtype=float) * 1e-6) / ng_vox_m[2]
    return x_vox, y_vox, z_vox


def df_for_role(df_src: pd.DataFrame, role: str) -> pd.DataFrame:
    if role == "presynapse":
        x_um, y_um, z_um = df_src["pre_x_um"], df_src["pre_y_um"], df_src["pre_z_um"]
        sid = df_src.get("id", df_src.index).astype(str)
        partner = df_src.get("post_pt_root_id", pd.Series([""] * len(df_src)))
    else:
        x_um, y_um, z_um = df_src["post_x_um"], df_src["post_y_um"], df_src["post_z_um"]
        sid = df_src.get("id", df_src.index).astype(str)
        partner = df_src.get("pre_pt_root_id", pd.Series([""] * len(df_src)))

    x_vox, y_vox, z_vox = um_to_ng_vox(x_um, y_um, z_um)
    return pd.DataFrame({
        "x": x_vox,
        "y": y_vox,
        "z": z_vox,
        "synapse_id": sid,
        "partner_segment": partner.astype(str),
    })


def rows_to_annotations(df: pd.DataFrame, role: str, soma_id: int) -> list:
    annotations = []
    tag = "pre" if role == "presynapse" else "post"
    for i, row in df.reset_index(drop=True).iterrows():
        annotations.append({
            "type": "point",
            "id": f"{tag}_{soma_id}_{i + 1}",
            "point": [float(row["x"]), float(row["y"]), float(row["z"])],
            "description": (
                f"cell={soma_id} | {'output' if tag == 'pre' else 'input'} | "
                f"synapse_id={row['synapse_id']} | partner={row['partner_segment']}"
            ),
        })
    return annotations


def build_layer(annotations: list, layer_name: str, color_hex: str | None = None) -> dict:
    layer = {
        "type": "annotation",
        "source": {
            "url": "local://annotations",
            "transform": {
                "outputDimensions": {
                    "x": [NG_VOX_M[0], "m"],
                    "y": [NG_VOX_M[1], "m"],
                    "z": [NG_VOX_M[2], "m"],
                }
            },
        },
        "tool": "annotatePoint",
        "tab": "annotations",
        "name": layer_name,
        "annotations": annotations,
    }
    if color_hex:
        layer["annotationColor"] = color_hex
    return layer


def build_midline_colored_segmentation_layer(
    soma_ids: List[int],
    nuclei_records: List[Tuple[int, Tuple[float, float, float]]] | None,
    layer_name: str = "fish1_v250915_midline_distance",
) -> dict | None:
    """Build a Neuroglancer segmentation layer with cells colored by left-right distance from the midline."""
    if not nuclei_records:
        return None

    nucleus_pos_by_soma = {
        int(soma_id): pos
        for soma_id, pos in nuclei_records
        if pos is not None
    }

    cells_with_pos = [int(soma_id) for soma_id in soma_ids if int(soma_id) in nucleus_pos_by_soma]
    if not cells_with_pos:
        return None

    # Color by distance to the midline, normalized independently on each side.
    # On each side of the brain:
    #   closest cell to the midline = blue
    #   farthest cell from the midline = red
    # This makes the left and right sides use the same full color extent.
    signed_distances_px = np.asarray(
        [float(nucleus_pos_by_soma[soma_id][1]) - float(MIDLINE_POSITION[1]) for soma_id in cells_with_pos],
        dtype=float,
    )
    distances_um = np.abs(signed_distances_px) * float(VOX_SIZE_UM_EM[1])

    if not np.any(np.isfinite(distances_um)):
        return None

    norm_distances = np.full_like(distances_um, 0.5, dtype=float)
    for side_mask in [signed_distances_px < 0, signed_distances_px >= 0]:
        valid = side_mask & np.isfinite(distances_um)
        if not np.any(valid):
            continue

        side_distances = distances_um[valid]
        side_min = float(np.nanmin(side_distances))
        side_max = float(np.nanmax(side_distances))

        if np.isclose(side_min, side_max):
            norm_distances[valid] = 0.0
        else:
            norm_distances[valid] = np.clip(
                (side_distances - side_min) / (side_max - side_min),
                0.0,
                1.0,
            )

    cmap = matplotlib.colormaps.get_cmap("coolwarm")
    segment_colors = {}
    for soma_id, value in zip(cells_with_pos, norm_distances):
        rgba = cmap(float(value))
        rgb = tuple(int(round(255 * channel)) for channel in rgba[:3])
        segment_colors[str(soma_id)] = "#{:02x}{:02x}{:02x}".format(*rgb)

    return {
        "type": "segmentation",
        "source": SEGMENTATION_SOURCE,
        "tab": "segments",
        "segments": [str(soma_id) for soma_id in cells_with_pos],
        "segmentColors": segment_colors,
        "name": layer_name,
    }


FUNCTIONAL_LAYER_COLORS = {
    "pag_contra": "#ff00ff",  # magenta
    "pag_ipsi": "#000000",    # black
    "ni": "#00C141",          # green
}


def classify_functional_projection(functional_id: object, projection: object) -> str | None:
    """Assign a metadata category for the PAG/NI Neuroglancer layer.

    NI takes precedence over PAG if both strings appear in ``functional_id``.
    PAG cells require a recognizable ipsilateral/contralateral projection.
    """
    functional = "" if pd.isna(functional_id) else str(functional_id).strip().lower()
    projection_text = "" if pd.isna(projection) else str(projection).strip().lower()

    # Keep this deliberately tolerant of IDs such as
    # fish1_v250915_184922_ni_rel3a_#2.
    if re.search(r"(?:^|[_\s-])ni(?:[_\s-]|$)", functional) or "_ni_" in functional:
        return "ni"
    if "pag" not in functional:
        return None
    if "contra" in projection_text:
        return "pag_contra"
    if "ipsi" in projection_text:
        return "pag_ipsi"
    return None


def is_pag_functional_id(functional_id: object) -> bool:
    """Return whether a metadata functional ID denotes a PAG neuron."""
    text = "" if pd.isna(functional_id) else str(functional_id).strip().lower()
    return "pag" in text


def build_functional_projection_segmentation_layer(
    soma_ids: List[int],
    functional_by_soma: Dict[int, str],
    projection_by_soma: Dict[int, str],
    layer_name: str = "PAG contra / ipsi / NI neurons",
) -> dict | None:
    """Build a colored segmentation layer from functional/projection metadata."""
    selected: list[int] = []
    segment_colors: dict[str, str] = {}
    category_counts = {category: 0 for category in FUNCTIONAL_LAYER_COLORS}
    skipped_pag_projection = 0

    for soma_id in soma_ids:
        category = classify_functional_projection(
            functional_by_soma.get(int(soma_id), ""),
            projection_by_soma.get(int(soma_id), ""),
        )
        if category is None:
            functional = str(functional_by_soma.get(int(soma_id), "")).lower()
            if "pag" in functional:
                skipped_pag_projection += 1
            continue
        selected.append(int(soma_id))
        segment_colors[str(int(soma_id))] = FUNCTIONAL_LAYER_COLORS[category]
        category_counts[category] += 1

    if not selected:
        return None

    print(
        "[functional layer] "
        f"PAG contra={category_counts['pag_contra']}, "
        f"PAG ipsi={category_counts['pag_ipsi']}, "
        f"NI={category_counts['ni']}"
    )
    if skipped_pag_projection:
        print(
            f"[functional layer] skipped {skipped_pag_projection} PAG row(s) "
            "without ipsi/contra projection metadata"
        )

    return {
        "type": "segmentation",
        "source": SEGMENTATION_SOURCE,
        "tab": "segments",
        "segments": [str(soma_id) for soma_id in selected],
        "segmentColors": segment_colors,
        "name": layer_name,
    }


def write_layers_state_snippet(path: Path, layers: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"layers": layers}, indent=2))

# Metadata helpers

def parse_nucleus_pixel_field(val) -> Tuple[float, float, float] | None:
    if pd.isna(val):
        return None
    s = str(val).strip().strip("[]()")
    parts = s.split(",")
    if len(parts) < 3:
        return None
    try:
        return float(parts[0]), float(parts[1]), float(parts[2])
    except Exception:
        return None


def read_metadata(args) -> tuple:
    """Load neuron annotations and retain exact integer segmentation IDs."""
    csv_path = Path(args.root_ids_csv or DEFAULT_CELLS_CSV)
    if not csv_path.is_file():
        raise SystemExit(f"Neuron metadata not found: {csv_path}")
    df = pd.read_csv(csv_path, dtype={"soma_id": "string"})
    required = {"soma_id", "functional_id"}
    if not required.issubset(df.columns):
        raise SystemExit(f"Metadata must contain columns: {sorted(required)}")
    df = df.dropna(subset=["soma_id"]).copy()
    if df["soma_id"].duplicated().any():
        raise SystemExit("Metadata must contain one row per soma_id.")

    if args.root_ids:
        if args.complete_only:
            print("[WARN] --complete-only is ignored when --root-ids are supplied directly.")
        requested = {int(value) for value in args.root_ids}
        available = {int(value) for value in df["soma_id"]}
        if requested - available:
            raise SystemExit(f"Root IDs missing from metadata: {sorted(requested - available)}")
        df = df[df["soma_id"].map(int).isin(requested)].copy()
    elif args.complete_only:
        if "reconstruction" not in df.columns:
            raise SystemExit("--complete-only requires a reconstruction column.")
        before = len(df)
        eligible = df["reconstruction"].astype(str).str.contains(
            r"\b(?:complete|sufficient)\b", case=False, na=False, regex=True)
        df = df[eligible].copy()
        print(f"[filter] complete-only: kept {len(df)} / {before} rows")

    soma_ids = [int(value) for value in df["soma_id"]]
    nuclei_records = None
    if "nucleus_position" in df.columns:
        positions = [parse_nucleus_pixel_field(value) for value in df["nucleus_position"]]
        nuclei_records = [(int(sid), pos) for sid, pos in zip(df["soma_id"], positions)
                          if pos is not None]

    def annotations(column):
        if column not in df.columns:
            return {}
        return {int(sid): str(value).strip()
                for sid, value in zip(df["soma_id"], df[column]) if pd.notna(value)}

    projection = {sid: value.lower() for sid, value in annotations("projection classifier").items()}
    functional = annotations("functional_id")
    incomplete = {int(sid): not bool(re.search(r"\b(?:complete|sufficient)\b", str(status),
                                              flags=re.IGNORECASE))
                  for sid, status in zip(df["soma_id"], df["reconstruction"])} if "reconstruction" in df else {}
    return soma_ids, nuclei_records, projection, incomplete, functional

# Plot helpers

def order_cells_by_left_right_position(
    cells: List[int],
    nuclei_records: List[Tuple[int, Tuple[float, float, float]]] | None,
) -> List[int]:
    if not nuclei_records:
        return list(range(len(cells)))

    nucleus_pos_by_soma = {int(soma_id): pos for soma_id, pos in nuclei_records if pos is not None}

    def sort_key(idx: int):
        soma_id = int(cells[idx])
        pos = nucleus_pos_by_soma.get(soma_id)
        if pos is None:
            return (1, 0.0, soma_id)
        # In this dataset, y - midline is positive on the anatomical left.
        # Sort from positive to negative to display the brain from left to right.
        left_right_position = float(pos[1]) - float(MIDLINE_POSITION[1])
        return (0, -left_right_position, soma_id)

    return sorted(range(len(cells)), key=sort_key)


def get_left_right_positions(
    cells: List[int],
    nuclei_records: List[Tuple[int, Tuple[float, float, float]]] | None,
) -> np.ndarray:
    nucleus_pos_by_soma = {int(soma_id): pos for soma_id, pos in (nuclei_records or []) if pos is not None}
    values = []
    for cell in cells:
        cell_id = int(cell)
        if cell_id in nucleus_pos_by_soma:
            values.append(float(nucleus_pos_by_soma[cell_id][1]) - float(MIDLINE_POSITION[1]))
        else:
            values.append(np.nan)
    return np.asarray(values, dtype=float)


def make_xticklabels(cells: List[int], left_right_positions: np.ndarray, mode: str) -> List[str]:
    """Return signed soma distance from the midline, always in µm.

    Cell/root IDs are deliberately never used as x-axis labels.  ``mode`` is
    retained in the signature for compatibility with older callers, but is
    ignored so all exported figures use the same anatomical coordinate.
    """
    labels = []
    for pos_px in left_right_positions:
        if np.isfinite(pos_px):
            pos_um = float(pos_px) * float(VOX_SIZE_UM_EM[1])
            labels.append(f"{pos_um:.1f}")
        else:
            labels.append("NA")
    return labels


def x_positions_with_midline_gap(
    left_right_positions: np.ndarray,
    gap_width: float = 1.0,
) -> Tuple[np.ndarray, float | None]:
    """Return bar positions with a gap between left (+) and right (-) brain sides."""
    x = np.arange(len(left_right_positions), dtype=float)
    finite = np.isfinite(left_right_positions)
    left_indices = np.flatnonzero(finite & (left_right_positions > 0))
    right_indices = np.flatnonzero(finite & (left_right_positions <= 0))

    if len(left_indices) == 0 or len(right_indices) == 0:
        return x, None

    left_end = int(left_indices[-1])
    right_start = int(right_indices[0])
    x[right_start:] += gap_width
    midline_x = float((x[left_end] + x[right_start]) / 2.0)
    return x, midline_x


def group_small_regions_for_plot(
    mat_in: np.ndarray,
    mat_out: np.ndarray,
    regions: List[str],
    threshold: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Group each cell's nonzero region proportions below threshold into Other."""
    grouped_in = mat_in.copy()
    grouped_out = mat_out.copy()

    small_in = (grouped_in > 0) & (grouped_in < threshold)
    small_out = (grouped_out > 0) & (grouped_out < threshold)
    other_in = np.where(small_in, grouped_in, 0.0).sum(axis=1)
    other_out = np.where(small_out, grouped_out, 0.0).sum(axis=1)
    grouped_in[small_in] = 0.0
    grouped_out[small_out] = 0.0

    keep = np.any(grouped_in > 0, axis=0) | np.any(grouped_out > 0, axis=0)
    plot_regions = [region for region, keep_region in zip(regions, keep) if keep_region]
    grouped_in = grouped_in[:, keep]
    grouped_out = grouped_out[:, keep]

    if np.any(other_in > 0) or np.any(other_out > 0):
        grouped_in = np.column_stack([grouped_in, other_in])
        grouped_out = np.column_stack([grouped_out, other_out])
        plot_regions.append("Other")

    return grouped_in, grouped_out, plot_regions


def reorder_region_stacks(
    mat_in: np.ndarray,
    mat_out: np.ndarray,
    regions: List[str],
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Order stack columns bottom-to-top as Other, hindbrain, midbrain, diencephalon."""
    if mat_in.ndim != 2 or mat_out.ndim != 2:
        raise ValueError("Region profile matrices must be two-dimensional.")
    if mat_in.shape[1] != len(regions) or mat_out.shape[1] != len(regions):
        raise ValueError("Region labels must match the matrix column counts.")

    def stack_rank(region: str) -> Tuple[int, int]:
        if region == "Other":
            return 0, 0
        # Keep any future Forebrain or unclassified labels next to Other while
        # preserving them as distinct columns and leaving their values intact.
        if region in HINDBRAIN_MECE_LEVEL_1_REGIONS:
            return 2, HINDBRAIN_STACK_RANK[region]
        if region in MIDBRAIN_STACK_REGIONS:
            return 3, MIDBRAIN_STACK_RANK[region]
        if region in DIENCEPHALON_STACK_REGIONS:
            return 4, DIENCEPHALON_STACK_RANK[region]
        return 1, 0

    column_order = sorted(
        range(len(regions)),
        key=lambda index: (*stack_rank(regions[index]), index),
    )
    return (
        mat_in[:, column_order],
        mat_out[:, column_order],
        [regions[index] for index in column_order],
    )




def projection_category(value: object) -> str:
    text = "" if pd.isna(value) else str(value).strip().lower()
    if "ipsi" in text:
        return "ipsilateral"
    if "contra" in text:
        return "contralateral"
    return "unassigned"


def order_cells_by_projection(
    cells: List[int],
    nuclei_records: List[Tuple[int, Tuple[float, float, float]]] | None,
    projection_by_soma: Dict[int, str],
    gap_width: float = 1.5,
) -> Tuple[List[int], np.ndarray, List[Tuple[str, float, float]]]:
    """Group bars by projection while preserving left-right order within groups."""
    ordered_indices: list[int] = []
    blocks: list[Tuple[str, float, float]] = []
    x_values: list[float] = []
    cursor = 0.0
    for category in ("ipsilateral", "contralateral", "unassigned"):
        group_indices = [
            index
            for index, soma_id in enumerate(cells)
            if projection_category(projection_by_soma.get(int(soma_id), "")) == category
        ]
        if not group_indices:
            continue
        group_cells = [cells[index] for index in group_indices]
        local_order = order_cells_by_left_right_position(group_cells, nuclei_records)
        group_indices = [group_indices[index] for index in local_order]
        start = cursor
        for index in group_indices:
            ordered_indices.append(index)
            x_values.append(cursor)
            cursor += 1.0
        end = cursor - 1.0
        blocks.append((category, start, end))
        cursor += gap_width
    return ordered_indices, np.asarray(x_values, dtype=float), blocks


def plot_input_output_stacked_side_by_side(
    mat_in: np.ndarray,
    mat_out: np.ndarray,
    cells: List[int],
    regions: List[str],
    region_colors: Dict[str, str],
    nuclei_records: List[Tuple[int, Tuple[float, float, float]]] | None,
    projection_by_soma: Dict[int, str],
    incomplete_axon_by_soma: Dict[int, bool],
    out_path: Path,
    x_label_mode: str,
) -> None:
    """Plot input/output stacks with ipsilateral and contralateral cells separated."""
    fig_width_in = 147.8 / 25.4
    fig_height_in = 3.2

    order, _, projection_blocks = order_cells_by_projection(
        cells,
        nuclei_records,
        projection_by_soma,
    )
    if nuclei_records is None:
        print("[WARN] No nucleus positions available; keeping input order within projection groups.")

    cells_plot = [cells[i] for i in order]
    lr_positions = get_left_right_positions(cells_plot, nuclei_records)
    # Preserve the anatomical midline gap used in the original combined plot.
    x, midline_x = x_positions_with_midline_gap(lr_positions)
    xticklabels = make_xticklabels(cells_plot, lr_positions, x_label_mode)

    mat_in_plot = mat_in[order, :]
    mat_out_plot = mat_out[order, :]
    incomplete_axon_indices = np.asarray(
        [
            i
            for i, soma_id in enumerate(cells_plot)
            if incomplete_axon_by_soma.get(int(soma_id), False)
        ],
        dtype=int,
    )

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(147.8/25.4, 3.2),
        sharey=True,
    )

    fig.subplots_adjust(
        left=0.08,
        right=0.73,   # reserve ~27% of the figure width for the legend
        bottom=0.20,
        top=0.90,
        wspace=0.15,
    )
    if len(projection_blocks) == 1:
        fig.suptitle(projection_blocks[0][0].capitalize(), fontsize=7, y=0.98)

    for ax, mat_plot, title in zip(
        axes,
        [mat_in_plot, mat_out_plot],
        ["Inputs", "Outputs"],
    ):
        base = np.zeros(len(cells_plot))
        for region in regions:
            j = regions.index(region)
            ax.bar(
                x,
                mat_plot[:, j],
                bottom=base,
                color=region_colors[region],
                align="center",
                edgecolor="black",
                linewidth=0.2,
                width=0.85,
            )
            base += mat_plot[:, j]

        # Mark incomplete axons above the output bars.  The input panel is
        # intentionally left unmarked because this annotation concerns the
        # completeness of the traced axon/output compartment.
        if title == "Outputs" and len(incomplete_axon_indices) > 0:
            ax.scatter(
                x[incomplete_axon_indices],
                np.full(len(incomplete_axon_indices), 1.025),
                marker="*",
                s=22,
                color="#D55E00",
                zorder=5,
            )

        for block_index, (category, start, end) in enumerate(projection_blocks):
            if block_index > 0:
                boundary = start - 0.75
                ax.axvline(
                    boundary,
                    color="#777777",
                    linestyle=(0, (2, 2)),
                    linewidth=0.6,
                    zorder=4,
                )
            if len(projection_blocks) > 1:
                label = category.capitalize()
                ax.text(
                    (start + end) / 2.0,
                    1.025,
                    label,
                    ha="center",
                    va="bottom",
                    transform=ax.get_xaxis_transform(),
                    fontsize=5,
                )

        if midline_x is not None:
            ax.axvline(
                midline_x,
                color="#555555",
                linestyle=(0, (3, 2)),
                linewidth=0.8,
                zorder=4,
            )

        ax.set_title(title, fontsize=7)
        ax.set_xticks(x)
        ax.set_xticklabels(
            xticklabels,
            rotation=90,
            ha="center",
            fontsize=5,
        )
        ax.tick_params(axis="y", labelsize=7)
        ax.set_ylim(0, 1.075)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("Proportion", fontsize=7)
    fig.supxlabel("Signed distance from midline (µm; + = left)", fontsize=7)
    handles, labels = grouped_region_legend(regions, region_colors)
    if len(incomplete_axon_indices) > 0:
        handles.append(
            Line2D(
                [0],
                [0],
                marker="*",
                linestyle="None",
                markerfacecolor="#D55E00",
                markeredgecolor="#D55E00",
                markersize=6,
            )
        )
        labels.append("Incomplete axon")
    fig.legend(
        handles,
        labels,
        title="Regions",
        loc="center left",
        bbox_to_anchor=(0.75, 0.55),
        bbox_transform=fig.transFigure,
        frameon=False,
        fontsize=5,
        title_fontsize=5,
        handlelength=1.0,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=PLOT_DPI, bbox_inches="tight", format="pdf")
    plt.close()
    print(f"[saved] {out_path}")

# Global analysis

def run_global_analysis(
    soma_ids: List[int],
    nuclei_records: List[Tuple[int, Tuple[float, float, float]]] | None,
    projection_by_soma: Dict[int, str],
    incomplete_axon_by_soma: Dict[int, bool],
    functional_by_soma: Dict[int, str],
    synaptome_x_labels: str = "distance",
) -> None:
    label, region_names, mask_vox_um = load_label_volume(H5_MASKS)
    print(f"[label_volume] shape={label.shape}, regions={len(region_names)}, voxel_size_um={mask_vox_um}")

    # The synaptome summaries and synapse layers are intentionally restricted
    # to PAG neurons.  NI cells remain available for the separate functional
    # Neuroglancer layer, but are not included in this synaptome analysis.
    pag_soma_ids = [
        int(soma_id)
        for soma_id in soma_ids
        if is_pag_functional_id(functional_by_soma.get(int(soma_id), ""))
    ]
    if not pag_soma_ids:
        print("[ERROR] No PAG neurons were identified from functional_id; synaptome analysis skipped.")
        return
    print(f"[synaptome] restricting analysis to {len(pag_soma_ids)} PAG neuron(s)")

    global_presyn_ann = []
    global_postsyn_ann = []
    inputs_profiles = []
    outputs_profiles = []
    counts_rows = []

    for soma_id in pag_soma_ids:
        in_path = cell_csv_inputs(soma_id)
        out_path = cell_csv_outputs(soma_id)
        if not in_path.exists() or not out_path.exists():
            print(f"[WARN] Missing per-cell CSVs for soma_id {soma_id}; skipping.", file=sys.stderr)
            continue

        df_inputs = pd.read_csv(in_path)
        df_outputs = pd.read_csv(out_path)

        mapped_inputs = map_df_to_regions(df_inputs, "post", label, region_names, mask_vox_um)
        mapped_outputs = map_df_to_regions(df_outputs, "pre", label, region_names, mask_vox_um)

        in_counts_df = counts_from_region_lists(mapped_inputs["post_regions"])
        out_counts_df = counts_from_region_lists(mapped_outputs["pre_regions"])

        in_probs = normalize(dict(zip(in_counts_df["region"], in_counts_df["count"]))) if len(in_counts_df) else {}
        out_probs = normalize(dict(zip(out_counts_df["region"], out_counts_df["count"]))) if len(out_counts_df) else {}

        inputs_profiles.append({"ROOT_ID": soma_id, **in_probs})
        outputs_profiles.append({"ROOT_ID": soma_id, **out_probs})
        counts_rows.append({"ROOT_ID": soma_id, "n_inputs": len(df_inputs), "n_outputs": len(df_outputs)})

        pre_for_ng = df_for_role(df_outputs, role="presynapse")
        post_for_ng = df_for_role(df_inputs, role="postsynapse")
        global_presyn_ann.extend(rows_to_annotations(pre_for_ng, role="presynapse", soma_id=soma_id))
        global_postsyn_ann.extend(rows_to_annotations(post_for_ng, role="postsynapse", soma_id=soma_id))

    if not inputs_profiles and not outputs_profiles:
        print("[ERROR] No per-cell synaptome profiles found.")
        return

    ng_all_dir = RESULTS_DIR / "ng_layers"
    ng_all_dir.mkdir(parents=True, exist_ok=True)
    global_pre_layer = build_layer(global_presyn_ann, layer_name="presyn_all_cells", color_hex="#ff2b00")
    global_post_layer = build_layer(global_postsyn_ann, layer_name="postsyn_all_cells", color_hex="#2b00ff")
    (ng_all_dir / "inputs_all_presyn.json").write_text(json.dumps(global_pre_layer, indent=2))
    (ng_all_dir / "outputs_all_postsyn.json").write_text(json.dumps(global_post_layer, indent=2))
    write_layers_state_snippet(ng_all_dir / "all_cells_inputs_outputs_layers.json", [global_post_layer, global_pre_layer])
    print(f"[OK] Global NG layers written: {ng_all_dir}")

    functional_layer = build_functional_projection_segmentation_layer(
        soma_ids,
        functional_by_soma,
        projection_by_soma,
    )
    if functional_layer is not None:
        functional_layer_path = ng_all_dir / "pag_ni_projection_classified_cells.json"
        functional_layer_path.write_text(json.dumps(functional_layer, indent=2))
        print(f"[saved] {functional_layer_path}")
    else:
        print(
            "[INFO] No PAG/NI functional metadata matched; "
            "skipping functional projection layer."
        )

    midline_seg_layer = build_midline_colored_segmentation_layer(soma_ids, nuclei_records)
    if midline_seg_layer is not None:
        midline_seg_path = ng_all_dir / "all_cells_colored_by_midline_distance.json"
        midline_seg_path.write_text(json.dumps(midline_seg_layer, indent=2))
        print(f"[saved] {midline_seg_path}")
    else:
        print("[INFO] No nucleus positions available; skipping midline-colored segmentation layer.")

    inputs_df = pd.DataFrame(inputs_profiles).fillna(0.0)
    outputs_df = pd.DataFrame(outputs_profiles).fillna(0.0)
    inputs_df.to_csv(RESULTS_DIR / "inputs_summary.csv", index=False)
    outputs_df.to_csv(RESULTS_DIR / "outputs_summary.csv", index=False)
    print(f"[saved] {RESULTS_DIR / 'inputs_summary.csv'}")
    print(f"[saved] {RESULTS_DIR / 'outputs_summary.csv'}")

    in_regions = [c for c in inputs_df.columns if c != "ROOT_ID"]
    out_regions = [c for c in outputs_df.columns if c != "ROOT_ID"]
    all_regions = sorted(set(in_regions) | set(out_regions))
    if not all_regions:
        print("[WARN] No mapped brain-region columns found; skipping stacked plots.")
        return

    def ensure_cols(df: pd.DataFrame) -> pd.DataFrame:
        for region in all_regions:
            if region not in df.columns:
                df[region] = 0.0
        return df[["ROOT_ID"] + all_regions]

    inputs_df = ensure_cols(inputs_df)
    outputs_df = ensure_cols(outputs_df)
    all_cells = sorted(set(inputs_df["ROOT_ID"].tolist()) | set(outputs_df["ROOT_ID"].tolist()))

    def matrix_from_summary(df: pd.DataFrame) -> np.ndarray:
        df_idx = (
            df.groupby("ROOT_ID", as_index=True)[all_regions]
            .sum()
            .reindex(all_cells)
            .fillna(0.0)
        )
        mat = df_idx[all_regions].to_numpy(dtype=float)
        row_sums = mat.sum(axis=1, keepdims=True)
        return np.divide(mat, row_sums, out=np.zeros_like(mat), where=row_sums > 0)

    mat_in = matrix_from_summary(inputs_df)
    mat_out = matrix_from_summary(outputs_df)
    # Keep individual regions for granularity, use related shades within each
    # broad anatomical class, and pool contributions below 5% into black Other.
    plot_mat_in, plot_mat_out, plot_regions = group_small_regions_for_plot(
        mat_in,
        mat_out,
        all_regions,
        threshold=0.05,
    )
    plot_mat_in, plot_mat_out, plot_regions = reorder_region_stacks(
        plot_mat_in,
        plot_mat_out,
        plot_regions,
    )
    region_colors = build_region_color_map(plot_regions)

    # Keep ipsilateral and contralateral populations in separate figures so
    # their connectivity profiles are not visually mixed.  Cells without a
    # projection annotation are exported separately when present.
    for category in ("ipsilateral", "contralateral", "unassigned"):
        keep = np.asarray(
            [
                projection_category(projection_by_soma.get(int(soma_id), "")) == category
                for soma_id in all_cells
            ],
            dtype=bool,
        )
        if not np.any(keep):
            continue

        category_cells = [soma_id for soma_id, include in zip(all_cells, keep) if include]
        category_mat_in = plot_mat_in[keep, :]
        category_mat_out = plot_mat_out[keep, :]
        category_label = category.capitalize()
        print(f"[plot] {category_label}: {len(category_cells)} cells")
        plot_input_output_stacked_side_by_side(
            category_mat_in,
            category_mat_out,
            category_cells,
            plot_regions,
            region_colors,
            nuclei_records,
            projection_by_soma,
            incomplete_axon_by_soma,
            figure_pdf_path(f"inputs_outputs_region_proportions_stacked_{category}"),
            "distance",
        )

    if counts_rows:
        syn_df = pd.DataFrame(counts_rows).sort_values("ROOT_ID")
        syn_df.to_csv(RESULTS_DIR / "synapse_counts_per_cell.csv", index=False)
        print(f"[saved] {RESULTS_DIR / 'synapse_counts_per_cell.csv'}")

        # All analyzed dendrites are complete, so summarize inputs across every
        # loaded cell. Restrict output summaries to cells without an incomplete
        # axon annotation.
        # Use 0.0 SD when a summary contains only one cell.
        output_summary_syn_df = syn_df[
            ~syn_df["ROOT_ID"].map(
                lambda soma_id: incomplete_axon_by_soma.get(int(soma_id), False)
            )
        ]
        input_mean = syn_df["n_inputs"].mean()
        input_sd = syn_df["n_inputs"].std(ddof=1) if len(syn_df) > 1 else 0.0
        if output_summary_syn_df.empty:
            print(
                "[WARN] No cells without an incomplete axon annotation; "
                "output synapse-count mean ± SD is unavailable."
            )
            print(
                "[synapse counts] "
                f"inputs: {input_mean:.1f} ± {input_sd:.1f} "
                f"(mean ± SD across all {len(syn_df)} cells)"
            )
        else:
            output_mean = output_summary_syn_df["n_outputs"].mean()
            output_sd = (
                output_summary_syn_df["n_outputs"].std(ddof=1)
                if len(output_summary_syn_df) > 1
                else 0.0
            )
            print(
                "[synapse counts] "
                f"inputs: {input_mean:.1f} ± {input_sd:.1f}; "
                f"outputs: {output_mean:.1f} ± {output_sd:.1f} "
                f"(inputs: all {len(syn_df)} cells; "
                f"outputs: {len(output_summary_syn_df)} cells without incomplete axons)"
            )

        syn_cells = syn_df["ROOT_ID"].astype(int).tolist()
        count_order = order_cells_by_left_right_position(syn_cells, nuclei_records)
        syn_plot = syn_df.iloc[count_order].reset_index(drop=True)
        cells_plot = syn_plot["ROOT_ID"].astype(int).tolist()
        lr_positions = get_left_right_positions(cells_plot, nuclei_records)
        idx, midline_x = x_positions_with_midline_gap(lr_positions)
        xticklabels = make_xticklabels(cells_plot, lr_positions, synaptome_x_labels)
        projection_marker_styles = {
            "ipsilateral": ("o", "#0072B2", "Ipsilateral"),
            "contralateral": ("^", "#CC79A7", "Contralateral"),
            "unassigned": ("s", "#777777", "Unassigned"),
        }
        projection_categories = [
            projection_category(projection_by_soma.get(int(soma_id), ""))
            for soma_id in cells_plot
        ]
        width = 0.78
        figure_width_mm = 47.8
        figure_height_mm = 48.6
        fig, ax = plt.subplots(
            figsize=(figure_width_mm / 25.4, figure_height_mm / 25.4)
        )
        fig.subplots_adjust(
            left=0.12,
            right=0.98,
            bottom=0.36,
            top=0.92,
        )
        ax.bar(
            idx,
            syn_plot["n_inputs"].to_numpy(int),
            width,
            label="Inputs",
            color="#4D4D4D",
            edgecolor="black",
            linewidth=0.3,
        )
        ax.bar(
            idx,
            syn_plot["n_outputs"].to_numpy(int),
            width,
            bottom=syn_plot["n_inputs"].to_numpy(int),
            label="Outputs",
            color="#BDBDBD",
            edgecolor="black",
            linewidth=0.3,
        )
        total_counts = (
            syn_plot["n_inputs"].to_numpy(int)
            + syn_plot["n_outputs"].to_numpy(int)
        )
        marker_offset = max(float(total_counts.max()) * 0.03, 1.0)
        # Use a symbol above each bar to identify its projection class rather
        # than adding projection text to the crowded x-axis labels.
        for category, (marker, color, legend_label) in projection_marker_styles.items():
            category_indices = np.asarray(
                [i for i, value in enumerate(projection_categories) if value == category],
                dtype=int,
            )
            if len(category_indices) == 0:
                continue
            ax.scatter(
                idx[category_indices],
                total_counts[category_indices] + marker_offset,
                marker=marker,
                s=22,
                color=color,
                label=legend_label,
                zorder=5,
            )
        incomplete_indices = np.array(
            [
                i
                for i, soma_id in enumerate(cells_plot)
                if incomplete_axon_by_soma.get(int(soma_id), False)
            ],
            dtype=int,
        )
        if len(incomplete_indices) > 0:
            ax.scatter(
                idx[incomplete_indices],
                total_counts[incomplete_indices] + 2 * marker_offset,
                marker="*",
                s=22,
                color="#D55E00",
                label="Incomplete axon",
                zorder=5,
            )
        ax.set_ylim(top=float(total_counts.max()) + 3 * marker_offset)
        if midline_x is not None:
            ax.axvline(
                midline_x,
                color="#555555",
                linestyle=(0, (3, 2)),
                linewidth=0.8,
                zorder=4,
            )
        ax.set_xticks(
            idx,
            xticklabels,
            rotation=90 if synaptome_x_labels == "distance" else 45,
            ha="center" if synaptome_x_labels == "distance" else "right",
            fontsize=5,
        )
        ax.set_xlabel("Signed distance from midline (µm; + = left)", fontsize=7)
        ax.set_ylabel("Synapses per cell", fontsize=7)
        ax.set_title("Number of input and output synapses per cell", fontsize=5)
        ax.tick_params(axis="y", labelsize=5)
        ax.legend(fontsize=5, frameon=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # Do not use bbox_inches="tight" here: it changes the physical PDF
        # page size after figsize has been set. The exported page remains
        # exactly 47.8 mm wide.
        fig.savefig(
            figure_pdf_path("synapse_counts_per_cell"),
            dpi=PLOT_DPI,
            format="pdf",
        )
        plt.close(fig)
        print(f"[saved] {figure_pdf_path('synapse_counts_per_cell')}")

    if nuclei_records:
        xs_px = np.array([p[1][0] for p in nuclei_records], dtype=float)
        ys_px = np.array([p[1][1] for p in nuclei_records], dtype=float)
        zs_px = np.array([p[1][2] for p in nuclei_records], dtype=float)
        xs_um, ys_um, zs_um = px_to_um_arrays(xs_px, ys_px, zs_px, VOX_SIZE_UM_EM)
        xi, yi, zi = um_to_voxel(xs_um, ys_um, zs_um, mask_vox_um)

        X, Y, Z = label.shape
        nuc_regions = []
        for a, b, c in zip(xi, yi, zi):
            a, b, c = int(a), int(b), int(c)
            if 0 <= a < X and 0 <= b < Y and 0 <= c < Z:
                lab = int(label[a, b, c])
                nuc_regions.append(region_names[lab - 1] if lab > 0 else None)
            else:
                nuc_regions.append(None)

        nuclei_assign_df = pd.DataFrame({
            "soma_id": [sid for sid, _ in nuclei_records],
            "nucleus_region": nuc_regions,
        })
        nuclei_assign_df.to_csv(RESULTS_DIR / "nuclei_region_assignments.csv", index=False)
        print(f"[saved] {RESULTS_DIR / 'nuclei_region_assignments.csv'}")
    else:
        print("[INFO] No nucleus positions provided; skipping nuclei-region assignments.")

# Command-line interface

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Run global synaptome analysis from per-cell outputs.")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("--root-ids", nargs="+", type=int, help="List of soma IDs to analyze.")
    group.add_argument("--root-ids-csv", type=str, default=None, help=f"Metadata CSV. Default: {DEFAULT_CELLS_CSV}")
    parser.add_argument("--atlas", type=Path, default=H5_MASKS,
                        help="MECE HDF5 file containing label_volume and region_names.")
    parser.add_argument("--cells-dir", type=Path, default=CELLS_BASE,
                        help="Input directory containing cell_<root_id> CSV folders.")
    parser.add_argument("--output-dir", type=Path, default=RESULTS_DIR,
                        help="Directory for summary tables, figures, and Neuroglancer layers.")
    parser.add_argument(
        "--complete-only",
        action="store_true",
        help="Only include rows classified as 'complete' or 'sufficient', including qualified sufficient entries.",
    )
    parser.add_argument(
        "--synaptome-x-labels",
        choices=["segment_id", "distance"],
        default="distance",
        help="Deprecated compatibility option; all plots use left-right distance from the midline in µm.",
    )
    return parser.parse_args(argv)


def main(argv=None):
    global H5_MASKS, CELLS_BASE, RESULTS_DIR, FIGURES_DIR
    args = parse_args(argv)
    H5_MASKS = args.atlas.expanduser().resolve()
    CELLS_BASE = args.cells_dir.expanduser().resolve()
    RESULTS_DIR = args.output_dir.expanduser().resolve()
    FIGURES_DIR = RESULTS_DIR / "figures"
    if not H5_MASKS.is_file():
        raise SystemExit(f"MECE atlas not found: {H5_MASKS}")
    if not CELLS_BASE.is_dir():
        raise SystemExit(f"Synapse input directory not found: {CELLS_BASE}")
    ensure_dirs()
    (
        soma_ids,
        nuclei_records,
        projection_by_soma,
        incomplete_axon_by_soma,
        functional_by_soma,
    ) = read_metadata(args)
    print(f"[CELLS_BASE] {CELLS_BASE}")
    print(f"[n_cells] {len(soma_ids)}")
    run_global_analysis(
        soma_ids,
        nuclei_records,
        projection_by_soma,
        incomplete_axon_by_soma,
        functional_by_soma,
        synaptome_x_labels=args.synaptome_x_labels,
    )


if __name__ == "__main__":
    main()

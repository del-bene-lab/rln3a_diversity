#!/usr/bin/env python3
"""Extract input and output synapses from CAVE and assign brain regions.

Query selected neurons, convert EM pixel coordinates to micrometers, and write
per-neuron synapse tables, region histograms, and Neuroglancer annotation
layers. This optional acquisition step requires CAVE access. The published
regional analysis uses frozen CSV inputs and does not require this step.

By default, newly acquired tables are written under results/per_neuron, leaving
the bundled data unchanged. Run ``python scripts/per_cell_synaptome.py --help``
from synapse_regions for path options. See ../README.md for usage.
"""

from __future__ import annotations

import sys
import json
import argparse
from pathlib import Path
from typing import List, Dict, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py


try:
    from .mece_region_assignment import coarsen_hindbrain_to_mece_level_1
except ImportError:  # Direct execution: python scripts/per_cell_synaptome.py
    from mece_region_assignment import coarsen_hindbrain_to_mece_level_1

# Configuration

# Neurodata / CAVE
CAVE_URL       = "https://global.brain-wire-test.org/"
DATASTACK      = "fish1_full"
SYNAPSE_TABLE  = "synapses_axde"


# EM voxel size (µm/px) used by the synapse positions in CAVE (16×16×30 nm)
VOX_SIZE_UM_EM = (0.016, 0.016, 0.030)

ANALYSIS_DIR = Path(__file__).resolve().parents[1]
H5_MASKS = ANALYSIS_DIR / "data" / "fish1.0_MECE_labels.hdf5"
CELLS_BASE = ANALYSIS_DIR / "results" / "per_neuron"
DEFAULT_CELLS_CSV = ANALYSIS_DIR / "data" / "neurons.csv"

# Plot settings (Arial everywhere)
PLOT_DPI = 300
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
})

# Neuroglancer voxel size in meters (16×16×30 nm)
NG_VOX_M = (16e-9, 16e-9, 30e-9)

# Path helpers

def ensure_dirs() -> None:
    CELLS_BASE.mkdir(parents=True, exist_ok=True)



def cell_dir(soma_id: int) -> Path:
    return CELLS_BASE / f"cell_{soma_id}"

def cell_csv_inputs(soma_id: int) -> Path:
    return cell_dir(soma_id) / f"cell_{soma_id}_inputs.csv"

def cell_csv_outputs(soma_id: int) -> Path:
    return cell_dir(soma_id) / f"cell_{soma_id}_outputs.csv"

def cell_plots_dir(soma_id: int) -> Path:
    return cell_dir(soma_id) / "plots"

def cell_ng_dir(soma_id: int) -> Path:
    return cell_dir(soma_id) / "ng"

def ensure_cell_dirs(soma_id: int) -> None:
    for p in (cell_dir(soma_id), cell_plots_dir(soma_id), cell_ng_dir(soma_id)):
        p.mkdir(parents=True, exist_ok=True)

# CAVE queries and coordinate conversion

def px_to_um_arrays(x_px, y_px, z_px, vox_um) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sx, sy, sz = vox_um
    return x_px.astype(float)*sx, y_px.astype(float)*sy, z_px.astype(float)*sz

def extract_from_array_col(series) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract xyz from a list/tuple/ndarray or a string like '[x,y,z]'. Returns float arrays."""
    import ast
    xs, ys, zs = [], [], []
    for v in series:
        try:
            arr = v if isinstance(v, (list, tuple, np.ndarray)) else ast.literal_eval(str(v))
            x, y, z = float(arr[0]), float(arr[1]), float(arr[2])
        except Exception:
            x = y = z = np.nan
        xs.append(x); ys.append(y); zs.append(z)
    return np.asarray(xs, float), np.asarray(ys, float), np.asarray(zs, float)

def extract_xyz_pixels(df: pd.DataFrame, prefix: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return xyz *in pixels* for 'pre' or 'post' from cave tables."""
    arr_col = f"{prefix}_pt_position"
    if arr_col in df.columns:
        return extract_from_array_col(df[arr_col])
    cols = [f"{prefix}_pt_x", f"{prefix}_pt_y", f"{prefix}_pt_z"]
    if all(c in df.columns for c in cols):
        return (
            pd.to_numeric(df[cols[0]], errors="coerce").to_numpy(float),
            pd.to_numeric(df[cols[1]], errors="coerce").to_numpy(float),
            pd.to_numeric(df[cols[2]], errors="coerce").to_numpy(float),
        )
    n = len(df)
    return (np.full(n, np.nan), np.full(n, np.nan), np.full(n, np.nan))

def fetch_synapses_for_root(root_id: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Query CAVE for all synapses where this cell is post (inputs) or pre (outputs)."""
    try:
        from caveclient import CAVEclient
    except ImportError:
        raise SystemExit("Install caveclient first: pip install caveclient cloud-volume")
    import datetime as dt

    client = CAVEclient(datastack_name=DATASTACK, server_address=CAVE_URL)
    _ = client.materialize.get_versions()
    ts = dt.datetime.now(dt.timezone.utc)

    df_inputs = client.materialize.query_table(
        table=SYNAPSE_TABLE, timestamp=ts, filter_equal_dict={"post_pt_root_id": int(root_id)}
    )
    df_outputs = client.materialize.query_table(
        table=SYNAPSE_TABLE, timestamp=ts, filter_equal_dict={"pre_pt_root_id": int(root_id)}
    )

    if not isinstance(df_inputs, pd.DataFrame):  df_inputs = pd.DataFrame(df_inputs)
    if not isinstance(df_outputs, pd.DataFrame): df_outputs = pd.DataFrame(df_outputs)
    return df_inputs, df_outputs

def convert_positions_um(df_inputs: pd.DataFrame, df_outputs: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Add µm columns (post_x/y/z_um and pre_x/y/z_um) based on VOX_SIZE_UM_EM."""
    post_x, post_y, post_z = extract_xyz_pixels(df_inputs, "post")
    pre_x,  pre_y,  pre_z  = extract_xyz_pixels(df_outputs, "pre")

    post_x_um, post_y_um, post_z_um = px_to_um_arrays(post_x, post_y, post_z, VOX_SIZE_UM_EM)
    pre_x_um,  pre_y_um,  pre_z_um  = px_to_um_arrays(pre_x,  pre_y,  pre_z,  VOX_SIZE_UM_EM)

    df_inputs  = df_inputs.copy()
    df_outputs = df_outputs.copy()
    df_inputs[["post_x_um","post_y_um","post_z_um"]] = np.c_[post_x_um,post_y_um,post_z_um]
    df_outputs[["pre_x_um","pre_y_um","pre_z_um"]]   = np.c_[pre_x_um,pre_y_um,pre_z_um]
    return df_inputs, df_outputs

# Hdf5 label volume — mapping µm→voxel→region

def load_label_volume(h5_path: str) -> Tuple[np.ndarray, List[str], Tuple[float,float,float]]:
    """Load MECE label volume and region names; return (label_volume, region_names, voxel_size_um)."""
    import numpy as np, h5py
    with h5py.File(h5_path, "r") as f:
        if "label_volume" not in f:
            raise RuntimeError("HDF5 does not contain /label_volume")
        dset = f["label_volume"]

        # Force a copy from HDF5 (works for chunked/compressed datasets)
        try:
            label = dset[()]  # = np.array(dset, copy=True)
        except Exception:
            # Very old h5py fallback
            label = np.array(dset, copy=True)

        # Ensure expected dtype & contiguous layout
        if label.dtype != np.uint16:
            label = label.astype(np.uint16, copy=False)
        label = np.ascontiguousarray(label)

        # Region names (attribute may be bytes or an HDF5 array)
        names_attr = f.attrs.get("region_names", None)
        if names_attr is None:
            n = int(label.max())
            region_names = [f"region_{i}" for i in range(1, n+1)]
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
    """µm→voxel indices, rounded to nearest integer voxel index."""
    vx, vy, vz = voxel_size_um
    return (
        np.round(np.asarray(x_um, float) / vx).astype(np.int64),
        np.round(np.asarray(y_um, float) / vy).astype(np.int64),
        np.round(np.asarray(z_um, float) / vz).astype(np.int64),
    )

def map_df_to_regions(df: pd.DataFrame, prefix: str,
                      label: np.ndarray, region_names: List[str],
                      voxel_size_um: Tuple[float,float,float]) -> pd.DataFrame:
    """For each row's µm position, attach `{prefix}_regions` (list) and `{prefix}_n_regions`."""
    X, Y, Z = label.shape
    xi, yi, zi = um_to_voxel(df[f"{prefix}_x_um"], df[f"{prefix}_y_um"], df[f"{prefix}_z_um"], voxel_size_um)
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
    out[f"{prefix}_regions"]    = regions
    out[f"{prefix}_n_regions"]  = out[f"{prefix}_regions"].apply(len)
    return out

def counts_from_region_lists(series_of_lists: pd.Series) -> pd.DataFrame:
    """Explode list-of-regions and count occurrences per region."""
    exploded = series_of_lists.explode()
    exploded = exploded[exploded.notna() & (exploded != "")]
    return (exploded.value_counts()
            .rename_axis("region")
            .reset_index(name="count")
            .sort_values("region")
            .reset_index(drop=True))

# Plot helpers

def normalize(counts: Dict[str, int]) -> Dict[str, float]:
    tot = float(sum(counts.values()) or 1.0)
    return {k: v/tot for k, v in counts.items()}

def plot_counts_simple(counts_df: pd.DataFrame, title: str, out_path: Path) -> None:
    """Simple bar of absolute counts (per-cell), saved as PDF."""
    plt.figure(figsize=(6, 6))
    ax = plt.gca()

    if counts_df.empty:
        ax.text(
            0.5,
            0.5,
            "No mapped synapses",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        counts_df.sort_values("count", ascending=False).plot.bar(
            x="region",
            y="count",
            legend=False,
            ax=ax,
            edgecolor="black",
            linewidth=0.5,
        )
        plt.xticks(rotation=60, ha="right")
        ax.set_ylabel("Count")

    ax.set_title(title)
    plt.tight_layout()

    out_path = out_path.with_suffix(".pdf")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=PLOT_DPI, bbox_inches="tight", format="pdf")
    plt.close()
    print(f"[saved] {out_path}")

# Neuroglancer layer helpers

def um_to_ng_vox(x_um, y_um, z_um, ng_vox_m=NG_VOX_M):
    """µm→NG voxel coordinates (meters / voxel size)."""
    x_vox = (np.asarray(x_um, dtype=float) * 1e-6) / ng_vox_m[0]
    y_vox = (np.asarray(y_um, dtype=float) * 1e-6) / ng_vox_m[1]
    z_vox = (np.asarray(z_um, dtype=float) * 1e-6) / ng_vox_m[2]
    return x_vox, y_vox, z_vox

def df_for_role(df_src: pd.DataFrame, role: str) -> pd.DataFrame:
    """Extract NG-ready coordinates and metadata from inputs/outputs DataFrame."""
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
        "x": x_vox, "y": y_vox, "z": z_vox,
        "synapse_id": sid, "partner_segment": partner.astype(str),
    })

def rows_to_annotations(df: pd.DataFrame, role: str, soma_id: int) -> list:
    """Build NG 'point' annotations with stable ids and descriptive text."""
    anns = []
    tag = "pre" if role == "presynapse" else "post"
    for i, r in df.reset_index(drop=True).iterrows():
        anns.append({
            "type": "point",
            "id": f"{tag}_{soma_id}_{i+1}",
            "point": [float(r["x"]), float(r["y"]), float(r["z"])],
            "description": f"cell={soma_id} | {'output' if tag=='pre' else 'input'} | "
                           f"synapse_id={r['synapse_id']} | partner={r['partner_segment']}",
        })
    return anns

def build_layer(annotations: list, layer_name: str, color_hex: str | None = None) -> dict:
    """Construct a stand-alone NG layer JSON object."""
    layer = {
        "type": "annotation",
        "source": {
            "url": "local://annotations",
            "transform": {
                "outputDimensions": {
                    "x": [NG_VOX_M[0], "m"],
                    "y": [NG_VOX_M[1], "m"],
                    "z": [NG_VOX_M[2], "m"]
                }
            }
        },
        "tool": "annotatePoint",
        "tab": "annotations",
        "name": layer_name,
        "annotations": annotations
    }
    if color_hex:
        layer["annotationColor"] = color_hex
    return layer

def write_layers_state_snippet(path: Path, layers: list) -> None:
    """
    Write a small 'state snippet' that matches Neuroglancer's expected structure:
      { "layers": [ <layer_obj>, <layer_obj>, ... ] }
    This can be merged into an existing state.json by appending the objects
    to its top-level 'layers' array.
    """
    snippet = {"layers": layers}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(snippet, indent=2))

# Pipeline (per-cell processing + summaries + ng layers + nuclei regions)

def process_root(soma_id: int,
                 label: np.ndarray,
                 region_names: List[str],
                 mask_vox_um: Tuple[float,float,float],
                 global_presyn_ann: list,
                 global_postsyn_ann: list) -> Tuple[dict, dict, int, int]:
    """
    Per-cell:
      1) fetch synapses (inputs/outputs), convert to µm
      2) map to regions with HDF5 label volume
      3) save CSVs & per-cell plots
      4) produce per-cell NG layers (+ combined state snippet); add annotations to global lists
      5) return (normalized inputs profile, normalized outputs profile, n_inputs, n_outputs)
    """
    ensure_cell_dirs(soma_id)

    # 1) Fetch + convert
    df_inputs_raw, df_outputs_raw = fetch_synapses_for_root(soma_id)
    df_inputs, df_outputs         = convert_positions_um(df_inputs_raw, df_outputs_raw)

    n_inputs  = int(len(df_inputs))     # absolute counts for side-by-side histogram
    n_outputs = int(len(df_outputs))

    # 2) Save raw-with-µm CSVs (per cell)
    in_csv  = cell_csv_inputs(soma_id)
    out_csv = cell_csv_outputs(soma_id)
    df_inputs.to_csv(in_csv, index=False)
    df_outputs.to_csv(out_csv, index=False)
    print(f"[saved] {in_csv}\n[saved] {out_csv}")

    # 3) Region mapping
    mapped_inputs  = map_df_to_regions(df_inputs,  "post", label, region_names, mask_vox_um)
    mapped_outputs = map_df_to_regions(df_outputs, "pre",  label, region_names, mask_vox_um)

    # 4) Per-cell region counts + plots
    in_counts_df  = counts_from_region_lists(mapped_inputs["post_regions"])
    out_counts_df = counts_from_region_lists(mapped_outputs["pre_regions"])

    plots_dir = cell_plots_dir(soma_id)
    plot_counts_simple(in_counts_df,  f"Inputs by region (cell {soma_id})",  plots_dir / f"root_{soma_id}_inputs_counts.pdf")
    plot_counts_simple(out_counts_df, f"Outputs by region (cell {soma_id})", plots_dir / f"root_{soma_id}_outputs_counts.pdf")

    # 5) Neuroglancer per-cell (+ combined) + accumulate global
    pre_for_ng  = df_for_role(df_outputs, role="presynapse")   # outputs (pre)
    post_for_ng = df_for_role(df_inputs,  role="postsynapse")  # inputs  (post)

    presyn_ann  = rows_to_annotations(pre_for_ng,  role="presynapse",  soma_id=soma_id)
    postsyn_ann = rows_to_annotations(post_for_ng, role="postsynapse", soma_id=soma_id)

    ng_dir = cell_ng_dir(soma_id)
    layer_pre_cell  = build_layer(presyn_ann,  layer_name=f"cell_{soma_id}_presyn",  color_hex="#ff2b00")
    layer_post_cell = build_layer(postsyn_ann, layer_name=f"cell_{soma_id}_postsyn", color_hex="#2b00ff")

    # Write separate per-cell layers
    (ng_dir / f"cell_{soma_id}_presyn.json").write_text(json.dumps(layer_pre_cell,  indent=2))
    (ng_dir / f"cell_{soma_id}_postsyn.json").write_text(json.dumps(layer_post_cell, indent=2))

    # Combined per-cell state snippet (both layers, matching NG "layers": [...] structure)
    combined_path = ng_dir / f"cell_{soma_id}_inputs_outputs_layers.json"
    # Layer order: inputs (postsynaptic, blue), then outputs (presynaptic, red).
    write_layers_state_snippet(combined_path, [layer_post_cell, layer_pre_cell])
    print(f"[OK] NG layers written for cell {soma_id}: {ng_dir}")

    # Accumulate for global layers
    global_presyn_ann.extend(presyn_ann)
    global_postsyn_ann.extend(postsyn_ann)

    # 6) Return normalized profiles and absolute counts
    in_probs  = normalize(dict(zip(in_counts_df["region"],  in_counts_df["count"])))  if len(in_counts_df)  else {}
    out_probs = normalize(dict(zip(out_counts_df["region"], out_counts_df["count"]))) if len(out_counts_df) else {}
    return in_probs, out_probs, n_inputs, n_outputs

def run_pipeline(soma_ids: List[int]) -> None:

    # Load label volume once
    label, region_names, mask_vox_um = load_label_volume(H5_MASKS)
    print(f"[label_volume] shape={label.shape}, regions={len(region_names)}, voxel_size_um={mask_vox_um}")

    for soma_id in soma_ids:
        try:
            process_root(
                soma_id, label, region_names, mask_vox_um,
                [], []
            )
        except Exception as e:
            print(f"[WARN] soma_id {soma_id} failed: {e}", file=sys.stderr)

    print(f"[DONE] Per-cell synaptome outputs generated for {len(soma_ids)} requested cells.")

# Command-line interface

def parse_nucleus_pixel_field(val) -> Tuple[float,float,float] | None:
    """Parse '[x,y,z]' | '(x,y,z)' | 'x,y,z' in pixels → (x,y,z) or None."""
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

def parse_args(argv=None):
    p = argparse.ArgumentParser(description="Compute region connectivity statistics for multiple cells.")
    g = p.add_mutually_exclusive_group(required=False)
    g.add_argument("--root-ids", nargs="+", type=int, help="List of soma IDs.")
    g.add_argument(
        "--root-ids-csv",
        type=str,
        default=DEFAULT_CELLS_CSV,
        help=f"CSV with columns 'soma_id' and optional 'nucleus_position'. Default: {DEFAULT_CELLS_CSV}",
    )
    p.add_argument(
        "--complete-only",
        action="store_true",
        help="Only include neurons whose reconstruction column contains 'complete'.",
    )
    p.add_argument("--atlas", type=Path, default=H5_MASKS,
                   help="MECE HDF5 label volume.")
    p.add_argument("--cells-dir", type=Path, default=CELLS_BASE,
                   help="Output directory for newly acquired per-neuron tables.")
    return p.parse_args(argv)

def main(argv=None):
    global H5_MASKS, CELLS_BASE
    args = parse_args(argv)
    H5_MASKS = args.atlas.expanduser().resolve()
    CELLS_BASE = args.cells_dir.expanduser().resolve()
    if not H5_MASKS.is_file():
        raise SystemExit(f"MECE atlas not found: {H5_MASKS}")
    ensure_dirs()

    if args.root_ids is not None and len(args.root_ids) > 0:
        if args.complete_only:
            print("[WARN] --complete-only is ignored when --root-ids are supplied directly.")
        soma_ids = [int(x) for x in args.root_ids]
        nuclei_records = None
    else:
        csv_path = args.root_ids_csv or DEFAULT_CELLS_CSV
        df = pd.read_csv(csv_path, dtype={"soma_id": "string"})
        if args.complete_only:
            if "reconstruction" not in df.columns:
                raise SystemExit(
                    "--complete-only requires a 'reconstruction' column in the CSV. "
                    f"Got columns: {list(df.columns)}"
                )
            n_before = len(df)
            df = df[
                df["reconstruction"]
                .astype(str)
                .str.contains("complete", case=False, na=False)
            ].copy()
            print(f"[filter] complete-only: kept {len(df)} / {n_before} rows")
        if "soma_id" not in df.columns:
            raise SystemExit(f"CSV must contain 'soma_id' column. Got columns: {list(df.columns)}")
        soma_ids = df["soma_id"].dropna().astype(int).tolist()

        nuclei_records = None
        nucleus_position_col = None
        if "nucleus_position" in df.columns:
            nucleus_position_col = "nucleus_position"
        elif "nucleus_position (um)" in df.columns:
            nucleus_position_col = "nucleus_position (um)"

        if nucleus_position_col is not None:
            parsed = [parse_nucleus_pixel_field(v) for v in df[nucleus_position_col]]
            nuclei_records = [(int(sid), pos) for sid, pos in zip(df["soma_id"], parsed)
                              if pd.notna(sid) and pos is not None]

    print(f"[CELLS_BASE] {CELLS_BASE}")
    print(f"[CELLS_CSV ] {args.root_ids_csv if (args.root_ids is None) else '(provided via --root-ids)'}")
    run_pipeline(soma_ids)

if __name__ == "__main__":
    main()

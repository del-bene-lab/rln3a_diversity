# Synapse region analysis

Regional distributions of input and output synapses for reconstructed rln3a PAG
neurons. Maintainer: Jonathan Boulanger-Weill.

The bundled inputs reproduce the regional analysis offline. CAVE access is
required only to acquire new synapse tables.

## Contents

```text
synapse_regions/
├── environment.yml
├── scripts/
│   ├── all_cells_synaptome.py      # Regional summaries and publication figures
│   ├── per_cell_synaptome.py       # Optional acquisition from CAVE
│   └── mece_region_assignment.py   # Hindbrain region hierarchy
├── data/
│   ├── README.md
│   ├── neurons.csv                # Neuron metadata and segmentation root IDs
│   ├── fish1.0_MECE_labels.hdf5    # Compact, lossless atlas label volume
│   ├── provenance.json            # Sources, row counts, and SHA-256 checksums
│   ├── per_neuron/                # Frozen input/output synapse CSVs
│   └── reference_summaries/        # Original summary tables for validation
└── results/
    ├── figures/
    ├── inputs_summary.csv
    ├── outputs_summary.csv
    ├── synapse_counts_per_cell.csv
    └── nuclei_region_assignments.csv
```

The metadata contain 30 neurons: 22 PAG neurons used in the synapse analysis
and eight NI neurons retained as metadata context. The 44 bundled synapse CSVs
cover the 22 PAG neurons, with 3,979 input and 1,988 output synapse records.

## Installation

Run from the repository root:

```bash
conda env create -f synapse_regions/environment.yml
conda activate rln3a-synapse-regions
```

The environment is adapted from the original project's `relaxin` environment
and pins the relevant packages to the versions used for local verification.
It excludes machine-specific paths and platform-specific build strings.
`caveclient` and `cloud-volume` support optional data acquisition; they are not
imported by the offline analysis. The scripts use a noninteractive Matplotlib
backend. Arial is used when available, with Matplotlib's default font fallback
on systems without Arial.

Loading the atlas requires approximately 2.3 GB of array memory, in addition
to the Python process and plotting overhead.

## Reproduce the analysis

From the repository root:

```bash
python synapse_regions/scripts/all_cells_synaptome.py
```

Paths are resolved relative to the script location, so the command also works
from another working directory when the script path is specified correctly.
All generated files go under `synapse_regions/results/`; frozen inputs in
`data/` are not modified.

The main figures are:

- `inputs_outputs_region_proportions_stacked_ipsilateral.pdf`
- `inputs_outputs_region_proportions_stacked_contralateral.pdf`
- `synapse_counts_per_cell.pdf`

Figures are written to `results/figures/`. Summary CSVs are written directly
to `results/`. Neuroglancer annotations are regenerated in `results/ng_layers/`.
Generated outputs are excluded from version control; the original reference
tables are retained under `data/reference_summaries/`.

Alternative inputs and output locations can be supplied explicitly:

```bash
python synapse_regions/scripts/all_cells_synaptome.py \
  --root-ids-csv synapse_regions/data/neurons.csv \
  --cells-dir synapse_regions/data/per_neuron \
  --atlas synapse_regions/data/fish1.0_MECE_labels.hdf5 \
  --output-dir /path/to/results
```

Use `--complete-only` to restrict the regional analysis to neurons annotated
as complete or sufficient, including qualified sufficient annotations. The
default includes all 22 PAG neurons and marks incomplete axons in the figures.
Individual neurons can be selected with `--root-ids`; their annotations are
read from the bundled metadata table.

## Analysis conventions

- Inputs are synapses where the neuron is postsynaptic; outputs are synapses
  where it is presynaptic. Mapping uses the corresponding post- or presynaptic
  coordinate, respectively.
- Original EM voxel dimensions are `(0.016, 0.016, 0.030)` micrometers in XYZ
  order. Atlas indices are obtained by rounding micrometer coordinates to the
  nearest atlas voxel. The atlas stores XYZ axes and its voxel dimensions in
  HDF5 attributes.
- Hindbrain assignments are pooled to Rhombomere 1 through Rhombomere 7 and
  Caudal Hindbrain, using `mece_region_assignment.py`. Residual `Hindbrain`
  voxels map to Caudal Hindbrain; the figure legend displays Caudal Hindbrain
  as Rhombomere 8, preserving the source analysis convention.
- Regional proportions are normalized over mapped synapses. Unmapped synapses
  are excluded from those proportions but included in total synapse counts.
- Plot segments below 5% are pooled into `Other` separately for each neuron
  and direction. The summary CSVs retain the regional values before this
  plotting-only pooling.
- Cells are displayed by signed distance from the midline, using
  `(y - 16358) * 0.016` micrometers; positive values indicate anatomical left.
  Ipsilateral and contralateral projections are plotted separately.
- Printed input-count statistics include all analyzed neurons. Printed
  output-count statistics exclude neurons whose reconstruction annotation
  contains neither the word `complete` nor `sufficient`. This preserves the
  source analysis's operational incomplete-axon classification.

## Acquire new synapse tables (optional)

The acquisition script requires an authenticated CAVE account with access to
the configured `fish1_full` datastack and `synapses_axde` table. It queries the
service at execution time; newly acquired data can differ from the frozen
exports used above.

```bash
python synapse_regions/scripts/per_cell_synaptome.py \
  --root-ids 864691128634696641 \
  --cells-dir synapse_regions/results/per_neuron
```

The default acquisition output is also `results/per_neuron/`, leaving the frozen
tables intact. Pass that directory to the global script with `--cells-dir`
when analyzing newly acquired data. Without explicit root IDs, acquisition
uses all neurons in the metadata; its legacy `--complete-only` option selects
annotations containing `complete`, whereas the global analysis also accepts
`sufficient`.

See [data/README.md](data/README.md) for the data schema, provenance, and snapshot
information. No credentials are stored in this analysis directory.


# rln3a diversity

Analysis code and supporting data for the rln3a diversity study.

The repository is organized into three analysis areas: calcium imaging and
behavior, maintained by Matthieu, and synapse region analysis, maintained by
Jonathan. Shared neuron metadata belong in `data/`.

## Repository structure

```text
rln3a_diversity/
├── README.md
├── data/
│   └── README.md                 # Shared neuron metadata and data documentation
├── calcium_imaging/              # Matthieu
│   ├── README.md
│   ├── notebooks/
│   │   └── Analysis_Figure3_part1.ipynb
│   └── scripts/
├── behavior/                     # Matthieu
│   ├── README.md
│   └── scripts/
└── synapse_regions/              # Jonathan
    ├── README.md
    ├── scripts/
    ├── data/                     # Synapse tables used in the analyses
    └── results/                  # Summary tables and final figures
```

## Analyses

| Analysis | Maintainer | Contents |
| --- | --- | --- |
| [Calcium imaging](calcium_imaging/README.md) | Matthieu | Calcium imaging analysis notebooks and supporting scripts. |
| [Behavior](behavior/README.md) | Matthieu | Behavioral analysis scripts and notebooks as they are added. |
| [Synapse regions](synapse_regions/README.md) | Jonathan | Synapse extraction, assignment to brain regions, and regional distribution analyses as they are added. |

The existing Figure 3 notebook has been copied without modification to
[`calcium_imaging/notebooks/Analysis_Figure3_part1.ipynb`](calcium_imaging/notebooks/Analysis_Figure3_part1.ipynb).
The original notebook at the repository root is retained. Use the organized
copy for subsequent work.

The notebook currently references local data paths and a separately maintained
`calcium_imaging.dFoverF` helper. Those dependencies need to be supplied and
documented before the analysis can run on another computer.

## Data and reproducibility

- Store the shared neuron metadata table, including segment/root IDs, in
  `data/neurons.csv` once the table has been curated. See
  [data documentation](data/README.md) for the intended fields.
- Keep the synapse tables used for the paper in `synapse_regions/data/` and
  derived summary tables and final plots in `synapse_regions/results/`.
- Document the input files, software environment, execution order, and paper
  figure panels in each analysis folder's README as the analyses are added.
- Use repository-relative paths or configurable data directories. Provide
  download links and dataset versions for large files hosted separately.

This is the initial repository structure. Neuron metadata, synapse analysis
code and tables, behavior analyses, and reproducible environment files will be
added as they are prepared for release. Paper citation and licensing information
will be added when finalized.

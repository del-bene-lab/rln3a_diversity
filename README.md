# rln3a diversity

Analysis code and supporting data for Tuffery et al., 2026. 

## Repository structure

```text
rln3a_diversity/
├── README.md
├── calcium_imaging/              
│   ├── README.md
│   ├── notebooks/
│   │   └── Analysis_Figure3_part1.ipynb
│   └── scripts/
├── behavior/                     
│   ├── README.md
│   └── scripts/
└── synapse_regions/              
    ├── README.md
    ├── environment.yml
    ├── scripts/
    ├── data/                     
    └── results/                  
```

## Synapse region analysis

[Synapse region analysis](synapse_regions/README.md) includes the scripts,
frozen synapse tables, neuron metadata with segment IDs, and compact MECE atlas
needed to reproduce regional distributions for 22 PAG neurons offline.

From the repository root:

```bash
conda env create -f synapse_regions/environment.yml
conda activate rln3a-synapse-regions
python synapse_regions/scripts/all_cells_synaptome.py
```

The analysis README documents the data sources, coordinate conventions,
filtering, generated figures, and optional CAVE acquisition workflow.


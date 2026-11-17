# Synapse region analysis

Maintainer: Jonathan.

This folder will contain the analysis of input and output synapses across brain
regions.

## Organization

- `scripts/`: synapse extraction, brain-region assignment, and summary plots.
- `data/`: the frozen input and output synapse tables used for the paper.
- `results/`: derived summary tables and final figures.
- `../data/neurons.csv`: the planned shared neuron metadata table, including
  segment/root IDs and analysis inclusion criteria.

The intended workflow is to extract synapses for the selected neurons, assign
synapses to brain regions, and summarize regional counts and proportions per
neuron and across neurons.

When importing the analysis, document the segmentation snapshot, synapse data
version, atlas and region hierarchy, coordinate units, software environment,
execution commands, and figure panels here. The analysis scripts and data have
not yet been copied into this repository.

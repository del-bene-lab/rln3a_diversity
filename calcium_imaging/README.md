# Calcium imaging

Maintainer: Matthieu.

## Contents

- [`notebooks/Analysis_Figure3_part1.ipynb`](notebooks/Analysis_Figure3_part1.ipynb):
  an unchanged copy of the notebook originally stored at the repository root.
- `scripts/`: supporting analysis functions and scripts as they are added.

The notebook imports NumPy, pandas, Matplotlib, seaborn, SciPy, and the external
`calcium_imaging.dFoverF` helper. It currently contains computer-specific paths
to that helper and the imaging data. The copy preserves those paths and the
notebook's existing outputs.

Document the software environment, data access, execution order, and figure
panels here when preparing the analysis for release. The helper and input data
are not yet included in this folder.

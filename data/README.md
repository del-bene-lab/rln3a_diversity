# Shared neuron metadata

This folder is reserved for the curated neuron metadata table used in the paper.
No neuron table has been added yet.

## Planned table: `neurons.csv`

Use one row per neuron and unique column names. Include the complete relevant
metadata, with segment/root IDs in the same table rather than a duplicate ID-only
list.

Suggested fields:

| Field | Description |
| --- | --- |
| `neuron_id` | Stable neuron label matching the paper. |
| `segment_id` | Segmentation root ID at the documented dataset snapshot. |
| `nucleus_id` | Associated nucleus ID, when available. |
| `hemisphere` | Left or right, using a documented convention. |
| `soma_x`, `soma_y`, `soma_z` | Soma coordinates in the documented coordinate system. |
| `reconstruction_status` | Reconstruction completeness or quality category. |
| `functional_class` | Functional annotation, when available. |
| `projection_class` | Projection annotation, when available. |
| `included_in_synapse_analysis` | Whether the neuron is included in the synapse analysis. |

Record the source table, segmentation dataset and snapshot date or version,
coordinate system, axis order, coordinate units, and missing-value conventions
when the table is added. If root IDs have changed over time, name historical
ID columns explicitly and document their dates or versions.

Preserve large integer IDs exactly; load them as strings or 64-bit integers,
without conversion through floating-point values. Remove unused empty columns
and resolve duplicate headers before adding a curated export.

Synapse-level tables belong in `../synapse_regions/data/`. Large imaging files
and atlas volumes can be hosted separately, with download links and versions
documented in the relevant analysis README.

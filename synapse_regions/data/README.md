# Data sources and schema

These files are the frozen inputs for the regional synapse analysis. The source
project is `rel3a_function_and_connectivity`; per-neuron exports and reference
summaries originate from `connectivity/synaptome/cells_synaptome_info/`.
File checksums and the metadata cleanup are recorded in
[`provenance.json`](provenance.json).

## Neuron metadata: `neurons.csv`

The source is `csv/rel3a_neurons_111126.csv`. All 30 nonempty neuron rows and all
nonempty field values are retained. Entirely empty columns were removed,
including two empty duplicate `soma_id` columns. Existing informative column
names were preserved for compatibility with the analysis.

| Column | Interpretation |
| --- | --- |
| `soma_id` | Segmentation root/segment ID used to identify the neuron, not the nucleus ID. |
| `functional_id` | Original neuron annotation; PAG/NI membership is determined from this string. |
| `nucleus_id` | Associated nucleus identifier. |
| `nucleus_position` | Comma-separated XYZ nucleus coordinates in original EM voxels. |
| `left/right` | Original hemisphere annotation. |
| `reconstruction` | Reconstruction status used for filtering and incomplete-axon markers. |
| `projection classifier` | Ipsilateral or contralateral projection annotation, when assigned. |
| `type`, `comment`, `neurotransmitter classifier` | Original biological annotations and notes. |
| `neuroglancer_link` | Original viewer link. |
| `imaging_modality`, `date_of_tracing`, `tracer names` | Original tracing provenance. |
| `# contra branches` | Original contralateral branch annotation. |

Empty fields represent missing annotations. Identifiers must be read as strings
or 64-bit integers, never converted through floating-point values.

## Synapse tables: `per_neuron/`

Each `cell_<soma_id>/` directory contains the neuron's unmodified
`cell_<soma_id>_inputs.csv` and `cell_<soma_id>_outputs.csv` exports. Only the
22 PAG neurons selected by the regional analysis are included.

Key columns are:

- `id`: synapse identifier.
- `pre_pt_root_id`, `post_pt_root_id`: presynaptic and postsynaptic segment IDs.
- `pre_pt_supervoxel_id`, `post_pt_supervoxel_id`: original supervoxel IDs.
- `pre_pt_position`, `post_pt_position`: original XYZ coordinates in EM voxels.
- `post_x_um`, `post_y_um`, `post_z_um`: postsynaptic coordinates in input tables.
- `pre_x_um`, `pre_y_um`, `pre_z_um`: presynaptic coordinates in output tables.
- `created`, `superceded_id`, `valid`, `size`, `ctr_pt_position`: fields retained
  unchanged from the source exports, including their original spelling.

Rows are neuron-specific synaptic records. A synapse between two listed neurons
can occur in one neuron's outputs and the other's inputs, so summing all table
rows does not represent a globally deduplicated synapse count.

The source configuration names the `fish1_full` CAVE datastack, the
`synapses_axde` synapse table, and a Neuroglancer segmentation source ending in
`fish1_v250915`. The original exports do not record their CAVE query timestamp
or materialization version. Those values are recorded as unknown in the
provenance manifest; the segmentation source name is not evidence of the
exports' exact snapshot. Reproduction therefore uses these frozen CSVs.

## Atlas: `fish1.0_MECE_labels.hdf5`

The source is `fish1.0-tools/brain_regions/fish1.0_MECE_masks.hdf5`.
This compact copy retains only `/label_volume` and the original file attributes
`region_names`, `shape_xyz`, and `voxel_size_um`. The individual binary region
masks are not used by the scripts and were omitted.

Every stored compressed label-volume chunk and all root attributes were
verified against the original. The compact file is approximately 13.7 MB;
the source file is approximately 74.6 MB. No resampling, relabeling, or changes
to the spatial extent were applied.

The volume is a `(2100, 1016, 544)` XYZ array of unsigned 16-bit label IDs.
Background is 0; label `i` uses `region_names[i - 1]`. Its voxel dimensions
are stored in `voxel_size_um` and read directly by the scripts.

## Reference tables

`reference_summaries/` contains unmodified copies of the source
`inputs_summary.csv`, `outputs_summary.csv`, `synapse_counts_per_cell.csv`, and
`nuclei_region_assignments.csv`. They are used to compare reproduced outputs,
not as inputs to the regional calculations.

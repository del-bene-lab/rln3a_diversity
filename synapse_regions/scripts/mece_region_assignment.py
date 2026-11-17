"""Shared MECE region-name normalization for synaptome analyses.

The dense atlas label volume stores one winning label per voxel, including
fine hindbrain nuclei.  Synaptome analyses intentionally stop hindbrain
assignments at MECE level 1: Rhombomere_1 through Rhombomere_7, plus
Caudal_Hindbrain.  Labels outside the hindbrain are returned unchanged.

Parentage below comes from the fish1.0 MECE level-2 and level-3 hierarchy
files in ``fish1.5-tools/regions_masks_upload/mece_lists_fish1.0``.
"""

from __future__ import annotations


HINDBRAIN_MECE_LEVEL_1_CHILDREN = {
    "Rhombomere_1": frozenset({
        "Interpeduncular_Nucleus",
        "Raphe_-_Superior",
        "Ventrolateral_Stripe_of_Serotonergic_Neurons",
        "Oculomotor_Nucleus_nIV",
        "RoL-R1",
        "Cerebellum",
        "Locus_Coeruleus",
        "Granular_Eminence",
        "Valvula_Cerebellum",
        "Cerebellar_Corpus",
        "Medial_Caudal_Lobe",
    }),
    "Rhombomere_2": frozenset({
        "Anterior_Cluster_of_nV_Trigeminal_Motorneurons",
        "Raphe_-_Superior1",
        "RoL2",
        "RoM2",
        "Locus_Coeruleus1",
        "RoM1",
        "Ventrolateral_Stripe_of_Serotonergic_Neurons1",
    }),
    "Rhombomere_3": frozenset({
        "Posterior_Cluster_of_nV_Trigeminal_Motorneurons",
        "RoM3",
        "RoV3",
        "Spiral_Fiber_Neuron_Posterior_cluster",
        "Raphe_-_Inferior",
        "RoL3",
        "Spiral_Fiber_Neuron_Anterior_cluster",
        "Ventrolateral_Stripe_of_Serotonergic_Neurons2",
    }),
    "Rhombomere_4": frozenset({
        "MiV1",
        "Mauthner",
        "MiR1",
        "Raphe_-_Inferior1",
        "Spinal_Backfill_Vestibular_Population",
        "Mauthner_Cell_Axon_Cap",
        "MiM1",
        "MiR2",
    }),
    "Rhombomere_5": frozenset({
        "MiV2",
        "Tangential_Vestibular_Nucleus",
        "VII_Facial_Motor_and_octavolateralis_efferent_neurons",
        "Medial_Vestibular_Nucleus",
        "MiD2",
        "Raphe_-_Inferior2",
    }),
    "Rhombomere_6": frozenset({
        "MiD3",
        "VII_Facial_Motor_and_octavolateralis_efferent_neurons1",
        "Medial_Vestibular_Nucleus_1",
        "Tangential_Vestibular_Nucleus1",
        "Raphe_-_Inferior3",
        "MiT",
    }),
    "Rhombomere_7": frozenset({
        "CaD",
        "CaV",
        "VII_Facial_Motor_and_octavolateralis_efferent_neurons2",
        "Raphe_-_Inferior4",
    }),
    "Caudal_Hindbrain": frozenset({
        "Area_Postrema",
        "Lateral_Reticular_Nucleus",
        "Small_cluster_of_TH_stained_neurons",
        "Inferior_Olive",
        "Noradrendergic_neurons_of_the_Interfascicular_and_Vagal_areas",
        "X_Vagus_motorneuron_cluster",
    }),
}

HINDBRAIN_MECE_LEVEL_1_REGIONS = tuple(HINDBRAIN_MECE_LEVEL_1_CHILDREN)

HINDBRAIN_REGION_TO_MECE_LEVEL_1 = {
    child: parent
    for parent, children in HINDBRAIN_MECE_LEVEL_1_CHILDREN.items()
    for child in children
}
HINDBRAIN_REGION_TO_MECE_LEVEL_1.update(
    {parent: parent for parent in HINDBRAIN_MECE_LEVEL_1_REGIONS}
)

# The atlas label volume was assembled alphabetically with later masks taking
# precedence.  The broad Hindbrain mask therefore overwrote Caudal_Hindbrain,
# while Rhombomere_1..7 (which sort later) remained explicit.  In the synapse
# dataset, residual ``Hindbrain`` voxels belong to the caudal-hindbrain mask.
HINDBRAIN_REGION_TO_MECE_LEVEL_1["Hindbrain"] = "Caudal_Hindbrain"


def coarsen_hindbrain_to_mece_level_1(region_name: str) -> str:
    """Return a hindbrain label's MECE-1 parent; preserve all other labels."""
    region_name = str(region_name)
    return HINDBRAIN_REGION_TO_MECE_LEVEL_1.get(region_name, region_name)

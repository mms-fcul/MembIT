# Example membrane-analysis workflows

MembIT is most useful when it is treated as one component of a larger membrane
analysis protocol. This document gives examples of how MembIT-style calculations
can be integrated with complementary GROMACS analyses.

The examples are based on scripts used for an ASIC1a + POPC membrane system, but
the concepts are general.

## 1. Global membrane/box convergence: projected area per lipid

A projected box-area area-per-lipid calculation can be used as a global
equilibration/convergence metric. It is not a local protein-corrected area per
lipid, because it uses the xy box area and the total number of lipids.

Conceptually:

```text
box_area_nm2 = Lx * Ly
projected_box_average_APL_A2 = box_area_nm2 * 2 / N_TOTAL_LIPIDS * 100
```

This metric is useful for monitoring global box/membrane behavior, but it should
not be overinterpreted as the local membrane response around the protein. It is
best interpreted together with local metrics such as thickness, deformation, or
first-shell deformation.

Typical tools:

- `gmx traj -ob` to extract box dimensions;
- a small refinement script to compute projected APL and smoothing windows.

## 2. Local membrane deformation with MembIT

Local deformation asks how far the membrane marker atoms near the protein are
displaced relative to a bulk membrane reference.

A typical command is:

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 6 1 0 6 25   -deformation   -o deformation
```

One useful interpretation is:

- `0-6 Å`: near-protein or first-shell region;
- `25 Å`: bulk cutoff used to define the membrane reference away from the
  protein.

The raw top/bottom signs are leaflet-dependent. Downstream plotting should
document the sign convention used.

## 3. Radial membrane-thickness profiles with MembIT

Radial thickness calculations use the same prepared MembIT input, but run without
`-deformation`. The output is absolute membrane half-thickness rather than
deformation relative to bulk.

Example:

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 1 1 0 60 25   -o radial_thickness
```

Useful downstream summaries include:

- bulk upper and lower half-thickness;
- near-protein upper and lower half-thickness;
- total near-protein thickness;
- near-minus-bulk total thickness;
- radial profile averaged over frames.

## 4. Region-specific insertion and tilt

For some proteins, whole-protein insertion is not the most meaningful metric.
For example, a large soluble or extracellular domain can dominate the whole
protein center of mass and hide the behavior of the transmembrane region.

A better strategy is to define analysis-specific groups such as:

```text
TM_core_CA
TM_low_CA
TM_high_CA
Upper_P31
Lower_P31
```

Then use tools such as `gmx distance -oxyz` to compute:

- the z-position of the TM core relative to upper and lower leaflet markers;
- the membrane-midplane-relative insertion of the TM core;
- a TM-axis vector and tilt angle.

This is complementary to MembIT's insertion/thickness logic and follows the same
general preparation principle: define groups that match the physical question,
then document the interpretation.

## Choosing the right analysis

A good membrane-analysis report often combines several metrics:

| Question | Useful analysis |
| --- | --- |
| Is the membrane/box globally equilibrated? | Projected box-area APL |
| Is the local membrane perturbed around the protein? | MembIT deformation |
| How does thickness vary with distance from the protein? | MembIT radial thickness |
| Is the transmembrane bundle moving vertically or tilting? | Region-specific insertion/tilt |
| Are top and bottom leaflets behaving differently? | Leaflet-specific deformation or thickness |

The goal is not to run every possible analysis. The goal is to choose the
metrics that match the biological question and to prepare the trajectory/index
accordingly.

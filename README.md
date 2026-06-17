# MembIT

MembIT analyzes membrane-protein trajectories to quantify membrane insertion,
local membrane thickness, and membrane deformation around a protein or another
center of interest.

MembIT is designed to be part of a broader membrane-analysis workflow. It can be
used together with GROMACS and downstream plotting/refinement scripts to study,
for example:

- local membrane deformation around a protein or selected center of interest;
- radial membrane-thickness profiles;
- near-protein versus bulk membrane thickness;
- insertion of a selected protein region relative to the membrane midplane;
- convergence/equilibration metrics when interpreted together with complementary
  analyses such as projected area per lipid.

The recommended workflow is **not** to run MembIT directly on a raw trajectory
from a molecular-dynamics engine. MembIT should be used on a trajectory prepared
for membrane analysis: molecules imaged consistently, the protein or chosen
solute centered, unnecessary atoms removed when appropriate, and a matching
structure/index generated for the processed atom set.

## What MembIT needs

A normal MembIT calculation needs:

- a trajectory supplied with `-f`;
- an index file supplied with `-n`;
- for binary trajectories such as XTC/TRR, a matching structure/topology supplied
  with `-s`;
- one or more analysis options, such as `-thickness`, `-deformation`, or
  `-insertion`.

The index must contain these groups:

- `Protein`
- `Center_of_Interest`
- `Monolayer1`
- `Monolayer2`

For thickness/deformation calculations, `Monolayer1` and `Monolayer2` should
normally contain leaflet marker atoms, such as phosphate/headgroup atoms. They
should not be complete lipid atom clouds.

## Recommended workflow

1. Prepare the raw MD trajectory with GROMACS or equivalent tools:
   - handle PBC and molecule imaging;
   - center the protein, solute, or complex of interest;
   - keep the membrane and protein in a consistent representation.
2. Remove atoms not needed by the analysis:
   - for MembIT thickness/deformation, this usually means keeping the protein or
     center of interest plus lipid leaflet marker atoms;
   - for complementary analyses, such as global projected area per lipid, the
     full system or full box information may still be useful.
3. Generate the structure and index against the same processed atom set used in
   the trajectory.
4. Run the MembIT calculation or complementary analysis.
5. Refine and plot the raw output using a documented sign convention.

For larger MembIT trajectories, XTC/TRR input is preferred over multi-frame PDB
because it is much smaller and faster to read. PDB input remains supported for
legacy workflows and validation.

## Example: local deformation

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 6 1 0 6 25   -deformation   -o analysis/deformation
```

## Example: radial thickness profile

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 1 1 0 60 25   -o analysis/radial_thickness
```

This runs the thickness calculation without `-deformation`, so the output is
absolute half-thickness rather than deformation relative to a bulk reference.

## Example: input diagnostic

Before running a long calculation, inspect input consistency:

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 6 1 0 6 25   --diagnose-index
```

## Documentation

See:

- [`docs/input_preparation.md`](docs/input_preparation.md)
- [`docs/index_files.md`](docs/index_files.md)
- [`docs/trajectory_formats.md`](docs/trajectory_formats.md)
- [`docs/usage_examples.md`](docs/usage_examples.md)
- [`docs/analysis_workflows.md`](docs/analysis_workflows.md)
- [`docs/troubleshooting.md`](docs/troubleshooting.md)
- [`docs/developer_validation.md`](docs/developer_validation.md)

## Citation

Please cite the original MembIT publication (https://doi.org/10.1142/S2737416523500254)
and the repository version used for your analysis.

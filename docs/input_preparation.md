# Preparing trajectories for MembIT

MembIT should be run on trajectories prepared for membrane analysis. A raw
trajectory written directly by the MD engine is usually not the best input.

## Why preparation matters

Membrane analyses depend on the relative geometry of the protein, center of
interest, and membrane leaflets. If the system drifts through the box, molecules
are split by periodic boundary conditions, or the protein and membrane are not
imaged consistently, local thickness and insertion calculations can become
physically meaningless even if the trajectory file is readable.

## Recommended preprocessing

A typical GROMACS-based workflow is:

1. Use `gmx trjconv` to treat PBC and center the relevant solute.
2. Concatenate treated trajectory segments if the simulation was run in parts.
3. Create a reduced analysis trajectory containing only:
   - the protein or center of interest;
   - the lipid marker atoms needed to define the leaflets.
4. Generate a matching structure file and MembIT index for the same processed
   atom set.

Waters and ions should usually be removed unless they are part of the scientific
question. This reduces I/O, avoids unnecessary memory use, and prevents users
from building misleading indexes.

## Example PBC/centering step

The exact groups depend on the system, but the common pattern is:

```bash
gmx trjconv \
  -f raw.xtc \
  -s run.tpr \
  -n index.ndx \
  -o centered.xtc \
  -center \
  -ur compact \
  -pbc mol
```

When prompted, choose the group used for centering, commonly `Solute` or the
protein/complex, and then choose the group written to the output trajectory.

## Reduced analysis trajectory

For membrane deformation/thickness workflows, a reduced trajectory often contains
`Protein + Phos`, where `Phos` means the chosen leaflet marker atoms. For POPC in
the workflow this documentation was based on, the marker atoms were:

```text
O31 P31 O32 O33 O34
```

The exact atom names are force-field dependent. Always verify them with GROMACS.

## Matching files

The trajectory, structure, and index must describe the same atom set in the same
atom-numbering scheme.

Good:

```text
treated_Protein_Phos.xtc
treated_Protein_Phos.gro
membit_index_for_treated_Protein_Phos.ndx
```

Bad:

```text
full_system.xtc
full_system.tpr
membit_index_generated_for_reduced_Protein_Phos.ndx
```

This mismatch is especially easy to create when moving from PDB-based workflows
to XTC/TRR workflows.

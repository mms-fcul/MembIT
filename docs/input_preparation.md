# Preparing trajectories for MembIT

MembIT should be run on trajectories prepared for membrane analysis. A raw
trajectory written directly by the MD engine is usually not the best input.

Preparation is not only a performance optimization. It is part of the scientific
definition of the analysis. The trajectory, structure, and index must represent
the same physical and atom-numbering model.

## Why preparation matters

Membrane analyses depend on the relative geometry of the protein, center of
interest, and membrane leaflets. If the system drifts through the box, molecules
are split by periodic boundary conditions, or the protein and membrane are not
imaged consistently, local thickness, deformation, and insertion calculations can
become physically meaningless even if the trajectory file is readable.

## Recommended preprocessing

A typical GROMACS-based workflow is:

1. Use `gmx trjconv` to treat PBC and center the relevant solute.
2. Concatenate treated trajectory segments if the simulation was run in parts.
3. Create one or more analysis-specific trajectories:
   - a reduced MembIT trajectory containing the protein or center of interest
     plus the lipid marker atoms needed to define the leaflets;
   - optionally, a centered full-system or solute trajectory for complementary
     GROMACS-based analyses that require box dimensions, full lipid molecules,
     waters, ions, or other groups.
4. Generate a matching structure file and MembIT index for the same processed
   atom set used by the MembIT trajectory.

Waters and ions should usually be removed from the MembIT trajectory unless they
are part of the scientific question. This reduces I/O, avoids unnecessary memory
use, and prevents users from building misleading indexes. Complementary analyses
may still use the full or solute trajectory when that is the correct scientific
input.

## Example PBC/centering step

The exact groups depend on the system, but the common pattern is:

```bash
gmx trjconv   -f raw.xtc   -s run.tpr   -n index.ndx   -o centered.xtc   -center   -ur compact   -pbc mol
```

When prompted, choose the group used for centering, commonly `Solute` or the
protein/complex, and then choose the group written to the output trajectory.

## Reduced MembIT trajectory

For membrane deformation/thickness workflows, a reduced trajectory often contains
`Protein + Phos`, where `Phos` means the chosen leaflet marker atoms. For POPC in
the workflow this documentation was based on, the marker atoms were:

```text
O31 P31 O32 O33 O34
```

The exact atom names are force-field dependent. Always verify them with GROMACS.
These are the names used in the documented AMBER14SB/POPC workflow.

A reduced trajectory is useful for:

- local deformation around a protein or center of interest;
- radial thickness profiles;
- near-protein versus bulk membrane-thickness comparisons;
- insertion calculations that use membrane marker atoms as a reference.

## Complementary analysis trajectories

Not every membrane analysis should use the same reduced MembIT trajectory.

For example:

- transmembrane-core insertion or tilt may need a custom index with specific
  protein-region groups and leaflet marker groups;
- water, ions, or full lipids should be retained only when they are directly
  relevant to the metric being calculated.

The important rule is not "always reduce everything"; the rule is:

```text
prepare an analysis-specific trajectory/index pair whose atom content matches
the scientific question.
```

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

Use `--diagnose-index` before long production runs.

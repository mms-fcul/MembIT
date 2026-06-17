# Usage examples

## Diagnose an input set

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 6 1 0 6 25   --diagnose-index
```

Use this before long XTC/TRR calculations to check that the trajectory,
structure, and index are consistent.

## Local thickness/deformation around a protein

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 6 1 0 6 25   -deformation   -o analysis/deformation
```

This is useful for first-shell deformation around a protein or another center of
interest.

## Radial membrane thickness

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 1 1 0 60 25   -o analysis/radial_thickness
```

This runs without `-deformation`, so the output is absolute membrane
half-thickness. Downstream scripts can summarize near-protein thickness, bulk
thickness, and radial profiles.

## Legacy PDB workflow

```bash
python membit.py   -f Phos_Prot_010.pdb   -n membit_index_010.ndx   -thickness 6 1 0 6 25   -deformation   -o def
```

PDB input is useful for legacy workflows and validation, but XTC/TRR is preferred
for larger production trajectories.

## Insertion

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -insertion closest   -o analysis/insertion
```

For proteins with large soluble domains, consider whether whole-protein insertion
is the right physical metric. A region-specific workflow, such as TM-core
insertion/tilt, may be more interpretable.

## Profile timing

```bash
python membit.py   -f treated_Protein_Phos.xtc   -s treated_Protein_Phos.gro   -n membit_index.ndx   -thickness 6 1 0 6 25   -deformation   --profile-timing   -o analysis/deformation
```

Use this to identify whether runtime is dominated by trajectory reading,
index/atom mapping, or the analysis calculation itself.

# Trajectory formats

MembIT supports both legacy PDB trajectories and binary trajectories such as XTC
through MDAnalysis.

## PDB mode

PDB mode is useful for legacy workflows and validation. The index atom numbers
refer to PDB atom serials.

Example:

```bash
python membit.py   -f Phos_Prot_010.pdb   -n membit_index_010.ndx   -thickness 6 1 0 6 25   -deformation   -o def
```

## XTC/TRR mode

For XTC/TRR, provide a matching structure with `-s`:

```bash
python membit.py   -f Phos_Prot_010.xtc   -s Phos_Prot_010.gro   -n membit_index_010.ndx   -thickness 6 1 0 6 25   -deformation   -o def
```

The structure file supplies atom names, residue names, atom ordering, and other
metadata required by MembIT.

## Format selection

The default `-format auto` detects the trajectory type from the file extension.
You can also set it explicitly:

```bash
-format xtc
-format pdb
-format trr
```

## Important constraint

For XTC/TRR workflows, the trajectory, structure, and index must match. Do not
mix:

- a reduced index with a full-system trajectory;
- a full-system structure with a reduced trajectory;
- a PDB-generated index with an XTC that has different atom numbering.

Use `--diagnose-index` before long production runs.

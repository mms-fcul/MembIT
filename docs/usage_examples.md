# Usage examples

## Diagnose an input set

```bash
python membit.py \
  -f treated_Protein_Phos.xtc \
  -s treated_Protein_Phos.gro \
  -n membit_index.ndx \
  -thickness 6 1 0 6 25 \
  --diagnose-index
```

## Thickness and deformation

```bash
python membit.py \
  -f treated_Protein_Phos.xtc \
  -s treated_Protein_Phos.gro \
  -n membit_index.ndx \
  -thickness 6 1 0 6 25 \
  -deformation \
  -o analysis/deformation
```

The common thickness argument pattern is:

```text
-thickness window_size window_step min_radius first_shell_cutoff bulk_cutoff
```

## Legacy PDB workflow

```bash
python membit.py \
  -f Phos_Prot_010.pdb \
  -n membit_index_010.ndx \
  -thickness 6 1 0 6 25 \
  -deformation \
  -o def
```

## Profile timing

```bash
python membit.py \
  -f treated_Protein_Phos.xtc \
  -s treated_Protein_Phos.gro \
  -n membit_index.ndx \
  -thickness 6 1 0 6 25 \
  -deformation \
  --profile-timing \
  -o analysis/deformation
```

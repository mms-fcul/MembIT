# MembIT

MembIT analyzes membrane-protein trajectories to quantify membrane insertion,
local membrane thickness, and membrane deformation around a protein or another
center of interest.

The recommended workflow is not to run MembIT directly on a raw trajectory from
a molecular-dynamics engine. MembIT should be used on a trajectory prepared for
membrane analysis: molecules imaged consistently, the protein or chosen solute
centered, unnecessary atoms removed, and a matching structure/index generated for
the processed atom set.

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
2. Remove atoms not needed by MembIT, usually water and ions.
3. Keep the protein or center of interest and the lipid marker atoms needed to
   define the two leaflets.
4. Generate the structure and index against the same processed atom set used in
   the trajectory.
5. Run MembIT.

For larger trajectories, XTC/TRR input is preferred over multi-frame PDB because
it is much smaller and faster to read. PDB input remains supported for legacy
workflows and validation.

## Example

```bash
python membit.py \
  -f treated_Protein_Phos.xtc \
  -s treated_Protein_Phos.gro \
  -n membit_index.ndx \
  -thickness 6 1 0 6 25 \
  -deformation \
  -o analysis/deformation
```

Before running a long calculation, inspect input consistency:

```bash
python membit.py \
  -f treated_Protein_Phos.xtc \
  -s treated_Protein_Phos.gro \
  -n membit_index.ndx \
  -thickness 6 1 0 6 25 \
  --diagnose-index
```

## Documentation

See:

- [`docs/input_preparation.md`](docs/input_preparation.md)
- [`docs/index_files.md`](docs/index_files.md)
- [`docs/trajectory_formats.md`](docs/trajectory_formats.md)
- [`docs/usage_examples.md`](docs/usage_examples.md)
- [`docs/troubleshooting.md`](docs/troubleshooting.md)
- [`docs/developer_validation.md`](docs/developer_validation.md)

## Citation

Please cite the original MembIT publication and the repository version used for
your analysis.

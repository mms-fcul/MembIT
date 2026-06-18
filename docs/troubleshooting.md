# Troubleshooting

## `Please decrease the cutoff input value`

This message can mean the cutoff is inappropriate, but it can also be triggered
when the membrane groups are not suitable for the calculation.

Check:

- Are `Monolayer1` and `Monolayer2` leaflet marker atoms?
- Were the trajectory, structure, and index generated for the same atom set?
- Are you accidentally using a reduced index with a full-system trajectory?
- Is the trajectory centered/imaged consistently?
- Does the protein remain in a meaningful position relative to the membrane?

Run:

```bash
python membit.py   -f trajectory.xtc   -s structure.gro   -n membit_index.ndx   -thickness 6 1 0 6 25   --diagnose-index
```

## Monolayer groups contain no marker atoms

This usually means the index was generated against a different structure or atom
numbering scheme. Generate the index again with GROMACS against the same
processed structure/trajectory that will be supplied to MembIT.

## Raw full-system trajectory is slow

This is expected and is not the recommended MembIT workflow. Prepare an
analysis-specific trajectory containing only the atoms MembIT needs.

Complementary GROMACS analyses may still use full-system or solute trajectories
when that is scientifically appropriate.

## XTC/TRR run fails without `-s`

Binary trajectories do not contain all metadata needed by MembIT. Supply a
matching `.gro` or `.tpr`:

```bash
-s structure.gro
```

or:

```bash
-s structure.tpr
```

## Whole-protein insertion is hard to interpret

For proteins with large soluble domains, the whole-protein center of mass may not
describe membrane insertion well. Consider defining a region-specific group, such
as a transmembrane-core group, and compare it to leaflet marker atoms.

## Top/bottom signs look confusing

MembIT raw deformation signs are leaflet-dependent. A common plotting convention
is:

```text
upper_outward = upper_raw
lower_outward = -lower_raw
mean_outward  = 0.5 * (upper_raw - lower_raw)
imbalance     = upper_raw + lower_raw
```

Document the convention used in downstream plots.

# MembIT index files

MembIT uses a GROMACS-style index file. The required groups are:

```text
[ Protein ]
[ Center_of_Interest ]
[ Monolayer1 ]
[ Monolayer2 ]
```

Additional groups can exist, but these four are the groups MembIT needs.

## Protein

The `Protein` group should contain the full protein atoms required to define the
protein/membrane relationship.

## Center_of_Interest

The `Center_of_Interest` group is the group around which local insertion,
thickness, or deformation is calculated. In many workflows this is the same as
`Protein`, but it can be a subset, such as one chain, one domain, one helix, or a
set of residues.

For some biological questions, a custom center of interest is more meaningful
than the whole protein. For example, a transmembrane core may be a better
insertion/tilt reporter than a full protein containing a large soluble domain.

## Monolayer1 and Monolayer2

For thickness/deformation calculations, the monolayer groups should normally
contain leaflet marker atoms. These are atoms that represent the leaflet surface
or headgroup region, such as phosphate/headgroup atoms.

They should not normally contain:

- all atoms from the lipids;
- lipid tails;
- waters or ions;
- atom groups generated against a different structure or trajectory.

A useful rule is:

```text
Monolayer1 + Monolayer2 = all leaflet marker atoms used by MembIT
```

## Analysis-specific indexes

Different analyses may require different indexes.

Examples:

- deformation/thickness:
  - `Protein`
  - `Center_of_Interest`
  - `Monolayer1`
  - `Monolayer2`

- radial thickness:
  - usually the same MembIT index as deformation/thickness;

- region-specific insertion/tilt:
  - custom groups such as `TM_core_CA`, `TM_low_CA`, `TM_high_CA`,
    `Upper_P31`, and `Lower_P31`.

The important requirement is that every index is generated against the same
structure/trajectory atom set used by the calculation.

## Leaflet assignment

One practical approach is:

1. Select one marker atom per lipid, such as `P31`.
2. Sort those atoms by z-coordinate in a representative frame.
3. Find the largest z-gap between consecutive marker atoms.
4. Assign atoms above the gap to one leaflet and below the gap to the other.
5. If using a multi-atom marker group, assign all marker atoms from each lipid to
   the same leaflet as its splitting atom.

## Diagnostics

Run:

```bash
python membit.py   -f trajectory.xtc   -s structure.gro   -n membit_index.ndx   -thickness 6 1 0 6 25   --diagnose-index
```

The diagnostic report checks group sizes, atom-number ranges, common atom names,
common residues, and the presence of common marker atoms such as `O31`, `P31`,
`O32`, `O33`, and `O34`.

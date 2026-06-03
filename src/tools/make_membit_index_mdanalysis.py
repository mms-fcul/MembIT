#!/usr/bin/env python3
"""Create a MembIT index from a structure/trajectory using MDAnalysis.

This helper is meant for the XTC workflow, where MembIT index numbers must be
1-based atom positions in the structure/topology supplied with ``membit.py -s``.
It avoids the reduced-PDB renumbering problem by writing atom.index + 1.

Typical full-system usage:

  python tools/make_membit_index_mdanalysis.py \
      -s min1.tpr \
      -f 1_010-full_5frames.xtc \
      -o membit_index_010_fullsystem.ndx \
      --protein-selection "protein" \
      --coi-selection "protein" \
      --phos-selection "name O31 P31 O32 O33 O34" \
      --split-name P31 \
      --mode Phos

The leaflet split follows the same basic idea as the existing preparation
script: sort phosphate split atoms by z, find the largest z gap, assign atoms
above the gap to Monolayer1 and atoms below the gap to Monolayer2.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def write_group(handle, name: str, atom_indices: list[int], per_line: int = 15) -> None:
    """Write one GROMACS/MembIT-style index group.

    MDAnalysis atom indices are zero-based; MembIT/GROMACS index files are
    one-based, so the caller must pass already-converted atom numbers.
    """
    handle.write(f"[ {name} ]\n")
    for i in range(0, len(atom_indices), per_line):
        handle.write(" ".join(str(v) for v in atom_indices[i:i + per_line]) + "\n")
    handle.write("\n")


def unique_in_order(values: list[int]) -> list[int]:
    """Return values without duplicates while preserving their first order."""
    seen = set()
    out = []
    for value in values:
        if value not in seen:
            out.append(value)
            seen.add(value)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a MembIT index using MDAnalysis atom numbering.")
    parser.add_argument("-s", "--structure", required=True, type=Path, help="Structure/topology file used with MembIT -s")
    parser.add_argument("-f", "--trajectory", default=None, type=Path, help="Optional trajectory; frame 0 is used by default")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output MembIT index")
    parser.add_argument("--frame", type=int, default=0, help="Trajectory frame index used for leaflet assignment")
    parser.add_argument("--protein-selection", default="protein", help="MDAnalysis selection for Protein group")
    parser.add_argument("--coi-selection", default=None, help="MDAnalysis selection for Center_of_Interest; defaults to protein selection")
    parser.add_argument("--phos-selection", default="name O31 P31 O32 O33 O34", help="MDAnalysis selection for phosphate-group atoms")
    parser.add_argument("--split-name", default="P31", help="Atom name used to split the two monolayers")
    parser.add_argument("--mode", choices=["P", "Phos"], default="Phos", help="Monolayer groups contain only split atoms or full phosphate groups")
    args = parser.parse_args()

    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise SystemExit("This tool requires MDAnalysis: python -m pip install MDAnalysis") from exc

    if args.trajectory:
        universe = mda.Universe(str(args.structure), str(args.trajectory))
        universe.trajectory[args.frame]
    else:
        universe = mda.Universe(str(args.structure))

    coi_selection = args.coi_selection or args.protein_selection

    protein = universe.select_atoms(args.protein_selection)
    coi = universe.select_atoms(coi_selection)
    phos = universe.select_atoms(args.phos_selection)
    split_atoms = phos.select_atoms(f"name {args.split_name}")

    if len(protein) == 0:
        raise SystemExit(f"Protein selection is empty: {args.protein_selection!r}")
    if len(coi) == 0:
        raise SystemExit(f"Center_of_Interest selection is empty: {coi_selection!r}")
    if len(phos) == 0:
        raise SystemExit(f"Phosphate selection is empty: {args.phos_selection!r}")
    if len(split_atoms) < 2:
        raise SystemExit(f"Need at least two split atoms named {args.split_name!r}; found {len(split_atoms)}")

    # Largest z-gap leaflet split.  This mirrors the current project script but
    # uses MDAnalysis atom.index values so the output is valid for XTC input.
    split_sorted = sorted(split_atoms, key=lambda atom: atom.position[2])
    largest_gap = -1.0
    gap_low = None
    gap_high = None
    for lower, upper in zip(split_sorted[:-1], split_sorted[1:]):
        gap = abs(float(upper.position[2] - lower.position[2]))
        if gap >= largest_gap:
            largest_gap = gap
            gap_low = float(lower.position[2])
            gap_high = float(upper.position[2])

    if gap_low is None or gap_high is None:
        raise SystemExit("Could not determine a leaflet split gap")

    ml1_p = [atom for atom in split_atoms if float(atom.position[2]) >= gap_high]
    ml2_p = [atom for atom in split_atoms if float(atom.position[2]) <= gap_low]

    if not ml1_p or not ml2_p:
        raise SystemExit(f"Leaflet assignment failed: ML1={len(ml1_p)} ML2={len(ml2_p)}")

    phos_names = set(atom.name for atom in phos)

    def leaflet_numbers(p_atoms):
        numbers = []
        for p_atom in p_atoms:
            if args.mode == "P":
                atoms = [p_atom]
            else:
                # In full-system topology, residue grouping should be reliable.
                # This is safer than relying on compact PDB block order.
                atoms = [atom for atom in p_atom.residue.atoms if atom.name in phos_names]
            numbers.extend(atom.index + 1 for atom in sorted(atoms, key=lambda atom: atom.index))
        return unique_in_order(numbers)

    protein_numbers = [atom.index + 1 for atom in protein]
    coi_numbers = [atom.index + 1 for atom in coi]
    ml1_numbers = leaflet_numbers(ml1_p)
    ml2_numbers = leaflet_numbers(ml2_p)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        write_group(handle, "Protein", protein_numbers)
        write_group(handle, "Center_of_Interest", coi_numbers)
        write_group(handle, "Monolayer1", ml1_numbers)
        write_group(handle, "Monolayer2", ml2_numbers)

    print(f"Wrote {args.output}")
    print(f"Protein atoms             : {len(protein_numbers)}")
    print(f"Center_of_Interest atoms  : {len(coi_numbers)}")
    print(f"Phosphate atoms selected  : {len(phos)}")
    print(f"Split atoms selected      : {len(split_atoms)}")
    print(f"Largest z gap             : {largest_gap:.6f} A")
    print(f"Monolayer1 index atoms    : {len(ml1_numbers)}")
    print(f"Monolayer2 index atoms    : {len(ml2_numbers)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# MembIT

Python script to perform insertion and thickness analysis on lipid bilayer systems

## Requirements
numpy module

## Input Files
- PDB file
A trajectory comprising all atoms included in the INDEX file
- INDEX file
Required groups in the index file:
Protein - all atoms from the inserting molecule.
Center_of_Interest - all atoms of the inserting molecule group
          (residue, motif, atom, etc) whose geometric center will
          be the reference for the insertion calculations.
Monolayer1 - all atoms from one of the monolayers.
Monolayer2 - all atoms from the other monolayer.

## Main Modes
### Insertion
Reports the insertion of the Center_of_Interest.
This insertion is measured as the 2D minimum distances
from the Membrane atoms to the Center_of_Interest.

One can choose between three insertion definitions:
- closest (insertion to closest membrane atom)
- average (insertion to average membrane z position)
- cutoff  (insertion to membrane atoms within a cutoff)

### Thickness
Reports the thickness of each Monolayer.
This thickness is reported according to the radial distance between
the protein and the membrane atoms.
The thickness is defined as the difference between the z coordinate
average of Monolayer1 and Monolayer2 atoms within a given xy window.

### Simplethickness
Reports a difference between the average z of both leaflets


##Example usage:
```
membit.py -f test.pdb -n index.ndx -thickness 10 0.1
```
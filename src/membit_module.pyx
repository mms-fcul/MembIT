#!python
#cython: language_level=2, boundscheck=False, wraparound=False, cdivision=True, profile=False

import numpy as np
cimport numpy as np

"""
TODO: further optimization can be done in saveOutput,
      in saveThickness and in and in the way getThickness handles output.
      (getThickness is creating a string that is splited in saveOutput, this is rather inefficient
       and getInsertion does the same...)
"""

cpdef float PYX2Ddist(float p1_x, float p1_y,
                     float p2_x, float p2_y,
                     float box_x, float box_y):
    cdef float dx, dy

    dx = abs(p1_x - p2_x)
    dy = abs(p1_y - p2_y)

    if not dx < box_x / 2:
        dx = box_x - dx

    if not dy < box_y / 2:
        dy = box_y - dy

    return dx ** 2 + dy ** 2

cpdef float PYX3Ddist(float p1_x, float p1_y, float p1_z,
                     float p2_x, float p2_y, float p2_z,
                     float box_x, float box_y, float box_z):
    cdef float dx, dy, dz

    dx = abs(p1_x - p2_x)
    dy = abs(p1_y - p2_y)
    dz = abs(p1_z - p2_z)

    if not dx < box_x / 2:
        dx = box_x - dx

    if not dy < box_y / 2:
        dy = box_y - dy

    if not dz < box_z / 2:
        dz = box_z - dz

    return dx ** 2 + dy ** 2 + dz ** 2


cpdef PYXgetAllAtomsMinDist2(float[:] molAtoms_x, float[:] molAtoms_y, float[:] molAtoms_z,
                             float[:] membAtoms_x, float[:] membAtoms_y, float[:] membAtoms_z,
                             box, int criteria, membraneAtoms, moleculeAtoms):
    cdef int membraneAtomsLen, moleculeAtomsLen, membatom_i, atom_i
    cdef float minDist, dist
    cdef float box_x, box_y, box_z
    cdef float atom_x, atom_y, atom_z
    cdef float membCoords_x, membCoords_y,membCoords_z
    membraneAtomsLen = len(membAtoms_x)
    moleculeAtomsLen = len(molAtoms_x)

    atoms_distances_arr = np.zeros(membraneAtomsLen, dtype="f")
    cdef float [:] atoms_distances = atoms_distances_arr

    box_x, box_y, box_z = box[0:3]

    for membatom_i in range(membraneAtomsLen):
        membCoords_x = membAtoms_x[membatom_i]
        membCoords_y = membAtoms_y[membatom_i]
        membCoords_z = membAtoms_z[membatom_i]
        minDist = 999999
        for atom_i in range(moleculeAtomsLen):
            atom_x = molAtoms_x[atom_i]
            atom_y = molAtoms_y[atom_i]
            atom_z = molAtoms_z[atom_i]
            if criteria == 2:
                dist = PYX2Ddist(atom_x, atom_y,
                                 membCoords_x, membCoords_y,
                                 box_x, box_y)
            elif criteria == 3:
                dist = PYX3Ddist(atom_x, atom_y, atom_z,
                                 membCoords_x, membCoords_y, membCoords_z,
                                 box_x, box_y, box_z)
            else:
                raise Exception('Wrong criteria it should be either 2 or 3')
            if dist < minDist:
                minDist = dist
                membatom = membraneAtoms[membatom_i]
                atom = moleculeAtoms[atom_i]
        #if membatom.getNumber() == '437':
        #    print 'MinDist2', dist ** 0.5, minDist ** 0.5, atom.getNumber(), \
        #        membatom.get3DPosition(), atom.get3DPosition()

        atoms_distances[membatom_i] = minDist ** 0.5

    return np.asarray(atoms_distances)

cpdef getAllAtomsMinDist2(moleculeAtoms, membraneAtoms, box, criteria):
    """Calculates the minimum distance
    between all membraneAtoms to the moleculeAtoms
    Ensures:
      atoms_distances is list of tuples
      with the membrane atom object and
      its minimum distance to the molecule
    """
    cdef int atom_i
    cdef float x, y, z
    nMolAtoms, nMembAtoms = len(moleculeAtoms), len(membraneAtoms)
    molAtoms_x_arr = np.zeros(nMolAtoms, dtype="f")
    cdef float [:] molAtoms_x = molAtoms_x_arr
    molAtoms_y_arr = np.zeros(nMolAtoms, dtype="f")
    cdef float [:] molAtoms_y = molAtoms_y_arr
    molAtoms_z_arr = np.zeros(nMolAtoms, dtype="f")
    cdef float [:] molAtoms_z = molAtoms_z_arr

    membAtoms_x_arr = np.zeros(nMembAtoms, dtype="f")
    cdef float [:] membAtoms_x = membAtoms_x_arr
    membAtoms_y_arr = np.zeros(nMembAtoms, dtype="f")
    cdef float [:] membAtoms_y = membAtoms_y_arr
    membAtoms_z_arr = np.zeros(nMembAtoms, dtype="f")
    cdef float [:] membAtoms_z = membAtoms_z_arr

    for atom_i in range(nMolAtoms):
        atom = moleculeAtoms[atom_i]
        x, y, z = atom.get3DPosition()
        molAtoms_x[atom_i] = x
        molAtoms_y[atom_i] = y
        molAtoms_z[atom_i] = z

    for atom_i in range(nMembAtoms):
        atom = membraneAtoms[atom_i]
        x, y, z = atom.get3DPosition()
        membAtoms_x[atom_i] = x
        membAtoms_y[atom_i] = y
        membAtoms_z[atom_i] = z

        return PYXgetAllAtomsMinDist2(molAtoms_x, molAtoms_y, molAtoms_z,
                                      membAtoms_x, membAtoms_y, membAtoms_z,
                                      box, criteria, membraneAtoms, moleculeAtoms)
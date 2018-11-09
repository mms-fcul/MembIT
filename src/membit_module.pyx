#!python
#cython: language_level=2, boundscheck=False, wraparound=False, cdivision=True, profile=False

import numpy as np
cimport numpy as np
import os

""" 
TODO: further optimization can be done in saveOutput,
      in saveThickness and in and in the way getThickness handles output.
      (getThickness is creating a string that is splited in saveOutput, this is rather inefficient
       and getInsertion does the same...)
"""

cdef float PYX2Ddist(float p1_x, float p1_y,
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

cdef float PYX3Ddist(float p1_x, float p1_y, float p1_z,
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
                #if membatom.getNumber() == '62':
                #    print 'MinDist2', dist ** 0.5, minDist ** 0.5, atom.getNumber(), \
                #        membatom.get3DPosition(), atom.get3DPosition()

        atoms_distances[membatom_i] = minDist ** 0.5

    return np.asarray(atoms_distances)
    

class Atom:
    def __init__(self, number):
        self._number = number

    def addProperties(self, atype, residue, x, y, z):
        self._atype   = atype
        self._residue = residue
        self._x       = x
        self._y       = y
        self._z       = z

    def addLeaflet(self, leaflet):
        self._leaflet = leaflet

    def getLeaflet(self):
        return self._leaflet

    def setDistance2Protein(self, distance):
        self._distance2protein = distance

    def setDistance2CoI(self, distance):
        self._distance2CoI = distance

    def getNumber(self):
        return self._number

    def getDistance2Protein(self):
        return self._distance2protein

    def getDistance2CoI(self):
        return self._distance2CoI

    def get3DPosition(self):
        return (self._x, self._y, self._z)

    def distTo2D(self, p2_x, p2_y, box_x, box_y):
        p1_x, p1_y = self._x, self._y
        return PYX2Ddist(p1_x, p1_y, p2_x, p2_y, box_x, box_y)

    def distTo3D(self, p2_x, p2_y, p2_z, box_x, box_y, box_z):
        p1_x, p1_y, p1_z = self._x, self._y, self._z
        return PYX3Ddist(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z,
                         box_x, box_y, box_z)


class AtomCollections:
    def __init__(self):
        self._atoms = {}

    def addAtom(self, number):
        newAtom = Atom(number)
        self._atoms[number] = newAtom

    def addProperties(self, number, atype, residue, x, y, z):
        self._atoms[number].addProperties(atype, residue, x, y, z)

    def getAtomsNumbers(self):
        return self._atoms

    def getAtoms(self):
        return self._atoms.values()

    def getSize(self):
        counter = 0
        for atom in self.getAtoms():
            counter += 1
        return counter

    def distToZ(self, reference_membrane, other_membrane,
                inserting_point, box):
        if reference_membrane > other_membrane:
            reference = 'top'
        else:
            reference = 'down'

        if reference == 'top':            
            if inserting_point > other_membrane:
                insertion = inserting_point - reference_membrane
            elif inserting_point < other_membrane:
                insertion = box - (reference_membrane - inserting_point)

        elif reference == 'down':
            if inserting_point < other_membrane:
                insertion = reference_membrane - inserting_point
            elif inserting_point > other_membrane:
                insertion = box - (inserting_point - reference_membrane)

        return insertion


class Membrane(AtomCollections):
    def __init__(self):
        AtomCollections.__init__(self)
        self._membraneOne = []
        self._membraneTwo = []

        self._closestLeaflet = None
        self._otherLeaflet = None
        self._closestAtom = None
        self._furthestAtom = None

        self._half_z = None

        self._topLeaflet = None
        self._bottomLeaflet = None

        self._thickness_top = {}
        self._thickness_bottom = {}

    def addAtom(self, number, leaflet):
        newAtom = Atom(number)
        self._atoms[number] = newAtom
        if leaflet == 'one':
            self._membraneOne.append(newAtom)
            newAtom.addLeaflet('one')
        elif leaflet == 'two':
            self._membraneTwo.append(newAtom)
            newAtom.addLeaflet('two')

    def getLeafletAtoms(self, leaflet):
        if leaflet == 'one':
            return self._membraneOne
        elif leaflet == 'two':
            return self._membraneTwo
        elif leaflet == 'top':
            return self.getLeafletAtoms(self._topLeaflet)
        elif leaflet == 'bottom':
            return self.getLeafletAtoms(self._bottomLeaflet)

    def getClosestLeaflet(self):
        return self._closestLeaflet

    def getClosestAtom(self):
        return self._closestAtom

    def getFurthestAtom(self):
        return self._furthestAtom

    def getOtherLeaflet(self):
        return self._otherLeaflet

    def setHalfZ(self, half_z):
        self._half_z = half_z

    def getHalfZ(self):
        return self._half_z

    def setLeafletsOrder(self, leaftlet, order):
        if order == 'top':
            self._topLeaflet = leaftlet
        if order == 'bottom':
            self._bottomLeaflet = leaftlet

    def chooseClosestLeaflet(self, protein, box, criteria):
        center = protein.getCenter()

        closest = 999999
        furthest = -1
        closest_atom = None
        furthest_atom = None
        for membraneAtom in self.getAtoms():
            if criteria == '2D':
                new_dist = membraneAtom.distTo2D(center[0], center[1], box[0], box[1])
                membraneAtom.setDistance2CoI(new_dist)
            elif criteria == '3D':
                new_dist = membraneAtom.distTo3D(center[0], center[1], center[2], box[0], box[1], box[2])
                membraneAtom.setDistance2CoI(new_dist)

            if new_dist < closest:
                closest = new_dist
                closest_atom = membraneAtom
            elif new_dist > furthest:
                furthest = new_dist
                furthest_atom = membraneAtom

        closest_monolayer = closest_atom.getLeaflet()

        self._closestLeaflet = closest_monolayer
        self._closestAtom = closest_atom
        self._furthestAtom = furthest_atom

        if closest_monolayer == 'one':
            self._otherLeaflet = 'two'
        else:
            self._otherLeaflet = 'one'

    def getAverageZ(self, atoms):
        total = 0.0
        counter = 0
        for atom in atoms:
            counter += 1
            total += atom.get3DPosition()[2]

        return total / counter

    def getReferenceML2CoIDistance(self):
        atoms_distances = []
        membAtoms = self.getLeafletAtoms(self._closestLeaflet)
        for membAtom in membAtoms:
            distance = membAtom.getDistance2CoI() ** 0.5
            atoms_distances.append(distance)
        return (atoms_distances, membAtoms)

    def getAllAtomsMinDist2(self, moleculeAtoms, membraneAtoms, box, criteria):
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

    def getSlices(self, window_size, window_step,
                  min_window, max_window, distances, membrane_atoms):
        cdef int natoms
        cdef float distance, window_begin, window_end, window_half
        window_begin = min_window
        window_end = window_begin + window_size

        natoms = len(membrane_atoms)

        membrane_atoms_2sort = []
        for i in range(natoms):
            membrane_atoms_2sort.append((membrane_atoms[i], distances[i]))
                       
        membrane_atoms_2sort.sort(key=lambda x: x[1])

        #with open('tmpnew', 'w') as f:
        #    text = ''
        #    for atom_i in range(natoms):
        #        atom = membrane_atoms_2sort[atom_i][0]
        #        distance = membrane_atoms_2sort[atom_i][1]
        #        #if atom.getNumber() == '62':
        #        text += '{} {} {}\n'.format('slices_o', atom.getNumber(), round(distance, 3))
        #    f.write(text)


        tmp_window = []
        atoms = []
        #print ('a', min([i[1] for i in membrane_atoms]),
        #       max([i[1] for i in membrane_atoms]))
        for atom_i in range(natoms):
            atom = membrane_atoms_2sort[atom_i][0]
            distance = membrane_atoms_2sort[atom_i][1]
            if distance >= max_window:
                continue

            while distance > window_end:
                # Save tmp_window in windows
                window_half = window_begin + window_size / 2.0
                yield window_half, atoms, tmp_window

                # Update the window
                window_begin += window_step
                window_end += window_step

                # Remove points from new window
                new_tmp_window = []
                new_atoms = []
                i = 0
                while i < len(tmp_window) and tmp_window[i] <= window_end:
                    if tmp_window[i] >= window_begin:
                        new_tmp_window.append(tmp_window[i])
                        new_atoms.append(atoms[i])
                    i += 1

                tmp_window = new_tmp_window
                atoms = new_atoms

            if distance >= window_begin:
                tmp_window.append(distance)
                atoms.append(atom)

        # Deal with the last windows
        if window_step == 0:
            window_step = 999999

        window_end_converted = int(round(window_end, 3) * 1000)
        max_window_converted = int(max_window * 1000)

        while window_end_converted <= max_window_converted:
            new_tmp_window = []
            new_atoms = []
            i = 0
            if len(tmp_window) != 0:
                while i < len(tmp_window) and tmp_window[i] < window_end:
                    if tmp_window[i] > window_begin:
                        new_tmp_window.append(tmp_window[i])
                        new_atoms.append(atoms[i])
                    i += 1

                tmp_window = new_tmp_window
                atoms = new_atoms
                window_half = window_begin + window_size / 2.0
                yield window_half, atoms, tmp_window
            else:
                window_half = window_begin + window_size / 2.0
                yield window_half, [], []

            window_begin += window_step
            window_end += window_step
            window_half = window_begin + window_size / 2.0
            #yield window_half, [], []
            window_end_converted = int(round(window_end, 3) * 1000)
            max_window_converted = int(max_window * 1000)

    def calcHalfMembraneZ(self, protein, parameters, box):
        """
        """
        cdef int atom_i, counter_one, counter_two
        cdef float distance, cutoff, total_one, total_two, oneML_z, twoML_z, half_membrane_z
        nargs_insertion = len(parameters)
        if nargs_insertion == 5:
            cutoff = float(parameters[4])
        else:
            cutoff = 0

        # Get lipids beyond a specified 2D cutoff
        protein = protein.getAtoms()
        membrane = self.getAtoms()
        dists2Protein = self.getAllAtomsMinDist2(protein, membrane, box, 2)
        natoms = len(membrane)

        #for atom_i in range(natoms):
        #    atom = membrane[atom_i]
        #    distance = dists2Protein[atom_i]
        #    if atom.getNumber() == '62':
        #        print "halfmemb", atom.getNumber(), distance
        
        
        counter_one = 0
        counter_two = 0
        total_one = 0
        total_two = 0
        for atom_i in range(natoms):
            atom = membrane[atom_i]
            distance = dists2Protein[atom_i]
            if distance >= cutoff:
                z = atom.get3DPosition()[2]
                if atom.getLeaflet() == 'one':
                    counter_one += 1
                    total_one += z
                else:
                    counter_two += 1
                    total_two += z

        # Get outter layer lipids half membrane z value
        oneML_z = total_one / counter_one
        twoML_z = total_two / counter_two

        if oneML_z < twoML_z:
            self.setLeafletsOrder('one', 'bottom')
            self.setLeafletsOrder('two', 'top')
        else:
            self.setLeafletsOrder('one', 'top')
            self.setLeafletsOrder('two', 'bottom')

        half_membrane_z = (oneML_z + twoML_z) / 2

        self.setHalfZ(half_membrane_z)

    def getSimpleThickness(self, outputnameThickness):
        cdef int ntop, nbottom
        cdef float topZ, bottomZ

        topAtoms = self.getLeafletAtoms('one')
        bottomAtoms = self.getLeafletAtoms('two')

        topZ = 0
        ntop = 0
        for atom in topAtoms:
            z = atom.get3DPosition()[2]
            topZ += z
            ntop += 1
        topZ /= ntop

        bottomZ = 0
        nbottom = 0
        for atom in bottomAtoms:
            z = atom.get3DPosition()[2]
            bottomZ += z
            nbottom += 1
        bottomZ /= nbottom

        output = ' {0} '.format(abs(topZ - bottomZ))

        return output

    def getThickness(self, CoI, leaflet, box,
                     parameters, outputfile, printnatoms):
        box_x = box[0]
        box_y = box[1]
        box_z = box[2]

        box_size = (box_x ** 2 + box_y ** 2) ** 0.5

        nargs_insertion = len(parameters)
        if nargs_insertion >= 4:
            max_window = float(parameters[3])
            min_window = float(parameters[2])
        if nargs_insertion >= 2:
            if nargs_insertion == 3:
                min_window = float(parameters[2])
                max_window = box_size
            elif nargs_insertion == 2:
                min_window = 0                
                max_window = box_size

            window_step = float(parameters[1])
            window_size = float(parameters[0])

            ouput_radius = ''
            if window_size > box_size:
                window_size = int(box_size)
                ouput_radius = 'box_size'

            # Get distances to coi            
            CoI = CoI.getAtomsInLeaflet(leaflet)
            membrane = self.getLeafletAtoms(leaflet)

            dists2Protein = self.getAllAtomsMinDist2(CoI, membrane, box, 2)

            slices = self.getSlices(window_size, window_step,
                                    min_window, max_window,
                                    dists2Protein, membrane)

            #with open('slices', 'w') as f:
            #    text = ''
            #    for window_half, atoms, distances in slices:
            #        if len(distances) > 0:
            #            text += '{} {}\n'.format(window_half, sum(distances)/len(distances))
            #        else:
            #            text += '{} {}\n'.format(window_half, 0)
            #    f.write(text)

            output = ''
            output_natoms = ''
            windows = []
            for window_half, atoms, distances in slices:
                counter = 0
                memb_total = 0.0
                for atom in atoms:
                    atom_z = atom.get3DPosition()[2]
                    atomML = atom.getLeaflet()

                    counter += 1
                    memb_total += atom_z
                if counter > 0:
                    memb_average = memb_total / counter
                    thickness_value = abs(memb_average - self.getHalfZ())
                    thickness_value = round(thickness_value, 3)
                    self.saveThickness(window_half, thickness_value, leaflet)
                    #print window_half, thickness_value, leaflet, memb_average, self.getHalfZ()
                else:
                    thickness_value = "NaN"
                windows.append(window_half)

                if printnatoms:
                    output += ' {0} {1} '.format(thickness_value, counter)
                else:
                    output += ' {0} '.format(thickness_value)

            if not os.path.isfile(outputfile):
                new_output = 'time '
                for i in windows:
                    new_output += '{0} '.format(round(i, 3))
                output = new_output + '\n' + output

            return output

    def saveThickness(self, window, thickness, leaflet):
        if leaflet == 'top':
            dictionary = self._thickness_top
        else:
            dictionary = self._thickness_bottom

        if window in dictionary.keys():
            dictionary[window].append(thickness)
        else:
            dictionary[window] = [thickness]

    def calcThicknessAvg(self):
        avgs_both = []
        windows_both = []

        for thickness_dict in (self._thickness_top, self._thickness_bottom):
            avgs = []
            windows = []
            
            window_list = thickness_dict.keys()
            window_list.sort()
            for window in window_list:
                total = 0
                npoints = 0
                windows.append(window)
                for point in thickness_dict[window]:
                    npoints += 1
                    total += point
                avgs.append(total / npoints)

            avgs_both.append(avgs)
            windows_both.append(windows)

        return avgs_both[0], windows_both[0], avgs_both[1], windows_both[1]

class Protein(AtomCollections):
    def __init__(self):
        AtomCollections.__init__(self)

        self._proteinTop = []
        self._proteinBottom = []
        
        self._center_x = None
        self._center_y = None
        self._center_z = None

    def addAtomToLeaflet(self, atom, leaflet):
        if leaflet == 'top':
            self._proteinTop.append(atom)
        elif leaflet == 'bottom':
            self._proteinBottom.append(atom)

    def clearLeafletAtoms(self):
        self._proteinTop = []
        self._proteinBottom = []

    def getAtomsInLeaflet(self, leaflet):
        if leaflet == 'top':
            atoms = self._proteinTop
        elif leaflet == 'bottom':
            atoms = self._proteinBottom

        return atoms

    def calcCenter(self):
        counter = 0
        x_total = 0
        y_total = 0
        z_total = 0
        for atom in self.getAtoms():
            x, y, z = atom.get3DPosition()
            x_total += x
            y_total += y
            z_total += z

            counter += 1

        center_x = x_total / counter
        center_y = y_total / counter
        center_z = z_total / counter

        self._center_x = center_x
        self._center_y = center_y
        self._center_z = center_z

    def getCenter(self):
        return (self._center_x, self._center_y, self._center_z)

    def getInsertion(self, membrane, parameters, box, outputfile, traj):
        center = self.getCenter()
        center_x = center[0]
        center_y = center[1]
        center_z = center[2]

        box_x = box[0]
        box_y = box[1]
        box_z = box[2]

        box_size = (box_x ** 2 + box_y ** 2) ** 0.5

        if parameters[0] == 'zero':
            membrane_center = membrane.getHalfZ()

            insertion_value = center_z - membrane_center

            return str(round(insertion_value, 3))

        other_ml = membrane.getOtherLeaflet()
        other_ml_atoms = membrane.getLeafletAtoms(other_ml)
        mean_other_memb = membrane.getAverageZ(other_ml_atoms)

        nargs_insertion = len(parameters)
        if nargs_insertion >= 4:
            max_window = float(parameters[3])
            min_window = float(parameters[2])
        if nargs_insertion >= 2:
            #closest_atom = membrane.getClosestAtom()
            #furthest_atom = membrane.getFurthestAtom()

            #closest_distance = closest_atom.getDistance2CoI() ** 0.5
            #furthest_distance = furthest_atom.getDistance2CoI()  ** 0.5

            if nargs_insertion == 3:
                min_window = float(parameters[2])
                max_window = box_size
            elif nargs_insertion == 2:
                min_window = 0                
                max_window = box_size

            window_step = float(parameters[1])
            window_size = float(parameters[0])

            # get reference membrane atoms
            refMLdists2CoI, membAtoms = membrane.getReferenceML2CoIDistance()
            slices = membrane.getSlices(window_size, window_step,
                                        min_window, max_window,
                                        refMLdists2CoI, membAtoms)
            output = ''
            windows = []
            for window_half, atoms, distances in slices:
                total = 0.0
                counter = 0
                for atom in atoms:
                    counter += 1
                    total += atom.get3DPosition()[2]
                if counter > 0:
                    average_z = total / counter
                    insertion_value = self.distToZ(average_z,
                                                   mean_other_memb,
                                                   center_z, box_z)
                    insertion_value = round(insertion_value, 3)
                    natoms = counter
                else:
                    natoms = 'closest'
                    if 'noNaN' in parameters:
                        closest_atom = membrane.getClosestAtom()
                        closest_atom_z = closest_atom.get3DPosition()[2]
                        insertion_value = self.distToZ(closest_atom_z, mean_other_memb,
                                                       center_z, box_z)
                        insertion_value = round(insertion_value, 3)
                    else:
                        average_z = 0
                        insertion_value = "NaN"
                windows.append(window_half)

                if traj._printnatoms:
                    output += ' {0} {1} '.format(insertion_value, natoms)
                else:
                    output += ' {0} '.format(insertion_value)

            if traj.getInsertionOutput() == '':
                new_output = 'time '
                for i in windows:
                    if traj._printnatoms:
                        new_output += '{0} natoms '.format(round(i, 3))
                    else:
                        new_output += '{0} '.format(round(i, 3))
                output = new_output + '\n' + output

            return output

        elif parameters[0] == 'closest':
            closest_atom = membrane.getClosestAtom()
            closest_atom_z = closest_atom.get3DPosition()[2]

            insertion_value = self.distToZ(closest_atom_z, mean_other_memb,
                                           center_z, box_z)
            return str(round(insertion_value, 3))

        elif parameters[0] == 'average':
            closest_ml = membrane.getClosestLeaflet()
            closest_ml_atoms = membrane.getLeafletAtoms(closest_ml)
            average_z = membrane.getAverageZ(closest_ml_atoms)

            insertion_value = self.distToZ(average_z, mean_other_memb,
                                           center_z, box_z)
            return str(round(insertion_value, 3))

    def IndexandTrajAtomsMatch(self):
        counter = 0
        for i in self.getAtomsNumbers():
            counter += 1

        positions = 0
        for i in self._atoms.values():
            try:
                i.get3DPosition()
                positions += 1
            except:
                pass

        if positions == counter:
            return True
        else:
            return False

    def calcAtomsClosestML(self, membrane):
        half_z = membrane.getHalfZ()
        for atom in self.getAtoms():
            atom_z = atom.get3DPosition()[2]
            if atom_z < half_z:
                self.addAtomToLeaflet(atom, 'bottom')
            else:
                self.addAtomToLeaflet(atom, 'top')


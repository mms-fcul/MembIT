import numpy as np
cimport numpy as np
import os
from atom import Atom
from membit_module import PYXgetAllAtomsMinDist2
from atomcollections import AtomCollections

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

    def getReferenceML2CoIDistanceSorted(self):
        atoms_distances = []
        membAtoms = self.getLeafletAtoms(self._closestLeaflet)
        for membAtom in membAtoms:
            distance = membAtom.getDistance2CoI() ** 0.5
            atoms_distances.append(distance)

        membrane_atoms_sorted = []
        for i, atom in enumerate(membAtoms):
            membrane_atoms_sorted.append((membAtoms[i], atoms_distances[i]))

        membrane_atoms_sorted.sort(key=lambda x: x[1])
        return membrane_atoms_sorted


    def getSlices(self, window_size, window_step,
                  min_window, max_window, membrane_atoms_sorted):
        cdef int natoms
        cdef float distance, window_begin, window_end, window_half
        window_begin = min_window
        window_end = window_begin + window_size

        natoms = len(membrane_atoms_sorted)

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
            atom = membrane_atoms_sorted[atom_i][0]
            distance = membrane_atoms_sorted[atom_i][1]

            #if atom.getNumber() == '437':
            #    print '437_slices', distance

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
                #if atom.getNumber() == '437':
                #    window_half = window_begin + window_size / 2.0
                #    print '437_slices_within', distance, window_begin, window_end, window_half


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

        if counter_one + counter_two == 0:
            raise IOError('Please decrease the cutoff input value.')

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

            with open("log.txt", 'a') as f:
                f.write('\n####')
            output = ''
            output_natoms = ''
            windows = []
            for window_half, atoms, distances in slices:
                counter = 0
                memb_total = 0.0

                text = ''

                for atom in atoms:
                    atom_z = atom.get3DPosition()[2]
                    atomML = atom.getLeaflet()

                    counter += 1
                    memb_total += atom_z


                    text += '{0} '.format(atom.getNumber())

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

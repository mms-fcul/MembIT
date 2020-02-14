import os

class Trajectory:
    def __init__(self, trajfile, indexfile, distance_criteria,
                 outputfile, thickness, insertion, printnatoms):
        """Instanciates a Trajectory object and checks some input the
        consistency of the input arguments

        Requires: 
        trajfile 
        indexfile 
        outputfile 
        thickness 
        insertion

        Ensures:
        The input arguments are correctly assigned to the attributes,
        considering the help messages provided to the user
        """

        self._trajfile   = trajfile
        self._indexfile  = indexfile
        self._distance_criteria = distance_criteria

        if outputfile:
            self._outputfile = outputfile
        else:
            self._outputfile = None

        self._printnatoms = printnatoms
        self._thickness = thickness
        if thickness:
            self._thicknessOutput1 = ''
            self._thicknessOutput2 = ''
            nargs_thickness = len(thickness)
            if nargs_thickness < 2:
                raise IOError('The thickness argument should have at least 2'
                              ' fields (the window size and step)')
            elif nargs_thickness > 5:
                raise IOError('The thickness argument should have at most 5 '
                              'fields (the window size, step, minimum and '
                              'and maximum values)')

        self._insertion = insertion
        if insertion:
            self._insertionOutput = ''            
            nargs_insertion = len(insertion)
            if nargs_insertion == 1:
                if insertion[0] == 'closest' or \
                   insertion[0] == 'average' or \
                   insertion[0] == 'zero':
                    self._insertion_window = insertion[0]
                else:
                    raise IOError('The insertion modes allowed are "closest" '
                                  'if the reference atom to calculate the '
                                  'insertion is the closest membrane atom and '
                                  '"average" if the reference to is the '
                                  'average z position of the membrane atoms.\n'
                                  'Usage: -insertion closest\n'
                                  '                  average\n'
                                  '                  1 0.1 0 10')
            elif nargs_insertion > 5:
                raise IOError('The insertion argument should have at most 5 '
                              'fields')

        self._curtime = None
        self._box = None
        self._protein = Protein()
        self._CoI = Protein()
        self._membrane = Membrane()

        self.loadIndex()

        proteinCounter = 0
        coiCounter = 0
        for i in self._protein.getAtomsNumbers():
            proteinCounter += 1

        for i in self._CoI.getAtoms():
            coiCounter += 1

        if self._insertion and coiCounter < 1:
            raise IOError('The provided index file should have at least one '
                          'atom belonging to the Center_of_Interest group')
        elif self._thickness and proteinCounter < 1:
            raise IOError('The provided index file should have at least one '
                          'atom belonging to the Protein group')

        top_memb_size = len(self._membrane.getLeafletAtoms('one'))
        bottom_memb_size = len(self._membrane.getLeafletAtoms('two'))

        if top_memb_size < 1 or bottom_memb_size < 1:
            raise IOError('The provided index file should have at least one '
                          'atom in both Monolayer1 and '
                          'Monolayer2 groups')

        if not insertion and not thickness:
            raise IOError('This script can calculate thickness and insertion '
                          'provided you use the -thickness or -insertion '
                          'arguments respectively')
    def getInsertionOutput(self):
        return self._insertionOutput

    def analyseTrajectory(self):
        traj = self.loadTrajectory()

        if self._insertion:
            outputnameInsertion = self.getOutputName("insertion")
            os.system('rm -f {0}'.format(outputnameInsertion))

        if self._thickness:
            outputnameThickness1 = self.getOutputName("thickness1")
            os.system('rm -f {0}'.format(outputnameThickness1))

            outputnameThicknessAvg1 = self.getOutputName("thickness1_avg")
            os.system('rm -f {0}'.format(outputnameThicknessAvg1))

            outputnameThickness2 = self.getOutputName("thickness2")
            os.system('rm -f {0}'.format(outputnameThickness2))

            outputnameThicknessAvg2 = self.getOutputName("thickness2_avg")
            os.system('rm -f {0}'.format(outputnameThicknessAvg2))

        for frame in traj:
            if self._insertion:
                # Calculate geometric center of Center_of_Interest
                self._CoI.calcCenter()

                if 'zero' == self._insertion[0]:
                    # Calculate the Membrane Half Z
                    self._membrane.calcHalfMembraneZ(self._protein,
                                                     (0, 0, 0, 0,
                                                      self._insertion[1]),
                                                     self._box)
                else:
                    # Choose leaflet
                    self._membrane.chooseClosestLeaflet(self._CoI,
                                                        self._box,
                                                        self._distance_criteria)


                # Calculate insertion
                insertion = self._CoI.getInsertion(self._membrane,
                                                   self._insertion,
                                                   self._box,
                                                   outputnameInsertion,
                                                   self)

                # Save to Output
                self.saveOutput(outputnameInsertion, insertion)

            if self._thickness:
                # Calculate the Membrane Half Z
                self._membrane.calcHalfMembraneZ(self._protein,
                                                self._thickness,
                                                self._box)

                # Attribution of the Protein atoms to membrane
                # leaflets ('bottom' and 'top')
                self._CoI.calcAtomsClosestML(self._membrane)

                # Calculate the Thickness for ML1
                thickness1 = self._membrane.getThickness(self._CoI,
                                                         'top',
                                                         self._box,
                                                         self._thickness,
                                                         outputnameThickness1,
                                                         self._printnatoms)

                # Calculate the Thickness for ML2
                thickness2 = self._membrane.getThickness(self._CoI,
                                                         'bottom',
                                                         self._box,
                                                         self._thickness,
                                                         outputnameThickness2,
                                                         self._printnatoms)
                self._CoI.clearLeafletAtoms()
                # Save the Outputs
                self.saveOutput(outputnameThickness1, thickness1)
                self.saveOutput(outputnameThickness2, thickness2)


        # Write to Output
        if self._thickness:
            self.writeOutput(outputnameThickness1)
            self.writeOutput(outputnameThickness2)

            avgs_top, windows_top,\
                avgs_bottom, windows_bottom = self._membrane.calcThicknessAvg()

            self.writeAvgOutput(outputnameThicknessAvg1, avgs_top,
                                windows_top)
            self.writeAvgOutput(outputnameThicknessAvg2, avgs_bottom,
                                windows_bottom)

        if self._insertion:
            self.writeOutput(outputnameInsertion)


    def loadIndex(self):
        with open(self._indexfile) as f:
            addTo = None
            for line in f:
                line = line.strip()
                if '[ ' in line and ' ]' in line:
                    indexName = line.replace('[', '').replace(']', '')
                    indexName = indexName.replace(' ', '').lower()
                    if 'protein' == indexName:
                        addTo = 'protein'

                    elif 'center_of_interest' == indexName:
                        addTo = 'center_of_interest'

                    elif 'monolayer1' == indexName:
                        addTo = 'monolayer1'

                    elif 'monolayer2' == indexName:
                        addTo = 'monolayer2'

                    else:
                        addTo = None

                elif addTo:
                    for atomNumber in line.split():
                        if addTo == 'protein':
                            self._protein.addAtom(atomNumber)

                        elif addTo == 'center_of_interest':
                            self._CoI.addAtom(atomNumber)

                        elif addTo == 'monolayer1':
                            self._membrane.addAtom(atomNumber, 'one')

                        elif addTo == 'monolayer2':
                            self._membrane.addAtom(atomNumber, 'two')

    def loadTrajectory(self):        
        with open(self._trajfile) as f:
            for line in f:
                if line[0:4] == 'ATOM':
                    fields = line.split()

                    number  = fields[1]
                    atype   = fields[2]
                    residue = fields[4]
                    x       = float(fields[5])
                    y       = float(fields[6])
                    z       = float(fields[7])

                    if number in self._protein.getAtomsNumbers():
                        self._protein.addProperties(number, atype,
                                                    residue, x, y, z)

                    if number in self._CoI.getAtomsNumbers():
                        self._CoI.addProperties(number, atype,
                                                residue, x, y, z)

                    elif number in self._membrane.getAtomsNumbers():
                        self._membrane.addProperties(number, atype,
                                                     residue, x, y, z)
                elif line[0:6] == 'CRYST1':
                    fields = line.split()
                    box_x = float(fields[1])
                    box_y = float(fields[2])
                    box_z = float(fields[3])
                    self._box = box_x, box_y, box_z

                elif line[0:5] == 'TITLE':
                    line = line.strip()
                    time = line.split('t=')[1].replace(' ', '')
                    self._curtime = int(float(time))

                elif line[0:3] == 'TER':
                    if not self._CoI.IndexandTrajAtomsMatch():
                        raise IOError('Index file not correct. CoI group atoms in the index do '
                                      'not match the trajectory file')
                    if not self._protein.IndexandTrajAtomsMatch():
                        raise IOError('Index file not correct. Protein group atoms in the index do '
                                      'not match the trajectory file')

                    yield

    def getOutputName(self, prefix):
        if self._outputfile:
            outputname = '{0}_{1}.xvg'.format(self._outputfile, prefix)
        else:
            outputname = '{0}.xvg'.format(prefix)

        return outputname

    def saveOutput(self, outputname, data):
        if data[:4] == 'time':
            line = ''
        else:
            line = '{0:9f} '.format(self._curtime)

        nNaNs = 0
        for value in data.split(' '):
            if value == '\n':
                line += '\n{0:9f}\t'.format(self._curtime)
            else:
                line += '{0:5s} '.format(value)
                if value == 'NaN':
                    nNaNs += 1

        data_type = outputname.split('_')[-1].replace('.xvg', '')
        if data_type == 'insertion':
            self._insertionOutput += line + '\n'
        elif data_type == 'thickness1':
            # If all NaNs don't save the data            
            if nNaNs != (len(data.split(' ')) - 2 )/ 3:
                self._thicknessOutput1 += line + '\n'
        elif data_type == 'thickness2':
            # If all NaNs don't save the data
            if nNaNs != (len(data.split(' ')) - 2 )/ 3:
                self._thicknessOutput2 += line + '\n'

    def writeOutput(self, outputname):
        data_type = outputname.split('_')[-1].replace('.xvg', '')
        if data_type == 'insertion':
            data = self._insertionOutput
        elif data_type == 'thickness1':
            data = self._thicknessOutput1
        elif data_type == 'thickness2':
            data = self._thicknessOutput2

        with open(outputname, 'w') as f:
            if len(data) == 0:
                f.write('No occurrences in this monolayer\n')
            else:
                f.write(data)

    def writeAvgOutput(self, outputname, avgs, windows):
        text = ''
        with open(outputname, 'w') as f:
            for i in range(len(windows)):
                text += '{0:9} {1:9}\n'.format(windows[i], avgs[i])
            if len(text) == 0:
                text = 'No occurrences in this monolayer\n'
            f.write(text)

class Atom:
    def __init__(self, number):
        self._number  = number

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

    def distTo2D(self, point, box):
        dx = abs(self._x - point[0])
        dy = abs(self._y - point[1])

        box_x = box[0]
        box_y = box[1]

        if not dx < box_x / 2:
            dx = box_x - dx

        if not dy < box_y / 2:
            dy = box_y - dy

        return dx ** 2 + dy ** 2

    def distTo3D(self, point, box):
        dx = abs(self._x - point[0])
        dy = abs(self._y - point[1])
        dz = abs(self._z - point[2])

        box_x = box[0]
        box_y = box[1]
        box_z = box[2]

        if not dx < box_x / 2:
            dx = box_x - dx

        if not dy < box_y / 2:
            dy = box_y - dy

        if not dz < box_z / 2:
            dz = box_z - dz

        return dx ** 2 + dy ** 2 + dz ** 2

class AtomCollections:
    def __init__(self):
        self._atoms = {}

    def addAtom(self, number):
        newAtom = Atom(number)
        self._atoms[number] = newAtom

    def addProperties(self, number, atype, residue, x, y, z):
        self._atoms[number].addProperties(atype, residue, x, y, z)

    def getAtomsNumbers(self):
        for atom in self._atoms:
            yield atom

    def getAtoms(self):
        output_list = []
        for atom in self._atoms.values():
            output_list.append(atom)
        return output_list

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
                new_dist = membraneAtom.distTo2D(center, box)                
                membraneAtom.setDistance2CoI(new_dist)
            elif criteria == '3D':
                new_dist = membraneAtom.distTo3D(center, box)                
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
        for membAtom in self.getLeafletAtoms(self._closestLeaflet):
            distance = membAtom.getDistance2CoI() ** 0.5
            atoms_distances.append((membAtom, distance))
        return atoms_distances

    def getAllAtomsMinDist2(self, moleculeAtoms, membraneAtoms, box, criteria):
        """Calculates the minimum distance
        between all membraneAtoms to the moleculeAtoms
        Ensures:
          atoms_distances is list of tuples
          with the membrane atom object and 
          its minimum distance to the molecule
        """
        atoms_distances = []
        for membatom in membraneAtoms:
            membCoords = membatom.get3DPosition()
            minDist = 999999
            for atom in moleculeAtoms:                
                if criteria == '2D':
                    dist = atom.distTo2D(membCoords, box)
                else:
                    dist = atom.distTo3D(membCoords, box)
                if dist < minDist:
                    minDist = dist
            atoms_distances.append((membatom, minDist ** 0.5))

        return atoms_distances

    def getSlices(self, window_size, window_step,
                  min_window, max_window, membrane_atoms):
        window_begin = min_window
        window_end = window_begin + window_size

        membrane_atoms.sort(key=lambda x: x[1])
        tmp_window = []
        atoms = []
        #print ('a', min([i[1] for i in membrane_atoms]),
        #       max([i[1] for i in membrane_atoms]))
        for atom, distance in membrane_atoms:
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
        nargs_insertion = len(parameters)
        if nargs_insertion == 5:
            cutoff = float(parameters[4])
        else:
            cutoff = 0

        # Get lipids beyond a  specified 2D cutoff
        protein = protein.getAtoms()
        membrane = self.getAtoms()
        dists2Protein = self.getAllAtomsMinDist2(protein, membrane, box, '2D')

        counter_one = 0
        counter_two = 0
        total_one = 0
        total_two = 0
        for atom, distance in dists2Protein:                
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

            dists2Protein = self.getAllAtomsMinDist2(CoI, membrane, box, '2D')

            slices = self.getSlices(window_size, window_step,
                                    min_window, max_window,
                                    dists2Protein)
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
                else:
                    thickness_value = "NaN"
                windows.append(window_half)

                if printnatoms:
                    output += ' {0} {1} '.format(thickness_value, counter)
                else:
                    output += ' {0} '.format(thickness_value)

            with open(outputfile, 'a') as f:
                if os.path.getsize(outputfile) == 0:
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
            refMLdists2CoI = membrane.getReferenceML2CoIDistance()
            slices = membrane.getSlices(window_size, window_step,
                                        min_window, max_window,
                                        refMLdists2CoI)
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

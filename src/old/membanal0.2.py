import argparse
import os

"""
Os slices tem de ser fixos para nao variarem com os frames
Nan quando o slice nao tem pontos

na documentacao falar no problema da precisao das floats 
que apenas foi corrigida na condicao das janelas de output

ver pq e que o step 0 na insercao nao printa a primeira janela sequer
provavelmente o bug persiste na thickness
"""


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='\n'
                                 'Script to perform insertion and thickness calculations on '
                                 ' lipid bilayer systems')
parser.add_argument('-f', help='PDB trajectory file',
                    required=True, metavar='traj.pdb')

# Protein distancia minima a todos os atomos - so para thickness
# Center_of_Interest centro geometrico - so para insertion
parser.add_argument('-n', help='Required groups in the index file: \n'
                    'Protein - includes all atoms from the inserting molecule\n'
                    '          it is used to determine which membrane atoms will\n'
                    '          be considered for the thickness calculations\n'
                    'Center_of_Interest - includes all atoms of the inserting molecule\n'
                    '          group (residue, motif, atom, etc) whose geometric center will \n'
                    '          be the reference for the insertion calculations\n'
                    'Membrane - Monolayer1 + Monolayer2 atoms\n'
                    'Monolayer1 - all atoms from one of the monolayers\n'
                    'Monolayer2 - all atoms from the other monolayer',
                    required=True, metavar='index.ndx')
parser.add_argument('-o', help='Ouput identifier name\n'
                    'Ex: -o analysis -> filaname if -insertion = analysis_insertion.xvg\n'
                    '                -> filaname if -thickness = analysis_thickness.xvg', required=False,
                    metavar='analysis', default='')

# min max are optional and by default it should use the minimum distance to P
#and the maximum distance to P, accordingly
parser.add_argument('-thickness', help='Thickness parameters:\n'
                    'All window related distances are 2D minimum distances\n'
                    'from the Membrane atoms to the Protein Atoms.\n'
                    'The thickness is defined as the difference between the z coordinate average of '
                    'Monolayer1 and Monolayer2 atoms within a given xy window.'
                    '<window_size> <window_step> <min> <max>\n'
                    'window_size - output window size in Angstrom\n'
                    'window_step - moving window step in Angstrom\n'
                    'min - minimum distance between Membrane and Protein to be considered\n'
                    '      (the default is 0)\n'
                    'max - maximum distance between Membrane and Protein to be considered\n'
                    '      (the default is the box size in xy)\n'
                    'min and max are optional', required=False,
                    metavar='window step', default=None, nargs='+')
parser.add_argument('-insertion', help='Insertion paramenters:\n'
                    '<type> or <window> <step> <min> <max>\n'
                    'type - closest (insertion to closest membrane atom)\n'
                    '       average (insertion to average membrane z position)\n'
                    'All window related distances are 2D minimum distances\n'
                    'from the Membrane atoms to the geometric center of Center_of_Interest atoms.\n'
                    'The insertion is defined as the difference between the z coordinates of said'
                    ' geometric center and the average of the closest Monolayer atoms within a '
                    'given xy window.\n'
                    'window_size - output window size in Angstrom\n'
                    'window_step - moving window step in Angstrom\n'
                    'min - minimum distance between Membrane and Center_of_Interest to be considered\n'
                    '      (the default is 0)\n'
                    'max - maximum distance between Membrane and Center_of_Interest to be considered\n'                    
                    '      (the default is box_size in xy)\n'
                    'min and max are optional', required=False,
                    metavar='closest', default=None, nargs='+')


args = parser.parse_args()

# insertion
# escolha da monolayer p mais proximo 3d ao centro geometrico do Center_of_Interest
# media z do slice da monolayer - centro geometrico do Center_of_Interest

# thickness
# escolha dos P por distancia minima 2D de cada P a todos os atomos em Protein
# thickness = media z monolayer cima - media z monolayer baixo

# output insertion
# time reference_layer slice1 slice2 ... slicen

# output thickness
# time slice1 slice2 slice3 ... slicen      FILE1
# and
# slice1 average_thickness                  FILE2


class Trajectory:
    def __init__(self, trajfile, indexfile,
                 outputfile, thickness, insertion):
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

        if outputfile:
            self._outputfile = outputfile
        else:
            self._outputfile = None

        self._thickness = thickness
        if thickness:
            self._thicknessOutput = ''
            nargs_thickness = len(thickness)
            if nargs_thickness < 2:
                raise IOError('The thickness argument should have at least 2'
                              ' fields (the window size and step)')
            elif nargs_thickness > 4:
                raise IOError('The thickness argument should have at most 4 '
                              'fields (the window size, step, minimum and '
                              'and maximum values)')

        self._insertion = insertion
        if insertion:
            self._insertionOutput = ''            
            nargs_insertion = len(insertion)
            if nargs_insertion == 1:
                if insertion[0] == 'closest' or \
                   insertion[0] == 'average':
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
            elif nargs_insertion > 4:
                raise IOError('The insertion argument should have at most 4 '
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
            outputnameThickness = self.getOutputName("thickness")
            os.system('rm -f {0}'.format(outputnameThickness))

            outputnameThicknessAvg = self.getOutputName("thickness_avg")
            os.system('rm -f {0}'.format(outputnameThicknessAvg))

        for frame in traj:
            if self._insertion:
                # Calculate geometric center of Center_of_Interest
                self._CoI.calcCenter()

                # Choose leaflet
                self._membrane.chooseClosestLeaflet(self._CoI,
                                                    self._box)

                # Calculate the insertion
                insertion = self._CoI.getInsertion(self._membrane,
                                                   self._insertion,
                                                   self._box,
                                                   outputnameInsertion,
                                                   self)

                # Save to Output
                self.saveOutput(outputnameInsertion, insertion)

            if self._thickness:
                # Calculate the thickness
                thickness = self._membrane.getThickness(self._protein,
                                                        self._thickness,
                                                        self._box,
                                                        outputnameThickness)

                # Save the Output
                self.saveOutput(outputnameThickness, thickness)

        # Write to Output
        if self._thickness:
            self.writeOutput(outputnameThickness)

            avgs, windows = self._membrane.calcThicknessAvg()
            self.writeAvgOutput(outputnameThicknessAvg, avgs, windows)

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

        for value in data.split(' '):
            if value == '\n':
                line += '\n{0:9f}\t'.format(self._curtime)
            else:
                line += '{0:5s} '.format(value)

        data_type = outputname.split('_')[-1].replace('.xvg', '')
        if data_type == 'insertion':
            self._insertionOutput += line + '\n'

        elif data_type == 'thickness':
            self._thicknessOutput += line + '\n'

    def writeOutput(self, outputname):
        data_type = outputname.split('_')[-1].replace('.xvg', '')
        if data_type == 'insertion':
            data = self._insertionOutput
        elif data_type == 'thickness':
            data = self._thicknessOutput

        with open(outputname, 'w') as f:
            f.write(data)

    def writeAvgOutput(self, outputname, avgs, windows):
        text = ''
        with open(outputname, 'w') as f:
            for i in range(len(windows)):
                text += '{0:9} {1:9}\n'.format(windows[i], avgs[i])
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

        #if self._x == 51.82 and self._y == 40.38 and self._z == 29.0 and\
        #   point[0] == 53.65 and point[1] == 98.55 and point[2] == 25.02:
        #    print 'DISTANCES'
        #    print dx, box_x / 2
        #    print dy, box_y / 2
        #    print dz, box_z / 2

        if not dx < box_x / 2:
            dx = box_x - dx

        if not dy < box_y / 2:
            dy = box_y - dy

        if not dz < box_z / 2:
            dz = box_z - dz

        #if self._x == 51.82 and self._y == 40.38 and self._z == 29.0 and\
        #   point[0] == 53.65 and point[1] == 98.55 and point[2] == 25.02:
        #    print dx, dy, dz
        #    print 'END DISTANCES'

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
        for atom in self._atoms.values():
            yield atom

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

        self._thickness = {}
    def addAtom(self, number, leaflet):
        newAtom = Atom(number)
        self._atoms[number] = newAtom
        if leaflet == 'one':
            self._membraneOne.append(newAtom)
        elif leaflet == 'two':
            self._membraneTwo.append(newAtom)

    def getLeafletAtoms(self, leaflet):
        if leaflet == 'one':
            return self._membraneOne
        elif leaflet == 'two':
            return self._membraneTwo

    def getLeafletOf(self, atom):
        if atom in self._membraneOne:
            return 'one'
        elif atom in self._membraneTwo:
            return 'two'

    def getClosestLeaflet(self):
        return self._closestLeaflet

    def getClosestAtom(self):
        return self._closestAtom

    def getFurthestAtom(self):
        return self._furthestAtom

    def getOtherLeaflet(self):
        return self._otherLeaflet

    def chooseClosestLeaflet(self, protein, box):
        center = protein.getCenter()

        closest = 999999
        furthest = -1
        closest_atom = None
        furthest_atom = None
        for membraneAtom in self.getAtoms():
            new_dist = membraneAtom.distTo3D(center, box)
            new_dist2D = membraneAtom.distTo2D(center, box)
            membraneAtom.setDistance2CoI(new_dist2D)

            if new_dist < closest:
                closest = new_dist
                closest_atom = membraneAtom
            elif new_dist > furthest:
                furthest = new_dist
                furthest_atom = membraneAtom

        closest_monolayer = self.getLeafletOf(closest_atom)

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

    def getAllAtomsMinDist2(self, molecule, box):
        atoms_distances = []
        for membatom in self.getAtoms():
            membCoords = membatom.get3DPosition()
            minDist = 999999
            for atom in molecule.getAtoms():
                # 3D distance criteria
                #dist = atom.distTo3D(membCoords, box)
                # 2D distance criteria
                dist = atom.distTo2D(membCoords, box)
                if dist < minDist:
                    minDist = dist
                    idatom = atom
                #if str(membatom.getNumber()) == '5973' and \
                #   str(atom.getNumber()) == '164':
                #    print '164 - 5973'
                #    print atom.get3DPosition()
                #    print membatom.get3DPosition()
                #    print dist
                #    print '########################'
            atoms_distances.append((membatom, minDist ** 0.5))
            #if str(membatom.getNumber()) == '470':
            #    print minDist
            #    print idatom.getNumber(), idatom.get3DPosition()
            #    exit()

        return atoms_distances

    def getSlices(self, window_size, window_step,
                  min_window, max_window, membrane_atoms):
        window_begin = min_window
        window_end = window_begin + window_size

        membrane_atoms.sort(key=lambda x: x[1])
        tmp_window = []
        atoms = []
        for atom, distance in membrane_atoms:
            #print distance, window_begin, window_end
            #print 'distance >= max_window', distance >= max_window
            if distance >= max_window:
                continue

            #print 'distance > window_end', distance > window_end
            while distance > window_end:
                # Save tmp_window in windows
                window_half = window_begin + window_size / 2
                yield window_half, atoms, tmp_window
                #print 'yield ', window_half, tmp_window

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
                #print 'appended'
                tmp_window.append(distance)
                atoms.append(atom)

        # Deal with the last windows
        if window_step == 0:
            window_step = 999999

        window_end_converted = int(round(window_end, 3) * 1000)
        max_window_converted = int(max_window * 1000)
        #print window_end_converted <= max_window_converted
        while window_end_converted <= max_window_converted:
            #print window_end_converted
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
                window_half = window_begin + window_size / 2
                yield window_half, atoms, tmp_window
            else:
                window_half = window_begin + window_size / 2
                yield window_half, [], []

            window_begin += window_step
            window_end += window_step
            window_half = window_begin + window_size / 2
            #yield window_half, [], []
            window_end_converted = int(round(window_end, 3) * 1000)
            max_window_converted = int(max_window * 1000)


    def getThickness(self, protein, parameters, box, outputfile):
        """
        """
        box_x = box[0]
        box_y = box[1]
        box_z = box[2]

        box_size = (box_x ** 2 + box_y ** 2) ** 0.5

        nargs_insertion = len(parameters)
        if nargs_insertion == 4:
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

            # Get lipids within specified window
            dists2Protein = self.getAllAtomsMinDist2(protein, box)
            slices = self.getSlices(window_size, window_step,
                                    min_window, max_window,
                                    dists2Protein)
            output = ''
            windows = []
            for window_half, atoms, distances in slices:
                #print '######## WINDOW ', window_half, '########'
                counter_one = 0
                counter_two = 0
                membOne_total = 0.0
                membTwo_total = 0.0
                for atom in atoms:
                    atom_z = atom.get3DPosition()[2]
                    atomML = self.getLeafletOf(atom)

                    if atomML == 'one':
                        #print 'ATOM LF1 ', atom.getNumber(), atom.get3DPosition(), distances[counter_one+counter_two]
                        counter_one += 1
                        membOne_total += atom_z
                    elif atomML == 'two':
                        #print 'ATOM LF2 ', atom.getNumber(), atom.get3DPosition(), distances[counter_one+counter_two]
                        counter_two += 1
                        membTwo_total += atom_z

                if counter_one > 0 and counter_two > 0:
                    membOne_average = membOne_total / counter_one
                    membTwo_average = membTwo_total / counter_two
                    thickness_value = abs(membOne_average - membTwo_average)
                    thickness_value = round(thickness_value, 3)
                    #print ''
                    #print 'Leaftlet One Average', 'Leaftlet Two Average'
                    #print membOne_average, membTwo_average
                    #print 'Thickness -> ', thickness_value


                    self.saveThickness(window_half, thickness_value)
                else:
                    thickness_value = "NaN"
                windows.append(window_half)

                output += ' {0} '.format(thickness_value)

            with open(outputfile, 'a') as f:
                if os.path.getsize(outputfile) == 0:
                    new_output = 'time '
                    for i in windows:
                        new_output += '{0} '.format(round(i, 3))
                    output = new_output + '\n' + output

            return output

    def saveThickness(self, window, thickness):
        if window in self._thickness.keys():
            self._thickness[window].append(thickness)
        else:
            self._thickness[window] = [thickness]

    def calcThicknessAvg(self):
        avgs = []
        windows = []

        window_list = self._thickness.keys()
        window_list.sort()
        for window in window_list:
            total = 0
            npoints = 0
            windows.append(window)
            for point in self._thickness[window]:
                npoints += 1
                total += point
            avgs.append(total / npoints)

        return avgs, windows

class Protein(AtomCollections):
    def __init__(self):
        AtomCollections.__init__(self)
        self._center_x = None
        self._center_y = None
        self._center_z = None

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

        other_ml = membrane.getOtherLeaflet()
        other_ml_atoms = membrane.getLeafletAtoms(other_ml)
        mean_other_memb = membrane.getAverageZ(other_ml_atoms)

        nargs_insertion = len(parameters)
        if nargs_insertion == 4:
            max_window = float(insertion[3])
            min_window = float(insertion[2])
        if nargs_insertion >= 2:
            #closest_atom = membrane.getClosestAtom()
            #furthest_atom = membrane.getFurthestAtom()

            #closest_distance = closest_atom.getDistance2CoI() ** 0.5
            #furthest_distance = furthest_atom.getDistance2CoI()  ** 0.5

            if nargs_insertion == 3:
                min_window = float(insertion[2])
                max_window = box_size
            elif nargs_insertion == 2:
                min_window = 0                
                max_window = box_size

            window_step = float(insertion[1])
            window_size = float(insertion[0])

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
                else:
                    average_z = 0
                    insertion_value = "NaN"
                windows.append(window_half)

                output += ' {0} '.format(insertion_value)

            with open(outputfile, 'a') as f:
                if traj.getInsertionOutput() == '':
                    new_output = 'time '
                    for i in windows:
                        new_output += '{0} '.format(round(i, 3))
                    output = new_output + '\n' + output

            return output

        elif insertion[0] == 'closest':
            closest_atom = membrane.getClosestAtom()
            closest_atom_z = closest_atom.get3DPosition()[2]

            insertion_value = self.distToZ(closest_atom_z, mean_other_memb,
                                           center_z, box_z)
            return str(round(insertion_value, 3))

        elif insertion[0] == 'average':
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

if __name__ == '__main__':
    trajfile = args.f
    indexfile = args.n
    outputfile = args.o

    thickness = args.thickness
    insertion = args.insertion

    membanal = Trajectory(trajfile, indexfile,
                          outputfile, thickness, insertion)
    membanal.analyseTrajectory()

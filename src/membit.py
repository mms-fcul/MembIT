#! /usr/bin/python

import argparse
from protein import Protein
from membrane import Membrane
from atom import Atom
import os

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='\n'
                                 'Script to perform insertion and thickness calculations on '
                                 ' lipid bilayer systems')
parser.add_argument('-f', help='PDB trajectory file',
                    required=True, metavar='traj.pdb')

# Protein distancia minima a todos os atomos - so para thickness
# Center_of_Interest centro geometrico - so para insertion
parser.add_argument('-n', help='Required groups in the index file: \n'
                    'Protein - includes all atoms from the inserting molecule.\n'
                    '          it is used to determine which membrane atoms will\n'
                    '          be considered bulk membrane for the thickness calculations\n'
                    '          and also insertion calculations relative to the center of the membrane("zero")\n'
                    'Center_of_Interest - includes all atoms of the inserting molecule\n'
                    '          group (residue, motif, atom, etc) whose geometric center will \n'
                    '          be the reference for the insertion calculations. \n'
                    '          In thickness calculations all atoms of this index group \n'
                    '          will be used. \n'
                    'Monolayer1 - all atoms from one of the monolayers\n'
                    'Monolayer2 - all atoms from the other monolayer',
                    required=True, metavar='index.ndx')
parser.add_argument('-o', help='Ouput identifier name\n'
                    'Ex: -o analysis -> filaname if -insertion = analysis_insertion.xvg\n'
                    '                -> filaname if -thickness = analysis_thickness.xvg',
                    required=False, metavar='analysis', default='')

parser.add_argument('-simplethickness',
                    help='Reports a difference between the average z of both leaflets\n',
                    required=False, action='store_true')

# min max are optional and by default it should use the minimum distance to P
#and the maximum distance to P, accordingly
parser.add_argument('-thickness', help='Thickness parameters:\n'
                    ''
                    'All window related distances are 2D minimum distances\n'
                    'from the Membrane atoms to the Center_of_Interest Atoms.\n'
                    'The thickness is defined as the difference between the z coordinate average of '
                    'Monolayer1 and Monolayer2 atoms within a given xy window.\n'
                    '<window_size> <window_step> <min> <max> <cutoff>\n'
                    'window_size - output window size in Angstrom\n'
                    'window_step - moving window step in Angstrom\n'
                    'min - minimum distance between Membrane and Center_of_Interest to be considered\n'
                    '      (the default is 0)\n'
                    'max - maximum distance between Membrane and Center_of_Interest to be considered\n'
                    '      (the default is the box size in xy)\n'
                    'cutoff - membrane lipids within this cutoff will be ignored from the calculation'
                    ' of the center of the membrane, since it should only include "bulk" membrane atoms.\n'
                    '      (the default is 0, thus including all membrane atoms)\n', required=False,
                    metavar='window step', default=None, nargs='+')
parser.add_argument('-insertion', help='Insertion paramenters:\n'
                    '<type> or <window> <step> <min> <max> <noNaN|min><nclosest>\n'
                    'type - closest (insertion to closest membrane atom)\n'
                    '       average (insertion to average membrane z position)\n'
                    '       zero    (insertion to the center of the "bulk" membrane)\n'
                    '               requires a cutoff from which a bulk membrane is considered\n'
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
                    'noNaN - replaces NaN output entries where there are no membranes atoms in the \n'
                    '        specified cut off with the insertion relative to the closest atom \n'
                    'min - instead of noNaN the "min" option may be chosen. in this case the number of atoms \n'
                    '      specified in nclosest defines the minimum number of membrane atoms to be used in \n'
                    '      the insertion calculation. While noNaN is only trigger when there are no atoms \n'
                    '      within a given radius, min is always used.'
                    'min, max and noNaN|min are optional\n'
                    'nclosest - this argument can only be used with noNaN or min. It specifies the number of \n'
                    '           the closest membrane atoms to include in the calculation \n'
                    '           (if there is no membrane atom within the specified cutoff). ', required=False,
                    metavar='closest', default=None, nargs='+')

parser.add_argument('-distance', help='All distances between the membrane and the center_of_interest \n'
                    'will be calculated using 2 or 3 dimensions. The default is 3. \n'
                    'This is only for insertion, thickness is always 2D.'
                    'The choice of the closest membrane leaflet is based solely on the 3D distances.',
                    choices=['3D', '2D'], required=False, default='3D')

parser.add_argument('-printnatoms', help='Adds a column to the insertion output with the '
                    'number of membrane atoms reported', required=False, default=False, action='store_true')

parser.add_argument('-printclosestleaflet', help='Adds a column to the insertion output with the '
                    'membrane leaflet chosen as reference', required=False, default=False, action='store_true')


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
    def __init__(self, trajfile, indexfile, distance_criteria,
                 outputfile, thickness, simplethickness,
                 insertion, printnatoms):
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
        self._simplethickness = simplethickness

        if thickness and simplethickness:
            raise IOError('Incompatible arguments: simplethickness and thickness.')

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
        elif simplethickness:
            self._thicknessOutput = ''

        self._insertion = insertion
        if insertion:
            self._insertionOutput = ''
            nargs_insertion = len(insertion)

            if insertion[0] == 'closest' or \
               insertion[0] == 'average':
                if nargs_insertion == 1:
                    self._insertion_window = insertion[0]
                else:
                    print 'Warning: Extra arguments have been '\
                        'submitted and will be ignored'

            elif insertion[0] == 'zero':
                if nargs_insertion == 2:
                    self._insertion_window = insertion[0]
                elif nargs_insertion == 1:
                    raise IOError('Cutoff missing. The center of the '
                                  'membrane requires the definition of a cutoff '
                                  'beyond which bulk properties are assumed.')
                else:
                    print 'Warning: Extra arguments have been '\
                        'submitted and will be ignored'

            else:
                if nargs_insertion < 2:
                    raise IOError('The insertion argument requires '
                                  'at least 2 fields (window_size and step)')

                elif nargs_insertion > 5:
                    raise IOError('The insertion argument should '
                                  'have at most 5 fields')

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

        if not insertion and not (thickness or simplethickness):
            raise IOError('This script can calculate thickness and insertion '
                          'provided you use the -thickness or -insertion '
                          'arguments respectively')

    def getInsertionOutput(self):
        return self._insertionOutput

    def analyseTrajectory(self):
        def createOutputFile(filename):
            outputname = self.getOutputName(filename)
            os.system('rm -f {0}'.format(outputname))
            return outputname

        traj = self.loadTrajectory()

        if self._insertion:
            outputnameInsertion = createOutputFile("insertion")

        if self._thickness:
            outputnameThicknessTop    = createOutputFile("thicknessTop")
            outputnameThicknessAvg1 = createOutputFile("thicknessTop_avg")
            outputnameThicknessBottom    = createOutputFile("thicknessBottom")
            outputnameThicknessAvg2 = createOutputFile("thicknessBottom_avg")

        if self._simplethickness:
            outputnameThickness = createOutputFile("thickness")

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

                if args.printclosestleaflet:
                    insertion = '{0} {1}'.format(insertion, self._membrane._closestLeaflet)

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
                thicknessTop = self._membrane.getThickness(self._CoI,
                                                         'top',
                                                         self._box,
                                                         self._thickness,
                                                         outputnameThicknessTop,
                                                         self._printnatoms)

                # Calculate the Thickness for ML2
                thicknessBottom = self._membrane.getThickness(self._CoI,
                                                         'bottom',
                                                         self._box,
                                                         self._thickness,
                                                         outputnameThicknessBottom,
                                                         self._printnatoms)
                self._CoI.clearLeafletAtoms()
                # Save the Outputs
                self.saveOutput(outputnameThicknessTop, thicknessTop)
                self.saveOutput(outputnameThicknessBottom, thicknessBottom)

            if self._simplethickness:
                # Calculate the Membrane Thickness
                thickness = self._membrane.getSimpleThickness(outputnameThickness)

                # Save the Outputs
                self.saveOutput(outputnameThickness, thickness)

        # Write to Output
        if self._insertion:
            self.writeOutput(outputnameInsertion)

        if self._thickness:
            self.writeOutput(outputnameThicknessTop)
            self.writeOutput(outputnameThicknessBottom)

            avgs_top, windows_top,\
                avgs_bottom, windows_bottom = self._membrane.calcThicknessAvg()

            self.writeAvgOutput(outputnameThicknessAvg1, avgs_top,
                                windows_top)
            self.writeAvgOutput(outputnameThicknessAvg2, avgs_bottom,
                                windows_bottom)

        if self._simplethickness:
            self.writeOutput(outputnameThickness)


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
        def readLine(line):
            atype    = line[12:16].strip()
            residue  = line[23:26]
            x        = float(line[30:38])
            y        = float(line[38:46])
            z        = float(line[46:54])
            return atype, residue, x, y, z

        proteinAtoms  = self._protein.getAtomsNumbers()
        CoIAtoms      = self._CoI.getAtomsNumbers()
        membraneAtoms = self._membrane.getAtomsNumbers()
        with open(self._trajfile) as f:
            for line in f:
                if line[0:4] == 'ATOM':
                    number = line[4:11].strip()

                    if number in proteinAtoms:
                        atype, residue, x, y, z = readLine(line)
                        self._protein.addProperties(number, atype,
                                                    residue, x, y, z)

                    if number in CoIAtoms:
                        atype, residue, x, y, z = readLine(line)
                        self._CoI.addProperties(number, atype,
                                                residue, x, y, z)

                    elif number in membraneAtoms:
                        atype, residue, x, y, z = readLine(line)
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
                    time = line.split('t=')[1].split()[0]
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
                line = '{0}\n{1:9f}\t'.format(line, self._curtime)
            else:
                line = '{0}{1:5s} '.format(line, value)
                if value == 'NaN':
                    nNaNs += 1

        data_type = outputname.split('_')[-1].replace('.xvg', '')
        if data_type == 'insertion':
            self._insertionOutput += line + '\n'
        elif data_type == 'thicknessTop':
            # If all NaNs don't save the data
            if nNaNs != (len(data.split(' ')) - 2 )/ 3:
                self._thicknessOutput1 += line + '\n'
        elif data_type == 'thicknessBottom':
            # If all NaNs don't save the data
            if nNaNs != (len(data.split(' ')) - 2 )/ 3:
                self._thicknessOutput2 += line + '\n'
        elif data_type == 'thickness':
            self._thicknessOutput += line + '\n'

    def writeOutput(self, outputname):
        data_type = outputname.split('_')[-1].replace('.xvg', '')
        if data_type == 'insertion':
            data = self._insertionOutput
        elif data_type == 'thicknessTop':
            data = self._thicknessOutput1
        elif data_type == 'thicknessBottom':
            data = self._thicknessOutput2
        elif data_type == 'thickness':
            data = self._thicknessOutput

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


if __name__ == '__main__':
    trajfile = args.f
    indexfile = args.n
    outputfile = args.o

    simplethickness = args.simplethickness
    thickness = args.thickness
    insertion = args.insertion

    distance_criteria = args.distance

    printnatoms = args.printnatoms

    traj = Trajectory(trajfile, indexfile, distance_criteria,
                      outputfile, thickness, simplethickness,
                      insertion, printnatoms)

    traj.analyseTrajectory()

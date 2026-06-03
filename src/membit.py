#! /usr/bin/python

import argparse
import time
from contextlib import contextmanager
from protein import Protein
from membrane import Membrane
from atom import Atom
import os
from pathlib import Path


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description='\n'
                                 'Script to perform insertion and thickness calculations on '
                                 ' lipid bilayer systems')
# ``-f`` remains the main trajectory argument so old PDB commands keep working.
# New binary trajectory formats such as XTC are detected from this file name
# unless the user overrides detection with ``--format``.
parser.add_argument('-f', help='Trajectory file. PDB is supported natively; XTC/TRR/DCD/NC require -s/--structure and MDAnalysis.',
                    required=True, metavar='traj.pdb|traj.xtc')

# GROMACS XTC files contain coordinates but do not contain atom names, residue
# names, or enough topology information for MembIT to map index numbers to
# atoms.  MDAnalysis therefore needs a matching structure/topology file.
# This option is intentionally optional so legacy PDB-only usage remains:
#     python membit.py -f traj.pdb -n index.ndx ...
parser.add_argument('-s', '--structure',
                    help='Structure/topology file required for non-PDB trajectories, e.g. GRO or TPR when -f is XTC.',
                    required=False, metavar='structure.gro|structure.tpr', default=None)

# Auto-detection is convenient for day-to-day use, while an explicit format is
# useful for debugging, unusual file extensions, or scripted regression tests.
parser.add_argument('-format', '--format',
                    help='Trajectory format. Default: auto-detect from -f extension.',
                    choices=['auto', 'pdb', 'xtc', 'trr', 'dcd', 'nc'],
                    required=False, default='auto')

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
                    '<window_size> <window_step> <window_min> <window_max> <cutoff>\n'
                    'window_size - output window size in Angstrom\n'
                    'window_step - moving window step in Angstrom\n'
                    'window_min - minimum distance between Membrane and Center_of_Interest to be considered\n'
                    '      (the default is 0)\n'
                    'window_max - maximum distance between Membrane and Center_of_Interest to be considered\n'
                    '      (the default is the box size in xy)\n'
                    'cutoff - membrane lipids within this cutoff will be ignored from the calculation'
                    ' of the center of the membrane, since it should only include "bulk" membrane atoms.\n'
                    '      (the default is 0, thus including all membrane atoms)\n', required=False,
                    metavar='window step', default=None, nargs='+')
parser.add_argument('-deformation', help='The deformation flag replaces the thickness command (required) output profile with the local deformation profile\n'
                    'for each monolayer. The performed calculation uses the bulk lipids (>cutoff radius) to define the bulk monolayer\n'
                    'half-thickness. The local deformation is defined by the difference between the bulk thickness and the annulus thickness\n'
                    'for each trajectory frame.\n'
                    'Example:\n'
                    'python membit.py -f example.pdb -n template.ndx -o out -thickness 1 1 0 40 15 -deformation\n',  required=False, default=False, action='store_true')
parser.add_argument('-insertion', help='Insertion paramenters:\n'
                    '<type> or <window> <window_step> <window_min> <window_max> <noNaN|min><nclosest>\n'
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
                    'window_min - minimum distance between Membrane and Center_of_Interest to be considered\n'
                    '      (the default is 0)\n'
                    'window_max - maximum distance between Membrane and Center_of_Interest to be considered\n'
                    '      (the default is box_size in xy)\n'
                    'noNaN - replaces NaN output entries where there are no membranes atoms in a \n'
                    '        specific window with the insertion relative to the closest atom \n'
                    'min - instead of noNaN the "min" option may be chosen. in this case the number of atoms \n'
                    '      specified in nclosest defines the minimum number of membrane atoms to be used in \n'
                    '      the insertion calculation. While noNaN is only trigger when there are no atoms \n'
                    '      within a given radius, min is always used.'
                    'window_min, window_max and noNaN|min are optional\n'
                    'nclosest - this argument can only be used with noNaN or min. It specifies the number of \n'
                    '           the closest membrane atoms to include in the calculation \n'
                    '           (if there is no membrane atom in a specific cutoff). ', required=False,
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

parser.add_argument('--profile-timing',
                    help='Print a timing summary at the end of the run. This is intended for performance debugging; it does not change calculated outputs.',
                    required=False, default=False, action='store_true')


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
    def __init__(self, trajfile, indexfile, structurefile, traj_format,
                 distance_criteria, outputfile, thickness, deformation,
                 simplethickness, insertion, printnatoms, profile_timing=False):
        """Instanciates a Trajectory object and checks some input the
        consistency of the input arguments

        Requires:
        trajfile
        indexfile
        outputfile
        thickness
        deformation
        insertion

        Ensures:
        The input arguments are correctly assigned to the attributes,
        considering the help messages provided to the user
        """

        # ------------------------------------------------------------------
        # Optional profiling support
        # ------------------------------------------------------------------
        # ``--profile-timing`` is intentionally read-only: it records timing
        # information but does not alter calculations, output files, or data
        # flow.  Keeping this instrumentation inside Trajectory makes it easy
        # to profile both the original PDB reader and the new MDAnalysis reader
        # with the same command-line flag.
        self._profileTimingEnabled = profile_timing

        # Accumulated wall-clock seconds per labelled code section.  A normal
        # dict preserves insertion order in modern Python, which keeps the
        # final report stable and readable.
        self._profileTimings = {}

        # Number of times each labelled section was entered.  For example,
        # ``trajectory_reader_frame`` should match the number of analysed
        # frames, while ``mda_universe_init`` should normally be called once.
        self._profileCounts = {}

        # Frame counter used only for the profiling report.  The actual MembIT
        # calculations still use the trajectory time stored in ``self._curtime``.
        self._profileFrameCounter = 0

        self._trajfile = trajfile
        self._indexfile = indexfile

        # ``structurefile`` is ignored by the native PDB reader.  It is only
        # required for MDAnalysis-backed formats because binary trajectories do
        # not carry atom/residue metadata by themselves.
        self._structurefile = structurefile

        # Store the resolved format once at construction time so the rest of
        # the code can simply dispatch to the correct reader.
        self._traj_format = self._detectTrajectoryFormat(trajfile, traj_format)
        self._distance_criteria = distance_criteria

        if outputfile:
            self._outputfile = outputfile
        else:
            self._outputfile = None

        self._printnatoms = printnatoms
        self._thickness = thickness
        self._deformation = deformation
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
                    print('Warning: Extra arguments have been '\
                        'submitted and will be ignored')

            elif insertion[0] == 'zero':
                if nargs_insertion == 2:
                    self._insertion_window = insertion[0]
                elif nargs_insertion == 1:
                    raise IOError('Cutoff missing. The center of the '
                                  'membrane requires the definition of a cutoff '
                                  'beyond which bulk properties are assumed.')
                else:
                    print('Warning: Extra arguments have been '\
                        'submitted and will be ignored')

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

        # Index loading is usually small compared with trajectory analysis, but
        # it is measured separately because future full-system indexes may be
        # substantially larger than the reduced PDB indexes used so far.
        with self._profileSection('load_index'):
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

    @contextmanager
    def _profileSection(self, label):
        """Measure wall time for a named code section when profiling is enabled.

        The method is deliberately lightweight and safe to leave around normal
        production code.  When ``--profile-timing`` is not used, the context
        manager simply yields without recording anything.  This avoids changing
        the behavior or output of existing MembIT workflows.
        """
        if not self._profileTimingEnabled:
            yield
            return

        start = time.perf_counter()
        try:
            yield
        finally:
            elapsed = time.perf_counter() - start
            self._profileTimings[label] = self._profileTimings.get(label, 0.0) + elapsed
            self._profileCounts[label] = self._profileCounts.get(label, 0) + 1

    def _printProfileTimingReport(self):
        """Print a compact timing table for performance debugging.

        The report is printed only when ``--profile-timing`` is enabled.  It is
        sent to standard output so it can be captured easily with shell
        redirection or ``tee`` in the local benchmarking scripts.
        """
        if not self._profileTimingEnabled:
            return

        total = self._profileTimings.get('total_wall_clock', 0.0)
        if total <= 0.0:
            total = sum(self._profileTimings.values())

        print('')
        print('MembIT timing profile')
        print('=====================')
        print('Trajectory file : {0}'.format(self._trajfile))
        print('Structure file  : {0}'.format(self._structurefile if self._structurefile else 'None'))
        print('Reader format   : {0}'.format(self._traj_format))
        print('Frames analysed : {0}'.format(self._profileFrameCounter))
        print('')
        print('{0:<36s} {1:>12s} {2:>10s} {3:>12s}'.format('section', 'seconds', 'calls', 'percent'))
        print('{0:<36s} {1:>12s} {2:>10s} {3:>12s}'.format('-' * 36, '-' * 12, '-' * 10, '-' * 12))

        for label, seconds in self._profileTimings.items():
            calls = self._profileCounts.get(label, 1)
            percent = (100.0 * seconds / total) if total > 0.0 else 0.0
            print('{0:<36s} {1:12.6f} {2:10d} {3:11.2f}%'.format(label, seconds, calls, percent))

        print('')
        print('Notes:')
        print('  trajectory_reader_frame measures the time spent obtaining and populating one frame before analysis.')
        print('  frame_analysis_total measures the insertion/thickness/deformation calculations after a frame is loaded.')
        print('  Some nested sections overlap by design, so percentages are diagnostic rather than additive.')

    def getInsertionOutput(self):
        return self._insertionOutput

    def analyseTrajectory(self):
        """Run the selected MembIT analyses over every trajectory frame.

        The original implementation used a simple ``for frame in traj`` loop.
        For profiling, the loop is written explicitly with ``next()`` so that
        the time spent by the trajectory reader can be measured separately from
        the time spent by the analysis calculations.  With profiling disabled,
        the behavior and generated XVG files remain unchanged.
        """
        with self._profileSection('total_wall_clock'):
            def createOutputFile(filename):
                outputname = self.getOutputName(filename)
                os.system('rm -f {0}'.format(outputname))
                return outputname

            with self._profileSection('trajectory_reader_create'):
                traj = self.loadTrajectory()
                traj_iter = iter(traj)

            if self._insertion:
                outputnameInsertion = createOutputFile("insertion")

            if self._thickness:
                outputnameThicknessTop = createOutputFile("thicknessTop")
                outputnameThicknessAvg1 = createOutputFile("thicknessTop_avg")
                outputnameThicknessBottom = createOutputFile("thicknessBottom")
                outputnameThicknessAvg2 = createOutputFile("thicknessBottom_avg")

            if self._simplethickness:
                outputnameThickness = createOutputFile("thickness")

            while True:
                # ``next(traj_iter)`` performs the reader-specific frame work:
                # PDB parsing up to the next TER record, or MDAnalysis frame
                # advance plus MembIT Atom population for XTC/TRR/DCD/NC.
                # Measuring it separately shows whether the bottleneck is input
                # handling or the actual membrane analysis.
                with self._profileSection('trajectory_reader_frame'):
                    try:
                        next(traj_iter)
                    except StopIteration:
                        break

                self._profileFrameCounter += 1

                with self._profileSection('frame_analysis_total'):
                    if self._insertion:
                        with self._profileSection('insertion_total'):
                            # Calculate geometric center of Center_of_Interest.
                            with self._profileSection('insertion_calc_center'):
                                self._CoI.calcCenter()

                            if 'zero' == self._insertion[0]:
                                # Calculate the Membrane Half Z.
                                with self._profileSection('insertion_calc_half_z'):
                                    self._membrane.calcHalfMembraneZ(self._protein,
                                                                     (0, 0, 0, 0,
                                                                      self._insertion[1]),
                                                                     self._box)
                            else:
                                # Choose the closest leaflet used as insertion reference.
                                with self._profileSection('insertion_choose_leaflet'):
                                    self._membrane.chooseClosestLeaflet(self._CoI,
                                                                        self._box,
                                                                        self._distance_criteria)

                            # Calculate insertion.
                            with self._profileSection('insertion_calculation'):
                                insertion = self._CoI.getInsertion(self._membrane,
                                                                   self._insertion,
                                                                   self._box,
                                                                   outputnameInsertion,
                                                                   self)

                            if args.printclosestleaflet:
                                insertion = '{0} {1}'.format(insertion, self._membrane._closestLeaflet)

                            # Save to in-memory output buffer.
                            with self._profileSection('save_insertion_output'):
                                self.saveOutput(outputnameInsertion, insertion)

                    if self._thickness:
                        with self._profileSection('thickness_total'):
                            # Calculate the membrane reference plane / bulk half-thickness.
                            with self._profileSection('thickness_calc_half_z'):
                                self._membrane.calcHalfMembraneZ(self._protein,
                                                                 self._thickness,
                                                                 self._box)

                            # Attribute Center_of_Interest atoms to membrane leaflets
                            # ('bottom' and 'top') before calculating local profiles.
                            with self._profileSection('thickness_assign_coi_leaflets'):
                                self._CoI.calcAtomsClosestML(self._membrane)

                            # Calculate the thickness/deformation profile for ML1.
                            with self._profileSection('thickness_top_calculation'):
                                thicknessTop = self._membrane.getThickness(self._CoI,
                                                                           'top',
                                                                           self._box,
                                                                           self._thickness,
                                                                           outputnameThicknessTop,
                                                                           self._printnatoms,
                                                                           self._deformation)

                            # Calculate the thickness/deformation profile for ML2.
                            with self._profileSection('thickness_bottom_calculation'):
                                thicknessBottom = self._membrane.getThickness(self._CoI,
                                                                              'bottom',
                                                                              self._box,
                                                                              self._thickness,
                                                                              outputnameThicknessBottom,
                                                                              self._printnatoms,
                                                                              self._deformation)

                            self._CoI.clearLeafletAtoms()

                            # Save the outputs to the in-memory buffers.  Files are
                            # written only after all frames are processed, matching
                            # the original MembIT behavior.
                            with self._profileSection('save_thickness_output'):
                                self.saveOutput(outputnameThicknessTop, thicknessTop)
                                self.saveOutput(outputnameThicknessBottom, thicknessBottom)

                    if self._simplethickness:
                        with self._profileSection('simplethickness_total'):
                            # Calculate the Membrane Thickness.
                            thickness = self._membrane.getSimpleThickness(outputnameThickness)

                            # Save the output to the in-memory buffer.
                            self.saveOutput(outputnameThickness, thickness)

            # Write accumulated outputs to disk.  This is measured separately
            # because very long trajectories may spend non-trivial time writing
            # large XVG tables.
            with self._profileSection('write_output_total'):
                if self._insertion:
                    self.writeOutput(outputnameInsertion)

                if self._thickness:
                    self.writeOutput(outputnameThicknessTop)
                    self.writeOutput(outputnameThicknessBottom)

                    with self._profileSection('thickness_average_calculation'):
                        avgs_top, windows_top,\
                            avgs_bottom, windows_bottom = self._membrane.calcThicknessAvg()

                    self.writeAvgOutput(outputnameThicknessAvg1, avgs_top,
                                        windows_top)
                    self.writeAvgOutput(outputnameThicknessAvg2, avgs_bottom,
                                        windows_bottom)

                if self._simplethickness:
                    self.writeOutput(outputnameThickness)

        self._printProfileTimingReport()

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


    @staticmethod
    def _detectTrajectoryFormat(trajfile, traj_format):
        """Return the trajectory reader to use.

        The historical MembIT workflow uses PDB trajectories, so PDB remains
        the default behavior whenever ``-f`` ends in ``.pdb``.  Other supported
        extensions are read through MDAnalysis.  This keeps the old parser and
        output behavior isolated from the new XTC support.
        """
        if traj_format != 'auto':
            return traj_format.lower()

        suffix = Path(trajfile).suffix.lower().lstrip('.')
        if suffix in ('pdb', 'ent'):
            return 'pdb'
        if suffix in ('xtc', 'trr', 'dcd', 'nc'):
            return suffix

        raise IOError(
            'Could not auto-detect trajectory format from extension {0!r}. '
            'Use --format pdb or --format xtc.'.format(Path(trajfile).suffix)
        )

    def loadTrajectory(self):
        """Yield frames from the selected trajectory reader.

        The analysis code below this method expects the Protein, CoI and
        Membrane atom containers to be populated for the current frame before a
        bare ``yield`` happens.  Both readers follow that contract, which keeps
        the insertion/thickness/deformation calculations unchanged.
        """
        if self._traj_format == 'pdb':
            return self.loadPDBTrajectory()
        return self.loadMDAnalysisTrajectory()

    def loadPDBTrajectory(self):
        """Native PDB trajectory reader used by the original MembIT workflow.

        This code is deliberately kept as close as possible to the legacy
        implementation.  That makes it easier to verify that adding XTC support
        has not changed existing PDB behavior.
        """
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


    def loadMDAnalysisTrajectory(self):
        """Read XTC/TRR/DCD/NC trajectories through MDAnalysis.

        Atom-number convention for this reader:
            index atom number N -> MDAnalysis atom with zero-based index N - 1

        In practice this means the MembIT index must use the same 1-based atom
        numbering as the structure/topology file supplied with ``-s``.  For a
        full-system XTC plus full-system TPR/GRO, use a full-system MembIT index.
        For a reduced XTC/GRO containing only Protein+Phos atoms, use an index
        generated for that reduced structure.  Do not mix reduced-PDB numbering
        with a full-system XTC unless the numbers have been remapped.
        """
        if not self._structurefile:
            raise IOError(
                'A structure/topology file is required for {0} trajectories. '
                'Use -s structure.gro or -s structure.tpr.'.format(self._traj_format.upper())
            )

        with self._profileSection('mda_import'):
            try:
                import MDAnalysis as mda
            except ImportError as exc:
                raise ImportError(
                    'Reading {0} trajectories requires MDAnalysis. Install it with: '
                    'python -m pip install MDAnalysis'.format(self._traj_format.upper())
                ) from exc

        # MDAnalysis combines the topology/structure file and the trajectory
        # into one Universe.  Coordinates are updated in-place as we iterate
        # through universe.trajectory, so selected Atom objects can be cached.
        with self._profileSection('mda_universe_init'):
            universe = mda.Universe(self._structurefile, self._trajfile)

        proteinAtoms = self._protein.getAtomsNumbers()
        CoIAtoms = self._CoI.getAtomsNumbers()
        membraneAtoms = self._membrane.getAtomsNumbers()
        required_numbers = set(proteinAtoms) | set(CoIAtoms) | set(membraneAtoms)

        # Cache only the atoms requested in the MembIT index.  This avoids a
        # Python-level scan over every atom in a full system for every frame.
        with self._profileSection('mda_atom_cache'):
            atoms_by_number = {}
            natoms = len(universe.atoms)
            for number in sorted(required_numbers, key=lambda value: int(value)):
                atom_index = int(number) - 1
                if atom_index < 0 or atom_index >= natoms:
                    raise IOError(
                        'Index atom number {0} is outside the structure atom range 1..{1}. '
                        'For MDAnalysis/XTC input, MembIT index files must use 1-based '
                        'structure atom numbers, matching GROMACS .ndx convention.'.format(number, natoms)
                    )
                atoms_by_number[number] = universe.atoms[atom_index]

        for ts in universe.trajectory:
            if ts.dimensions is None or len(ts.dimensions) < 3:
                raise IOError('Trajectory frame has no unit-cell dimensions; MembIT requires box vectors.')

            # MDAnalysis reports GROMACS-like coordinate files in Angstrom, which
            # is the unit expected by the existing MembIT calculations.
            self._box = float(ts.dimensions[0]), float(ts.dimensions[1]), float(ts.dimensions[2])

            # Existing XVG output stores time as an integer.  Preserve that
            # behavior so PDB and XTC paths can be compared directly.
            self._curtime = int(float(ts.time))

            # This block is the current compatibility bridge between MDAnalysis
            # and the legacy MembIT object model.  It is intentionally profiled
            # as its own section because it is a prime candidate for future
            # optimization: the current implementation updates one Atom object
            # at a time, preserving behavior before we attempt bulk updates.
            with self._profileSection('mda_populate_membit_atoms'):
                for number, atom in atoms_by_number.items():
                    x, y, z = atom.position
                    atype = atom.name
                    residue = getattr(atom.residue, 'resname', '')

                    if number in proteinAtoms:
                        self._protein.addProperties(number, atype, residue, float(x), float(y), float(z))

                    if number in CoIAtoms:
                        self._CoI.addProperties(number, atype, residue, float(x), float(y), float(z))

                    elif number in membraneAtoms:
                        self._membrane.addProperties(number, atype, residue, float(x), float(y), float(z))

            with self._profileSection('trajectory_index_match_checks'):
                if not self._CoI.IndexandTrajAtomsMatch():
                    raise IOError('Index file not correct. CoI group atoms in the index do not match the trajectory file')
                if not self._protein.IndexandTrajAtomsMatch():
                    raise IOError('Index file not correct. Protein group atoms in the index do not match the trajectory file')

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
    structurefile = args.structure
    traj_format = args.format
    outputfile = args.o

    simplethickness = args.simplethickness
    thickness = args.thickness
    deformation = args.deformation
    insertion = args.insertion

    distance_criteria = args.distance

    printnatoms = args.printnatoms
    profile_timing = args.profile_timing

    traj = Trajectory(trajfile, indexfile, structurefile, traj_format,
                      distance_criteria, outputfile, thickness, deformation,
                      simplethickness, insertion, printnatoms, profile_timing)

    traj.analyseTrajectory()

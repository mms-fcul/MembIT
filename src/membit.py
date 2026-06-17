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

parser.add_argument('--diagnose-index',
                    help='Print an index/topology diagnostic report and exit before trajectory analysis. '
                         'This is useful for checking whether a MembIT index, structure, and prepared trajectory match, especially for XTC/TRR inputs.',
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
                 simplethickness, insertion, printnatoms, profile_timing=False, diagnose_index=False):
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

        # Index/topology diagnostic support.  These reports are intentionally
        # read-only: they do not change selections or calculated outputs.
        self._diagnoseIndexOnly = diagnose_index
        self._indexDiagnostics = None
        self._indexDiagnosticWarnings = []
        self._recommendedMembraneMarkerAtoms = set(['O31', 'P31', 'O32', 'O33', 'O34'])

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

    def _coerceAtomNumber(self, value):
        """Return an integer atom number from legacy MembIT atom containers.

        Most collection-level methods return atom numbers directly, but some
        leaflet methods return Atom objects.  Diagnostics must normalize both
        representations before comparing index groups against topology metadata.
        """
        try:
            return int(value)
        except Exception:
            pass

        for attr_name in ['getNumber', 'getAtomNumber', 'getNum',
                          'number', 'atom_number', '_number',
                          'atomNumber', 'id']:
            attr = getattr(value, attr_name, None)
            if attr is None:
                continue
            try:
                candidate = attr() if callable(attr) else attr
                return int(candidate)
            except Exception:
                continue

        raise TypeError('Could not determine atom number from object: {0!r}'.format(value))

    def _normalizeAtomNumbers(self, values):
        """Normalize an iterable of atom numbers or Atom objects."""
        return [self._coerceAtomNumber(value) for value in values]

    def _getGroupAtomNumbers(self):
        """Return MembIT index groups as plain lists of integer atom numbers.

        The legacy Protein/Membrane classes own the actual data structures.  This
        helper exposes the group membership in one place so diagnostics can be
        produced without changing calculation logic.
        """
        return {
            'Protein': self._normalizeAtomNumbers(self._protein.getAtomsNumbers()),
            'Center_of_Interest': self._normalizeAtomNumbers(self._CoI.getAtomsNumbers()),
            'Monolayer1': self._normalizeAtomNumbers(self._membrane.getLeafletAtoms('one')),
            'Monolayer2': self._normalizeAtomNumbers(self._membrane.getLeafletAtoms('two')),
        }

    def _topCounts(self, values, limit=12):
        counts = {}
        for value in values:
            counts[value] = counts.get(value, 0) + 1
        return sorted(counts.items(), key=lambda item: (-item[1], str(item[0])))[:limit]

    def _formatTopCounts(self, counts):
        if not counts:
            return 'none'
        return ', '.join('{0}:{1}'.format(name, count) for name, count in counts)

    def _guessIndexNumberingIssue(self, diagnostics):
        """Generate human-readable warnings for common index/topology mistakes."""
        warnings = []
        groups = diagnostics.get('groups', {})
        structure_natoms = diagnostics.get('structure_natoms')
        traj_format = diagnostics.get('trajectory_format')

        protein = groups.get('Protein', {})
        mono1 = groups.get('Monolayer1', {})
        mono2 = groups.get('Monolayer2', {})

        marker_atoms = self._recommendedMembraneMarkerAtoms
        marker_text = ', '.join(sorted(marker_atoms))

        for name in ['Monolayer1', 'Monolayer2']:
            group = groups.get(name, {})
            count = group.get('count', 0)
            marker_count = group.get('recommended_marker_count', 0)
            if count > 0:
                marker_fraction = float(marker_count) / float(count)
                if marker_count == 0:
                    warnings.append(
                        '{0} contains no common phosphate/headgroup marker atoms ({1}). '
                        'For thickness/deformation calculations, monolayer groups should usually contain leaflet marker atoms, not all lipid/tail atoms.'.format(name, marker_text)
                    )
                elif marker_fraction < 0.50:
                    warnings.append(
                        '{0} contains only {1}/{2} common phosphate/headgroup marker atoms ({3:.1f}%). '
                        'This is suspicious for thickness/deformation calculations unless you intentionally use a different membrane marker definition.'.format(
                            name, marker_count, count, 100.0 * marker_fraction)
                    )

        # A common mistake is to mix atom-numbering schemes between the
        # trajectory, structure, and index.  For example, an index made for a
        # reduced Protein+Phos trajectory often looks like Protein = 1..N and
        # membrane marker atoms immediately after N.  That is valid only for the
        # matching reduced structure/trajectory, not for a larger structure.
        protein_max = protein.get('max')
        mono1_min = mono1.get('min')
        mono2_min = mono2.get('min')
        if (traj_format != 'pdb' and structure_natoms and protein_max and
                mono1_min and mono2_min and min(mono1_min, mono2_min) >= protein_max + 1):
            warnings.append(
                'The membrane atom numbers start immediately after the Protein group, while the structure contains {0} atoms. '
                'This pattern often means an index/trajectory atom-numbering mismatch, such as a reduced Protein+membrane-marker index being used with a larger structure. '
                'Prepare the trajectory for membrane analysis first, then generate the structure and GROMACS .ndx against the same processed atom set supplied to MembIT.'.format(structure_natoms)
            )

        return warnings

    def _buildIndexDiagnosticsFromMetadata(self, metadata_by_number, structure_natoms=None):
        """Build a diagnostic summary using structure/topology atom metadata."""
        diagnostics = {
            'trajectory_file': self._trajfile,
            'structure_file': self._structurefile,
            'index_file': self._indexfile,
            'trajectory_format': self._traj_format,
            'structure_natoms': structure_natoms,
            'recommended_marker_atoms': sorted(self._recommendedMembraneMarkerAtoms),
            'groups': {},
        }

        for group_name, numbers in self._getGroupAtomNumbers().items():
            atom_names = []
            residue_names = []
            missing = []
            int_numbers = []

            for number in numbers:
                try:
                    int_numbers.append(int(number))
                except Exception:
                    pass

                metadata = metadata_by_number.get(number)
                if metadata is None:
                    metadata = metadata_by_number.get(str(number))

                if metadata is None:
                    missing.append(number)
                    continue

                # metadata tuple: (atom_number, atom_type, residue_name, mda_index)
                atom_names.append(metadata[1])
                residue_names.append(metadata[2])

            marker_count = sum(1 for name in atom_names if name in self._recommendedMembraneMarkerAtoms)
            diagnostics['groups'][group_name] = {
                'count': len(numbers),
                'min': min(int_numbers) if int_numbers else None,
                'max': max(int_numbers) if int_numbers else None,
                'missing_metadata_count': len(missing),
                'recommended_marker_count': marker_count,
                'top_atom_names': self._topCounts(atom_names),
                'top_residue_names': self._topCounts(residue_names),
            }

        diagnostics['warnings'] = self._guessIndexNumberingIssue(diagnostics)
        return diagnostics

    def _formatIndexDiagnosticReport(self, title='MembIT index/topology diagnostic report'):
        diagnostics = self._indexDiagnostics
        if diagnostics is None:
            return '{0}\nNo topology-aware diagnostic information is available yet. For PDB input, use a matching PDB/index. For XTC/TRR/DCD/NC input, provide -s structure.gro or -s structure.tpr.'.format(title)

        lines = []
        lines.append(title)
        lines.append('=' * len(title))
        lines.append('Trajectory file : {0}'.format(diagnostics.get('trajectory_file')))
        lines.append('Structure file  : {0}'.format(diagnostics.get('structure_file')))
        lines.append('Index file      : {0}'.format(diagnostics.get('index_file')))
        lines.append('Reader format   : {0}'.format(diagnostics.get('trajectory_format')))
        if diagnostics.get('structure_natoms') is not None:
            lines.append('Structure atoms : {0}'.format(diagnostics.get('structure_natoms')))
        lines.append('Marker atoms checked for membrane sanity: {0}'.format(', '.join(diagnostics.get('recommended_marker_atoms', []))))
        lines.append('')
        lines.append('{0:<20s} {1:>8s} {2:>10s} {3:>10s} {4:>12s}'.format('group', 'count', 'min', 'max', 'marker_atoms'))
        lines.append('{0:<20s} {1:>8s} {2:>10s} {3:>10s} {4:>12s}'.format('-' * 20, '-' * 8, '-' * 10, '-' * 10, '-' * 12))

        for group_name in ['Protein', 'Center_of_Interest', 'Monolayer1', 'Monolayer2']:
            group = diagnostics.get('groups', {}).get(group_name, {})
            lines.append('{0:<20s} {1:8d} {2:>10s} {3:>10s} {4:12d}'.format(
                group_name,
                group.get('count', 0),
                str(group.get('min')),
                str(group.get('max')),
                group.get('recommended_marker_count', 0)))
            lines.append('  atom names : {0}'.format(self._formatTopCounts(group.get('top_atom_names', []))))
            lines.append('  residues   : {0}'.format(self._formatTopCounts(group.get('top_residue_names', []))))

        warnings = diagnostics.get('warnings', [])
        if warnings:
            lines.append('')
            lines.append('Potential problems detected:')
            for warning in warnings:
                lines.append('  - {0}'.format(warning))

        lines.append('')
        lines.append('Recommended input workflow:')
        lines.append('  - Prepare the trajectory before running MembIT: handle PBC, center/image the protein or solute, and keep the membrane/protein geometry consistent.')
        lines.append('  - Run MembIT on a treated analysis trajectory, usually containing the protein or center of interest plus leaflet marker atoms, not a raw full-system trajectory straight from mdrun.')
        lines.append('  - Generate the structure and .ndx against the same processed atom set supplied to MembIT; the atom numbering in -f, -s, and -n must match.')
        lines.append('  - For thickness/deformation, Monolayer1 and Monolayer2 should normally contain leaflet marker atoms such as phosphate/headgroup atoms, not complete lipid atom clouds.')
        lines.append('  - If waters, ions, or full lipid atoms are scientifically required, verify the assumptions carefully and use --diagnose-index before long runs.')
        lines.append('  - If your force field uses different marker atom names, verify them with gmx make_ndx/select and generate Monolayer1/Monolayer2 accordingly.')

        return '\n'.join(lines)

    def _raiseWithIndexDiagnostics(self, exc, context):
        message = []
        message.append(str(exc))
        message.append('')
        message.append('MembIT failed while calculating {0}.'.format(context))
        message.append('This can happen when the membrane groups in the index are not suitable for the requested calculation, or when the index atom numbering does not match the trajectory/structure.')
        message.append('')
        message.append(self._formatIndexDiagnosticReport('Index/topology diagnostic at failure'))
        raise IOError('\n'.join(message)) from exc

    def diagnoseIndexAndExit(self):
        """Build and print topology-aware index diagnostics, then exit."""
        if self._traj_format == 'pdb':
            # The native PDB reader can validate exact atoms only while reading a
            # frame.  Keep the message explicit rather than pretending to have
            # topology metadata.
            print(self._formatIndexDiagnosticReport())
            return

        if not self._structurefile:
            raise IOError('--diagnose-index for {0} input requires -s structure.gro or -s structure.tpr'.format(self._traj_format.upper()))

        try:
            import MDAnalysis as mda
        except ImportError as exc:
            raise ImportError('Index diagnostics for {0} trajectories require MDAnalysis.'.format(self._traj_format.upper())) from exc

        universe = mda.Universe(self._structurefile, self._trajfile)
        group_numbers = self._getGroupAtomNumbers()
        required_numbers = set()
        for numbers in group_numbers.values():
            required_numbers.update(numbers)

        metadata_by_number = {}
        natoms = len(universe.atoms)
        for number in sorted(required_numbers, key=lambda value: int(value)):
            atom_index = int(number) - 1
            if atom_index < 0 or atom_index >= natoms:
                raise IOError(
                    'Index atom number {0} is outside the structure atom range 1..{1}. '
                    'Generate the index against the same structure/topology supplied with -s.'.format(number, natoms)
                )
            atom = universe.atoms[atom_index]
            metadata_by_number[number] = (number, atom.name, getattr(atom.residue, 'resname', ''), atom_index)

        self._indexDiagnostics = self._buildIndexDiagnosticsFromMetadata(metadata_by_number, structure_natoms=natoms)
        print(self._formatIndexDiagnosticReport())

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
                                    try:
                                        self._membrane.calcHalfMembraneZ(self._protein,
                                                                         (0, 0, 0, 0,
                                                                          self._insertion[1]),
                                                                         self._box)
                                    except (IOError, OSError) as exc:
                                        self._raiseWithIndexDiagnostics(exc, 'insertion zero/bulk membrane reference')
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
                                try:
                                    self._membrane.calcHalfMembraneZ(self._protein,
                                                                     self._thickness,
                                                                     self._box)
                                except (IOError, OSError) as exc:
                                    self._raiseWithIndexDiagnostics(exc, 'thickness/deformation bulk membrane reference')

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

        Performance note
        ----------------
        The first profiled XTC implementation spent most of its extra runtime
        in ``mda_populate_membit_atoms``.  The slow path was not the membrane
        calculation itself; it was the compatibility bridge that copied
        coordinates from MDAnalysis Atom objects into the legacy MembIT Atom
        containers every frame.

        This optimized reader therefore performs all static work once before the
        trajectory loop:

            * validate index atom numbers against the structure atom count;
            * convert 1-based GROMACS/MembIT numbers to 0-based MDAnalysis
              indices;
            * cache atom names and residue names;
            * build one update-record list per MembIT collection.

        Inside the frame loop, the code only reads the frame coordinate array
        and updates MembIT containers from cached records.  This keeps the old
        object model and output behavior intact, but avoids repeated MDAnalysis
        Atom property access, repeated integer conversion, and repeated group
        membership checks for every frame.
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
        # through universe.trajectory.
        with self._profileSection('mda_universe_init'):
            universe = mda.Universe(self._structurefile, self._trajfile)

        # Convert frequently used index containers to sets.  The legacy PDB
        # reader still uses the historical containers directly, but the XTC path
        # performs many membership checks while preparing update records.  Sets
        # make that one-time preparation explicit and cheap.
        proteinAtoms = set(self._protein.getAtomsNumbers())
        CoIAtoms = set(self._CoI.getAtomsNumbers())
        membraneAtoms = set(self._membrane.getAtomsNumbers())
        required_numbers = proteinAtoms | CoIAtoms | membraneAtoms

        # ``update_records`` are the main optimization introduced here.
        #
        # Each record is a tuple:
        #     (atom_number, atom_type, residue_name, mda_index)
        #
        # All fields except the coordinates are static over a normal trajectory,
        # so they are read once here instead of being pulled from MDAnalysis Atom
        # objects for every frame.  The separate lists preserve the previous
        # update logic:
        #
        #     if number in Protein: update Protein
        #     if number in CoI:     update CoI
        #     elif number in Membrane: update Membrane
        #
        # In particular, this preserves the old ``CoI`` versus ``Membrane``
        # precedence for any atom number that might appear in both groups.
        with self._profileSection('mda_atom_cache'):
            natoms = len(universe.atoms)
            metadata_by_number = {}

            for number in sorted(required_numbers, key=lambda value: int(value)):
                atom_index = int(number) - 1
                if atom_index < 0 or atom_index >= natoms:
                    raise IOError(
                        'Index atom number {0} is outside the structure atom range 1..{1}. '
                        'For MDAnalysis/XTC input, MembIT index files must use 1-based '
                        'structure atom numbers, matching GROMACS .ndx convention.'.format(number, natoms)
                    )

                atom = universe.atoms[atom_index]
                atype = atom.name
                residue = getattr(atom.residue, 'resname', '')

                metadata_by_number[number] = (number, atype, residue, atom_index)

            protein_update_records = [metadata_by_number[number]
                                      for number in sorted(proteinAtoms, key=lambda value: int(value))]
            coi_update_records = [metadata_by_number[number]
                                  for number in sorted(CoIAtoms, key=lambda value: int(value))]
            membrane_update_records = [metadata_by_number[number]
                                       for number in sorted(membraneAtoms - CoIAtoms,
                                                            key=lambda value: int(value))]

            self._indexDiagnostics = self._buildIndexDiagnosticsFromMetadata(metadata_by_number, structure_natoms=natoms)

        for ts in universe.trajectory:
            if ts.dimensions is None or len(ts.dimensions) < 3:
                raise IOError('Trajectory frame has no unit-cell dimensions; MembIT requires box vectors.')

            # MDAnalysis reports GROMACS-like coordinate files in Angstrom, which
            # is the unit expected by the existing MembIT calculations.
            self._box = float(ts.dimensions[0]), float(ts.dimensions[1]), float(ts.dimensions[2])

            # Existing XVG output stores time as an integer.  Preserve that
            # behavior so PDB and XTC paths can be compared directly.
            self._curtime = int(float(ts.time))

            # Fast coordinate bridge:
            #
            # ``ts.positions`` is the coordinate array for the current frame.
            # Using it directly avoids the costly per-atom ``atom.position``
            # property access that the baseline XTC reader used.  The records
            # built above already contain the MDAnalysis array row and the static
            # metadata needed by MembIT's legacy ``addProperties`` methods.
            with self._profileSection('mda_populate_membit_atoms'):
                positions = ts.positions

                for number, atype, residue, atom_index in protein_update_records:
                    x, y, z = positions[atom_index]
                    self._protein.addProperties(number, atype, residue, float(x), float(y), float(z))

                for number, atype, residue, atom_index in coi_update_records:
                    x, y, z = positions[atom_index]
                    self._CoI.addProperties(number, atype, residue, float(x), float(y), float(z))

                for number, atype, residue, atom_index in membrane_update_records:
                    x, y, z = positions[atom_index]
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
    diagnose_index = args.diagnose_index

    traj = Trajectory(trajfile, indexfile, structurefile, traj_format,
                      distance_criteria, outputfile, thickness, deformation,
                      simplethickness, insertion, printnatoms, profile_timing, diagnose_index)

    if diagnose_index:
        traj.diagnoseIndexAndExit()
    else:
        traj.analyseTrajectory()

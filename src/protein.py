from atomcollections import AtomCollections

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

        if parameters[0] == 'closest':
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

        if nargs_insertion >= 4:
            max_window = float(parameters[3])
            min_window = float(parameters[2])
        if nargs_insertion >= 2:
            #closest_atom = membrane.getClosestAtom()
            #furthest_atom = membrane.getFurthestAtom()

            #closest_distance = closest_atom.getDistance2CoI() ** 0.5
            #furthest_dist2410ance = furthest_atom.getDistance2CoI()  ** 0.5

            if nargs_insertion == 3:
                min_window = float(parameters[2])
                max_window = box_size
            elif nargs_insertion == 2:
                min_window = 0
                max_window = box_size

            window_step = float(parameters[1])
            window_size = float(parameters[0])

            # get reference membrane atoms
            membrane_atoms_sorted = membrane.getReferenceML2CoIDistanceSorted()
            slices = membrane.getSlices(window_size, window_step,
                                        min_window, max_window,
                                        membrane_atoms_sorted)
            atom_sorted_bydist = [i[0] for i in membrane_atoms_sorted]

            output = ''
            windows = []
            for window_half, atoms, distances in slices:
                natoms = len(atoms)


                if natoms == 0 and 'noNaN' in parameters[-1]:
                    nclosest_atoms = parameters[-1].strip('noNaN')
                    if len(nclosest_atoms) == 0:
                        nclosest_atoms = 1
                    else:
                        nclosest_atoms = int(nclosest_atoms)
                    atoms = atom_sorted_bydist[:nclosest_atoms]
                    natoms = len(atoms)
                elif 'min' in parameters[-1]:
                    nclosest_atoms = parameters[-1].strip('min')
                    if len(nclosest_atoms) == 0:
                        nclosest_atoms = 1
                    else:
                        nclosest_atoms = int(nclosest_atoms)
                    if natoms < nclosest_atoms:
                        nclosest_atoms = int(nclosest_atoms)
                        atoms = atom_sorted_bydist[:nclosest_atoms]
                        natoms = len(atoms)

                if natoms > 0:
                    total = 0.0
                    for atom in atoms:
                        total += atom.get3DPosition()[2]

                    average_z = total / natoms
                    insertion_value = self.distToZ(average_z,
                                                   mean_other_memb,
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


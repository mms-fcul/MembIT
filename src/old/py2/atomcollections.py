from atom import Atom

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
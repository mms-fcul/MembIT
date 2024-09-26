from membit_module import PYX2Ddist, PYX3Ddist

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

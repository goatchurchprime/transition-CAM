from basicgeo import P3
import numpy

class NTriangleBarMesh:
    def __init__(self, trianglebarmesh):
        self.nnodes = numpy.zeros((len(trianglebarmesh.nodes), 3), numpy.float)
        self.nbars = numpy.zeros((len(trianglebarmesh.bars), 4), numpy.int)
        for i, tnode in enumerate(trianglebarmesh.nodes):
            self.nnodes[i,] = tnode.p
        for j, tbar in enumerate(trianglebarmesh.bars):
            self.nbars[j,0] = tbar.nodeback.i
            self.nbars[j,1] = tbar.nodefore.i
            self.nbars[j,2] = tbar.barforeright.i if tbar.barforeright else -1
            self.nbars[j,3] = tbar.barbackleft.i if tbar.barbackleft else -1

    def GetNodePoint(self, i):
        return P3(*self.nnodes[i,])
    def GetBarPoints(self, j):
        return (self.GetNodePoint(self.nbars[j,0]), self.GetNodePoint(self.nbars[j,1]))
    def GetTriPoints(self, j):
        jfr = self.nbars[j,2]
        k = self.nbars[jfr,1] if self.nbars[j,1] == self.nbars[jfr,0] else self.nbars[jfr,0]
        assert k != self.nbars[j,0] and k != self.nbars[j,1]
        return (self.GetNodePoint(self.nbars[j,0]), self.GetNodePoint(self.nbars[j,1]), self.GetNodePoint(k))


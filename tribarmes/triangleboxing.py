import math
from basicgeo import P3, Partition1
import basicgeo

# This is where 99% of the time is spent in two functions that are intended to be GPU-friendly 
# functions DistP() embedded in MakePointZoneRF() and CutbarRF() embedded in CutbarRF() of 
# which hundreds of calls should act in parallel per iteration.  
# Both of these functions are similar to 3D graphics depth buffer plotting (finding the closest point, 
# and finding first point of intersection along a ray), with the main difference being that the 
# parallel calls do not correspond to a regular grid of pixels intersected by interpolatable triangles, 
# which is the underlying assumption within GLSL.

class TriangleBox:
    def __init__(self):
        self.pointis = [ ]
        self.edgeis = [ ]
        self.triangleis = [ ]
        

def CYSupdate(cys, x0, y0, x1, y1, xlo, xhi):
    assert x0 <= x1 and xlo < xhi
    if xlo <= x0 <= xhi:
        cys.append(y0)
    if xlo <= x1 <= xhi:
        cys.append(y1)
    if x0 <= xlo < x1:
        lam = (xlo - x0) / (x1 - x0)
        cys.append(y0 * (1 - lam) + y1 * lam)
    if x0 <= xhi < x1:
        lam = (xhi - x0) / (x1 - x0)
        cys.append(y0 * (1 - lam) + y1 * lam)


class TriangleBoxing:
    def __init__(self, tbarmesh, xlo, xhi, ylo, yhi, boxwidth):
        self.tbarmesh = tbarmesh
        self.xpart = Partition1(xlo, xhi, int(math.ceil((xhi - xlo)/boxwidth)))
        self.ypart = Partition1(ylo, yhi, int(math.ceil((yhi - ylo)/boxwidth)))
        self.boxes = [ [ TriangleBox()  for iy in range(self.ypart.nparts) ]  for ix in range(self.xpart.nparts) ]

    def GetNodePoint(self, i):
        return self.tbarmesh.GetNodePoint(i)
    def GetBarPoints(self, i):
        return self.tbarmesh.GetBarPoints(i)
    def GetTriPoints(self, i):
        return self.tbarmesh.GetTriPoints(i)

    def AddPoint(self, x, y, i):
        ix = self.xpart.GetPart(x)
        assert ix < len(self.boxes)
        iy = self.ypart.GetPart(y)
        tb = self.boxes[ix][iy]
        tb.pointis.append(i)
        
    def SortPointis(self):
        for ix in range(self.xpart.nparts):
            for iy in range(self.ypart.nparts):
                self.boxes[ix][iy].pointis.sort(key=lambda i: self.tbarmesh.GetNodePoint(i).z)
                
    def SlicePointisZ(self, pointis, zlo, zhi):
        #return pointis   # to disable
        assert zlo <= zhi
        if not pointis:
            return ()
        if self.GetNodePoint(pointis[-1]).z < zlo:
            return ()
        if self.GetNodePoint(pointis[0]).z > zhi:
            return ()

        j0, j1 = -1, len(pointis)
        if j1 <= 3:
            return pointis
        while True:
            jo = (j1 - j0)//2
            if jo == 0:
                return (pointis[j]  for j in range(j0+1,j1))
            jmid = j0 + jo
            zmid = self.GetNodePoint(pointis[jmid]).z
            if zmid < zlo:
                j0 = jmid
            elif zmid > zhi:
                j1 = jmid
            else:
                break
        assert j0 == -1 or self.GetNodePoint(pointis[j0]).z < zlo
        assert j1 == len(pointis) or self.GetNodePoint(pointis[j1]).z > zhi
        assert zlo <= self.GetNodePoint(pointis[jmid]).z <= zhi
        j0u = jmid
        j1l = jmid
        del jmid
        while True:
            j1o = (j1 - j1l)//2
            if j1o == 0:
                break
            j1mid = j1l + j1o
            z1mid = self.GetNodePoint(pointis[j1mid]).z
            if z1mid > zhi:
                j1 = j1mid
            else:
                j1l = j1mid
        while True:
            j0o = (j0u - j0)//2
            if j0o == 0:
                break
            j0mid = j0 + j0o
            z0mid = self.GetNodePoint(pointis[j0mid]).z
            if z0mid < zlo:
                j0 = j0mid
            else:
                j0u = j0mid
        assert j0 == -1 or self.GetNodePoint(pointis[j0]).z < zlo
        assert zlo <= self.GetNodePoint(pointis[j0+1]).z <= zhi
        assert j1 == len(pointis) or self.GetNodePoint(pointis[j1]).z > zhi
        assert zlo <= self.GetNodePoint(pointis[j1-1]).z <= zhi
        return (pointis[j]  for j in range(j0+1,j1))
        
        
    def AddEdge(self, x0, y0, x1, y1, i):
        assert x0 <= x1
        ixl, ixh = self.xpart.GetPartRange(x0, x1)
        for ix in range(ixl, ixh + 1):
            cys = [ ]
            CYSupdate(cys, x0, y0, x1, y1, self.xpart.vs[ix], self.xpart.vs[ix+1])
            if cys:
                iyl, iyh = self.ypart.GetPartRange(min(cys), max(cys))
                for iy in range(iyl, iyh + 1):
                    tb = self.boxes[ix][iy]
                    tb.edgeis.append(i)

    def AddEdgeR(self, x0, y0, x1, y1, i):
        if x0 <= x1:
            self.AddEdge(x0, y0, x1, y1, i)
        else:
            self.AddEdge(x1, y1, x0, y0, i)

    def AddTriangle(self, x0, y0, x1, y1, x2, y2, i):
        assert x0 <= x1 <= x2
        ixl, ixh = self.xpart.GetPartRange(x0, x2)
        for ix in range(ixl, ixh + 1):
            cys = [ ]
            CYSupdate(cys, x0, y0, x1, y1, self.xpart.vs[ix], self.xpart.vs[ix+1])
            CYSupdate(cys, x0, y0, x2, y2, self.xpart.vs[ix], self.xpart.vs[ix+1])
            CYSupdate(cys, x1, y1, x2, y2, self.xpart.vs[ix], self.xpart.vs[ix+1])
            if cys:
                iyl, iyh = self.ypart.GetPartRange(min(cys), max(cys))
                for iy in range(iyl, iyh + 1):
                    tb = self.boxes[ix][iy]
                    tb.triangleis.append(i)

    def AddTriangleRR(self, x0, y0, x1, y1, x2, y2, i):
        assert x0 <= x1 and x0 <= x2
        if x1 < x2:
            self.AddTriangle(x0, y0, x1, y1, x2, y2, i)
        else:
            self.AddTriangle(x0, y0, x2, y2, x1, y1, i)
            
    def AddTriangleR(self, x0, y0, x1, y1, x2, y2, i):
        if x0 <= x1 and x0 <= x2:
            self.AddTriangleRR(x0, y0, x1, y1, x2, y2, i)
        elif x1 <= x0 and x1 <= x2:
            self.AddTriangleRR(x1, y1, x0, y0, x2, y2, i)
        else:
            self.AddTriangleRR(x2, y2, x1, y1, x0, y0, i)
        
    def GetTriangleBox(self, ixy):
        return self.boxes[ixy[0]][ixy[1]] 
        
    def CloseBoxeGenerator(self, xlo, xhi, ylo, yhi, r):
        ixlo, ixhi = self.xpart.GetPartRange(xlo - r, xhi + r)
        #basicgeo.Dplotrect(xlo, xhi, ylo, yhi, materialnumber=1)
        for ix in range(ixlo, ixhi + 1):
            assert ix + 1 < len(self.xpart.vs), (ixlo, ixhi, self.xpart.nparts, len(self.xpart.vs))
            if self.xpart.vs[ix] <= xhi and self.xpart.vs[ix + 1] >= xlo:
                ry = r
            else:
                if self.xpart.vs[ix] <= xlo:
                    drx = xlo - self.xpart.vs[ix + 1]
                else:
                    drx = self.xpart.vs[ix] - xhi
                assert drx >= 0
                ry = math.sqrt(max(r*r - drx*drx, 0))
            iylo, iyhi = self.ypart.GetPartRange(ylo - ry, yhi + ry)
            for iy in range(iylo, iyhi + 1):
                #basicgeo.Dplotrect(self.xpart.vs[ix], self.xpart.vs[ix+1], self.ypart.vs[iy], self.ypart.vs[iy+1], materialnumber=2)
                yield (ix, iy)
            
        
    
def MakeTriangleBoxing(tbarmesh, boxwidth=-1):
    if boxwidth == -1:
        boxwidth = min(10.1, (tbarmesh.xhi - tbarmesh.xlo)*0.2, (tbarmesh.yhi - tbarmesh.ylo)*0.2)
        boxwidth *= 0.05
    d = boxwidth + 0.5
    print("make triangle boxing at", boxwidth)
    tboxing = TriangleBoxing(tbarmesh, tbarmesh.xlo - d, tbarmesh.xhi + d, tbarmesh.ylo - d, tbarmesh.yhi + d, boxwidth)
    
    for node in tbarmesh.nodes:
        tboxing.AddPoint(node.p.x, node.p.y, node.i)
    tboxing.SortPointis()
        
    for iedge, bar in enumerate(tbarmesh.bars):
        tboxing.AddEdgeR(bar.nodeback.p.x, bar.nodeback.p.y, bar.nodefore.p.x, bar.nodefore.p.y, iedge)
        
        if bar.barforeright:
            node2 = bar.barforeright.GetNodeFore(bar.nodefore == bar.barforeright.nodeback)
            if node2.i > bar.nodeback.i:
                tboxing.AddTriangleR(bar.nodeback.p.x, bar.nodeback.p.y, bar.nodefore.p.x, bar.nodefore.p.y, node2.p.x, node2.p.y, iedge)

    return tboxing
    

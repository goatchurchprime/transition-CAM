from basicgeo import P3, AlongAcc

PZ_WITHIN_R = 1 # values for PointZone.izone
PZ_BEYOND_R = 2

class PointZone:
    def __init__(self, izone, r, v):
        self.izone = izone
        self.r = r
        self.v = v

# This might best require a normal at each node, which defines the plane of projection 
class Node:
    def __init__(self, p, i):
        self.p = p
        self.i = i  # index
        self.pointzone = None
        
    # will in future be flattened by the local plane for the node
    def cperpbardot(self, bar, v):
        assert self == bar.nodeback or self == bar.nodefore
        vb = (bar.nodefore.p - bar.nodeback.p)
        vbs = (self == bar.nodeback and +1 or -1)
        cperpvb = P3(vb.y*vbs, -vb.x*vbs, 0)
        return P3.Dot(cperpvb, v)
        
    # will in future be flattened by the local plane for the node
    def cperpbardotN(self, bar, v):
        assert self == bar.nodeback or self == bar.nodefore
        vb = (bar.nodefore.p - bar.nodeback.p)
        vbs = (self == bar.nodeback and +1 or -1)
        cperpvb = P3(vb.y*vbs, -vb.x*vbs, 0)
        return P3.Dot(cperpvb, v)/vb.Len()


# an object that identifies with a cell and so is referenced by all the bars surrounding on the approrpiate side
class CellMark:
    def __init__(self, cellcolour):
        # 'colour' here refers to the four-colour theorem partitioning the set of cells into subsets that have no shared edges
        self.cellcolour = cellcolour
        self.cellhitcounter = -1

class Bar:
    def __init__(self, nodeback, nodefore):
        self.nodeback = nodeback
        self.nodefore = nodefore
        self.barforeright = None
        self.barbackleft = None
        self.bbardeleted = False
        assert nodefore.i > nodeback.i
        self.nodemid = None   # used to specify a contour cut or a voronoi polygon boundary cut
        self.midcontournumber = -1  # will be set to -2 to denote connecting to out of tolerance/unchecked tolerance trailing contour segment
        self.cellmarkright = None
        self.cellmarkleft = None
        self.barvecN = P3.ZNorm(self.nodefore.p - self.nodeback.p)  # used to detect colinearity as it's preserved on splitting
        
    def SetForeRightBL(self, bforeright, bar):
        if bforeright:
            self.barforeright = bar
        else:
            self.barbackleft = bar

    def GetForeRightBL(self, bforeright):
        return self.barforeright  if bforeright  else self.barbackleft
    def GetNodeFore(self, bfore):
        return self.nodefore  if bfore  else self.nodeback
    def GetCellMarkRightL(self, bright):
        return self.cellmarkright  if bright  else self.cellmarkleft
            
    def GetBarForeLeft(self):
        assert not self.bbardeleted
        if not self.barbackleft:  
            return None
        barleft, barleftnodeback = self.barbackleft, self.nodeback
        Dcounter = 0
        while True:
            bbarleftnodefore = (barleft.nodeback == barleftnodeback)
            barleftnodeback = barleft.GetNodeFore(bbarleftnodefore)
            if barleftnodeback == self.nodefore:
                break
            barleft = barleft.GetForeRightBL(bbarleftnodefore)
            assert barleft
            assert not barleft.bbardeleted
            Dcounter += 1
            assert Dcounter < 1000, "Infinite loop in GetBarForeLeft"
        assert barleft.GetForeRightBL(barleft.nodefore == self.nodefore) == self, (barleft == self, Dcounter, self.nodeback.p, self.nodefore.p)
        return barleft
        
    def GetBarBackRight(self):
        assert not self.bbardeleted
        if not self.barforeright:  
            return None
        barright, barrightnodeback = self.barforeright, self.nodefore
        Dcounter = 0
        while True:
            bbarrightnodefore = (barright.nodeback == barrightnodeback)
            barrightnodeback = barright.GetNodeFore(bbarrightnodefore)
            if barrightnodeback == self.nodeback:
                break
            barright = barright.GetForeRightBL(bbarrightnodefore)
            assert barright, barrightnodeback.p
            assert not barright.bbardeleted
            Dcounter += 1
            assert Dcounter < 1000, "Infinite loop in GetBarBackRight"
        assert barright.GetForeRightBL(barright.nodefore == self.nodeback) == self
        return barright

    # to be used sparingly
    def GetForeLeftBR(self, bforeleft):
        return self.GetBarForeLeft()  if bforeleft  else self.GetBarBackRight()

def setcellmark(bar, node, newcellmark, Dcellmark):
    assert not bar.bbardeleted
    lbar, lnode = bar, node
    Dcount = 0
    while True:
        lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
        lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
        assert Dcellmark == None or Dcellmark == lbar.GetCellMarkRightL(lnode == lbar.nodeback)
        if lnode == lbar.nodeback:
            lbar.cellmarkright = newcellmark
        else:
            lbar.cellmarkleft = newcellmark
        if lbar == bar:
            break
        Dcount += 1
        assert Dcount < 1000


#     /       
#    /         
# --n  - v - > 
#    \       
#     b       


#########################
# Main object
#########################
class BarMesh:
    def __init__(self):
        self.nodes = [ ]
        self.bars = [ ]
        #self.xlo, self.xhi, self.ylo, self.yhi  # set in NewNode()
        self.startcontournumber = 0
        self.endcontournumber = 0
        self.cellhitcounterend = 0
        self.trianglebarmeshtype = False
        
    def NewNode(self, p):
        if self.nodes:
            if p.x < self.xlo:  self.xlo = p.x
            if p.x > self.xhi:  self.xhi = p.x
            if p.y < self.ylo:  self.ylo = p.y
            if p.y > self.yhi:  self.yhi = p.y
            if p.z < self.zlo:  self.zlo = p.z
            if p.z > self.zhi:  self.zhi = p.z
        else:
            self.xlo = self.xhi = p.x
            self.ylo = self.yhi = p.y
            self.zlo = self.zhi = p.z
        self.nodes.append(Node(p, len(self.nodes)))
        return self.nodes[-1]

        
    # +-y->+-y->+
    # ^    ^    ^
    # |    |    |
    # x    x    x
    # |    |    |
    # +-y->+-y->+
    def BuildRectBarMesh(self, xpart, ypart, z):
        nxs = xpart.nparts + 1
        nodes = self.nodes
        xbars = [ ]  # multiple of nxs
        ybars = [ ]  # multiple of nxs-1
        self.bars = [ ]
        for y in ypart.vs:
            bfirstrow = (len(nodes) == 0)
            for i in range(nxs):
                nnode = self.NewNode(P3(xpart.vs[i], y, z))
                assert nnode == nodes[-1]
                if not bfirstrow:
                    xbars.append(Bar(nodes[-nxs - 1], nodes[-1]))
                    self.bars.append(xbars[-1])
                    if i != 0:
                        xbars[-1].barbackleft = ybars[1-nxs]
                        ybars[1-nxs].barbackleft = xbars[-2]
                if i != 0:
                    ybars.append(Bar(nodes[-2], nodes[-1]))
                    self.bars.append(ybars[-1])
                    if not bfirstrow:
                        ybars[-1].barforeright = xbars[-1]
                        xbars[-2].barforeright = ybars[-1]
        assert len(self.bars) == len(xbars)+len(ybars)
        
        
            
    #   |       |             #   |       |
    # --+------>+--  becomes  # --+-->o<--+--
    #   |       |             #   |       |
    def InsertNodeIntoBarF(self, bar, newnode, bcolinear):
        assert newnode.p != bar.nodeback.p and newnode.p != bar.nodefore.p, newnode.p
        assert newnode in self.nodes
        
        if __debug__:
            if bcolinear:
                Dv = bar.nodefore.p - bar.nodeback.p
                Dv0 = newnode.p - bar.nodeback.p
                Dlam = P3.Dot(Dv0, Dv) / Dv.Lensq()
                assert 0 < Dlam < 1, Dlam
                Dvd = Dv0 - Dv * Dlam
                assert abs(P3.Dot(Dvd, Dv)) < 0.001
                assert Dvd.Len() < 0.01, (Dv, Dv0)
        
        barforeleft = bar.GetBarForeLeft()
        barbackright = bar.GetBarBackRight()

        barback = Bar(bar.nodeback, newnode)
        barfore = Bar(bar.nodefore, newnode)
        
        if bcolinear:  # copy over so we can identify colinearity from creation and not attempt to split along it
            assert (bar.barvecN - barback.barvecN).Len() < 1e-5, (bar.nodeback.p, newnode.p, bar.nodefore.p)
            assert (bar.barvecN + barfore.barvecN).Len() < 1e-5, (bar.nodeback.p, newnode.p, bar.nodefore.p) 
            barback.barvecN = bar.barvecN
            barfore.barvecN = -bar.barvecN
        
        barback.cellmarkleft, barback.cellmarkright = bar.cellmarkleft, bar.cellmarkright
        barfore.cellmarkright, barfore.cellmarkleft = bar.cellmarkleft, bar.cellmarkright

        if barbackright:
            barback.barforeright = barfore
            barfore.barbackleft = bar.barforeright
            barbackright.SetForeRightBL(barbackright.nodefore == bar.nodeback, barback)

        if barforeleft:
            barback.barbackleft = bar.barbackleft
            barfore.barforeright = barback
            barforeleft.SetForeRightBL(barforeleft.nodefore == bar.nodefore, barfore)

        assert not barforeleft or barback.GetBarForeLeft() == barfore, barback.GetBarForeLeft()
        assert not barbackright or barfore.GetBarForeLeft() == barback, barfore.GetBarForeLeft()

        self.bars.append(barback)
        self.bars.append(barfore)
        
        bar.barforeright = barfore
        bar.barbackleft = barback 
        bar.bbardeleted = True
        
        return newnode
        
    def DTestColinearityF(self, node1, bar1, node2, bar2):
        bar1a = bar1.GetForeRightBL(bar1.nodefore == node1)
        lbar, lnode = bar1a, bar1a.GetNodeFore(bar1a.nodeback == node1)
        if lnode == node2:
            return True
        Dcount = 0
        while lnode != node2:
            assert lnode != node1
            lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
            lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
            if lbar.barvecN != bar1a.barvecN and lbar.barvecN != -bar1a.barvecN:
                return False
            Dcount += 1
            assert Dcount < 1000
        return True


    #    |       b2   
    #    |       |   
    #    V       V   
    # --n1 - - > n2--  
    #    ^       ^    
    #    |       |    
    #   b1       |   
    def MakeBarBetweenNodesF(self, node1, bar1, node2, bar2):
        assert not bar1.bbardeleted and not bar2.bbardeleted
        assert node1.i < node2.i
        assert not self.DTestColinearityF(node1, bar1, node2, bar2), "applying MakeBarBetweenNodesF colinear with side"
        assert not self.DTestColinearityF(node2, bar2, node1, bar1), "applying MakeBarBetweenNodesF colinear with side"

        bar1a = bar1.GetForeRightBL(bar1.nodefore == node1)
        bar2a = bar2.GetForeRightBL(bar2.nodefore == node2)

        # check is within the angles
        assert node1.cperpbardot(bar1, node2.p - node1.p) < 0.0, node1.cperpbardotN(bar1, node2.p - node1.p)
        assert node1.cperpbardot(bar1a, node2.p - node1.p) > 0.0, ("barnotinangle", node1.cperpbardotN(bar1a, node2.p - node1.p), node2.p, node1.p, (bar1a.nodeback.p, bar1a.nodefore.p))
        assert node2.cperpbardot(bar2, node2.p - node1.p) > 0.0, (node1.p, node2.p)
        assert node2.cperpbardot(bar2a, node2.p - node1.p) < 0.0, ("barnotinangle", node2.cperpbardotN(bar2a, node2.p - node1.p))

        # check we can track around the polygon to the other nodes
        if __debug__:
            lbar, lnode = bar1a, bar1a.GetNodeFore(bar1a.nodeback == node1)
            Dcount = 0
            while lnode != node2:
                assert lnode != node1
                lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
                lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
                Dcount += 1
                assert Dcount < 1000
                
            lbar, lnode = bar2a, bar2a.GetNodeFore(bar2a.nodeback == node2)
            Dcount = 0
            while lnode != node1:
                assert lnode != node2
                lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
                lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
                Dcount += 1
                assert Dcount < 1000
            
        # actually make the new bar and link it in
        newbar = Bar(node1, node2)
        newbar.SetForeRightBL(False, bar1a)
        newbar.SetForeRightBL(True, bar2a)
        bar1.SetForeRightBL(bar1.nodefore == node1, newbar)
        bar2.SetForeRightBL(bar2.nodefore == node2, newbar)
        self.bars.append(newbar)
        return newbar
       
    def PlotCellLinks(self, sendactivity):
        # these are crossing edges across the face to illustrate the foreright and backleft links
        barconnsfr = [ ]
        barconnsbl = [ ]
        barconns = [ ]
        for bar in self.bars:
            if bar.bbardeleted:
                continue
            barconns.append((bar.nodeback.p, bar.nodefore.p))
            if bar.barforeright:
                barconnsfr.append((AlongAcc(0.3, bar.nodeback.p, bar.nodefore.p), AlongAcc(0.6, bar.barforeright.nodeback.p, bar.barforeright.nodefore.p)))
            if bar.barbackleft:
                barconnsbl.append((AlongAcc(0.3, bar.nodeback.p, bar.nodefore.p), AlongAcc(0.65, bar.barbackleft.nodeback.p, bar.barbackleft.nodefore.p)))
        sendactivity("contours", contours=barconnsfr, materialnumber=1)
        sendactivity("contours", contours=barconnsbl, materialnumber=2)
        sendactivity("contours", contours=barconns, materialnumber=0)


def PolyCutBar(ktbar, ktnode1, vc, vch, btop):
    assert not ktbar.bbardeleted, "polycutbar hits deletion"
    Dktbarstart = ktbar
    ktnoded1 = P3.Dot(vc, ktnode1.p)
    
    assert btop == (P3.Dot(vc, ktbar.nodemid.p) < vch), ("fg", btop, P3.Dot(vc, ktbar.nodemid.p), vch)
    #assert ktnoded0 < vch  # not necessarily true if bar points back with mid on end
    
    Dcount = 0
    while True:
        ktnode0 = ktnode1
        ktnode1 = ktbar.GetNodeFore(ktbar.nodeback == ktnode0)
        ktnoded0 = ktnoded1
        ktnoded1 = P3.Dot(vc, ktnode1.p)
        if btop == (ktnoded1 > vch):
            break
        ktbar = ktbar.GetForeRightBL(ktbar.nodefore == ktnode1)
        assert not ktbar.bbardeleted
        assert Dktbarstart != ktbar
        Dcount += 1
        assert Dcount < 1000
    ktbarlam = (vch - ktnoded0) / (ktnoded1 - ktnoded0)      
    return ktnode0, ktbar, ktbarlam


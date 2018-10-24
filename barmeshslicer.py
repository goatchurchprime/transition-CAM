import math, time
from basicgeo import P3, P2, AlongAcc, Along
import barmesh
from barmesh import PolyCutBar

# the main batch calculating functions are within MakePointZoneRFS and CutbarRFS
# these are is structured to make it possible to be done in parallel
    
def getcrosscolours(bar, node, Dcellmark):
    assert not bar.bbardeleted, "getcrosscolours hits deletion"
    lbar, lnode = bar, node
    Dcount = 0
    crosscolourset = set()
    while True:
        lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
        lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
        assert lbar, Dcount
        assert Dcellmark == None or Dcellmark == lbar.GetCellMarkRightL(lnode == lbar.nodeback)
        cellmarkcross = lbar.GetCellMarkRightL(lnode == lbar.nodefore)
        if cellmarkcross:
            crosscolourset.add(cellmarkcross.cellcolour)
        if lbar == bar:
            break
        Dcount += 1
        assert Dcount < 1000
    return crosscolourset
    


def GetCellMarkBar(bar, cizone):  # derived from BarMeshCell.GetCellMark
    node = bar.GetNodeFore(bar.nodeback.pointzone.izone == cizone)
    assert node == bar.GetNodeFore(bar.nodeback.pointzone.izone == cizone)
    assert bar.GetForeRightBL(bar.nodeback == node) != None   # not RefersToOuterCell
    return bar.GetCellMarkRightL(node == bar.nodeback)
    
def IsCutBar(bar, cizone):
    assert not bar.bbardeleted
    return (bar.nodeback.pointzone.izone == cizone) != (bar.nodefore.pointzone.izone == cizone)
    
def RefersToOuterCellBar(bar, cizone):
    node = bar.GetNodeFore(bar.nodeback.pointzone.izone == cizone)
    return bar.GetForeRightBL(bar.nodeback == node) == None

def setmincellcolourBar(bar, node):  # ick!
    crosscolourset = getcrosscolours(bar, node, None)
    newcellcolour = 0
    while newcellcolour in crosscolourset:
        newcellcolour += 1
    newcellmark = barmesh.CellMark(newcellcolour)
    barmesh.setcellmark(bar, node, newcellmark, None)
    return newcellmark


# actually should be called a HalfCell as it only goes round the part on the outside of the contour defined by the pointzone
class BarMeshCell:  # not used outside of this module
    def __init__(self, bar, cizone):
        assert not bar.bbardeleted
        assert IsCutBar(bar, cizone)
        self.bar = bar
        self.cizone = cizone
        self.node = bar.GetNodeFore(bar.nodeback.pointzone.izone == cizone)  # this is the node before entering the izone
        assert self.node.pointzone.izone != cizone
        self.cbar = None
        self.cnode = None
        self.vc = None
        #self.leadsplitbartop, self.splitnodetop, self.bsplitnnodetopnew
        #self.leadsplitbarbot, self.splitnodebot, self.bsplitnnodebotnew

    def RefersToOuterCell(self):
        return self.bar.GetForeRightBL(self.bar.nodeback == self.node) == None

    def GetCellMark(self):
        assert self.node == self.bar.GetNodeFore(self.bar.nodeback.pointzone.izone == self.cizone)
        return self.bar.GetCellMarkRightL(self.node == self.bar.nodeback)
        
    def MakeCutSide(self):  # recalculates cbar, cnode in case it's been changed by parts subdividing
        assert self.node == self.bar.GetNodeFore(self.bar.nodeback.pointzone.izone == self.cizone)
        assert self.cbar == None
        cbar, cnode = self.bar, self.node
        cellmark = self.GetCellMark()
        Dcount = 0
        while True:
            cnode = cbar.GetNodeFore(cbar.nodeback == cnode)
            if not (cnode.pointzone.izone == self.cizone):
                break
            cbar = cbar.GetForeRightBL(cbar.nodefore == cnode)
            assert cellmark == cbar.GetCellMarkRightL(cnode == cbar.nodeback)
            assert cbar != self.bar
            Dcount += 1
            assert Dcount < 1000
            
        self.cbar = cbar
        self.cnode = cnode
        self.vc = self.cbar.nodemid.p - self.bar.nodemid.p

        # check rest colour is consistent
        if __debug__:
            lbar = cbar.GetForeRightBL(cbar.nodefore == cnode)
            lnode = cnode
            while lbar != self.bar:  
                lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
                lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
                assert cellmark == lbar.GetCellMarkRightL(lnode == lbar.nodeback)
                Dcount += 1
                assert Dcount < 1000

    def shoulddivide(self, rd, contourdelta, contourdotdiff, b2dcontournormals):
        if self.vc.Len() <= contourdelta:
            return False
        if b2dcontournormals:
            nd = P2.Dot(P2.ZNorm(P2(self.bar.nodemid.pointzone.v.x, self.bar.nodemid.pointzone.v.y)), P2.ZNorm(P2(self.cbar.nodemid.pointzone.v.x, self.cbar.nodemid.pointzone.v.y)))
        else:
            nd = P3.Dot(self.bar.nodemid.pointzone.v, self.cbar.nodemid.pointzone.v)/(rd*rd)
        if nd >= contourdotdiff:
            return False
        return True
            
    def calcsplithalfpos(self):
        # midpoint, but could be informed by multiple subdivisions on the rate of curvature
        # should also measure the bc.calccellcuttangency() and find a spot where it is properly wide and doesn't need to be tested
        vch = P3.Dot(self.vc, self.bar.nodemid.p + self.vc*0.5)
        
        # find the splitting points 
        self.ktnodetop, self.ktbartop, self.ktbartoplam = PolyCutBar(self.bar, self.node, self.vc, vch, True)
        self.ktnodebot, self.ktbarbot, self.ktbarbotlam = PolyCutBar(self.cbar, self.cbar.GetNodeFore(self.cbar.nodeback == self.cnode), self.vc, vch, False)

    def setmincellcolour(self):
        assert self.GetCellMark() == None
        crosscolourset = getcrosscolours(self.bar, self.node, None)
        newcellcolour = 0
        while newcellcolour in crosscolourset:
            newcellcolour += 1
        newcellmark = barmesh.CellMark(newcellcolour)
        barmesh.setcellmark(self.bar, self.node, newcellmark, None)
        assert self.GetCellMark() == newcellmark
        
    # either splits the bar or selects one of its end-nodes that are close enough to the end
    # (could break system by setting a large lamendgap)
    # sets self.leadsplitbartop, self.splitnodetop, self.bsplitnnodetopnew, self.leadsplitbarbot, self.splitnodebot, self.bsplitnnodebotnew
    def makesplitnodes(self, btop, bm, lamendgap):
        if btop:
            ktnode, ktbar, ktbarlam = self.ktnodetop, self.ktbartop, self.ktbartoplam
        else:
            ktnode, ktbar, ktbarlam = self.ktnodebot, self.ktbarbot, self.ktbarbotlam
        assert not ktbar.bbardeleted
        ktnodefore = ktbar.GetNodeFore(ktbar.nodeback == ktnode)
        
        if lamendgap < ktbarlam < 1 - lamendgap:
            splitnode = bm.InsertNodeIntoBarF(ktbar, bm.NewNode(Along(ktbarlam, ktnode.p, ktnodefore.p)), True)
            assert ktbar.bbardeleted
            leadsplitbar = ktbar.GetForeRightBL(ktnode == ktbar.nodefore)
            assert not leadsplitbar.bbardeleted
            bnewnode = True
            
        else:
            if ktbarlam < 0.5:
                splitnode = ktnode
                if ktnode == ktbar.nodeback:
                    leadsplitbar = ktbar.GetBarBackRight()
                else:
                    leadsplitbar = ktbar.GetBarForeLeft()
                assert splitnode == leadsplitbar.nodeback or splitnode == leadsplitbar.nodefore
            else:
                splitnode = ktnodefore
                leadsplitbar = ktbar
                assert splitnode == leadsplitbar.nodeback or splitnode == leadsplitbar.nodefore
            bnewnode = False

        if btop:
            self.leadsplitbartop, self.splitnodetop = leadsplitbar, splitnode
            assert self.splitnodetop == self.leadsplitbartop.nodeback or self.splitnodetop == self.leadsplitbartop.nodefore
            self.bsplitnnodetopnew = bnewnode
        else:
            self.leadsplitbarbot, self.splitnodebot = leadsplitbar, splitnode            
            assert self.splitnodebot == self.leadsplitbarbot.nodeback or self.splitnodebot == self.leadsplitbarbot.nodefore
            self.bsplitnnodebotnew = bnewnode
        return bnewnode
        
    def calccellcuttangency(self):  # could be checked before splitting on the lines to give an option of picking a less tangential position
        assert not self.leadsplitbartop.bbardeleted and not self.leadsplitbarbot.bbardeleted
        node1, bar1, node2, bar2 = self.splitnodetop, self.leadsplitbartop, self.splitnodebot, self.leadsplitbarbot
        bar1a = bar1.GetForeRightBL(bar1.nodefore == node1)
        bar2a = bar2.GetForeRightBL(bar2.nodefore == node2)
        vtopbot = P3.ZNorm(self.splitnodebot.p - self.splitnodetop.p)
        
        db1 = -node1.cperpbardotN(bar1, vtopbot)
        db1a = node1.cperpbardotN(bar1a, vtopbot)
        db2 = node2.cperpbardotN(bar2, vtopbot)
        db2a = -node2.cperpbardotN(bar2a, vtopbot)
        
        res = min(db1, db1a, db2, db2a)
        assert res > -0.0001
        return res
        
        
    def makecellsplittingbar(self, bm):
        if self.splitnodetop.i < self.splitnodebot.i:
            splitbar = bm.MakeBarBetweenNodesF(self.splitnodetop, self.leadsplitbartop, self.splitnodebot, self.leadsplitbarbot)
        else:
            splitbar = bm.MakeBarBetweenNodesF(self.splitnodebot, self.leadsplitbarbot, self.splitnodetop, self.leadsplitbartop)
            
        # allocate a new colour here
        crosscoloursforrightcell = getcrosscolours(splitbar, splitbar.nodeback, None)
        crosscoloursforleftcell = getcrosscolours(splitbar, splitbar.nodefore, None)
        newcellcolour = -1
        cellmark = self.GetCellMark()
        assert cellmark.cellcolour != -1
        while True:   # this assigns first free colour to one or other side, then breaks
            newcellcolour += 1
            if newcellcolour == cellmark.cellcolour:
                continue
            if newcellcolour not in crosscoloursforrightcell:
                newcellmark = barmesh.CellMark(newcellcolour)
                barmesh.setcellmark(splitbar, splitbar.nodeback, newcellmark, None)
                assert splitbar.cellmarkleft == None
                splitbar.cellmarkleft = cellmark
                assert cellmark.cellcolour not in getcrosscolours(splitbar, splitbar.nodefore, cellmark)
                break
            if newcellcolour not in crosscoloursforleftcell:
                newcellmark = barmesh.CellMark(newcellcolour)
                barmesh.setcellmark(splitbar, splitbar.nodefore, newcellmark, None)
                assert splitbar.cellmarkright == None
                splitbar.cellmarkright = cellmark
                assert cellmark.cellcolour not in getcrosscolours(splitbar, splitbar.nodeback, cellmark)
                break
        if newcellcolour > bm.maxcellcolour:
            bm.maxcellcolour = newcellcolour
        return splitbar
        


class BarMeshSlicer:
    def __init__(self, bm, tgf, rd, rd2, contourdotdiff, contourdelta, lamendgap):
        self.bm = bm      # BarMesh
        self.tgf = tgf    # TriangleGroupFuncs
        self.rd = rd
        self.rd2 = rd2
        assert rd2 > rd
        self.rdsq = self.rd*self.rd

        self.contourdotdiff = contourdotdiff
        self.contourdelta = contourdelta
        self.barsplitdelta = contourdelta
        self.lamendgap = lamendgap
        
        self.bm.maxcellcolour = -1
        self.barpolycuts = [ ]

        self.totalproctime = 0.0
        
        self.pztime = 0.0
        self.pzcalls = 0
        self.cbtime = 0.0
        self.cbztime = 0.0
        self.cbcalls = 0
        
    def MakePointZoneRFS(self, nodes):
        ctime1 = time.clock()
        for node in nodes:  # should be executed in parallel as a batch
            self.pzcalls += 1
            pz = node.pointzone = barmesh.PointZone(0, self.rd2, None)
            self.tgf.DistP(pz, node.p)
            if pz.r > self.rd:
                pz.izone = barmesh.PZ_BEYOND_R
            else:
                pz.izone = barmesh.PZ_WITHIN_R
        self.pztime += time.clock() - ctime1
    
    def CutbarRFS(self, bars):
        lambars = [ ]
        
        ctime1 = time.clock()
        for bar in bars:  # should be executed in parallel as a batch
            self.cbcalls += 1
            assert not bar.bbardeleted
            assert IsCutBar(bar, barmesh.PZ_BEYOND_R)
            if bar.nodefore.pointzone.izone == barmesh.PZ_BEYOND_R:
                assert bar.nodeback.pointzone.v is not None
                lam = 1 - self.tgf.Cutpos(bar.nodefore.p, bar.nodeback.p - bar.nodefore.p, bar.nodeback.p + bar.nodeback.pointzone.v, self.rd)
            else:
                assert bar.nodefore.pointzone.v is not None
                lam = self.tgf.Cutpos(bar.nodeback.p, bar.nodefore.p - bar.nodeback.p, bar.nodefore.p + bar.nodefore.pointzone.v, self.rd)
            assert -0.001 <= lam <= 1.001, lam
            lambars.append((lam, bar))
        self.cbtime += time.clock() - ctime1
        
        ctime2 = time.clock()
        for lam, bar in lambars:
            if 0 < lam < 1:
                bar.nodemid = barmesh.Node(bar.nodeback.p * (1 - lam) + bar.nodefore.p * lam, -1)
                bar.nodemid.pointzone = barmesh.PointZone(0, self.rd2, None)
                self.tgf.DistP(bar.nodemid.pointzone, bar.nodemid.p)
            
            # cases where we have gone off the end, so push back to one of the ends
            elif lam < 0.5:
                bar.nodemid = bar.nodeback
            else:
                bar.nodemid = bar.nodefore
                
            bar.midcontournumber = -1
        self.cbztime += time.clock() - ctime2
        
    def initializecutsanddistances(self):
        ctime0 = time.clock()
        
        # set the closest approach value for each node
        self.MakePointZoneRFS(self.bm.nodes)
    
        # subdivide some of the bars where distance 
        dcbars = self.bm.bars
        while dcbars:
            dcbars = self.splitbarsdirectionchangesR(dcbars)

        barcontourcells = [ ]
        for bar in self.bm.bars:
            if bar.bbardeleted:
                continue
            if IsCutBar(bar, barmesh.PZ_BEYOND_R):
                barcontourcells.append(bar)
        self.CutbarRFS(barcontourcells)

        # set the cell colours and marks on all
        for bar in barcontourcells:
            bc = BarMeshCell(bar, barmesh.PZ_BEYOND_R)
            if not bc.RefersToOuterCell():
                if bc.GetCellMark() == None:  # can be not None if we have two cuts in same cell
                    bc.setmincellcolour()
                    if bc.GetCellMark().cellcolour > self.bm.maxcellcolour:
                        self.bm.maxcellcolour = bc.GetCellMark().cellcolour
                self.barpolycuts.append(bar)
                bar.midcontournumber = -2
        self.totalproctime += time.clock() - ctime0
        
    def getbarendclos(self, barstosplit):  # for plotting closest approaches
        pcs = [ ]
        for bar in barstosplit:
            if bar.bbardeleted:
                continue
            if not ((bar.nodeback.pointzone.izone == barmesh.PZ_BEYOND_R) and (bar.nodefore.pointzone.izone == barmesh.PZ_BEYOND_R)):
                continue
            vb = bar.nodefore.p - bar.nodeback.p
            vbsq = vb.Lensq()
            if bar.nodeback.pointzone.v is not None:
                lampcb = P3.Dot(bar.nodeback.pointzone.v, vb) / vbsq
                if 0 < lampcb < 1:
                    dvb = bar.nodeback.pointzone.v - vb * lampcb
                    assert abs(P3.Dot(dvb, vb)) < 0.001
                    dvbsq = dvb.Lensq()
                    if dvbsq < self.rdsq:
                        pcs.append(bar.nodeback.p + bar.nodeback.pointzone.v)
                        
            if bar.nodefore.pointzone.v is not None:
                lampcf = -P3.Dot(bar.nodefore.pointzone.v, vb) / vbsq
                if 0 < lampcf < 1:
                    dvf = bar.nodefore.pointzone.v + vb * lampcf
                    assert abs(P3.Dot(dvf, vb)) < 0.001
                    dvfsq = dvf.Lensq()
                    if dvfsq < self.rdsq:
                        pcs.append(bar.nodeback.p + bar.nodefore.pointzone.v)
        return pcs
        
    # applies where sampling is very coarse and seems to find places where the line between goes too close according to the normal that there must be a point within distance rd of it
    def splitbarsdirectionchangesR(self, barstosplit):
        barsplitsdc = [ ]
        bsdsq = self.barsplitdelta*self.barsplitdelta
        for bar in barstosplit:
            if bar.bbardeleted:
                continue
            if not ((bar.nodeback.pointzone.izone == barmesh.PZ_BEYOND_R) and (bar.nodefore.pointzone.izone == barmesh.PZ_BEYOND_R)):
                continue
            vb = bar.nodefore.p - bar.nodeback.p
            vbsq = vb.Lensq()
            if vbsq < bsdsq:
                continue
            lamclossplit = -1
            dvbsq = -1
            if bar.nodeback.pointzone.v is not None:
                lampcb = P3.Dot(bar.nodeback.pointzone.v, vb) / vbsq
                if 0 < lampcb < 1:
                    dvb = bar.nodeback.pointzone.v - vb * lampcb
                    assert abs(P3.Dot(dvb, vb)) < 0.001
                    dvbsq = dvb.Lensq()
                    if dvbsq < self.rdsq:
                        lamclossplit = lampcb
                        
            if bar.nodefore.pointzone.v is not None:
                lampcf = -P3.Dot(bar.nodefore.pointzone.v, vb) / vbsq
                if 0 < lampcf < 1:
                    dvf = bar.nodefore.pointzone.v + vb * lampcf
                    assert abs(P3.Dot(dvf, vb)) < 0.001
                    dvfsq = dvf.Lensq()
                    if dvfsq < self.rdsq:
                        if lamclossplit == -1 or dvfsq < dvbsq:
                            lamclossplit = 1 - lampcf
                            
            if lamclossplit != -1:
                barsplitsdc.append((bar, lamclossplit))
            
        # split the bars
        nowsplitbars = [ ]
        splitnodes = [ ]
        for bar, lamc in barsplitsdc:
            splitnode = self.bm.InsertNodeIntoBarF(bar, self.bm.NewNode(Along(lamc, bar.nodeback.p, bar.nodefore.p)), True)
            assert bar.bbardeleted
            splitnodes.append(splitnode)
            nowsplitbars.append(bar.barforeright)
            nowsplitbars.append(bar.barbackleft)
            
        # make the measurements on all the new splitnodes
        self.MakePointZoneRFS(splitnodes)
            
        # make the measurements on all the new bars
        newbarpolycuts = [ ]
        for bar in nowsplitbars:
            if IsCutBar(bar, barmesh.PZ_BEYOND_R):
                newbarpolycuts.append(bar)
                
        self.CutbarRFS(newbarpolycuts)
        #print("splitbarsdirectionchanges", len(barstosplit), len(nowsplitbars))
        return nowsplitbars
            
                
    def splitbarpolyscolour(self, currentcolour):
        ctime0 = time.clock()
        cellmarksdone = set()

        Dsbarsouttol = set(bar  for bar in self.barpolycuts  if not bar.bbardeleted)
        DsbarsouttolW = set(bar  for bar in self.bm.bars  if not bar.bbardeleted and IsCutBar(bar, barmesh.PZ_BEYOND_R) and not RefersToOuterCellBar(bar, barmesh.PZ_BEYOND_R) and bar.midcontournumber == -2)
        assert Dsbarsouttol == DsbarsouttolW, (len(Dsbarsouttol), len(DsbarsouttolW))

        nextbarpolycuts = [ ]  # replaces self.barpolycuts after we've passed through
        barpolycutsthiscolour = [ ]
        # separate out the current colour from everything else 
        for bar in self.barpolycuts:
            if bar.bbardeleted:  # replacements are already in the list
                continue
            assert not RefersToOuterCellBar(bar, barmesh.PZ_BEYOND_R)
            assert IsCutBar(bar, barmesh.PZ_BEYOND_R)
            assert bar.midcontournumber == -2
            bccellmark = GetCellMarkBar(bar, barmesh.PZ_BEYOND_R)  # can change after bc init due to rescaling
            if bccellmark.cellcolour != currentcolour:
                nextbarpolycuts.append(bar) # fall through
            elif bccellmark in cellmarksdone:  # can't do two cuts in same cell
                nextbarpolycuts.append(bar) # fall through
            else:
                bc = BarMeshCell(bar, barmesh.PZ_BEYOND_R)
                bc.MakeCutSide()
                if bc.shoulddivide(self.rd, self.contourdelta, self.contourdotdiff, self.tgf.Isb2dcontournormals()):
                    assert currentcolour not in getcrosscolours(bc.bar, bc.node, bccellmark)
                    cellmarksdone.add(bccellmark)
                    barpolycutsthiscolour.append(bc)
                else:
                    bar.midcontournumber = -1
            
        # find the places we'll want to split in the cells (read-only operation)
        currpolycuts = [ ]
        for bc in barpolycutsthiscolour:
            bccellmark = bc.GetCellMark()
            bc.calcsplithalfpos()   # follows round and finds the top and bottom cuts to the polygon
            bc.makesplitnodes(True, self.bm, self.lamendgap)
            bc.makesplitnodes(False, self.bm, self.lamendgap)
            assert bccellmark == bc.leadsplitbartop.GetCellMarkRightL(bc.leadsplitbartop.nodeback != bc.splitnodetop), "hhheeee"
            assert bccellmark == bc.leadsplitbarbot.GetCellMarkRightL(bc.leadsplitbarbot.nodeback != bc.splitnodebot), "hhhoooo"
            assert bc.bar.midcontournumber == -2
            currpolycuts.append(bc)
        
        # gather together all new nodes and bars created so we can apply the calculations on them
        newsplitnodes = [ ]
        newsplitbars = [ ]
        for bc in currpolycuts:
            if bc.bsplitnnodetopnew:
                newsplitnodes.append(bc.splitnodetop)
                newsplitbars.append(bc.leadsplitbartop)
                newsplitbars.append(bc.leadsplitbartop.GetForeRightBL(bc.leadsplitbartop.nodefore == bc.splitnodetop))
            if bc.bsplitnnodebotnew:
                newsplitnodes.append(bc.splitnodebot)
                newsplitbars.append(bc.leadsplitbarbot)
                newsplitbars.append(bc.leadsplitbarbot.GetForeRightBL(bc.leadsplitbarbot.nodefore == bc.splitnodebot))
            if not bc.bar.bbardeleted:
                nextbarpolycuts.append(bc.bar)  # copy old split thing back in

        # make the measurements on all the new nodes
        self.MakePointZoneRFS(newsplitnodes)
        
        # join the new nodes with the new bar
        for bc in currpolycuts:
            if bc.calccellcuttangency() > 0.001:  # should be established in calcsplithalfpos, so this should be an assert
                newsplitbar = bc.makecellsplittingbar(self.bm)
                newsplitbars.append(newsplitbar)
            else:
                print("skipping split tang", bc.node.p)
                
        # resplit any that appear on either side of the bound case
        dcsplitbars = newsplitbars
        while dcsplitbars:
            dcsplitbars = self.splitbarsdirectionchangesR(dcsplitbars)
            newsplitbars.extend(dcsplitbars)

        # make the measurements on all the new bars
        newbarpolycuts = [ ]   # we add more to this list later on from adjacent cells
        for bar in newsplitbars:
            if bar.bbardeleted:
                continue
            if IsCutBar(bar, barmesh.PZ_BEYOND_R):
                bar.midcontournumber = -2
                newbarpolycuts.append(bar)

        # do the batch cut mid points of all the new bars
        self.CutbarRFS(newbarpolycuts)

        # make the new barmeshcells
        cellmarkset = [ ]
        cellhitcounterG = self.bm.cellhitcounterend
        self.bm.cellhitcounterend += 1
        for bar in newsplitbars:
            if bar.bbardeleted:
                continue
            
            if bar.GetForeRightBL(True):
                cellmarkright = bar.GetCellMarkRightL(True)
                if cellmarkright == None:
                    cellmarkright = setmincellcolourBar(bar, bar.nodeback)
                    if cellmarkright.cellcolour > self.bm.maxcellcolour:
                        self.bm.maxcellcolour = cellmarkright.cellcolour
                if cellmarkright.cellhitcounter != cellhitcounterG:
                    cellmarkright.cellhitcounter = cellhitcounterG
                    cellmarkset.append((bar, bar.nodeback, cellmarkright))

            if bar.GetForeRightBL(False):
                cellmarkleft = bar.GetCellMarkRightL(False)
                if cellmarkleft == None:
                    cellmarkleft = setmincellcolourBar(bar, bar.nodefore)
                    if cellmarkleft.cellcolour > self.bm.maxcellcolour:
                        self.bm.maxcellcolour = cellmarkleft.cellcolour
                if cellmarkleft.cellhitcounter != cellhitcounterG:
                    cellmarkset.append((bar, bar.nodefore, cellmarkleft))

        for (bar, node, Dcellmark) in cellmarkset:
            if bar.bbardeleted:
                continue   # the pieces will have been put in as well
            assert Dcellmark == bar.GetCellMarkRightL(bar.nodeback == node)
            lbar, lnode = bar, node
            Dcount = 0
            while True:
                lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
                lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
                assert Dcellmark == lbar.GetCellMarkRightL(lnode == lbar.nodeback)
                if not RefersToOuterCellBar(bar, barmesh.PZ_BEYOND_R):
                    if lbar.midcontournumber != -2:  # already included
                        if IsCutBar(lbar, barmesh.PZ_BEYOND_R):
                            if lbar == bar or GetCellMarkBar(lbar, barmesh.PZ_BEYOND_R) == Dcellmark:  # only keep those on the marked cell side
                                nextbarpolycuts.append(lbar)
                                lbar.midcontournumber = -2
                if lbar == bar:
                    break
                Dcount += 1
                assert Dcount < 1000
                
        # copy over for next round
        self.barpolycuts = nextbarpolycuts
        self.totalproctime += time.clock() - ctime0

        return len(newsplitbars)


    # this would be necessary to use when the main loop is leaking (not happening now)
    def PutBackInvalidIntolValues(self):
        ninvalidintolvalues = 0
        for bar in self.bm.bars:
            if bar.bbardeleted or not IsCutBar(bar, barmesh.PZ_BEYOND_R) or RefersToOuterCellBar(bar, barmesh.PZ_BEYOND_R):
                continue
            lbc = BarMeshCell(bar, barmesh.PZ_BEYOND_R)
            lbc.MakeCutSide()
            if lbc.shoulddivide(self.rd, self.contourdelta, self.contourdotdiff, self.tgf.Isb2dcontournormals()):
                self.barpolycuts.append(bc.bar)  # brand new cell object without the cutside
                ninvalidintolvalues += 1
        print("Found %d not intolerance values" % ninvalidintolvalues)
        return ninvalidintolvalues
        
    def CheckTolerance(self):
        ncontsegs, ntolbad, ntolworking = 0, 0, 0
        conts = [ ]
        for bar in self.bm.bars:
            if not bar.bbardeleted and IsCutBar(bar, barmesh.PZ_BEYOND_R) and not RefersToOuterCellBar(bar, barmesh.PZ_BEYOND_R):
                ncontsegs += 1
                bc = BarMeshCell(bar, barmesh.PZ_BEYOND_R)
                bc.MakeCutSide()
                if bc.shoulddivide(self.rd, self.contourdelta, self.contourdotdiff, self.tgf.Isb2dcontournormals()):
                    if bar not in self.barpolycuts:
                        if ntolbad == 0:
                            print("out of tolerance", bc.node.p)
                        assert bar.midcontournumber == -2, bar.midcontournumber
                        ntolbad += 1
                    else:
                        ntolworking += 1
                    conts.append([bc.bar.nodemid.p, bc.cbar.nodemid.p])
        print("nsegs:", ncontsegs, "bad:", ntolbad, "working:", ntolworking)
        return conts
        
    # repeatedly pass over each colour (which defines a subset of non-adjacent cells) and split them if out of tolerance
    def fullmakeslice(self):
        self.initializecutsanddistances()
        if self.bm.maxcellcolour == -1:
            return    # nothing to cut
        currentcolour = -1
        ncount = 0
        Dnsplitslist = [ ]
        while self.barpolycuts:   # the stack of remaining BarMeshCells we are working down through
            assert self.bm.maxcellcolour != -1
            currentcolour += 1
            if currentcolour > self.bm.maxcellcolour:
                currentcolour = 0
            Dnsplits = self.splitbarpolyscolour(currentcolour)
            Dnsplitslist.append(Dnsplits)
            ncount += 1
            if ncount == 180:
                print("breaking after", ncount, "iterations with", len(self.barpolycuts), "barpolycuts remaining")
                break
        self.CheckTolerance()
        # if you think the result is out of tolerance, try:
        #   sendactivity(contours=bms.CheckTolerance(), materialnumber=3)
        # and bms.PutBackInvalidIntolValues() before rerunning the while-loop
        print(sum(Dnsplitslist), Dnsplitslist)

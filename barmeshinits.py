import barmesh
from basicgeo import P3
from math import sqrt

def BuildCubeBarMesh(bm, xlo, xhi, ylo, yhi, zlo, zhi):
    w = 1/sqrt(3)
    def CubeNode(ncode):
        bxlo, bylo, bzlo = (ncode&0b100, ncode&0b010, ncode&0b001)
        node = bm.NewNode(P3(xlo if bxlo else xhi, ylo if bylo else yhi, zlo if bzlo else zhi))
        node.pointzone = barmesh.PointZone(0, -1, P3(w if bxlo else -w, w if bylo else -w, w if bzlo else -w))
        return node

    nodenames = [0b000,  0b001, 0b010, 0b100,  0b110, 0b101, 0b011,  0b111]
    noded = [ CubeNode(ncode)  for ncode in nodenames ]
    def thnext(i):  return i+1 if i != 2 else 0
    def thprev(i):  return i-1 if i != 0 else 2
    bars0 = [ barmesh.Bar(noded[0], noded[i+1])  for i in range(3) ]
    barsM = [ ]
    for i in range(3):
        barsM.append(barmesh.Bar(noded[i+1], noded[thnext(i)+4]))
        barsM.append(barmesh.Bar(noded[i+1], noded[thprev(i)+4]))
    bars1 = [ barmesh.Bar(noded[i+4], noded[7])  for i in range(3) ]
    for i in range(3):
        bars0[i].barbackleft = bars0[thprev(i)]
        bars0[i].barforeright = barsM[i*2+1]
        barsM[i*2].barbackleft = bars0[i]
        barsM[i*2].barforeright = bars1[thnext(i)]
        barsM[i*2+1].barbackleft = barsM[i*2]
        barsM[i*2+1].barforeright = barsM[thnext(i)*2]
        bars1[i].barbackleft = barsM[thnext(i)*2+1]
        bars1[i].barforeright = bars1[thnext(i)]
    bm.bars.extend(bars0)
    bm.bars.extend(barsM)
    bm.bars.extend(bars1)

def GetAllBarNodeFaces(bm):
    barnodeset = set()
    for bar in bm.bars:
        if bar.barforeright:
            barnodeset.add((bar, bar.nodeback))
        if bar.barbackleft:
            barnodeset.add((bar, bar.nodefore))
    barnodefaces = [ ]
    while barnodeset:
        bar, node = barnodeset.pop()
        barnodefaces.append((bar, node))
        lbar, lnode = bar, node
        Dcount = 0
        while True:
            lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
            if lnode == node:
                break
            lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
            barnodeset.remove((lbar, lnode))
            Dcount += 1
            assert Dcount < 1000
    return barnodefaces

    
def StarSplitFace(bm, bar, node):
    lbar, lnode = bar, node
    Dcounts = 0
    pnodebars = [ ]
    while True:
        pnodebars.append((lnode, lbar))
        lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
        if lnode == node:
            break
        lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
        Dcounts += 1
        assert Dcounts < 1000
    cpt = sum((lnode.p  for (lnode, lbar) in pnodebars), P3(0,0,0))*(1/len(pnodebars))
    cnode = bm.NewNode(cpt)
    cnode.pointzone = barmesh.PointZone(0, -1, P3(0,0,0))
    np = len(pnodebars)
    cbars = [ barmesh.Bar(lnode, cnode)  for lnode, lbar in pnodebars ]
    for i in range(np):
        cbars[i].barforeright = cbars[i-1  if i!=0  else np-1]
        cbars[i].barbackleft = pnodebars[i][1] # cbars[i-1  if i!=0  else np-1]
        lnode, lbar = pnodebars[i]
        lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
        lbar.SetForeRightBL(lbar.nodefore == lnode, cbars[i+1  if i!=np-1  else 0])
    bm.bars.extend(cbars)

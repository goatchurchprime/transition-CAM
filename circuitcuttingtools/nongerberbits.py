from basicgeo import Partition1


def InvAlong(v, a, b):
    return (v - a)/(b - a)
class ProbingZ:
    def __init__(self, fname):
        fin = open(fname)
        probepts = [ list(map(float, l.split()[:3]))  for l in fin ]
        pxs = sorted(list(set(p[0]  for p in probepts)))
        pys = sorted(list(set(p[1]  for p in probepts)))
        parx = Partition1(pxs[0], pxs[-1], len(pxs)-1)
        pary = Partition1(pys[0], pys[-1], len(pys)-1)
        perrx = max(a-b  for a, b in zip(parx.vs, pxs))
        perry = max(a-b  for a, b in zip(pary.vs, pys))
        minz, maxz = min(p[2]  for p in probepts), max(p[2]  for p in probepts)
        print("Minmax z", minz, maxz, "Error on probe partition evenness XY", perrx, perry)
        parx.vs = pxs
        pary.vs = pys
        self.probeptsmap = dict(((p[0], p[1]), p[2])  for p in probepts)
        self.parx, self.pary = parx, pary
        
    def InterpZ(self, x, y):
        parx, pary = self.parx, self.pary
        x = max(self.parx.vs[0], min(self.parx.vs[-1], x))
        y = max(self.pary.vs[0], min(self.pary.vs[-1], y))
        ix = parx.GetPart(x)
        iy = pary.GetPart(y)
        lamx = InvAlong(x, parx.vs[ix], parx.vs[ix+1])
        lamy = InvAlong(y, pary.vs[iy], pary.vs[iy+1])
        z00 = self.probeptsmap[(parx.vs[ix], pary.vs[iy])]
        z01 = self.probeptsmap[(parx.vs[ix], pary.vs[iy+1])]
        z10 = self.probeptsmap[(parx.vs[ix+1], pary.vs[iy])]
        z11 = self.probeptsmap[(parx.vs[ix+1], pary.vs[iy+1])]
        z = (1-lamx)*(z00*(1-lamy) + z01*lamy) + lamx*(z10*(1-lamy) + z11*lamy)
        return z

def post(conts, zlo, zhi, ngcF, prbz):
    fout = open(ngcF, "w")
    fout.write("G1Z%.3fF800\n" % zhi)
    for cont in conts:
        p = cont[0]
        fout.write("G0X%.3fY%.3f\n" % p)
        fout.write("G1Z%.3fF800\n" % (prbz.InterpZ(p[0], p[1])+zlo))
        for p in cont[1:]:
            fout.write("X%.3fY%.3fZ%.3f\n" % (p[0], p[1], prbz.InterpZ(p[0], p[1])+zlo))
        fout.write("G0Z%.3f\n" % zhi)
    fout.write("M2\n")    
    fout.close()


def ReorderConts(lconts):
    x, y = -132, 300.1
    cconts = [ ]
    while lconts:
        m, i, j = min((abs(p[0]-x)+abs(p[1]-y), i, j)  for i, cont in enumerate(lconts)  for j, p in enumerate(cont))
        cont = lconts.pop(i)
        if len(cont) != 1:
            cont = cont[j:-1]+cont[:j+1]
            x, y = cont[-1]
            for k in range(1, min(10, len(cont))):  # run off section
                p = cont[k]
                cont.append(p)
                if abs(p[0]-x)+abs(p[1]-y)>1.8:
                    break
        x, y = cont[-1]
        cconts.append(cont)
    return cconts

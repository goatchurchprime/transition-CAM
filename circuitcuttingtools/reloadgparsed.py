import sys, os
import pickle
from basicgeo import P2

dr = "/home/goatchurch/geom3d/maghullmachinecalibrations/circuitboards"
geofile = os.path.join(dr, "pcb2.pckl")
(ocont, conts, dpts) = pickle.load(open(geofile, "rb"))

xmin, ymin = min(p[0]  for p in ocont), min(p[1]  for p in ocont)

def ThinTooshort(cont):
    res = [ cont[0] ]
    for p in cont[1:]:
        if (res[-1]-p).Len() < 1e-5:
            continue
        else:
            res.append(p)
    if res[-1] != res[0]:
        res[-1] = res[0]
    return res

conts = [ThinTooshort([ P2(p[0]-xmin, p[1]-ymin)  for p in cont])  for cont in conts]
ocont = ThinTooshort([P2(p[0]-xmin, p[1]-ymin)  for p in ocont])
dpts = [P2(p[0]-xmin, p[1]-ymin)  for p in dpts]
sendactivity(contours=conts)
sendactivity(contours=[ocont])
sendactivity(points=dpts)


from vor2d import DFullBuildVB
conts.append(ocont)
vbs, vbm = DFullBuildVB(conts)


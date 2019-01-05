import sys, os  # has to be Python27
sys.path.append("/home/goatchurch/geom3d/circuitcuttingtools")
sys.path.append("/home/goatchurch/geom3d/flatcam")
sys.path.append("/home/goatchurch/geom3d/barmesh")
from gerbergetbits import GetIsolationCuts, GetHoleCuts, GetEdgeCuts
from nongerberbits import ProbingZ, post, ReorderConts
from basicgeo import I1

dr = "/home/goatchurch/geom3d/maghullmachinecalibrations/circuitboards/bbcape2"
f = os.path.join(dr, "beaglebone_cape_without gnds-F_Cu.gtl")
f1 = os.path.join(dr, "beaglebone_cape_without gnds-B_Cu.gbl")
d = os.path.join(dr, "beaglebone_cape.drl")
ed = os.path.join(dr, "beaglebone_cape_without gnds-Edge_Cuts.gbr")
prb = os.path.join(dr, "probing8.txt")
ngcF = os.path.join(dr, "bbcape.ngc")

prbz = ProbingZ(prb)
sendactivity(points=[(x,y,z*50)  for (x,y), z in prbz.probeptsmap.items()])
xrg, yrg = I1(prbz.parx.vs[0], prbz.parx.vs[-1]), I1(prbz.pary.vs[0], prbz.pary.vs[-1])
pp = [(xrg.Along(i/50.0), yrg.Along(j/50.0))  for i in range(51)  for j in range(51)]
sendactivity(points=[(x,y,prbz.InterpZ(x,y)*50)  for x,y in pp], materialnumber=1)

copperorg = (-204, 260)
circuitorg = (175, -160)
rgs = (139, 97)
dx, dy = copperorg[0]-circuitorg[0], copperorg[1]-circuitorg[1]


conts, contoffsets = GetIsolationCuts(f, 0.12, dx, dy)
sendactivity(contours=conts, materialnumber=0)
sendactivity(contours=contoffsets, materialnumber=1)

conts1, contoffsets1 = GetIsolationCuts(f1, 0.12, dx, dy)
sendactivity(contours=conts1, materialnumber=0)
sendactivity(contours=contoffsets1, materialnumber=1)

conts, contoffsets = GetEdgeCuts(ed, 0.4, dx, dy)
econt, econtoffset = conts[0], contoffsets[0]
sendactivity(contours=[econt], materialnumber=2)
sendactivity(contours=[econtoffset], materialnumber=3)

holediampts = GetHoleCuts(d, dx, dy)
sendactivity(points=sum(holediampts.values(), []))

xlo, xhi = min(p[0]  for p in econt), max(p[0]  for p in econt)
ylo, yhi = min(p[1]  for p in econt), max(p[1]  for p in econt)
print("xyrange", (xlo, ylo), (xhi, yhi))

# Location pins work
#sendactivity(points=holediampts[3.0479999999999996], materialnumber=1)
#ymid = (yhi+ylo)*0.5
#locationpts = holediampts[3.0479999999999996][:]
#del locationpts[2:4]
locationpts = [(-199.93, 352.69), (-199.93, 263.79), (-69.12, 352.69), (-69.12, 263.79)]
ymid = 308.24
print(locationpts[0][1]-ymid, locationpts[1][1]-ymid)
midloc = (-122, ymid)
leftloc = (-186, 328)
locpins = locationpts+[midloc, leftloc, (leftloc[0],2*ymid-leftloc[1]), (-199,ymid)]
sendactivity(points=locpins, materialnumber=2)
post([[c]  for c in locpins], -0.1, 2, ngcF, prbz)


# back face isolation tracks
lconts = [[(p[0],ymid*2-p[1])  for p in cont]  for cont in contoffsets1]
cconts = ReorderConts(lconts)
post(cconts, -0.1, 2, ngcF, prbz)

# back face drill pecks
dconts = [[(p[0],ymid*2-p[1])]  for p in sum(holediampts.values(), [])]
cconts = ReorderConts(dconts)
post(cconts, -0.14, 2, ngcF, prbz)

# front face isolation tracks
lconts = [[(p[0],p[1])  for p in cont]  for cont in contoffsets]
cconts = ReorderConts(lconts)
post(cconts, -0.1, 2, ngcF, prbz)

# front face isolation tracks redo
lconts = [[(p[0],p[1])  for p in cont]  for cont in contoffsets  if cont[0][1]<268 and -158<cont[0][0]<-120]
cconts = ReorderConts(lconts)
sendactivity(contours=cconts, materialnumber=2)
post(cconts, -0.15, 2, ngcF, prbz)

# front face drill pecks
dconts = [[(p[0],p[1])]  for p in sum(holediampts.values(), [])]
cconts = ReorderConts(dconts)
post(cconts, -0.166, 2, ngcF, prbz)

# front face full drills
post(cconts, -2.7, 2, ngcF, prbz)

# front face edge scoring (should have done a back face one)
post([econtoffset], -0.1, 2, ngcF, prbz)





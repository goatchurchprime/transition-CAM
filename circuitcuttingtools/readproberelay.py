# need to copy over my probing file off the BB on the machine

import re, os

fdir = "/home/goatchurch/geom3d/flatcam/adrianscircuitdoorbot/"

fn = "adrianedgecutsfront.ngc"
fns = "adrianedgecutsfront_sh.ngc"

fin = open(os.path.join(fdir, fn))
xs = [float(mx.group(1))  for mx in [ re.search("X([\-\d\.]+)", l)  for l in fin ]  if mx]
xlo, xhi = min(xs), max(xs)
fin = open(os.path.join(fdir, fn))
ys = [float(my.group(1))  for my in [ re.search("Y([\-\d\.]+)", l)  for l in fin ]  if my]
ylo, yhi = min(ys), max(ys)
print((xlo,xhi), (ylo,yhi))
zlo, zhi = -0.15, 2

#G92X195Y-185
#setp kins-rotrad 1.993

fin = open(os.path.join(fdir, "probing5.txt"))
probepts = [ list(map(float, l.split()[:3]))  for l in fin ]
sendactivity(points=[(p[0],p[1],p[2]*50)  for p in probepts])
max(p[2]  for p in probepts)


def probez(x, y):
    return min(probepts, key=lambda X:abs(X[0]-x)+abs(X[1]-y))[2]


# do probing
fin = open(os.path.join(fdir, fn))
fout = open(os.path.join(fdir, fns), "w")
pos = { "X":xlo, "Y":ylo, "Z":zhi }
conts = [ [] ]
for l in fin:
    m = re.search("(([XYZ][\-\d\.]+)+)", l)
    if not m:
        fout.write(l)
        continue
    for c in "XYZ":
        mc = re.search("%s([\-\d\.]+)"%c, m.group(0))
        if mc:
            pos[c] = float(mc.group(1))
    if pos["Z"] == -0.1:
        ml = "X%.4fY%.4fZ%.4f" % (pos["X"], pos["Y"], probez(pos["X"], pos["Y"])+zlo)
        fout.write(l[:m.start(0)]+ml+l[m.end(0):])
        conts[-1].append((pos["X"], pos["Y"]))
    else:
        fout.write(l)
        if conts[-1]:
            conts.append([])
fout.close()
fin.close()
sendactivity(contours=conts)
sendactivity(points=[cont[0]  for cont in conts  if cont])
sendactivity(points=[cont[-1]  for cont in conts  if cont])
        


print("(PROBEOPEN probing5.txt)")

xlo, xhi = 10, 110
ylo, yhi = -50, 0
zlo, zhi = -0.3,0.6 

nsubsx, nsubsy = 22, 12 
print("G1Z%.3fF800" % zhi)
for i in range(0, nsubsx+1):
    lamx = i*1.0/nsubsx
    x = xlo * (1-lamx) + xhi * lamx
    for j in range(0, nsubsy+1):
        lamy = j*1.0/max(1,nsubsy)
        if (i%2)==1:
            lamy = 1-lamy
        y = ylo * (1-lamy) + yhi * lamy
        print("G1X%.3fY%.3fF800" % (x, y))
        print("G38.2Z%.3fF200" % zlo)
        print("G1Z%.3fF800" % zhi)


"""
z = 0
y = -8
for i in range(6, 18):
    print("G1Z%.3fF800" % (z-i*0.5))
    print("G1Y%.3fF800" % y)
    y = -25 if y == -8 else -8
"""
print("(PROBECLOSE)")
print("M2")




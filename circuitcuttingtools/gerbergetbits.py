import sys, os
sys.path.append("/home/goatchurch/geom3d/flatcam")
import camlib  # has to be Python27
import shapely

def ExtractContoursFromPoly(poly, dx, dy):
    conts = [ ]
    bdy = poly.boundary
    if type(bdy) == shapely.geometry.multilinestring.MultiLineString:
        for lin in list(bdy):
            conts.append([(p[0]+dx, p[1]+dy)  for p in list(lin.coords)])
    else:
        conts.append([(p[0]+dx, p[1]+dy)  for p in list(bdy.coords)])
    return conts

def GetIsolationCuts(fname, rad, dx, dy):
    g = camlib.Gerber()
    g.parse_file(fname)
    g.convert_units("MM")
    g.create_geometry()

    if type(g.solid_geometry) == shapely.geometry.polygon.Polygon:
        polys = [g.solid_geometry]
        polyoffsets = [g.solid_geometry.buffer(rad)]
    else:
        polys = g.solid_geometry.geoms
        polyoffsets = g.solid_geometry.buffer(rad)
    conts = sum((ExtractContoursFromPoly(poly, dx, dy)  for poly in polys), [])
    contoffsets = sum((ExtractContoursFromPoly(poly, dx, dy)  for poly in polyoffsets), [])
    return conts, contoffsets
    
def GetEdgeCuts(fname, rad, dx, dy):
    g = camlib.Gerber()
    g.parse_file(fname)
    g.convert_units("MM")
    g.create_geometry()

    poly = g.solid_geometry
    polyoffset = poly.buffer(rad)
    conts = ExtractContoursFromPoly(poly, dx, dy)
    contoffsets = ExtractContoursFromPoly(polyoffset, dx, dy)
    return conts, contoffsets


def GetHoleCuts(fname, dx, dy):
    ge = camlib.Excellon()
    ge.parse_file(fname)
    ge.convert_units("MM")
    holediampts = { }
    for p in ge.drills:
        diam = ge.tools[p["tool"]]["C"]
        if diam not in holediampts:
            holediampts[diam] = [ ]
        q = (p["point"].coords)[0]
        holediampts[diam].append((q[0]+dx, q[1]+dy))
    return holediampts


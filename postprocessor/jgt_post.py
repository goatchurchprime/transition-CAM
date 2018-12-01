# ***************************************************************************
# *   (c) goatchurch (shopinthewoods@gmail.com) 2018                        *
# *                                                                         *
# *   This file is part of the FreeCAD CAx development system.              *
# *                                                                         *
# *   This program is free software; you can redistribute it and/or modify  *
# *   it under the terms of the GNU Lesser General Public License (LGPL)    *
# *   as published by the Free Software Foundation; either version 2 of     *
# *   the License, or (at your option) any later version.                   *
# *   for detail see the LICENCE text file.                                 *
# *                                                                         *
# *   FreeCAD is distributed in the hope that it will be useful,            *
# *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *   GNU Lesser General Public License for more details.                   *
# *                                                                         *
# *   You should have received a copy of the GNU Library General Public     *
# *   License along with FreeCAD; if not, write to the Free Software        *
# *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  *
# *   USA                                                                   *
# *                                                                         *
# ***************************************************************************/
from __future__ import print_function
import FreeCAD
import Path
from PathScripts import PostUtils
from PathScripts import PathUtils

TOOLTIP = '''
Structured post-processor for development for KinetiC-NC CNC machine 
https://github.com/DoESLiverpool/somebody-should/wiki/CNC-Router
'''


#
# Functions to preprocess the incoming FreeCAD object into a 
# list of [ ({tooldef}, [motion cmds]) ].
#    (If these become standardized then they can be moved into a separate module or embedded into /Mod/Paths/PathScripts/PastPost.py.)

import re

# iterates down the groupd and concatenates the path objects to create a list of cmd objects
def flattenobjectlisttocommands(objectslist):
    cmds = [ ]
    robjectstack = list(reversed(objectslist))
    while robjectstack:
        pathobj = robjectstack.pop()
        if hasattr(pathobj, "Group"):
            robjectstack.extend(reversed(pathobj.Group))
        elif hasattr(pathobj, "Path"):
            cmds.extend(pathobj.Path.Commands)
    return cmds

# regroup sequences of pure motion in between tool change definition commands
motionsequencename = set(["G0", "G1", "G2", "G3"])
def findmotionsequenceindexes(cmds):
    imotionsequences = [ ]
    for i, cmd in enumerate(cmds):
        if cmd.Name in motionsequencename:
            if not imotionsequences or imotionsequences[-1][1] != i-1:
                imotionsequences.append([i, i])
            else:
                imotionsequences[-1][1] = i
    return imotionsequences

# merge a series of tool definition commands into a single object
def extracttooldef(cmds):
    tooldef = { }
    for cmd in cmds:
        if cmd.Name in ["M6", "M3"]:
            tooldef.update(cmd.Parameters)
        else:
            commentparams = re.findall("(Diameter): ([\d\.]+)", cmd.Name)
            tooldef.update(dict(commentparams))
    return tooldef

def extractfirstpos(cmds):
    pos = { }
    for cmd in cmds:
        for w in "XYZ":
            if w not in pos and w in cmd.Parameters:
                pos[w] = cmd.Parameters[w]
        if len(pos) == 3:
            break
    return pos

# returns [ ({tooldef}, [motion cmds]) ]
#     where tooldef = { i, T, S, Diameter, prevpos:{XYZ}, firstpos:{XYZ}, firstpos:{XYZ} }
def flattenandgroup(postlist):
    cmds = flattenobjectlisttocommands(postlist)
    imotionsequences = findmotionsequenceindexes(cmds)
    tooldefmotions = [ ]
    prevb = -1
    for i, (a, b) in enumerate(imotionsequences):
        motioncmds = cmds[a:b+1]
        tooldef = extracttooldef(cmds[prevb+1:a])
        tooldef["i"] = i
        tooldef["firstpos"] = extractfirstpos(motioncmds)
        tooldef["lastpos"] = extractfirstpos(reversed(motioncmds))
        if tooldefmotions:
            tooldef["prevpos"] = tooldefmotions[-1][0]["lastpos"]
        tooldefmotions.append((tooldef, motioncmds))
    return tooldefmotions
    

#
# Main three functions that make up the core of this post processor
#
import datetime
def writetooldefheader(fout, tooldef, currpos):
    if tooldef["i"] == 0:
        fout.write("%\n")
        fout.write("(Start %s)\n" % tooldef["filename"])
        fout.write("(Exported by FreeCAD)\n")
        fout.write("(Post Processor: %s)\n" % __name__)
        fout.write("(Output Time: %s)\n" % str(datetime.datetime.now()))
        fout.write("G90 G94\nG17\nG21\nG28 G91 Z0\nG90\n\n")
    fout.write("(Diameter: %s)\n" % tooldef.get("Diameter", "unknown"))
    fout.write("M9\n")
    fout.write("T%d M6\n" % tooldef.get("T", 0))
    fout.write("S%d M3\n" % tooldef.get("S", 0))
    fout.write("G54\nM8\n")
    currpos.update(tooldef["firstpos"])
    fout.write("G0 X%.3f Y%.3f\n" % (currpos["X"], currpos["Y"]))
    fout.write("G43 Z%.3f H1\n" % currpos["Z"])
    currpos["Name"] = "G0"

def writemotioncmds(fout, motioncmds, currpos):
    for cmd in motioncmds:
        fline = [ ]
        if cmd.Name != currpos.get("Name") or cmd.Name in ["G2", "G3"]:
            fline.append("%s " % cmd.Name)
            currpos["Name"] = cmd.Name
        for w in "XYZ":
            if w in cmd.Parameters and cmd.Parameters[w] != currpos.get(w):
                fline.append("%s%.3f " % (w, cmd.Parameters[w]))
                currpos[w] = cmd.Parameters[w]
        if cmd.Name in ["G2", "G3"]:
            fline.append("I%.3f J%.3f " % (cmd.Parameters["I"], cmd.Parameters["J"]))
        if fline:
            fout.write("%s\n" % "".join(fline))

def writefilefooter(fout, currpos):
    fout.write("M9\n")
    fout.write("G28 G91 Z0\n")
    fout.write("G90\n")
    fout.write("G28 G91 X0 Y0\n")
    fout.write("G90\n")
    fout.write("M30\n")
    fout.write("%")
    

#
# The entry point function that:
#    (1) Calls the regrouping code, 
#    (2) Calls the three core functions, and 
#    (3) saves the file after popping up an editor
#
from StringIO import StringIO
def export(objectslist, filename, argstring):
    print("postprocessing...")
    tooldefmotions = flattenandgroup(objectslist)
    if tooldefmotions:
        tooldefmotions[0][0]["filename"] = filename

    fout = StringIO()
    currpos = {}
    for tooldef, motioncmds in tooldefmotions:
        writetooldefheader(fout, tooldef, currpos)
        writemotioncmds(fout, motioncmds, currpos)
    writefilefooter(fout, currpos)

    gcode = fout.getvalue()
    if FreeCAD.GuiUp: 
        dia = PostUtils.GCodeEditorDialog()
        dia.editor.setText(fout.getvalue())
        result = dia.exec_()
        if result:
            gcode = dia.editor.toPlainText()
    print("done postprocessing.")

    if filename != '-':
        gfile = pythonopen(filename, "w")
        gfile.write(gcode)
        gfile.close()
        

print(__name__ + " gcode postprocessor loaded.")



#
# Useful functions to keep safe for when running posts outside of FreeCAD
#
"""import sys
freecadpath = "/home/julian/extrepositories/FreeCAD/freecad-build/lib"
sys.path.append(freecadpath)
import FreeCAD
import PathScripts.PathJob as PathJob
import PathScripts.PathLog as PathLog
import PathScripts.PathToolController as PathToolController
import PathScripts.PathUtil as PathUtil
"""

def extractpostlistfromfile(fname):
    doc = FreeCAD.open(fname)
    #jobs = [o  for o in FreeCAD.ActiveDocument.Objects  if hasattr(o, "Proxy") and isinstance(o.Proxy, PathJob.ObjectJob)]
    job = doc.getObject("Job")
    postlist = []
    currTool = None
    for obj in job.Operations.Group:
        PathLog.debug("obj: {}".format(obj.Name))
        tc = PathUtil.toolControllerForOp(obj)
        if tc is not None:
            if tc.ToolNumber != currTool:
                postlist.append(tc)
                currTool = tc.ToolNumber
        postlist.append(obj)
    return postlist

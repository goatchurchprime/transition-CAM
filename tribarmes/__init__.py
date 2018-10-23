from .trianglebarmesh import TriangleBarMesh
from .triangleboxing import MakeTriangleBoxing


bnumpyexists = True
try:  import numpy
except ImportError:  bnumpyexists = False

if bnumpyexists:
    from .ntrianglebarmesh import NTriangleBarMesh
else:
    NTriangleBarMesh = None
    
    

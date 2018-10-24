import math
from basicgeo import P3
import barmesh

# this module handles balls touching triangles along a vector and closest approach of the triangles to a point

class DistPZ:
    def __init__(self, p, r):
        self.p = p
        self.r = r
        self.v = None
        
    def DistPpointPZ(self, p0):
        lv = p0 - self.p
        if abs(lv.z) <= self.r:
            lvlen = lv.Len()
            if lvlen < self.r:
                self.r = lvlen
                self.v = lv
        """# this is slower!:
        vz = p0.z - self.p.z
        if abs(vz) > self.r:
            return
        vx = p0.x - self.p.x
        vy = p0.y - self.p.y
        lvsq = vx*vx + vy*vy + vz*vz
        if lvsq >= self.rsq:
            return
        lvlen = math.sqrt(lvsq)
        if lvlen < self.r:
            self.r = lvlen
            self.rsq = lvsq
            self.v = P3(vx, vy, vz)
        """
    
    def DistPedgePZ(self, p0, p1):
        lv = self.p - p0
        v = p1 - p0
        vsq = v.Lensq()
        lam = P3.Dot(v, lv) / vsq
        if 0.0 < lam < 1.0:
            vd = lv - v * lam
            assert abs(P3.Dot(vd, v)) < 0.0001
            vdlen = vd.Len()
            if vdlen < self.r:
                self.r = vdlen
                self.v = -vd

    def DistPtrianglePZ(self, p0, p1, p2):
        v1 = p1 - p0
        v2 = p2 - p0
        v1sq = v1.Lensq()
        v2sq = v2.Lensq()
        v1dv2 = P3.Dot(v1, v2)
        det = v1sq*v2sq - v1dv2**2
        if det == 0.0:
            return   # near zero width triangle has no face to touch
        invdet = 1.0 / det
        
        # should do distance by the dot on the crossproduct normal
        
        lv = self.p - p0
        v1dlv = P3.Dot(v1, lv)
        v2dlv = P3.Dot(v2, lv)
        
        # solve vd = lv - v1 * lam1 - v2 * lam2, where vd.v1 = vd.v2 = 0
        # (v1sq   v1dv2)   ( lam1 )   ( v1dlv )
        # (v1dv2   v2sq) . ( lam2 ) = ( v2dlv )
        lam1 = (v2sq * v1dlv - v1dv2 * v2dlv) * invdet
        lam2 = (-v1dv2 * v1dlv + v1sq * v2dlv) * invdet
        if 0 < lam1 and 0 < lam2 and lam1 + lam2 < 1:
            vd = lv - v1 * lam1 - v2 * lam2
            assert abs(P3.Dot(vd, v1)) < 0.001
            assert abs(P3.Dot(vd, v2)) < 0.001
            vdlen = vd.Len()
            if vdlen < self.r:
                self.r = vdlen
                self.v = -vd


class DistLamPZ:
    def __init__(self, p, vp, r):
        self.p = p
        self.vp = vp
        self.vpsq = vp.Lensq()
        self.r = r
        self.rsq = r*r
        self.lam = 2.0
        
    def DistLamPpointPZ(self, p0):
        lv = p0 - self.p
        if lv.z < min(self.vp.z, 0) - self.r or lv.z > max(self.vp.z, 0) + self.r:
            return
        # |lv - vp * lam| = r
        # qa = vpsq
        qb2 = -P3.Dot(self.vp, lv)
        qc = lv.Lensq() - self.rsq
        
        qdq = qb2*qb2 - self.vpsq*qc
        if qdq < 0.0:
            return
        qs = math.sqrt(qdq) / self.vpsq
        qm = -qb2 / self.vpsq
        assert abs(qc + (qm + qs)*(2*qb2 + (qm + qs)*self.vpsq)) < 0.002
        if qm + qs <= 0.0:
            return
        laml = qm - qs
        if laml < 0.0:
            self.lam = 0.0  # shouldn't happen
        elif laml < self.lam:
            self.lam = laml

    def DistLamPedgePZ(self, p0, p1):
        v = p1 - p0
        vsq = v.Lensq()
        
        lv = self.p - p0
        # solve |lv + vp * lam - v * mu| == r, where (lv + vp * lam - v * mu) . vp == 0
        mu0 = P3.Dot(lv, v)/vsq
        lvf = lv - v*mu0
        vpdv = P3.Dot(self.vp, v)
        muvp = vpdv/vsq
        
        vpf = self.vp - v*muvp
        vpfsq = vpf.Lensq()
        if vpfsq == 0.0:
            return
        assert abs(vpfsq - (self.vpsq - muvp * vpdv)) < 0.001
        
        lvfdvpf = P3.Dot(lvf, vpf)
        lamc = -lvfdvpf / vpfsq
        cp = lvf + vpf * lamc
        cpsq = cp.Lensq()
        lvfsq = lvf.Lensq()
        assert abs(cpsq - (lvfsq + 2 * lvfdvpf * lamc + vpfsq * lamc * lamc)) < 0.001
        assert abs(P3.Dot(cp, vpf)) < 0.001
        llamdsq = self.rsq - cp.Lensq()
        if llamdsq < 0.0:
            return
        lamd = math.sqrt(llamdsq / vpfsq)
        if lamc + lamd < 0.0:
            return
        lam = lamc - lamd
        if lam < 0.0:
            return  # check closer stuff
        if lam > self.lam:
            return
        mu = mu0 + muvp * lam
        if mu < 0 or mu > 1:
            return
        dv = lv + self.vp * lam - v * mu
        assert abs(dv.Len() - self.r) < 0.001
        assert abs(P3.Dot(dv, v)) < 0.001
        self.lam = lam
        
    def DistLamPtrianglePZ(self, p0, p1, p2):
        # solve vd = lv + vp * lam - v1 * lam1 - v2 * lam2, where |vd| = r and vd.v1 = vd.v2 = 0
        # solve +-r = (lv + vp * lam) . vnorm = lv . vnorm + vp . vnorm lam
        v1 = p1 - p0
        v2 = p2 - p0
        vcross = P3.Cross(v1, v2)
        assert abs(P3.Dot(vcross, v1)) < 0.001
        assert abs(P3.Dot(vcross, v2)) < 0.001
        vnorm = P3.ZNorm(vcross)
        
        lv = self.p - p0
        lvdvnorm = P3.Dot(lv, vnorm)
        vpdvnorm = P3.Dot(self.vp, vnorm)
        if vpdvnorm == 0.0:
            return
        # lam = (+-r - lvdvnorm)/vpdvnorm
        if vpdvnorm > 0.0:
            lam = (-self.r - lvdvnorm)/vpdvnorm
        else:
            lam = (self.r - lvdvnorm)/vpdvnorm
        if lam < 0 or lam > self.lam:
            return

        lvl = lv + self.vp * lam
        v1dlv = P3.Dot(v1, lvl)
        v2dlv = P3.Dot(v2, lvl)
        v1sq = v1.Lensq()
        v2sq = v2.Lensq()
        v1dv2 = P3.Dot(v1, v2)
        det = v1sq*v2sq - v1dv2**2
        if det == 0.0:
            return  # no face size
        invdet = 1.0 / det
        
        # solve vd = lv - v1 * lam1 - v2 * lam2, where vd.v1 = vd.v2 = 0
        # (v1sq   v1dv2)   ( lam1 )   ( v1dlv )
        # (v1dv2   v2sq) . ( lam2 ) = ( v2dlv )
        lam1 = (v2sq * v1dlv - v1dv2 * v2dlv) * invdet
        lam2 = (-v1dv2 * v1dlv + v1sq * v2dlv) * invdet
        if 0 < lam1 and 0 < lam2 and lam1 + lam2 < 1:
            vd = lvl - v1 * lam1 - v2 * lam2
            assert abs(P3.Dot(vd, v1)) < 0.001
            assert abs(P3.Dot(vd, v2)) < 0.001
            assert abs(vd.Len() - self.r) < 0.001
            self.lam = lam
            
    
class ImplicitAreaBallOffset:
    def __init__(self, tboxing):
        self.tbarmesh = tboxing.tbarmesh
        self.tboxing = tboxing
        self.hitreg = [0]*len(self.tbarmesh.bars)
        self.nhitreg = 0

    def Isb2dcontournormals(self):
        return False
    
    def DistP(self, pz, p):    # pz=PointZone
        dpz = DistPZ(p, pz.r)  # temporarily convert to a different implementation of pointzone with some functions on it (could make the functions more global)
        
        for ix, iy in self.tboxing.CloseBoxeGenerator(p.x, p.x, p.y, p.y, pz.r):
            tbox = self.tboxing.boxes[ix][iy]
            for i in self.tboxing.SlicePointisZ(tbox.pointis, p.z-pz.r, p.z+pz.r):
            #for i in tbox.pointis:
                dpz.DistPpointPZ(self.tboxing.GetNodePoint(i))
                
        self.nhitreg += 1
        for ix, iy in self.tboxing.CloseBoxeGenerator(p.x, p.x, p.y, p.y, pz.r):
            tbox = self.tboxing.boxes[ix][iy]
            for i in tbox.edgeis:
                if self.hitreg[i] != self.nhitreg:
                    dpz.DistPedgePZ(*self.tboxing.GetBarPoints(i))
                    self.hitreg[i] = self.nhitreg
                    
        # then do the triangles
        self.nhitreg += 1
        for ix, iy in self.tboxing.CloseBoxeGenerator(p.x, p.x, p.y, p.y, pz.r):
            tbox = self.tboxing.boxes[ix][iy]
            for i in tbox.triangleis:
                if self.hitreg[i] != self.nhitreg:
                    dpz.DistPtrianglePZ(*self.tboxing.GetTriPoints(i))
                    self.hitreg[i] = self.nhitreg
                    
        pz.r = dpz.r
        pz.v = dpz.v
        

    def Cutpos(self, p, vp, cp, r):  # point, vector, cp=known close point to narrow down the cutoff search
        dlpz = DistLamPZ(p, vp, r)  
        if cp is not None:
            dlpz.DistLamPpointPZ(cp)
            assert dlpz.lam != 2.0
        
        # solve |p0 + vp * lam - p| == r where Dist(p0) >= r >= Dist(p0 + vp)
        rexp = r + 0.01
        for ix, iy in self.tboxing.CloseBoxeGenerator(min(p.x, p.x+vp.x), max(p.x, p.x+vp.x), min(p.y, p.y+vp.y), max(p.y, p.y+vp.y), rexp):
            tbox = self.tboxing.boxes[ix][iy]
            for i in self.tboxing.SlicePointisZ(tbox.pointis, min(p.z,p.z+vp.z)-rexp, max(p.z,p.z+vp.z)+rexp):
            #for i in tbox.pointis:
                dlpz.DistLamPpointPZ(self.tboxing.GetNodePoint(i))
                
        self.nhitreg += 1
        for ix, iy in self.tboxing.CloseBoxeGenerator(min(p.x, p.x+vp.x*dlpz.lam), max(p.x, p.x+vp.x*dlpz.lam), min(p.y, p.y+vp.y*dlpz.lam), max(p.y, p.y+vp.y*dlpz.lam), r + 0.01):
            tbox = self.tboxing.boxes[ix][iy]
            for i in tbox.edgeis:
                if self.hitreg[i] != self.nhitreg:
                    dlpz.DistLamPedgePZ(*self.tboxing.GetBarPoints(i))
                    self.hitreg[i] = self.nhitreg
        
        self.nhitreg += 1
        for ix, iy in self.tboxing.CloseBoxeGenerator(min(p.x, p.x+vp.x*dlpz.lam), max(p.x, p.x+vp.x*dlpz.lam), min(p.y, p.y+vp.y*dlpz.lam), max(p.y, p.y+vp.y*dlpz.lam), r + 0.01):
            tbox = self.tboxing.boxes[ix][iy]
            for i in tbox.triangleis:
                if self.hitreg[i] != self.nhitreg:
                    dlpz.DistLamPtrianglePZ(*self.tboxing.GetTriPoints(i))
                    self.hitreg[i] = self.nhitreg
        
        if __debug__:
            if not (dlpz.lam == 0.0 or dlpz.lam == 2.0):
                pz = barmesh.PointZone(0, r + 1.1, None)
                self.DistP(pz, dlpz.p + dlpz.vp*dlpz.lam)
                assert abs(pz.r - r) < 0.002, ("cutposbad", pz.r, dlpz.lam, dlpz.p, dlpz.vp)
            
        return dlpz.lam



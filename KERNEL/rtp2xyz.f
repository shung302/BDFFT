c convert a point at spherical coordinate to cartesian
c
      subroutine rtp2xyz(r,t,p,x,y,z)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
c      real*8 r,t,p,x,y,z,st,ct,sp,cp

      if (t.ne.pi2) go to 10
      st=1.0d0
      ct=0.0d0
10    st=dsin(t)
      ct=dcos(t)
      sp=dsin(p)
      cp=dcos(p)
      x=r*st*cp
      y=r*st*sp
      z=r*ct

      return
      end

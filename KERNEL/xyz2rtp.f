c convert a point at cartesian to spherical coordinate
c
      subroutine xyz2rtp(x,y,z,r,t,p)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

      r=dsqrt(x*x+y*y+z*z)
      t=dacos(z/r)
      if (x.eq.0.0d0) then
         p=pi2
         return
      endif
      p=datan2(y,x)

      return
      end

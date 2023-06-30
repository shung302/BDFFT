c transform the path of the ray in (r,theta,phi) to (x,y,z)
c
      subroutine ray2xyz(ni)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      COMMON /xyzpath/ xr(NSTPMX),yr(NSTPMX),zr(NSTPMX)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      INCLUDE 'paths.inc'

      t=pi2
      st=1.0d0
      ct=0.0d0
      do i=1,ni
         r=yp(1,i)
         p=yp(3,i)
         sp=dsin(p)
         cp=dcos(p)
         xr(i)=r*cp
         yr(i)=r*sp
         zr(i)=0.0d0
c         write(31,31)r,p*dpr,t*dpr,xr(i),yr(i),zr(i)
31    format(6f12.4)
      enddo

      return
      end


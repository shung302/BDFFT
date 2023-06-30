c subroutine to do rotation using euler angles
c input: vector (x,y,z), cosine of three euler angles 
c
      subroutine euler_rot(rs,ts,ps,rr,tr,pr,g2e,e2g)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      dimension g2e(3,3),e2g(3,3)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

      call cart(ts,ps,xi,yi,zi)
      call cart(tr,pr,xr,yr,zr)
c      call rtp2xyz(rs,ts,ps,xi,yi,zi)
c      call rtp2xyz(rr.0d0,tr,pr,xr,yr,zr)
c      write(6,*)'source and receiver in x,y,z'
c      write(6,*)xi,yi,zi
c      write(6,*)xr,yr,zr
c pole (xk,yk,zk) = (xi,yi,zi) x (xr,yr,zr) of great circle plane between s and r
c  
      xk=yi*zr-yr*zi
      yk=xr*zi-xi*zr
      zk=xi*yr-xr*yi
      call cross(xi,yi,zi,xr,yr,zr,xk,yk,zk)

      dk=dsqrt(xk*xk+yk*yk+zk*zk)
cc       write(6,*)'dk=',dk
      xk=xk/dk
      yk=yk/dk
      zk=zk/dk
c (xj,yj,zj) = (xk,yk,zk) x (xi,yi,zi)
      xj=yk*zi-yi*zk
      yj=xi*zk-xk*zi
      zj=xk*yi-xi*yk
      dj=dsqrt(xj*xj+yj*yj+zj*zj)
      xj=xj/dj
      yj=yj/dj
      zj=zj/dj

c rotation matrix cos(x',x)
      g2e(1,1)=xi
      g2e(2,1)=xj
      g2e(3,1)=xk
      g2e(1,2)=yi
      g2e(2,2)=yj
      g2e(3,2)=yk
      g2e(1,3)=zi
      g2e(2,3)=zj
      g2e(3,3)=zk

c      write(6,*)'a matrix'
c      do i=1,3
c         write(6,*)(g2e(i,j),j=1,3)
c      enddo

      e2g(1,1)=xi
      e2g(2,1)=yi
      e2g(3,1)=zi
      e2g(1,2)=xj
      e2g(2,2)=yj
      e2g(3,2)=zj
      e2g(1,3)=xk
      e2g(2,3)=yk
      e2g(3,3)=zk


c      xse=g2e(1,1)*xi+g2e(1,2)*yi+g2e(1,3)*zi
c      yse=g2e(2,1)*xi+g2e(2,2)*yi+g2e(2,3)*zi
c      zse=g2e(3,1)*xi+g2e(3,2)*yi+g2e(3,3)*zi
c      xre=g2e(1,1)*xr+g2e(1,2)*yr+g2e(1,3)*zr
c      yre=g2e(2,1)*xr+g2e(2,2)*yr+g2e(2,3)*zr
c      zre=g2e(3,1)*xr+g2e(3,2)*yr+g2e(3,3)*zr
c      write(6,*)'rotating'
c      write(6,*)xse,yse,zse,xre,yre,zre
c      xso=e2g(1,1)*xse+e2g(1,2)*yse+e2g(1,3)*zse
c      yso=e2g(2,1)*xse+e2g(2,2)*yse+e2g(2,3)*zse
c      zso=e2g(3,1)*xse+e2g(3,2)*yse+e2g(3,3)*zse
c      xro=e2g(1,1)*xre+e2g(1,2)*yre+e2g(1,3)*zre
c      yro=e2g(2,1)*xre+e2g(2,2)*yre+e2g(2,3)*zre
c      zro=e2g(3,1)*xre+e2g(3,2)*yre+e2g(3,3)*zre

c      delta=dcos(ts)*dcos(tr)+dsin(ts)*dsin(tr)*dcos(pr-ps)
c      xre1=delta
c      delta=dacos(delta)
c      write(6,*)'delta=',delta*dpr
c      yre1=dsin(delta)
c      write(6,*)xre1,yre1
c      write(6,*)'source=',xi,yi,zi,xso,yso,zso
c      write(6,*)'receiver=',xr,yr,zr,xro,yro,zro

      return
      end


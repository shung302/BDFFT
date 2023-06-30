c______________________________________________________________________
      subroutine tran(on,at,onp,atp,x,y,z,id)
c     (on,at) is lat,lon of Geographic Earth
c     (onp,atp) is lat,lon where center of the region is (0,0)
c     id = 1 --> (on,at) to (onp,atp)
c        =-1 --> (onp,atp) to (on,at)
      double precision trans(3,3)
      common/meshb2/trans
      double precision pi,rad
      parameter(pi=3.14159265d0,rad=180.d0/pi)
      double precision on,at,onp,atp,x,y,z,xp,yp,zp,ontmp,attmp

      if(id.eq.1)then
       ontmp=on/rad
       attmp=at/rad
       x=dcos(attmp)*dcos(ontmp)
       y=dcos(attmp)*dsin(ontmp)
       z=dsin(attmp)
       xp=trans(1,1)*x+trans(1,2)*y+trans(1,3)*z
       yp=trans(2,1)*x+trans(2,2)*y+trans(2,3)*z
       zp=trans(3,1)*x+trans(3,2)*y+trans(3,3)*z
       r=dsqrt(xp*xp+yp*yp+zp*zp)
       atp=dasin(zp/r)*rad
       onp=datan2(yp,xp)*rad
c       if(onp.lt.0.d0)onp=onp+360.d0
       if(onp.gt.180.d0)onp=onp-360.d0
       x=xp/r
       y=yp/r
       z=zp/r
      elseif(id.eq.-1)then
       ontmp=onp/rad
       attmp=atp/rad
       xp=dcos(attmp)*dcos(ontmp)
       yp=dcos(attmp)*dsin(ontmp)
       zp=dsin(attmp)
       x=trans(1,1)*xp+trans(2,1)*yp+trans(3,1)*zp
       y=trans(1,2)*xp+trans(2,2)*yp+trans(3,2)*zp
       z=trans(1,3)*xp+trans(2,3)*yp+trans(3,3)*zp
       r=dsqrt(x*x+y*y+z*z)
       at=dasin(z/r)*rad
       on=datan2(y,x)*rad
       if(on.lt.0.d0)on=on+360.d0
       x=x/r
       y=y/r
       z=z/r
      else
       stop' not permitted in tran !'
      endif
      return
      end


c______________________________________________________________________
      subroutine linit
c ... initialize for conjugate gradient least squares
      parameter(mlx=7,mly=7,mlz=7,
     ,          mnx=2**(mlx-1)+1,mny=2**(mly-1)+1,mnz=2**(mlz-1)+1,
     ,          mvolu=mnx*mny*mnz,mdata=60000,mpath=240)
c      include 'lsqrcom1'
      n = 0
      m = 0
      nel=0
      return
      end


      SUBROUTINE PIKSRT(N,ARR)
      double precision ARR(N)
      DO 12 J=2,N
        A=ARR(J)
        DO 11 I=J-1,1,-1
          IF(ARR(I).LE.A)GO TO 10
          ARR(I+1)=ARR(I)
11      CONTINUE
        I=0
10      ARR(I+1)=A
12    CONTINUE
      RETURN
      END


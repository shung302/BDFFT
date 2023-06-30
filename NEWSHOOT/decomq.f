c subroutine to decompose q vector to q1 and q2
c
      subroutine decomq(nq,ie,xe,ye,ze,q,q1,q2)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      DIMENSION xe(3),ye(3),ze(3),q(3,3),q1(3),q2(3),ie(3)
      DIMENSIOn pk(3),sv(3),sh(3)
      INCLUDE 'paths.inc'
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

      ct=0.0d0
      st=1.0d0
      do iq=1,nq
         i=ie(iq)
         ra=yp(1,i)
         ai=yp(2,i)
         ph=yp(3,i) 
         cp=dcos(ph)
         sp=dsin(ph)
         ci=dcos(ai)
         si=dsin(ai)
         pk(1)=st*cp*ci-sp*si
         pk(2)=st*sp*ci+cp*si
         pk(3)=ct*ci         
         sv(1)=st*cp*si+sp*ci
         sv(2)=st*sp*si-cp*ci
         sv(3)=ct*si
         sh(1)=ct*cp
         sh(2)=ct*sp
         sh(3)=-st
c         write(6,*)'q=',(q(ii,iq),ii=1,3)
c         write(6,*)'pv=',(pk(ii),ii=1,3),(sv(ii),ii=1,3)
         call projectv(q(1,iq),sv,q1(iq))
         call projectv(q(1,iq),sh,q2(iq))
      enddo
      return
      end

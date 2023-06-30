c find wave vector for a given incident angle
c asuume ray at equatorial plane
      subroutine coordl(ai,ph,rotm)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      DIMENSION rotm(3,3)

      ct=0.0d0
      st=1.0d0
      cp=dcos(ph)
      sp=dsin(ph)
      ci=dcos(ai)
      si=dsin(ai)
      p1=st*cp*ci-sp*si
      p2=st*sp*ci+cp*si
      p3=ct*ci
      sv1=st*cp*si+sp*ci
      sv2=st*sp*si-cp*ci
      sv3=ct*si
c      write(6,*)'p and sv=',p1,p2,p3,sv1,sv2,sv3
      call cross(p1,p2,p3,sv1,sv2,sv3,sh1,sh2,sh3)
c rotation matrix
c  x'=rotm * x; where x' is the coordinates in the x original frame
c  rotm= | p1 sv1 sh1 |
c        | p2 sv2 sh2 |
c        | p3 sv3 sh3 |
c
      rotm(1,1)=p1
      rotm(2,1)=p2
      rotm(3,1)=p3
      rotm(1,2)=sv1
      rotm(2,2)=sv2
      rotm(3,2)=sv3
      rotm(3,1)=sh1
      rotm(3,2)=sh2
      rotm(3,3)=sh3

      return
      end


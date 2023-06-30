c find wave vector for a given incident angle
c asuume ray at equatorial plane
      subroutine pkv(ai,ph,pk)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      DIMENSION pk(3)

      ct=0.0d0
      st=1.0d0
      cp=dcos(ph)
      sp=dsin(ph)
      ci=dcos(ai)
      si=dsin(ai)
      pk(1)=st*cp*ci-sp*si
      pk(2)=st*sp*ci+cp*si
      pk(3)=ct*ci
      return
      end


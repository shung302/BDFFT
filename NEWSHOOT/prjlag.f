c	
      subroutine prjlag(radistrs,rs,ts,ps,nlags,lagb,lage,lag)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
c      INTEGER nrbp,irbp(10)
      INTEGER nlags,lagb(10),lage(10)
      INCLUDE 'paths.inc'

      if (ps.le.0.0d0) then
         lag=1
         return
      endif
      if (ps.ge.radistrs) then
         lag=nlags
         return
      endif
      do i=1,nlags
         phi1=yp(3,lagb(i))
         phi2=yp(3,lage(i))
         if (ps.ge.phi1.and.ps.le.phi2) then
               lag=i
               return
         endif
      enddo

      return
      end

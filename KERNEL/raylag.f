c find index for ray bounce points
      subroutine raylag(iraytype,ni,nlags,lagb,lage,rayp,ramp,thmp,phmp)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INTEGER nlags,lagb(10),lage(10)
      INCLUDE 'paths.inc'
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

      thmp=pi2
      if (iraytype.eq.1.or.iraytype.eq.21) then
         nlags=1
         lagb(1)=1
         lage(1)=ni
         do i=1,ni
            if ( (yp(2,i)-pi2)*(yp(2,i+1)-pi2).lt.0.0d0) then
               cb=(vp(i)+vp(i+1))*0.5d0
               ramp=cb*rayp
               phmp=yp(3,i)+dtan(yp(2,i))*(ramp/yp(1,i)-1.0d0)
               return
            endif
         enddo
      endif

      nlags=0
      i=1
      lagb(1)=1
      do 100 while ( i .lt. ni ) 
         i=i+1
         if (yp(1,i).eq.rsurf.or.yp(1,i).eq.rcmb) then
            nlags=nlags+1
            lage(nlags)=i-1
            lagb(nlags+1)=i
            ramp=yp(1,i)
            phmp=yp(3,i)
         endif
100   enddo
      nlags=nlags+1
      lage(nlags)=ni

      return
      end

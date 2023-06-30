c find index for ray bounce points
      subroutine raylag(iraytype,ni,nlags,lagb,lage,rayp,ramp,thmp,phmp)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INTEGER nlags,lagb(10),lage(10)
      INCLUDE 'paths.inc'
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)

      thmp=pi2
      if (iraytype.eq.1.or.iraytype.eq.21) then
         nlags=1
         lagb(1)=1
         lage(1)=ni
         do i=1,ni
            if ( (yp(2,i)-pi2)*(yp(2,i+1)-pi2).lt.0.0d0) then
c               cb=(vp(i)+vp(i+1))*0.5d0
c             ramp=yp(1,i)+(pi2-yp(2,i))
c     &           *(dvp(i)/vp(i)-1./yp(1,i))/dtan(yp(2,i))
c fixed by Zhang, 2012/11/19, Tian et al (2007)
             ramp=yp(1,i)+(pi2-yp(2,i))
     &           /(dvp(i)/vp(i)-1./yp(1,i))/dtan(yp(2,i))
c             ramp=cb*rayp
             if (ramp.lt.rcmb) ramp=vdis(2,2,iwp(i))*rayp
             phmp=yp(3,i)+dtan(yp(2,i))*(ramp/yp(1,i)-1.0d0)
             return
            endif
         enddo
      endif

c for sks or pkp phases
      if ( (iraytype.ge.16.and.iraytype.le.20) .or.
     &     (iraytype.eq.36.and.iraytype.le.40) ) then
         nlags=1
         lagb(1)=1
         lage(1)=ni
         do i=1,ni
            if ( (yp(2,i)-pi2)*(yp(2,i+1)-pi2).lt.0.0d0) then
c               cb=(vp(i)+vp(i+1))*0.5d0
c             ramp=yp(1,i)+(pi2-yp(2,i))
c     &           *(dvp(i)/vp(i)-1./yp(1,i))/dtan(yp(2,i))
             ramp=yp(1,i)+(pi2-yp(2,i))
     &           /(dvp(i)/vp(i)-1./yp(1,i))/dtan(yp(2,i))
c             ramp=cb*rayp
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

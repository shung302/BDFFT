c subroutine to calcuate instrument response for a given T* at
c 1-pole seismometer 
c
      subroutine calc_stf1p(tstar,npts,wb,we,wx,ww,tper,amp2,sdown)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      parameter (nfmx=500)
      dimension wx(*),ww(*),amp2(*),sdown(*)
      dimension resp(nfmx)
      COMPLEX*16 pl,plc,cmu0
      character frq*80
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /poles/ pl,plc,cmu0,omega0,alpha0

c c generate tabulated abscissas and weights for gaussian-legendre
c integration between wb and we
      flow=twopi*wb
      fhigh=twopi*we
      call dgauleg(flow,fhigh,wx,ww,npts)
      
      omega0=twopi/tper
      alpha0=pi/tper
      omega2=omega0*omega0
      alpha2=alpha0*alpha0
      pl=dcmplx(omega0,alpha0)
      plc=dcmplx(omega0,-alpha0)
      do i=1,npts
         w2=wx(i)*wx(i)
         down=(w2-omega2-alpha2)**2+4.0d0*w2*alpha2
         amp2(i)=dexp(-wx(i)*tstar)/down
         sdown(i)=(w2+omega2+alpha2)**2-4.0d0*w2*alpha2
         write(2,'(4e15.6)')wx(i)/twopi,sqrt(amp2(i)),sdown(i)
      enddo

      if (int(tstar).lt.10) then
      write(frq,'(''resp1p.ts'',i1,''.sac'')')int(tstar)
      else
      write(frq,'(''resp1p.ts'',i2,''.sac'')')int(tstar)
      endif
      df=(we-wb)/dble(nfmx)
      do i=1,nfmx
         w=twopi*(wb+df*dble(i-1))
         w2=w*w
         resp(i)=dexp(-w*tstar)/((w2-omega2-alpha2)**2+4.0d0*w2*alpha2)
      enddo
c      call wsac1(frq,resp,nfmx,wb,df,ierr)

      return
      end


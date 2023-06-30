c calculate spectram at regularated frequencies for a 
c gauss-type source time function 

      subroutine calc_stfgs(alpha,npts,wb,we,wx,ww,amp2)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      dimension wx(*),ww(*),amp2(*)
      real*4 tt,tdiff,expo,respi(1000),tdelay,dt,tp
      character*80 frq
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

      open(51,file='stfw.out')
      dt=0.5/real(we)
      tp=sqrt(real(twopi*twopi/alpha))
      tdelay=tp*0.5
      write(6,*)dt,we

      flow=twopi*wb
      fhigh=twopi*we
      call dgauleg(flow,fhigh,wx,ww,npts)

      do i=1,npts
        w2=wx(i)*wx(i)
        w4=w2*w2
        amp2(i)=w4*dexp(-w2/alpha*0.50d0)
        write(51,111)wx(i)/twopi,wx(i),amp2(i)
      enddo
      close(51)

c     output stf waveforms

      if (int(tp).lt.10) then
         write(frq,'(''respG.ts'',i1,''.sac'')')int(tp)
      else
         write(frq,'(''respG.ts'',i2,''.sac'')')int(tp)
      endif
      
      do i=1,npts
         tt=dt*(i-1)
         tdiff=tt-tdelay
         expo=exp(-real(alpha)*tdiff**2)
         respi(i)=-2.0*real(alpha)*tdiff*expo
      enddo
c      call wsac1(frq,respi,1000,0.,dt,ierr)
111   format(3e15.6)
      return
      end


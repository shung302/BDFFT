c subroutine to read RESP.STN.NET.CHN (IRIS station response files)
c to get instrument spectrum
c
      subroutine calc_stfra(tstar,npts,wb,we,wx,ww,amp2)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      parameter (nfmx=50000,nfmx2=nfmx*2)
c      include 'parameters.inc'
      dimension wx(*),ww(*),amp2(*)
      dimension freq(nfmx),resp(nfmx2)
      character file*80,datime*20,fresp*80,fq*80,frq*80,fsac*80
      character*5  unts,sta,cha,net,locid,rtyp
      character*10 vbs
      integer evresp
      integer npts, start_stage, stop_stage, stdio_flag
      real*4 respi(1000),q(1000),f,df,dt2,dt,tstr,rr,ri,
     &   omega,wsave(2015)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

      open(2,file='resp.in') 
      read(2,*)
      read(2,'(a5)')sta
      read(2,'(a5)')cha
      read(2,'(a5)')net
      read(2,'(a20)')datime
c iuntflag=0 for velocity (defaults) and =1 for for acceleration, =2 for displacement'
      read(2,*)iuntflag
      close(2)
      locid='*'
      unts='VEL'
c set start_stage to -1, stop_stage to 0 (this will cause all of the
c response stages to be included in the response that is calculated)

      start_stage = -1
      stop_stage = 0

c set the stdio_flag to zero (the response is returned to the calling routine

      stdio_flag = 0
c if above line commented out and following line not, must
c define the SEEDRESP environment variable
      file = ''
      vbs = '-v'
      rtyp = 'CS'

c      write(6,*)'open resp.in'
      twopi=2.0d0*pi
      flow=twopi*wb
      fhigh=twopi*we
      call dgauleg(flow,fhigh,wx,ww,npts)
      do i=1,npts
         freq(i)=wx(i)/twopi
         resp(i)=0.0d0
      enddo
c
c now call evresp, assume resp() will contain the multiplexed output
c

       iflag = evresp(sta,cha,net,locid,datime,unts,file,freq,npts,resp,
     .                rtyp,vbs,start_stage,stop_stage,stdio_flag)
       if (iflag.ne.0) then
          write(*,31) 'ERROR processing: ',sta,cha,datime
 31       format(a,1x,a5,1x,a5,1x,a20)
          return
       endif
c
c  write out file of frequency, amplitude, phase
c    
      open(8,file='evtest.out')
      j = 1
      if (iuntflag-1) 3,4,5
 3    do 30 i=1,npts
        amp2(i) = resp(j)**2 + resp(j+1)**2
        amp = dsqrt(amp2(i))
        pha = datan2(resp(j+1),resp(j)) * 180. / pi
        write(8,'(5e15.6)') freq(i),amp,pha,resp(j),resp(j+1)
        j = j + 2
 30   continue
      go to 10
c if acceleration, divided by (-iw)
 4    do 40 i=1,npts
        rimag=resp(j)/wx(i)
        rreal=-resp(j+1)/wx(i)
        amp2(i) = rreal**2 + rimag**2
        amp = dsqrt(amp2(i))
        pha = datan2(rimag,rreal) * 180. / pi
        write(8,'(5e15.6)') freq(i),amp,pha,rreal,rimag
        resp(j)=rreal
        resp(j+1)=rimag
        j = j + 2
 40   continue
      go to 10
c if displacement,multiply by (iw)
 5    do 50 i=1,npts
        rimag=resp(j)*wx(i)
        rreal=-resp(j+1)*wx(i)
        amp2(i) = rreal**2 + rimag**2
        amp = dsqrt(amp2(i))
        pha = datan2(rimag,rreal) * 180. / pi
        write(8,'(5e15.6)') freq(i),amp,pha,rreal,rimag
        resp(j)=rreal
        resp(j+1)=rimag
        j = j + 2
 50   continue
 10   close(8)

c  interpolate spectrum to get a regular-spaced array for fft
c  time series in 120 s long
      dt2=0.2
      nn=600
      nn2=nn/2
      df=1./(nn*dt2)
      respi(1)=0.
c      respi(2)=0.
c      j=3
      j=2
      do 110 i=1,nn2
         f=i*df
         omega=f*twopi
         call intresp(f,npts,freq,resp,rr,ri)
         respi(j)=rr
         respi(j+1)=ri
110      j=j+2
      respi(j+1)=0.
      call rffti(nn,wsave)
      call rfftb(nn,respi,wsave)

c      call fftl(respi,nn,-2,ierr)

      write(fresp,'(''resp.'',a4)')sta
      write(fsac,'(''resp.'',a4,''.sac'')')sta
      open(11,file=fresp)
c      write(11,'(f8.1)')dt2
      do i=1,nn
         ti=dt2*(i-1)
         write(11,'(f8.1,e13.4)')ti,respi(i)
      enddo
      close(11)

c      call wsac1(fsac,respi,nn,0.,.2,ierr)

      j=0
      do 120 i=1,nn,5
         j=j+1
120      respi(j)=respi(i)
       nn=nn/5
      
c      call wsac1('resp.sac',respi,nn,0.,1.,ierr)

c multiply t* operator
      dt=1.
      fct=1.0d-10
      if (tstar.ne.0.0d0) then
         tstr=tstar
c         write(6,*)'tstar=',tstr
         call qualfc(tstr,dt,q,nq)
c         write(6,*)'tstar and nq=',tstr,nq
         call conv(respi,q,nn,nq)
c  convolve instrument response with t* operator
         do i=1,npts
            amp2(i)=(dexp(-wx(i)*tstr)*amp2(i)*wx(i)*wx(i))*fct
         enddo
      endif

      if (int(tstar).lt.10) then
      write(fq,'(''Q.ts'',i1,''.sac'')')int(tstar)
      write(frq,'(''respQ.ts'',i1,''.sac'')')int(tstar)
      else
      write(fq,'(''Q.ts'',i2,''.sac'')')int(tstar)
      write(frq,'(''respQ.ts'',i2,''.sac'')')int(tstar)
      endif

c      call wsac1(fq,q,nq,0.,dt,ierr)
c      call wsac1(frq,respi,nn,0.,dt,ierr)

      return
      end
c interpolate intrument response for a given w
        subroutine intresp(f,npts,tfreq,resp,rr,ri)
        real*8 resp(*),tfreq(*)
C for linear interpolation
        real*8 tmp0
        real rr,ri,f
        tmp=dble(f)
        if(tmp.le.tfreq(1)) then
          rr=resp(1)
          ri=resp(2)
          return
        endif
        j=3
        do k = 2, npts
        if(tmp.lt.tfreq(k)) then
        tmp0=(tmp-tfreq(k-1))/(tfreq(k)-tfreq(k-1))
        rr=tmp0*(resp(j)-resp(j-2))+resp(j-2)
        ri=tmp0*(resp(j+1)-resp(j-1))+resp(j-1)
        return
        endif
        j=j+2
        end do
c       write(6,*) 'Freq=',tmp,' > Nyquist Freq'
        rr=resp(npts*2-1)
        ri=resp(npts*2)

        return
        end


      subroutine qualfc(tstr,dt,q,n)
c*** compute t* operator based on absorption band model (see e.g. doornbos
c*** 1983).normalisation is such that q tends to a delta funcion as t* tends
c*** to zero(which is not the same as a fourier integral).slight ringing is
c*** possible for small values of t* but it is assumed that this will be
c*** wiped out by the instrument response.note that the q achieves about
c*** 1/1000 of its peak value in a time of about 30t* secs.
c*** guy's fft code seems not working for t* >=4
c*** replaced by netlib fft subroutines, SH Hung

      dimension q(*),wsave(2015)
      data pi/3.1415927/
      p2=pi*2
c*** tau=.3 secs seems appropriate
      tau=1./pi
      nfilt=30.*tstr/dt
      if (nfilt.gt.1800) nfilt=1800
      n=(int(nfilt/100.)+1)*100.
c get an even number acceptable to FFTL
c      call factor(n)
      nn=n/2
      fnn=nn
      df=p2/(n*dt)
      q(1)=1.
c      q(2)=0.
      nshft=tstr/dt
      nshft=2*nshft+2
      tshft=nshft*dt
c     kk=2 for fftl subroutines
      kk=1
      do 10 i=1,nn
      w=i*df
      wt=w*tau
      amp=exp(-w*tstr/2.)
      phas=-w*tstr*alog(1.+1./(wt*wt))/p2-w*tshft
      kk=kk+1
      q(kk)=amp*cos(phas)
      kk=kk+1
   10 q(kk)=amp*sin(phas)
c      nold=kk-1
c      q(nold)=sqrt(q(nold)*q(nold)+q(kk)*q(kk))
      q(kk+1)=0.
      call rffti(n,wsave)
      call rfftb(n,q,wsave)
c      call fftl(q,n,-2,ierr)
c      if (ierr.ne.0) print *,'******ierr in fftl = ',ierr
      i1=nshft+1
      do 20 i=i1,n
   20 q(i-nshft)=q(i)/fnn
      n=n-nshft
c
c      qmax=0
c      do 40 i=1,n
c         if (q(i).gt.qmax) qmax=q(i)
c         if (q(i).lt.qmax/100.) go to 50
c40    continue
c      print *,'***WARNING: q(n).gt.qmax/100'
c      return
c50    n=i
      return
      end

c---------------------------------------------------------------
      subroutine conv(x,b,lx,lb)
c   this subroutine convolves x and b.  lx and lb are the
c   lengths of x and b respectively.  ly is the length of the
c   resultant time series and  should equal lx+lb-1.  the first
c   lx points of the result are overwritten onto x.
      dimension x(1),b(1)
      real*8 y
      dimension y(8000)
      ly=lx+lb-1
      do 10 i=1,ly
   10 y(i)=0.d0
      do 20 i=1,lx
      do 20 j=1,lb
   20 y(i+j-1)=y(i+j-1)+x(i)*b(j)
      do 30 i=1,lx
   30 x(i)=y(i)
      return
      end
c---------------------------------------------------------------------


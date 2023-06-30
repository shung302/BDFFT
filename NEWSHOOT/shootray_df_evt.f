c use two-point boundary-value problem solver (Bisection shooting method)
c to do two-point ray tracing
c ======================================================================
c three first-order equations are solved simultaneously for (r,i,phi) that 
c integrate w.r.t. l (l=arc length until phi=Delta) where Delta=epicentral distance
c d(r)/d(l)=cos(i)
c d(i)/d(l)=p(dv/dr-v/r)/r
c d(phi)/d(l)=sin(i)/r
c change independent variable for integration from phi to l (arc length
c ======================================================================
c right now, the phases that work for PREM model include P(S),PP(SS),
c PPP(SSS),PcP,ScS,SKS
c 8/20/99 S.-H. Hung; revised at 5/3/00
c note: small dh (dphi, integration interval) might be changed depending on
c the takeoff angle (steeper angle requires smaller dh and of course take 
c more time) in order to avoid the jump of discontinuity at shallower 
c crustal layers
c given the radius, theta, arclength, the kernels are calculated at points 
c surrounding the ray  ;
c subroutine called is kernel_tube
c PKPcd (PKiKP)
      program shootray_df_evt
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INCLUDE 'phases.inc'        
      data pha /'P','PP','PPP','PcP','PcP2','PcP3','PcP4',
     &     'PcP5','PcP6','pP',
     &     'pPP','Pd','PS','PcS','PKS','PKPab','PKPbc',
     &     'PKPdf','PKKPdf','PKPcd',
     &     'S','SS','SSS','ScS','ScS2','ScS3','ScS4',
     &     'ScS5','ScS6','sS',
     &     'sSS','Sd','SP','ScP','SKP','SKSac','SKKSac',
     &     'SKSdf','SKKSdf','SKiKP'/

      INTEGER n2,nvar,kmax,kount
      REAL*4 rtbeg,rtend,rta(2)
      REAL*8 x1,x2,dvdr(2),vr(2),ystart(NMX)
      REAL*8 vb1,vb2,fvec1,fvec2,aibd(4)
      REAL*8 wx(NTMX),ww(NTMX),stfa2(NTMX,20),
     &       tstar(20),tpstf(20),sdown(NTMX)
      REAL*8 e2g(3,3),g2e(3,3)
      REAL*8 tk(20),tf(20)
      INTEGER irp,nrbp,irbp(10)
      INTEGER nlags,lagb(10),lage(10)
      REAL*8 h1,hmin,y(NMX),rayprm(2),grang(2)
      EXTERNAL derivs,dfint1,dfint2
      CHARACTER*120 fid,fout,fevst,fnm,fnmray,fnmbp
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /rcvsrc/ rar,thrdeg,phrdeg,ras,thsdeg,phsdeg
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /rayint/ ips,idistrs,idepsrc,nsurf,ncmb,nocb,niob 
      COMMON /ptau/ srd(NPTMX),stt(NPTMX),srp(NPTMX),stau(NPTMX),
     &              dtdp2(NPTMX)
      COMMON /xyzpath/ xr(NSTPMX),yr(NSTPMX),zr(NSTPMX)
      COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),rkn(NKMX,MSMX,20)
      COMMON /diffr/ phip,phiq,tlp,tlq,tpq,cb,rb,rmu0c,crbt

      character*5 stnnm,qual*1,header*100
      real blat,blon,bdep,del,azs,ctim,ptim,tsatt

      DIMENSION slat(NOBSMX),slon(NOBSMX),sdep(NOBSMX),
     & rlat(NOBSMX),rlon(NOBSMX),rdep(NOBSMX),tdobs(NOBSMX)
      DIMENSION stnnm(NOBSMX),s2n(NOBSMX),nobs(NOBSMX),header(NOBSMX),
     & blat(NOBSMX),blon(NOBSMX),bdep(NOBSMX),
     & del(NOBSMX),azs(NOBSMX),ctim(NOBSMX),ptim(NOBSMX),
     & tsatt(NOBSMX),qual(NOBSMX)
      DIMENSION iblk(NBLKMX),rblk(NBLKMX),tblk(NBLKMX),pblk(NBLKMX),
     & skblk(20,NBLKMX),rbke(NBLKMX),tbke(NBLKMX),pbke(NBLKMX)

      REAL*8 rptbl(51,71),ditbl(51),dstbl(71)
      data dtt/0.20d0/, ttmx/100.0d0/, ms/4/
      data ndi /51/, nds /71/
      data ddi /5.0d0/, dds /10.0d0/
      data di0 /20.0d0/, ds0 /0.0d0/
      data irpflag /0/
      data rptbl /3621*0.0d0/
      data tol /1.0d-2/

      INCLUDE 'paths.inc'
      INCLUDE 'PQpaths.inc'
      INCLUDE 'models.inc'
      INCLUDE 'prem.inc'
      INCLUDE 'iasp.inc'
      INCLUDE 'grad.inc'
      INCLUDE 'flat.inc'
      INCLUDE 'homo.inc'
      INCLUDE 'lyrs.inc'
      INCLUDE 'ak135.inc'
c set up constants
      pi=dasin(1.0d0)*2.0d0
      pi2=pi/2.0d0
      twopi=2.0d0*pi
      rpd=pi/180.0d0
      dpr=180.0d0/pi
      vp0=8.0d0
      vs0=4.50d0
      dpdr=-0.0025d0
      dsdr=-0.0025d0
c 
c  input 
c
      write(6,*)'1: interactive; 2:read in from kern.in(2) ?'
      read(5,*) inflag
      if (inflag.eq.1) then
         write(6,*)'1-D background velocity structure:'
	 write(6,*)'-3: constant-velocity model;
     &              -2: homogeneous-layered model;
     &              -1: flattening model; 
     &               0: linear gradient model; 
     &               1: compute PREM at T=ivelin;
     &              -4: compute IASP91 model '
	 read(5,*) ivelin
         write(6,9)
         read(5,*) iraytype
         write(6,*)'radius and beginning theta of the ray-tube scatterer space'
         read(5,*)rtb,tbd0
         write(6,*)'integ. intervals dl (=h), radius, and theta in (km,km,deg) ?'
         read(5,*)dl,drtb,dtbd
         write(6,*)'# of points for gaussian-legendre integration'
         read(5,*)nn
c input info for source time function
         write(6,*)'source time function (0) Gaussian (1) instrument response'
         read(5,*) insflag, ntstf

c if interactive input, only one event-station pair will calculate
         if (insflag.eq.0) then
            write(6,*)'input source time function parameters:'
            write(6,*)'read in dt='
            read(5,*)dtstf
            write(6,*)'read in tau='
            read(5,*)(tpstf(i),i=1,ntstf)
c source time function
c lowest and highest frequencies
c skip source time function calculation when doing rayshooting
c            wb=0.0d0
c            we=0.5d0/dtstf
c           tp=dsqrt(2.0d0)*tp
c            do i=1,ntstf
c            alpha=4.0d0*pi*pi/(tpstf(i)*tpstf(i))
c            call calc_stfgs(alpha,nn,wb,we,wx,ww,stfa2(1,i))
c            enddo
         else
            write(6,*) 'read in lowest and highest frequencies'
            read(5,*)wb,we
c tpstf=tstar
            read(5,*)(tpstf(i),i=1,ntstf)
c            do i=1,ntstf
c            call calc_stfra(tpstf(i),nn,wb,we,wx,ww,stfa2(1,i))
c            enddo
         endif

         write(6,*)'event-station info and points at which kernels to be calculated:'
         write(6,*)'(0) interactive(nobs;slat,slon,sdep,rlat,rlon,rdep);'
         write(6,*)'(1) input from file(filename)'
         read(5,*)ievstflag
         if (ievstflag.eq.0) then
c         read in 'source position (latitude,longitude,depth):'
c         read in receiver position (latitude,longitude,depth):'
            read(5,*)nobs0
            do i=1,nobs0
            read(5,*)slat(i),slon(i),sdep(i),
     &      rlat(i),rlon(i),rdep(i)
            enddo
         else
            read(5,'(a80)')fevst
            open(1,file=fevst)
            i=0
            do ie=1,nobsmx
               read(1,'(a100)',end=7)header(ie)
               read(1,*)nobs(ie)
               do iobs=1,nobs(ie)
                  i=i+1
                  read(1,14,end=7)stnnm(i),rlat(i),rlon(i),rdep(i),
     &            slat(i),slon(i),sdep(i),ctim(i),s2n(i),ptim(i)
c                  write(6,*)i,rlat(i),stnnm(i)
                  rdep(i)=0.
               enddo
            enddo
7           nevt=ie-1
            close(1)
        endif
      else
c for inflag=2 input from file kern.in
         open(1,file='kern.in')
         read(1,*)
         read(1,*)ivelin
         read(1,*)iraytype
         read(1,*)rtb,tbd0
         read(1,*)dl,drtb,dtbd
         read(1,*)nn
         read(1,*)insflag,ntstf
         if (insflag.eq.0) then
            read(1,*)dtstf
            read(1,*)(tpstf(i),i=1,ntstf)
c            do i=1,ntstf
c            wb=0.0d0
c            we=0.5d0/dtstf
c            alpha=4.0d0*pi*pi/(tpstf(i)*tpstf(i))
c            call calc_stfgs(alpha,nn,wb,we,wx,ww,stfa2(1,i))
c            enddo
         else
            read(1,*)wb,we
            read(1,*)(tpstf(i),i=1,ntstf)
c            do i=1,ntstf
c            call calc_stfra(tpstf(i),nn,wb,we,wx,ww,stfa2(1,i))
c            enddo
         endif

         read(1,*)ievstflag
         if (ievstflag.eq.0) then
            read(1,*)nobs0
            do i=1,nobs0
            read(1,*)slat(i),slon(i),sdep(i),
     &      rlat(i),rlon(i),rdep(i),tdobs(i)
            enddo
            close(1)
         else
            read(1,'(a80)')fevst
            open(1,file=fevst)
            i=0
            do ie=1,nobsmx
               read(1,'(a100)',end=8)header(ie)
               read(1,*)nobs(ie)
               do iobs=1,nobs(ie)
                  i=i+1
                  read(1,14,end=8)stnnm(i),rlat(i),rlon(i),rdep(i),
     &            slat(i),slon(i),sdep(i),ctim(i),s2n(i),ptim(i)
c                  write(6,*)i,rlat(i),ctim(i),stnnm(i)
                  rdep(i)=0.
               enddo
            enddo
8           nevt=ie-1
            close(1)
        endif
6        format(i4)
      endif
12    format('> ',i4,2x,i3,2x,2(i2,2x),f6.3,f9.4,2x,f9.4,2x,f5.1)
c14    format(a5,1x,4(f9.4,1x),2f6.1)
14    format(a5,2(f11.6,1x),f9.4,1x,2(f9.4,1x),f6.1,f9.4,1x,f9.5,f10.3)
16    format(a5,1x,i4,2(2x,i2),1x,f6.3,4f8.2,2f6.1)
26    format('> ',a5,1x,i4,4f8.2,f6.1,2f8.2,f7.1,f7.2,f8.2,f9.3,1x,f6.1)
c25    format(a5,1x,2(f11.6,1x),2(f9.4,1x),f6.1,2f8.2,f7.1,f7.2,f8.2,f9.4,1x,f9.5,f10.3)
25    format(a5,1x,2(f11.6,1x),f8.4,1x,2(f9.4,1x),f6.1,2f8.2,f7.1,f7.2,f8.2,f9.4,f9.5,f9.3,f9.4)

c36    format(f10.3,2f10.4,f10.4,f10.3,f10.4,f8.2)
36    format(f10.3,2f10.4,f10.4,f10.3,f10.4,f8.2,f8.2)
46    format(2(f10.2,2f8.2))

9     format('Phases:',/2x,
     &         '1: P',/2x,
     &         '2: PP (or PKP2)',/2x,
     &         '3: PPP',/2x,
     &         '4: PcP',/2x,
     &         '5: PcP2',/2x,
     &         '6: PcP3',/2x,
     &         '7: PcP4',/2x,
     &         '8: PcP5',/2x,
     &         '9: PcP6',/2x,
     &         '10: pP',/2x,
     &         '11: pPP',/2x,
     &         '12: Pd',/2x,
     &         '13: PS',/2x,
     &         '14: PcS',/2x,
     &         '15: PKS',/2x,
     &         '16: PKPab',/2x,
     &         '17: PKPbc',/2x,
     &         '18: PKPdf',/2x,
     &         '19: PKKPdf',/2x,
     &         '20: PKPcd',/2x,
     &         '21: S',/2x,
     &         '22: SS',/2x,
     &         '23: SSS',/2x,
     &         '24: ScS',/2x,
     &         '25: ScS2',/2x,
     &         '26: ScS3',/2x,
     &         '27: ScS4',/2x,
     &         '28: ScS5',/2x,
     &         '29: ScS6',/2x,
     &         '30: sS',/2x,
     &         '31: sSS',/2x,
     &         '32: Sd',/2x,
     &         '33: SP',/2x,
     &         '34: ScP',/2x,
     &         '35: SKP',/2x,
     &         '36: SKSac',/2x,
     &         '37: SKKSac',/2x,
     &         '38: SKSdf',/2x,
     &         '39: SKKSdf',/2x,
     &         '40: SKiKP')

      call bpt(iraytype,ips,iupdown,nsurf,ncmb,nocb,niob)
      ips0=ips
      call rayptbl(iraytype,ndi,nds,ddi,dds,di0,ds0,
     &        irpflag,ditbl,dstbl,rptbl)

c      write(6,*)'ncmb=',ncmb
c       write(6,*)pha(iraytype)

c initialize array vdis
c      pi=datan(1.0d0)*4.0d0
c      pi2=pi/2.0d0
c      rpd=pi/180.d0

      rsurf=6371.0d0
      rn=1.0d0/rsurf
      r0=rsurf
      if (ivelin.eq.-4) then
      rcmb=3482.0d0
      riob=1217.1d0
      elseif (ivelin.ge.2) then
      rcmb=3479.5d0
      riob=1217.1d0
      else
      rcmb=3480.0d0
      riob=1221.5d0
      endif
      r660= 5711.0d0
      r670= 5701.0d0
      r400= 5971.0d0
      ntt=int(ttmx/dtt)+1
      do it=1,ntt
         ttdf(it)=dble(it-1)*dtt
      enddo

      write(6,21)

      write(fnmray,'(''shootrays.'',a5)')pha(iraytype)
      open(17,file=fnmray)
      write(fnmbp,'(''shootrays.bp.'',a5)')pha(iraytype)
      open(18,file=fnmbp)
      ibeg=1
      do 2100 ievt=1,nevt
         if (nobs(ievt).ne.0) then
         write(17,'(a100)') header(ievt)
         write(17,'(i3)') nobs(ievt)
         endif
         iend=nobs(ievt)+ibeg-1
c         write(6,*)ievt,nobs(ievt),header(ievt)
      do 2000 iobs=ibeg,iend
c calculate azimuth and delta (distrs) in radius and degree)
c      slat(iobs)=90.0d0-slat(iobs)
      ts=90.0d0-slat(iobs)
      ps=slon(iobs)
c      ds=sdep(iobs)
c      rlat(iobs)=90.0d0-rlat(iobs)
      tr=90.0d0-rlat(iobs)
      pr=rlon(iobs)
      call stadis_sph(ts,ps,tr,pr,distrs,azsrdeg)
      radistrs=distrs*rpd
c        call azimdist(90.-slat(iobs),90.0d0-rlat(iobs),slon(iobs),rlon(iobs),
c     &    azsr,azrs,azsrdeg,azrsdeg,distrs,radistrs)
c        write(6,*)'epicentral distance= ',ievt,iobs,distrs,
c     &  slat(iobs),slon(iobs),sdep(iobs),rlat(iobs),rlon(iobs),rdep(iobs)
       del(iobs)=distrs
       azs(iobs)=azsrdeg
       radsrc=rsurf-sdep(iobs)
       radrcv=rsurf-rdep(iobs)
       idepsrc=int(sdep(iobs))
       idistrs=int(distrs)

       sth=ts*rpd
       sph=ps*rpd
       rth=tr*rpd
       rph=pr*rpd

c euler_rot: calculate rotation matrix from finding cos(x',x)
         call euler_rot(radsrc,sth,sph,rsurf,rth,rph,g2e,e2g)
c       write(6,*)'backazimuth and epicentral distance=',azsrdeg,distrs
	
c get file id that composesd by rsdist,source depth and phase info.
         call getfid(ivelin,iraytype,idistrs,idepsrc,
     &    distrs,deps,nid,fid)
c determine velocity at discontinuities
c
         call velodisc(ivelin,rayprm,grang)

c start ray tracing
      nvar=4
      x1=0.0d0
      x2=radistrs
      dh=dl
c tracing diffracted waves
      call raydiff(ivelin,iraytype,rayprm,grang,nvar,dh,
     &  x1,x2,ystart,ni,nid,fid,nlags,lagb,lage,nidf,isdflag)

      if (isdflag.eq.1) then
         ramp=rcmb
         thmp=pi2
         phmp=(phip+phiq)*0.5
         ni=nidf
         go to 1800
      endif

c find the index for bounce points which will be used to determine projection point
c find the geographical coordinate (rpg,tpg,ppg) of p and (rqg,tqg,pqg) of q
c
      call rtp2xyz(yp(1,lage(1)),pi2,phip,xmp,ymp,zmp)
      call rotate(xmp,ymp,zmp,e2g,xmpg,ympg,zmpg)
      call xyz2rtp(xmpg,ympg,zmpg,rpg,tpg,ppg)
      call rtp2xyz(yp(1,lagb(3)),pi2,phiq,xmp,ymp,zmp)
      call rotate(xmp,ymp,zmp,e2g,xmpg,ympg,zmpg)
      call xyz2rtp(xmpg,ympg,zmpg,rqg,tqg,pqg)

c find the index for bounce points which will be used to determine projection
c points
       call raylag(iraytype,ni,nlags,lagb,lage,rayp,ramp,thmp,phmp)
1800   continue
       call rtp2xyz(ramp,thmp,phmp,xmp,ymp,zmp)
       call rotate(xmp,ymp,zmp,e2g,xmpg,ympg,zmpg)
       call xyz2rtp(xmpg,ympg,zmpg,rmpg,tmpg,pmpg)
       write(18,46)rsurf-ramp,90.-tmpg*dpr,pmpg*dpr,rsurf-rmpg,90.-tmpg*dpr,pmpg*dpr
       blat(iobs)=90.-tmpg*dpr
       blon(iobs)=pmpg*dpr
       bdep(iobs)=rsurf-rmpg
c tc= topo correction
       tc=rdep(iobs)*dcos(yp(2,ni))/vp(ni)

       write(17,25)stnnm(iobs),rlat(iobs),rlon(iobs),rdep(iobs),slat(iobs),slon(iobs),
     &  sdep(iobs),blat(iobs),blon(iobs),bdep(iobs),
     &  del(iobs),azs(iobs),ctim(iobs),s2n(iobs),ptim(iobs),tc
       do i=1,ni
c transform from equatorial plane to great-circle plane between s and r
         ramp=yp(1,i)
         phmp=yp(3,i)
         call rtp2xyz(ramp,pi2,phmp,xmp,ymp,zmp)
         call rotate(xmp,ymp,zmp,e2g,xmpg,ympg,zmpg)
         call xyz2rtp(xmpg,ympg,zmpg,rmpg,tmpg,pmpg)
          write(17,36)rsurf-rmpg,90.-tmpg*dpr,pmpg*dpr,
     &     xp(i),yp(1,i),yp(3,i)*dpr,yp(2,i)*dpr,vp(i)
       enddo

c       do i=1,nlags
c       write(6,*)i,lagb(i),lage(i)
c       enddo
c       write(6,*)'finish calculate hessian matrix !'
c39    format(3f12.4)

c calculate hessian matrix using rk4 method
19     format(f10.3,3e16.8,f8.3,f10.3)

2000   continue
       ibeg=iend+1
2100   continue

21    format('distance, azs2r, src depth, iangle, rayprm, geom. spreading, midpt (d,t,p)')
31    format(3f10.2,f8.1,f12.4,e13.4,f8.1,2f8.2)
c       write(2,'(''time spent in kernel='',f14.6)')rtend-rtbeg
       stop
       end

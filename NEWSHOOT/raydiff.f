c raytracing diffracted P or S waves (iraytype=12 or 32)
c
      subroutine raydiff(ivelin,iraytype,rayprm,grang,nvar,dl,
     &  x1,x2,ystart,ni,nid,fid,nlags,lagb,lage,nidf,isdflag)

      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      PARAMETER (EPS=1.0d-8,TOL=1.0d-6)
      INTEGER nlags,lagb(10),lage(10),nidf
      REAL*8 rayprm(2),grang(2)
      REAL*8 x1,x2,dvdr(2),vr(2),ystart(NMX)
      INTEGER iraytype,idistrs,idepsrc,nsurf,ncmb
      CHARACTER*80 fid,fout
      EXTERNAL derivs_df
      COMMON /rayint/ ips,idistrs,idepsrc,nsurf,ncmb,nocb,niob,ndlb
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn,rdlb
      COMMON /diffr/ phip,phiq,tlp,tlq,tpq,cb,rb,rmu0c,crbt
      INCLUDE 'paths.inc'

      isdflag=1
      nlags=3
      dh=dl
c for diffracted P or S waves
      call findr(ivelin,rsurf,irrec)
      call velo(ivelin,irrec,rsurf,dvdr,vr)
      if (iraytype.eq.12) then
         rayp=rayprm(1)
         aincs=grang(1)
         aincr=pi-dasin(rayp*vr(1)/rsurf)
      elseif (iraytype.eq.32) then
         rayp=rayprm(2)
         aincs=grang(2)
         aincr=pi-dasin(rayp*vr(2)/rsurf)
      endif
      rayp2=rayp*rayp
c      fout='shoot.'//fid(1:nid)
c      open(1,file=fout)
      ystart(1)=radsrc
      ystart(3)=0.0d0
      ystart(2)=aincs
      ystart(4)=0.0d0
      call tracer_df(ivelin,ips,nsurf,ncmb,nocb,nvar,ystart,
     &  x1,x2,dh,fvec1,ni,isurf,icmb,iocb,iiob,idlb)
      xlp=xpdf(ni)
      phip=ypdf(3,ni)
      tlp=ypdf(4,ni)
      do i=1,ni
         xp(i)=xpdf(i)
         do j=1,nvar
            yp(j,i)=ypdf(j,i)
         enddo
         vp(i)=vpdf(i)
         dvp(i)=dvpdf(i)
         d2vp(i)=d2vpdf(i)
         irp(i)=irpdf(i)
         iwp(i)=iwpdf(i)
         write(1,39)xp(i),yp(1,i),yp(2,i)*dpr,yp(3,i)*dpr,yp(4,i),
     &   vp(i),dvp(i),d2vp(i),irp(i),iwp(i)
c         write(6,39)xp(i),yp(1,i),yp(2,i)*dpr,yp(3,i)*dpr,yp(4,i),
c     &   vp(i),dvp(i),d2vp(i),irp(i),iwp(i)
      enddo
      ii=ni
      lagb(1)=1
      lage(1)=ni
      lagb(2)=ni+1
c      write(6,*)'phi at grazing=',yp(1,ni),yp(2,ni)*dpr,yp(3,ni)*dpr
c from receiver to q
      dh=dl
      ystart(1)=rsurf
      ystart(3)=0.0d0
      ystart(2)=aincr
      ystart(4)=0.0d0
      call tracer_df(ivelin,ips,nsurf,ncmb,nocb,nvar,ystart,
     &  x1,x2,dh,fvec1,ni,isurf,icmb,iocb,iiob,idlb)
      phiq=radistrs-ypdf(3,ni)
c      if (phiq.lt.phip) then
c         write(6,*)'no grazing ray at epicentral distance=',radistrs*dpr
c         isdflag=0
c         return
c      endif
c      do k=1,ni
c         write(6,39)xp(k),ypdf(1,k),ypdf(2,k)*dpr,ypdf(3,k)*dpr,ypdf(4,k)
c      enddo
      tq=ypdf(4,ni)
      rq=ypdf(1,ni)
      aincq=ypdf(2,ni)
      phipq=phiq-phip
      nphi=int(phipq*rcmb/dl)+1
      dphi=phipq/dble(nphi)

      call velo1(ivelin,ips,3,rcmb,c,dc,d2c)
      cb=c
      dcb=dc
      rb=1.0d0/rcmb-dcb/cb
      rb=1.0d0/rb
c      write(6,*)'c=',c,dc,d2c,dphi,dl
c along cmb
      do i=1,nphi-1
         ii=ii+1
         xlb=dble(i)*dphi*rcmb
         xl=xlp+xlb
         phi=phip+dble(i)*dphi
         tl=tlp+xlb/c
         xp(ii)=xl
         yp(1,ii)=rcmb
         yp(2,ii)=pi2
         yp(3,ii)=phi
         yp(4,ii)=tl
         vp(ii)=c
         dvp(ii)=dc
         d2vp(ii)=d2c
         irp(ii)=3
         iwp(ii)=ips
c         write(1,39)xp(ii),yp(1,ii),yp(2,ii)*dpr,yp(3,ii)*dpr,yp(4,ii),
c     &   vp(ii),dvp(ii),d2vp(ii),irp(ii),iwp(ii)
      enddo
      lage(2)=ii
      
c along segment III
      dh=dl
      xlf=(phiq-phi)*rcmb
      xl=xl+xlf
      xlq=xl
      phi=phiq
      tl=tl+xlf/c
      tlq=tl
c      x1=xl
c      x2=radistrs
c      ystart(1)=rq
c      ystart(2)=pi-aincq
c      ystart(3)=phiq
c      ystart(4)=tlq
c      call tracer_df(ivelin,ips,nsurf,ncmb,nocb,nvar,ystart,
c     &  x1,x2,dh,fvec1,ni,isurf,icmb,iocb,iiob)
c      write(6,*)'ypdf=',ypdf(1,ni),yp(2,ni)*dpr
      lagb(3)=ii+1
      do i=ni,1,-1
         ii=ii+1
         xp(ii)=xl+(xpdf(ni)-xpdf(i))
         yp(1,ii)=ypdf(1,i)
         yp(2,ii)=pi-ypdf(2,i)
         yp(3,ii)=phi+ypdf(3,ni)-ypdf(3,i)
         yp(4,ii)=tl+ypdf(4,ni)-ypdf(4,i)
         vp(ii)=vpdf(i)
         dvp(ii)=dvpdf(i)
         d2vp(ii)=d2vpdf(i)
         irp(ii)=irpdf(i)
         iwp(ii)=iwpdf(i)
c         write(1,39)xp(ii),yp(1,ii),yp(2,ii)*dpr,yp(3,ii)*dpr,yp(4,ii),
c     &   vp(ii),dvp(ii),d2vp(ii),irp(ii),iwp(ii)
c         write(6,39)xp(ii),yp(1,ii),yp(2,ii)*dpr,yp(3,ii)*dpr,yp(4,ii),
c     &   vp(ii),dvp(ii),d2vp(ii),irp(ii),iwp(ii)
      enddo
      lage(3)=ii
      nidf=ii
      tpq=tlq-tlp
c      write(6,*)'nidf=',nidf
c      close(1)
c      write(6,*)'cb,rb=',cb,rb
39    format(4f14.6,f8.2,3e14.6,2i6)
c      do i=1,ni
c         xp(i)=xpdf(i)
c         yp(1,i)=ypdf(1,i)
c         yp(2,i)=ypdf(2,i)
c         yp(3,i)=ypdf(3,i)
c         yp(4,i)=ypdf(4,i)
c      enddo
      return
      end

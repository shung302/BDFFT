c rotate any given two points to form an equator plane
c
        IMPLICIT DOUBLEPRECISION (a-h,o-z)
        COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
        real*8 mth,mph,mlat,mlon
        dimension g2e(3,3),e2g(3,3)

	write(6,*)'(lon,lat) of two points to form the equator plane'
	read(5,*)alon,alat,blon,blat
        write(6,*)'the mesh assumed dx=dy'
        write(6,*)'read in grids of nx, ny, nz (33,33,33) or (65,65,65)'
        read(5,*)nx,ny,nz 
        write(6,*)'the depth extent between d1 and d2 (d1<d2), enter d1, d2'
        read(5,*)d1,d2

        
        pi=dasin(1.0d0)*2.0d0
        twopi=pi*2.0d0       
        pi2=pi*0.5d0
        rpd=pi/180.0d0
        dpr=180.0d0/pi
        rsurf=6371.0d0
        rnorm=1.0d0
        rd1=rsurf-d1
        rd2=rsurf-d2
        ts=90.0d0-alat
        ps=alon
        tr=90.0d0-blat
        pr=blon
        sth=ts*rpd
        sph=ps*rpd
        rth=tr*rpd
        rph=pr*rpd
c  ts,ps and tr,pr are colatitude and longitude
        call stadis_sph(ts,ps,tr,pr,distrs,azsrdeg)
        write(6,*)'distrs=',distrs
c sth,sph are latitude and longitude
        call euler_rot(rnorm,sth,sph,rnorm,rth,rph,g2e,e2g,plat,plon)
        write(6,'(4f15.5)')plon,plat,mlon,mlat
c find mid point
        rdist2=distrs*rpd*0.5d0
	write(6,*)'rdist2=',rdist2,pi2
        call rtp2xyz(rnorm,pi2,rdist2,xtmp,ytmp,ztmp)
        call rotate(xtmp,ytmp,ztmp,e2g,xtmpg,ytmpg,ztmpg)
        call xyz2rtp(xtmpg,ytmpg,ztmpg,rm,tm,pm)
        mth=tm
        mph=pm
        tm=tm*dpr
        pm=pm*dpr
        mlat=90.0d0-tm
        mlon=pm
        write(6,*)'mid point=',pm,90.0d0-tm
        call stadis_sph(tm,pm,tr,pr,distrs2,azsrdeg)
	write(6,*)'distrs2=',distrs2
        distrs=distrs2*2.0d0
c sth,sph are latitude and longitude
        call euler_rot(rnorm,mth,mph,rnorm,rth,rph,g2e,e2g,plat,plon)

	dx=distrs/dble(nx-1)
        dy=dx
        if (ny.eq.nx) then
          yfrac=1.0d0
        else
          yfrac=dble(ny-1)/dble(nx-1)
        endif

        ylendeg=distrs*yfrac
        xlendeg=distrs
        ybdeg=-ylendeg*0.5d0
        yedeg=ylendeg*0.5d0
        xbdeg=-xlendeg*0.5d0
        xedeg=xlendeg*0.5d0
        dz=(rd1-rd2)/dble(nz-1)
        zb=rd2
        ze=rd1
        open(1,file='mesh.xy')
        open(2,file='mesh.xy1')
        open(11,file='oblique_mesh.config')

        levelx=nint(log(real(nx))/log(2.)+0.5)
        levely=nint(log(real(ny))/log(2.)+0.5)
        levelz=nint(log(real(nz))/log(2.)+0.5)

        write(11,'(6i5,3f15.5)')levelx,levely,levelz,nx,ny,nz
        write(11,'(6f15.5)')xbdeg,ybdeg,rsurf-d2,dx,dy,dz
        write(11,'(3f15.5)')((g2e(i,j),j=1,3),i=1,3)
c        write(11,'(3f15.5)')((e2g(i,j),j=1,3),i=1,3)
c        write(11,'(4f15.5,2f11.1)')alon,alat,blon,blat,d1,d2
c        write(11,'(4f15.5)')plon,plat,mlon,mlat
c        write(11,'(3f15.5)')((e2g(i,j),j=1,3),i=1,3)
        write(6,'(4f15.5)')plon,plat,mlon,mlat
        
        yrad=(90.-29.23399d0)*rpd
        xrad=88.850998d0*rpd
        call rtp2xyz(rnorm,yrad,xrad,xtmp,ytmp,ztmp)
        write(6,*)'xtmp=',xtmp,ytmp,ztmp
        call rotate(xtmp,ytmp,ztmp,g2e,xtmpg,ytmpg,ztmpg)
        write(6,*)'xtmpg=',xtmpg,ytmpg,ztmpg
        call xyz2rtp(xtmpg,ytmpg,ztmpg,rg,tg,pg)
        write(6,*)'project',90.-tg*dpr,pg*dpr
 
        do j=1,ny
           y=ybdeg+(j-1)*dy
           yrad=(90.0d0-y)*rpd
           write(1,'(''> y='',f10.3)')y
           write(2,'(''> y='',f10.3)')y
           do i=1,nx
              x=xbdeg+(i-1)*dx
              xrad=x*rpd
              call rtp2xyz(rnorm,yrad,xrad,xtmp,ytmp,ztmp)
              call rotate(xtmp,ytmp,ztmp,e2g,xtmpg,ytmpg,ztmpg)
              call xyz2rtp(xtmpg,ytmpg,ztmpg,rg,tg,pg)
              tg=90.0d0-tg*dpr
              pg=pg*dpr
              write(1,'(4f10.3)')sngl(pg),sngl(tg),x,y
              if (i.eq.1.or.i.eq.nx) write(2,'(4f10.3)')sngl(pg),sngl(tg),x,y
           enddo
        enddo
         
        do i=1,nx
           x=xbdeg+(i-1)*dx
           xrad=x*rpd
           write(1,'(''> x='',f10.3)')x
           write(2,'(''> x='',f10.3)')x
           do j=1,ny
              y=ybdeg+(j-1)*dy
              yrad=(90.0d0-y)*rpd
              call rtp2xyz(rnorm,yrad,xrad,xtmp,ytmp,ztmp)
              call rotate(xtmp,ytmp,ztmp,e2g,xtmpg,ytmpg,ztmpg)
              call xyz2rtp(xtmpg,ytmpg,ztmpg,rg,tg,pg)
              tg=90.0d0-tg*dpr
              pg=pg*dpr
              write(1,'(4f10.3)')sngl(pg),sngl(tg),x,y
              if (j.eq.1.or.j.eq.ny) write(2,'(4f10.3)')sngl(pg),sngl(tg),x,y
           enddo
        enddo

        stop
	end

c-------------------------------------------------------
c subroutine to do rotation using euler angles
c input: vector (x,y,z), cosine of three euler angles 
c
      subroutine euler_rot(rs,ts,ps,rr,tr,pr,g2e,e2g,tpole,ppole)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      dimension g2e(3,3),e2g(3,3)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

      call cart(ts,ps,xi,yi,zi)
      call cart(tr,pr,xr,yr,zr)
c      call rtp2xyz(rs,ts,ps,xi,yi,zi)
c      call rtp2xyz(rr.0d0,tr,pr,xr,yr,zr)
c      write(6,*)'source and receiver in x,y,z'
c      write(6,*)xi,yi,zi
c      write(6,*)xr,yr,zr
c pole (xk,yk,zk) = (xi,yi,zi) x (xr,yr,zr) of great circle plane between s and r
c  
      xk=yi*zr-yr*zi
      yk=xr*zi-xi*zr
      zk=xi*yr-xr*yi
      call xyz2rtp(xk,yk,zk,rpole,tpole,ppole)
      tpole=90.0d0-tpole*dpr
      ppole=ppole*dpr
      call cross(xi,yi,zi,xr,yr,zr,xk,yk,zk)

      dk=dsqrt(xk*xk+yk*yk+zk*zk)
cc       write(6,*)'dk=',dk
      xk=xk/dk
      yk=yk/dk
      zk=zk/dk
c (xj,yj,zj) = (xk,yk,zk) x (xi,yi,zi)
      xj=yk*zi-yi*zk
      yj=xi*zk-xk*zi
      zj=xk*yi-xi*yk
      dj=dsqrt(xj*xj+yj*yj+zj*zj)
      xj=xj/dj
      yj=yj/dj
      zj=zj/dj

c rotation matrix cos(x',x)
      g2e(1,1)=xi
      g2e(2,1)=xj
      g2e(3,1)=xk
      g2e(1,2)=yi
      g2e(2,2)=yj
      g2e(3,2)=yk
      g2e(1,3)=zi
      g2e(2,3)=zj
      g2e(3,3)=zk
c      write(6,*)'a matrix'
c      do i=1,3
c         write(6,*)(g2e(i,j),j=1,3)
c      enddo

      e2g(1,1)=xi
      e2g(2,1)=yi
      e2g(3,1)=zi
      e2g(1,2)=xj
      e2g(2,2)=yj
      e2g(3,2)=zj
      e2g(1,3)=xk
      e2g(2,3)=yk
      e2g(3,3)=zk

      return
      end

c--------------------------------------------------------------------------
c subroutine packages to calculate delta, azimuth, 
c and bounce points, etc...
c from Guy Masters (geographic coords are used when calculating
c epicentral distance and azimuth;ie., assuming geographic coords.)
c
      subroutine stadis_sph(colat,colon,scolat,scolon,del,az)

c Computes the epicentral distance and azimuth from source to receiver.
c Latitudes are converted to geocentric latitudes prior to performing
c the computations (it is assumed that input latitudes are geographic).
c  input:
c    colat = source colatitude, (degrees)
c    colon =   "    colongitude,  "      
c    scolat    = station colatidue, 
c    scolon    =    "    colongitude, 

c  output:
c    del   = epicentral distance (degrees)
c    az    = azimuth from source to receiver, measure from North (degrees)

      implicit real*8 (a-h,o-z)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

c  first do eq coords.
c      t0=colat/rad
c      t0=geocen(colat*rpd)    
      t0=colat*rpd
      p0=colon*rpd
      c0=dcos(t0)
      s0=dsin(t0)
c  now do station coords.
c     t2=scolat*rpd
c      t2=geocen(scolat*rpd)           
      t2=scolat*rpd
      c1=dcos(t2)
      s1=dsin(t2)
      p1=scolon*rpd
c  now calculate distance
      dp=p1-p0
      co=c0*c1+s0*s1*dcos(dp)
      si=dsqrt(1.d0-co*co)
      del=datan2(si,co)*dpr
c  now calculate azimuth
      caz=(c1-c0*co)/(si*s0)
      dp2=-dp
      saz=-s1*dsin(dp2)/si
      az=datan2(saz,caz)*dpr
c	write(6,*)'del=',del,colat,colon,rpd,dpr,rcmb,rsurf
      if(az.lt.0.0d0) az=360.0d0 + az  
      return
      end

c-----------------------------------------------------------------
c convert a point at spherical coordinate to cartesian
c
      subroutine rtp2xyz(r,t,p,x,y,z)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
c      real*8 r,t,p,x,y,z,st,ct,sp,cp

      if (t.ne.pi2) go to 10
      st=1.0d0
      ct=0.0d0
10    st=dsin(t)
      ct=dcos(t)
      sp=dsin(p)
      cp=dcos(p)
      x=r*st*cp
      y=r*st*sp
      z=r*ct

      return
      end

c--------------------------------------------------------------
c rotate coordinate in one to another given a rotation matrix
	subroutine rotate(xo,yo,zo,rotm,xn,yn,zn)
	implicit REAL*8 (a-h,o-z)
        dimension rotm(3,3)

	xn=xo*rotm(1,1)+yo*rotm(1,2)+zo*rotm(1,3)
	yn=xo*rotm(2,1)+yo*rotm(2,2)+zo*rotm(2,3)
	zn=xo*rotm(3,1)+yo*rotm(3,2)+zo*rotm(3,3)
c	vlen=sqrt(xn**2+yn**2+zn**2)
c	xn=xn/vlen
c	yn=yn/vlen
c	zn=zn/vlen

	return
	end
	
c-----------------------------------------------------------------
c convert a point at cartesian to spherical coordinate
c
      subroutine xyz2rtp(x,y,z,r,t,p)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn


      r=dsqrt(x*x+y*y+z*z)
      t=dacos(z/r)
      if (x.eq.0.0d0) then
         p=pi2
         return
      endif
      p=datan2(y,x)

      return
      end

c
c--------------------------------------------------------------------------
c
      subroutine cart(thet,phi,x,y,z)
      implicit real*8(a-h,o-z)
      s=dsin(thet)
      x=s*dcos(phi)
      y=s*dsin(phi)
      z=dcos(thet)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine cross(sx,sy,sz,rx,ry,rz,px,py,pz)
      implicit real*8(a-h,o-z)
      px = sy*rz - sz*ry
      py = sz*rx - sx*rz
      pz = sx*ry - sy*rx
      return
      end
c
c----------------------------------------------------------------------



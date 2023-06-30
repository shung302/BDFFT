      subroutine tran(on,at,onp,atp,x,y,z,id)
c     (on,at) is lat,lon of Geographic Earth
c     (onp,atp) is lat,lon where center of the region is (0,0)
c     id = 1 --> (on,at) to (onp,atp)
c        =-1 --> (onp,atp) to (on,at)
      double precision trans(3,3)
      common/meshb2/trans
      double precision pi,rad
      parameter(pi=3.14159265d0,rad=180.d0/pi)
      double precision on,at,onp,atp,x,y,z,xp,yp,zp,ontmp,attmp

      if(id.eq.1)then
       ontmp=on/rad
       attmp=at/rad
       x=dcos(attmp)*dcos(ontmp)
       y=dcos(attmp)*dsin(ontmp)
       z=dsin(attmp)
       xp=trans(1,1)*x+trans(1,2)*y+trans(1,3)*z
       yp=trans(2,1)*x+trans(2,2)*y+trans(2,3)*z
       zp=trans(3,1)*x+trans(3,2)*y+trans(3,3)*z
       r=dsqrt(xp*xp+yp*yp+zp*zp)
       atp=dasin(zp/r)*rad
       onp=datan2(yp,xp)*rad
c       if(onp.lt.0.d0)onp=onp+360.d0
       if(onp.gt.180.d0)onp=onp-360.d0
       x=xp/r
       y=yp/r
       z=zp/r
      elseif(id.eq.-1)then
       ontmp=onp/rad
       attmp=atp/rad
       xp=dcos(attmp)*dcos(ontmp)
       yp=dcos(attmp)*dsin(ontmp)
       zp=dsin(attmp)
       x=trans(1,1)*xp+trans(2,1)*yp+trans(3,1)*zp
       y=trans(1,2)*xp+trans(2,2)*yp+trans(3,2)*zp
       z=trans(1,3)*xp+trans(2,3)*yp+trans(3,3)*zp
       r=dsqrt(x*x+y*y+z*z)
       at=dasin(z/r)*rad
       on=datan2(y,x)*rad
       if(on.lt.0.d0)on=on+360.d0
       x=x/r
       y=y/r
       z=z/r
      else
       stop' not permitted in tran !'
      endif
      return
      end


c______________________________________________________________________
      subroutine inside(on1,at1,dp1,on2,at2,dp2,id,npath,npaths,path,
     ,                  nlocat,locat)
      parameter(mpath=2000)
      dimension npaths(mpath),locat(mpath)
      double precision path(mpath),spath(mpath)
      double precision on1,at1,dp1,onp1,atp1,on2,at2,dp2,onp2,atp2
      double precision x0,y0,z0,dx,dy,dz,x1,y1,z1,x2,y2,z2
      double precision trans(3,3),rad,pi,
     ,                 xt1,yt1,zt1,xt2,yt2,zt2,xt3,yt3,zt3
      parameter(pi=3.14159265d0,rad=pi/180.d0)
      parameter(eps=1e-5)
      double precision rr,a,b,c,d,e,ss,xx,yy,zz
      common/meshb1/x0,y0,z0,dx,dy,dz,nx,ny,nz,nface,nvolu
      common/meshb2/trans
      common/storex/onp2,atp2,x2,y2,z2,ix21,iy21,iz21,ix22,iy22,iz22

      npath=0
      if(id.eq.0)then

       call tran(on1,at1,onp1,atp1,x1,y1,z1,1)
       x1=x1*dp1
       y1=y1*dp1
       z1=z1*dp1
       iz1=idint((dp1-z0)/dz)+1
       ix1=idint((onp1-x0)/dx)+1
       iy1=idint((atp1-y0)/dy)+1
       if(iz1.lt.1.or.ix1.lt.1.or.iy1.lt.1.or.
     ,    iz1.gt.nz.or.ix1.gt.nx.or.iy1.gt.ny)then
        return
       else
        id=1
       endif

      elseif(id.eq.1)then

       onp1=onp2
       atp1=atp2
       ix1=ix2
       iy1=iy2
       iz1=iz2
       x1=x2
       y1=y2
       z1=z2

      else
       stop' not permitted in inside !'
      endif

      call tran(on2,at2,onp2,atp2,x2,y2,z2,1)
      x2=x2*dp2
      y2=y2*dp2
      z2=z2*dp2
      iz2=idint((dp2-z0)/dz)+1
      ix2=idint((onp2-x0)/dx)+1
      iy2=idint((atp2-y0)/dy)+1
      if(iz2.lt.1.or.ix2.lt.1.or.iy2.lt.1.or.
     ,   iz2.gt.nz.or.ix2.gt.nx.or.iy2.gt.ny)then
       id=0
       return
      endif

      rr=dsqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))
c if r1 and r2 are the same points, no need to compute ray segments
      if (rr.lt.eps) return

      npath=2
      spath(1)=0.d0
      spath(2)=1.d0
      iz31=min0(iz1,iz2)-1
      if(iz31.lt.1)iz31=1
      iz32=max0(iz1,iz2)+1
      if(iz32.gt.nz)iz32=nz
      ix31=min0(ix1,ix2)-1
      if(ix31.lt.1)ix31=1
      ix32=max0(ix1,ix2)+1
      if(ix32.gt.nx)ix32=nx
      iy31=min0(iy1,iy2)-1
      if(iy31.lt.1)iy31=1
      iy32=max0(iy1,iy2)+1
      if(iy32.gt.ny)iy32=ny

c     crossing r-surface
      a=rr*rr
      b=2.d0*(x1*(x2-x1)+y1*(y2-y1)+z1*(z2-z1))/a
      do i=iz31,iz32
       zz=z0+dble(i-1)*dz
       c=(x1*x1+y1*y1+z1*z1-zz*zz)/a
       d=b*b-4.d0*c
       if (d.gt.0.)then
          d=sqrt(d)
          ss=(-b+d)/2.d0
          if((ss.lt.1.d0).and.(ss.gt.0.d0))then
             npath=npath+1
             spath(npath)=ss
          endif
          ss=(-b-d)/2.d0
          if((ss.lt.1.d0).and.(ss.gt.0.d0))then
             npath=npath+1
             spath(npath)=ss
          endif
       endif
      end do

c     crossing lat-surface
      do i=iy31,iy32
       yy=y0+dble(i-1)*dy
       yy=yy*rad
       b=dsin(yy)*6371.d0
       ss=(b-z1)/(z2-z1)
       if((ss.lt.1.d0).and.(ss.gt.0.d0))then
        npath=npath+1
        spath(npath)=ss
       endif
      end do

c     crossing lon-surface
      do i=ix31,ix32
       xx=x0+dble(i-1)*dx
       xx=xx*rad
       a=dcos(xx)*6371.d0
       b=dsin(xx)*6371.d0
       c=0.d0
       ss=(a*y1-b*x1)/(b*(x2-x1)-a*(y2-y1))
       if((ss.lt.1.d0).and.(ss.gt.0.d0))then
        npath=npath+1
        spath(npath)=ss
       endif
      end do

      call piksrt(npath,spath)
      xt1=x1
      yt1=y1
      zt1=z1
      a=x2-x1
      b=y2-y1
      c=z2-z1
      nseg=npath
      npath=0
      nlocat=0
      do i=2,nseg
       ss=spath(i)
       xt3=x1+ss*a
       yt3=y1+ss*b
       zt3=z1+ss*c
       xt2=(xt1+xt3)*0.5d0
       yt2=(yt1+yt3)*0.5d0
       zt2=(zt1+zt3)*0.5d0
       zz=dsqrt(xt2*xt2+yt2*yt2+zt2*zt2)
       iz=idint((zz-z0)/dz)+1
       yy=dasin(zt2/zz)/rad
       iy=idint((yy-y0)/dy)+1
       xx=datan2(yt2,xt2)/rad
       ix=idint((xx-x0)/dx)+1
       nlocat=nlocat+1
       locat(nlocat)=(iz-1)*(nx-1)*(ny-1)+(iy-1)*(nx-1)+ix
       rr=dsqrt((xt3-xt1)*(xt3-xt1)+
     +          (yt3-yt1)*(yt3-yt1)+(zt3-zt1)*(zt3-zt1))/8.d0
       xt1=xt3
       yt1=yt3
       zt1=zt3
       ibase0=(iz-1)*nx*ny
       npath=npath+1
       if(npath.gt.mpath)stop' mpath too small!'
       npaths(npath)=ibase0+(iy-1)*nx+ix
       path(npath)=rr
       npath=npath+1
       if(npath.gt.mpath)stop' mpath too small!'
       npaths(npath)=ibase0+(iy-1)*nx+ix+1
       path(npath)=rr
       npath=npath+1
       if(npath.gt.mpath)stop' mpath too small!'
       npaths(npath)=ibase0+(iy)*nx+ix
       path(npath)=rr
       npath=npath+1
       if(npath.gt.mpath)stop' mpath too small!'
       npaths(npath)=ibase0+(iy)*nx+ix+1
       path(npath)=rr
       ibase0=(iz)*nx*ny
       npath=npath+1
       if(npath.gt.mpath)stop' mpath too small!'
       npaths(npath)=ibase0+(iy-1)*nx+ix
       path(npath)=rr
       npath=npath+1
       if(npath.gt.mpath)stop' mpath too small!'
       npaths(npath)=ibase0+(iy-1)*nx+ix+1
       path(npath)=rr
       npath=npath+1
       if(npath.gt.mpath)stop' mpath too small!'
       npaths(npath)=ibase0+(iy)*nx+ix
       path(npath)=rr
       npath=npath+1
       if(npath.gt.mpath)stop' mpath too small!'
       npaths(npath)=ibase0+(iy)*nx+ix+1
       path(npath)=rr
c	write(*,*)i,npath,xx,yy,zz,ix,iy,iz,rr
      end do
      return
      end

c______________________________________________________________________
      subroutine linit
c ... initialize for conjugate gradient least squares
      parameter(mlx=7,mly=7,mlz=7,
     ,          mnx=2**(mlx-1)+1,mny=2**(mly-1)+1,
     ,          mnz=2**(mlz-1)+1,
     ,          mvolu=mnx*mny*mnz,mdata=200000,mpath=2000,mst=1000,mevt=10000,
     ,          mvar=mvolu+mst+mevt)
      include 'lsqrcom1'
      n = 0
      m = 0
      nel=0
      return
      end


c$$$      subroutine lldrow(coef,jdx,ncoef)
c$$$c ... add one equation of data into row of matrix
c$$$c     for conjugate gradient method (slsqr) of solution
c$$$      dimension jdx(ncoef)
c$$$      double precision  coef(ncoef)
c$$$      parameter(mlx=7,mly=7,mlz=7,
c$$$     ,          mnx=2**(mlx-1)+1,mny=2**(mly-1)+1,
c$$$     ,          mnz=2**(mlz-1)+1,
c$$$     ,          mvolu=mnx*mny*mnz,mdata=200000,mpath=2000,mst=1000,mevt=3000,
c$$$     ,          mvar=mvolu+mst+mevt)
c$$$      include 'lsqrcom1'
c$$$
c$$$      if (ncoef.gt.0) then
c$$$       nc_temp=0
c$$$       do i=1,ncoef
c$$$        if(coef(i).ne.0.d0)then
c$$$         n=max0(n,jdx(i))
c$$$        else
c$$$         nc_temp=nc_temp+1
c$$$        endif
c$$$       enddo
c$$$       m = m+1
c$$$       b(m) = rhs
c$$$       na(m) = ncoef-nc_temp
c$$$       do i=1,ncoef
c$$$        if(coef(i).ne.0.d0)then
c$$$         nel = nel+1
c$$$         if (nel.gt.big) then
c$$$          write (*,*) 'sorry, out of room in coeff. array for lldrow'
c$$$          write(*,*)big,nel
c$$$          stop
c$$$         end if
c$$$         ra(nel) = coef(i)
c$$$         ja(nel) = jdx(i)
c$$$        endif
c$$$       end do
c$$$      endif
c$$$      return
c$$$      end

      SUBROUTINE PIKSRT(N,ARR)
      double precision ARR(N)
      DO 12 J=2,N
        A=ARR(J)
        DO 11 I=J-1,1,-1
          IF(ARR(I).LE.A)GO TO 10
          ARR(I+1)=ARR(I)
11      CONTINUE
        I=0
10      ARR(I+1)=A
12    CONTINUE
      RETURN
      END


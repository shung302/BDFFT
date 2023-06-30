c     preplot_ray_map.f; plot velocity perturbations retrieved portrayed by ray G-file in map view
c
      parameter(mlx=8,mly=8,mlz=7,
     ,          mnx=2**(mlx-1)+1,mny=2**(mly-1)+1,mnz=2**(mlz-1)+1,
     ,          mface=mnx*mny,mvolu=mface*mnz,mdata=100000,mpath=1000)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      character*250 ia
      double precision lat0,lon0,lat1,lon1,dp0,dp1,thres
      double precision gtg(mvolu),x,y,z,
     ,                 image(mface),ima_max,ima_min
      double precision x0,y0,z0,dx,dy,dz
      double precision trans(3,3)
      common/meshb1/x0,y0,z0,dx,dy,dz,nx,ny,nz,nface,nvolu
      common/meshb2/trans
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
      include 'vmodels.inc'

      call mesh_car(lx,ly,lz,ndata)
      write(*,*)lx,ly,lz,ndata
      write(*,*)' GTG or ABSG file to be displayed:'
      read(*,'(a)')ia
      open(1,file=ia,err=73)
      read(1,*,err=73,end=73)(gtg(i),i=1,nvolu)
      close(1)
c      gtgmx=-1.d10
c      do i=1,nvolu
c         gtgmx=dmax1(dabs(gtg(i)),gtgmx)
c      enddo
c      do i=1,nvolu
c         gtg(i)=gtg(i)/gtgmx
c      enddo

      open(1,file='gtg.slice')
      write(*,*)' there are ',nz,' layers:'
      write(*,*)' which layer you want to inspect?'
      read(*,*)ilayer
      if(ilayer.gt.nz.or.ilayer.lt.1)ilayer=nz/2
      
      ima_max=-1.d6
      ima_min=1.d6
      iz0=(ilayer-1)*nface
c      zrad=z0+(ilayer-1)*dz-dz/2.0d0
c      zrad=z0+(ilayer-1)*dz+dz/2.0d0
      zrad = z0+(ilayer-1)*dz
      vmean=0.
      do i=1,ny
       do j=1,nx
        iface=(i-1)*nx+j
        image(iface)=gtg(iface+iz0)
        ima_min=dmin1(ima_min,image(iface))
        ima_max=dmax1(ima_max,image(iface))
        vmean=vmean+image(iface)
       enddo
      enddo
      vmean=vmean/dble(nface)
      write(6,*)'mean velocity perturbation at layer=',ilayer,vmean
      do i=1,ny
        lat0=y0+dble(i-1)*dy
       do j=1,nx
        lon0=x0+dble(j-1)*dx
        iface=(i-1)*nx+j
        call tran(lon1,lat1,lon0,lat0,x,y,z,-1)
        write(1,'(2f10.4,e12.4)')lon1,lat1,
     &  image(iface)
       end do
      end do

c      write(*,*)ima_min,ima_max
      stop
73    write(*,*)' ERROR reading '
      write(*,'(a)')ia
      end


c______________________________________________________________________
      subroutine mesh_car(lx,ly,lz,ndata)
      double precision x0,y0,z0,dx,dy,dz
      double precision trans(3,3)
      common/meshb1/x0,y0,z0,dx,dy,dz,nx,ny,nz,nface,nvolu
      common/meshb2/trans

      open(2,file='mesh.config',err=1)
      read(2,'(6i5,3f15.5)',err=1,end=1)lx,ly,lz,nx,ny,nz
      read(2,'(6f15.5)',err=1,end=1)x0,y0,z0,dx,dy,dz
      read(2,'(3f15.5)',err=1,end=1)((trans(i,j),j=1,3),i=1,3)
c      read(2,'(i8)',err=1,end=1)ndata
      nface=nx*ny
      nvolu=nface*nz
      close(2)
      return
1     stop' ERROR opening or reading mesh.config !'
      end


c______________________________________________________________________
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

c subroutine to find the range of r
c
      subroutine findr(ivin,r,ir)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn

c for flat and linear gradient model
      if (ivin.le.0.and.ivin.ge.-2) then
         if (r.le.rdis(1)) then
             ir=1
         elseif (r.gt.rdis(1).and.r.le.rdis(2)) then
             ir=2
         elseif (r.gt.rdis(2).and.r.le.rdis(11)) then
             ir=3
         else
             ir=4
         endif
c         if (r.eq.rcmb) ir=3
c PREM or 11 homogeneous layered models or IASP91
      elseif (ivin.eq.1.or.ivin.le.-3) then
c	 write(6,*)ivin,r,rdis(1),rdis(11)
         if (r.le.rdis(1)) then
         ir=1
         return
         endif

         if (r.gt.rdis(11)) then
         ir=12
         return
         endif
         do i=2,11
         if (r.gt.rdis(i-1).and.r.le.rdis(i)) then
            ir=i
            return
         endif
         enddo
c         if (r.eq.rcmb) ir=3
      elseif (ivin.eq.2) then
c for ak135 continent model
         if (r.le.rdis(1)) then
            ir=1
            return
         endif
         if (r.gt.rsurf) then
            ir=10
            return
         endif
         do i=2,9
         if (r.gt.rdis(i-1).and.r.le.rdis(i)) then
            ir=i
            return
         endif
         enddo
c         if (r.eq.rcmb) ir=3
      elseif (ivin.eq.3) then
c for ak135 spherical average model
         if (r.le.rdis(1)) then
            ir=1
            return
         endif
         if (r.gt.rsurf) then
            ir=13
            return
         endif
         do i=2,12
         if (r.gt.rdis(i-1).and.r.le.rdis(i)) then
            ir=i
            return
         endif
         enddo
c         if (r.eq.rcmb) ir=3
      endif
      return
      end



c    subroutine to find velocity at a given depth
c
      subroutine velo(ivin,ir,r,dvdr,vr)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,dvdr(2),vr(2)

      if (ivin.eq.-3) then
         call vlyrs(ir,r,dvdr,vr)
      elseif (ivin.eq.-2) then
         call vhomo(ir,r,dvdr,vr)
      elseif (ivin.eq.-1) then
         call vflat(ir,r,dvdr,vr)
      elseif (ivin.eq.0) then
         call vgrad(ir,r,dvdr,vr)
      elseif (ivin.eq.1) then
         call vprem(ir,r,dvdr,vr)
      elseif (ivin.eq.-4) then
         call viasp(ir,r,dvdr,vr)
      elseif (ivin.ge.2) then
         call vak135(ir,r,dvdr,vr)
      endif

      return
      end

c subroutine to get P and S wave velocity for a constant homo model
c
      subroutine vhomo(ir,r,dvdr,vr)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,dvdr(2),vr(2)
      COMMON /homo_model/ hc(3,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
      
      dvdr(1)=0.0d0
      dvdr(2)=0.0d0
      if (ir.gt.3.and.r.ge.rsurf) then
         vr(1)=hc(3,1)
         vr(2)=hc(3,2)
      else
         vr(1)=hc(ir,1)
         vr(2)=hc(ir,2)
      endif

      return
      end

c subroutine to get P and S wave velocity for a linear gradient model
c
      subroutine vflat(ir,r,dvdr,vr)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,dvdr(2),vr(2)
      COMMON /flat_model/ fc(2,3,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
      
      if (ir.gt.3.and.r.ge.rsurf) then
         vr(1)=fc(1,3,1)
         vr(2)=fc(1,3,2)
         dvdr(1)=fc(1,3,1)/rsurf-fc(2,3,1)
         dvdr(2)=fc(1,3,2)/rsurf-fc(2,3,2)
      else
         x=r/rsurf
         xlg=dlog(x)
         vr(1)=x*(fc(1,ir,1)-fc(2,ir,1)*rsurf*xlg)
         vr(2)=x*(fc(1,ir,2)-fc(2,ir,2)*rsurf*xlg)
         dvdr(1)=fc(1,ir,1)/rsurf-fc(2,ir,1)-fc(2,ir,1)*xlg
         dvdr(2)=fc(1,ir,2)/rsurf-fc(2,ir,2)-fc(2,ir,2)*xlg
      endif

      return
      end

c subroutine to get P and S wave velocity for a linear gradient model
c
      subroutine vgrad(ir,r,dvdr,vr)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,dvdr(2),vr(2)
      COMMON /grad_model/ gc(2,3,2)
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
      
      if (ir.gt.3.and.r.ge.rsurf) then
         x=0.0d0
         vr(1)=gc(1,3,1)+gc(2,3,1)*x
         vr(2)=gc(1,3,2)+gc(2,3,2)*x
         dvdr(1)=gc(2,3,1)
         dvdr(2)=gc(2,3,2)
      else
         x=r-rdis(ir)
         vr(1)=gc(1,ir,1)+gc(2,ir,1)*x
         vr(2)=gc(1,ir,2)+gc(2,ir,2)*x
         dvdr(1)=gc(2,ir,1)
         dvdr(2)=gc(2,ir,2)
      endif

      return
      end

c subroutine to get P and S wave velocity, and derivatives
c at a given radius r
c
      subroutine vlyrs(ir,r,dvdr,vr)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,dvdr(2),vr(2)
      INTEGER ir
      COMMON /lyrs_model/ yc(11,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
c isotropic model between 24.4 (6346.6) and 220 (6151) km
c
        dvdr(1)=0.0d0
        dvdr(2)=0.0d0
        if (ir.gt.11.and.r.ge.rsurf) then
        vr(1)=yc(11,1)
        vr(2)=yc(11,2)
        else
        vr(1)=yc(ir,1)
        vr(2)=yc(ir,2)
        endif

        return
        end


c subroutine to get P and S wave velocity, and derivatives
c at a given radius r
c
      subroutine vprem(ir,r,dvdr,vr)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,dvdr(2),vr(2)
      INTEGER ir
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
c isotropic model between 24.4 (6346.6) and 220 (6151) km
c
        x=r*rn
        if (ir.gt.11.and.r.ge.rsurf) then
        vr(1)=pc(1,11,1)+pc(2,11,1)+pc(3,11,1)+pc(4,11,1)
        dvdr(1)=rn*(pc(2,11,1)+2.0d0*pc(3,11,1)+3.0d0*pc(4,11,1))
        vr(2)=pc(1,11,2)+pc(2,11,2)+pc(3,11,2)+pc(4,11,2)
        dvdr(2)=rn*(pc(2,11,2)+2.0d0*pc(3,11,2)+3.0d0*pc(4,11,2))
        else
        vr(1)=pc(1,ir,1)+x*(pc(2,ir,1)+x*(pc(3,ir,1)+x*pc(4,ir,1)))
        dvdr(1)=rn*(pc(2,ir,1)+x*(2.0d0*pc(3,ir,1)+x*3.0d0*pc(4,ir,1)))
        vr(2)=pc(1,ir,2)+x*(pc(2,ir,2)+x*(pc(3,ir,2)+x*pc(4,ir,2)))
        dvdr(2)=rn*(pc(2,ir,2)+x*(2.0d0*pc(3,ir,2)+x*3.0d0*pc(4,ir,2)))
        endif

        return
        end

c subroutine to get P and S wave velocity, and derivatives
c at a given radius r
c
      subroutine viasp(ir,r,dvdr,vr)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,dvdr(2),vr(2)
      INTEGER ir
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn

c   Emiask returns model parameters for the IASPEI working model
c   (September 1990.1). 
c   Given non-dimensionalized radius x0, emiasp returns
c   non-dimensionalized density, ro, compressional velocity, vp, and
c   shear velocity, vs.  Non-dimensionalization is according to the
c   scheme of Gilbert in program EOS:  x0 by a (the radius of the
c   Earth), ro by robar (the mean density of the Earth), and velocity
c   by a*sqrt(pi*G*robar) (where G is the universal gravitational
c   constant.
c
      if ((ir.gt.11).and.(r.ge.rsurf)) then
      x=1.0d0
      vr(1)=apc(1,11,1)+apc(2,11,1)+apc(3,11,1)+apc(4,11,1)
      dvdr(1)=rn*(apc(2,11,1)+2.0d0*apc(3,11,1)+3.0d0*apc(4,11,1))
      vr(2)=apc(1,11,2)+apc(2,11,2)+apc(3,11,2)+apc(4,11,2)
      dvdr(2)=rn*(apc(2,11,2)+2.0d0*apc(3,11,2)+3.0d0*apc(4,11,2))
      else
      x=r*rn
      vr(1)=apc(1,ir,1)+x*(apc(2,ir,1)+x*(apc(3,ir,1)+x*apc(4,ir,1)))
      dvdr(1)=rn*(apc(2,ir,1)+x*(2.0d0*apc(3,ir,1)+x*3.0d0*apc(4,ir,1)))
      vr(2)=apc(1,ir,2)+x*(apc(2,ir,2)+x*(apc(3,ir,2)+x*apc(4,ir,2)))
      dvdr(2)=rn*(apc(2,ir,2)+x*(2.0d0*apc(3,ir,2)+x*3.0d0*apc(4,ir,2)))
      endif
      return
      end

c subroutine to get P and S wave velocity, and derivatives
c at a given radius r
c for ak135 model (136 layers)

      subroutine vak135(ir,r,dvdr,vr)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,dvdr(2),vr(2)
      INTEGER ir
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
c isotropic model between 24.4 (6346.6) and 220 (6151) km
c
      if (r.ge.rsurf) then
         vr(1)=akc(nlak,1)
         vr(2)=akc(nlak,2)
         dvdr(1)=akdc(nlak,1)
         dvdr(2)=akdc(nlak,2)
         return
      endif
      do i=1,nlak-1
         if (r.ge.akr(i).and.r.lt.akr(i+1)) then
            dvdr(1)=akdc(i,1)
            dvdr(2)=akdc(i,2)
            x=r-akr(i)
            vr(1)=akc(i,1)+x*akdc(i,1)
            vr(2)=akc(i,2)+x*akdc(i,2)
         endif
      enddo
      return
      end


c subroutine to get c,dc/dr,dc^2/dr^2
c
      subroutine velo1(ivin,iw,ir,r,c,dc,d2c)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,c,dc,d2c
      INTEGER iw,ir
      COMMON /homo_model/ hc(3,2)
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn

      iw0=iw
      if (ivin.eq.-3) then
         if (ir.eq.2) iw=1
         call vlyrs1(iw,ir,r,c,dc,d2c)
      elseif (ivin.eq.-2) then
         if (ir.eq.2) iw=1
         dc=0.0d0
         d2c=0.0d0
         if (r.gt.rsurf) ir=3
         c=hc(ir,iw)
      elseif (ivin.eq.-1) then
         if (ir.eq.2) iw=1
         call vflat1(iw,ir,r,c,dc,d2c)
      elseif (ivin.eq.0) then
         if (ir.eq.2) iw=1
         call vgrad1(iw,ir,r,c,dc,d2c)
      elseif (ivin.eq.1) then
         if (ir.eq.2) iw=1
         call vprem1(iw,ir,r,c,dc,d2c)
      elseif (ivin.eq.-4) then
         if (ir.eq.2) iw=1
         call viasp1(iw,ir,r,c,dc,d2c)
      elseif (ivin.eq.2) then
c ak135 model
         if (ir.eq.2) iw=1
         call vak1351(iw,ir,r,c,dc,d2c)
      elseif (ivin.eq.3) then
         if (ir.eq.2.or.ir.eq.11) iw=1
         call vak1351(iw,ir,r,c,dc,d2c)
      endif
      iw=iw0
      return
      end
         

c subroutine to get P and S wave velocity for a linear gradient model
c
      subroutine vflat1(iw,ir,r,c,dc,d2c)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 c,dc,d2c,r
      COMMON /flat_model/ fc(2,3,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
      
      if (ir.gt.3.and.r.ge.rsurf) then
         x=1.0d0
         c=fc(1,3,iw)
         dc=fc(1,3,iw)/rsurf-fc(2,3,iw)
         d2c=-fc(2,ir,iw)/r
      else
         x=r/rsurf
         xlg=dlog(x)
         c=x*(fc(1,ir,iw)-fc(2,ir,iw)*rsurf*xlg)
         dc=fc(1,ir,iw)/rsurf-fc(2,3,iw)-fc(2,ir,iw)*xlg
         d2c=-fc(2,ir,iw)/r
      endif

      return
      end

c
c subroutine to get P or S wave velocity, and
c 1st and 2nd derivatives at a given radius r,
c
      subroutine vgrad1(iw,ir,r,c,dc,d2c)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,c,dc,d2c
      INTEGER iw,ir
      COMMON /grad_model/ gc(2,3,2)
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn

      d2c=0.0d0
      if (ir.gt.3.and.r.ge.rsurf) then
         x=0.0d0
         c=gc(1,3,iw)+gc(2,3,iw)*x
         dc=gc(2,3,iw)
      else
         x=r-rdis(ir)
         c=gc(1,ir,iw)+gc(2,ir,iw)*x
         dc=gc(2,ir,iw)
      endif

      return
      end

c subroutine to get P or S wave velocity, and 
c 1st and 2nd derivatives at a given radius r,
c
      subroutine vlyrs1(iw,ir,r,c,dc,d2c)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,c,dc,d2c
      INTEGER iw,ir
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /lyrs_model/ yc(11,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
c isotropic model between 24.4 (6346.6) and 220 (6151) km
c

        dc=0.0d0
        d2c=0.0d0
        if (ir.gt.11.or.r.gt.rsurf) then
           c=yc(11,iw)
        else
           c=yc(ir,iw)
        endif

        return
        end


c subroutine to get P or S wave velocity, and 
c 1st and 2nd derivatives at a given radius r,
c
      subroutine vprem1(iw,ir,r,c,dc,d2c)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,c,dc,d2c
      INTEGER iw,ir
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
c isotropic model between 24.4 (6346.6) and 220 (6151) km
c

        rn2=rn*rn
        if (ir.gt.11.or.r.gt.rsurf) then
           c=pc(1,11,iw)+pc(2,11,iw)+pc(3,11,iw)+pc(4,11,iw)
           dc=rn*(pc(2,11,iw)+2.0d0*pc(3,11,iw)+3.0d0*pc(4,11,iw))
           d2c=rn2*(2.0d0*pc(3,11,iw)+6.0d0*pc(4,11,iw))
        else
           x=r*rn
           c=pc(1,ir,iw)+x*(pc(2,ir,iw)+x*(pc(3,ir,iw)+x*pc(4,ir,iw)))
           dc=rn*(pc(2,ir,iw)+x*(2.0d0*pc(3,ir,iw)+x*3.0d0*pc(4,ir,iw)))
           d2c=rn2*(2.0d0*pc(3,ir,iw)+x*6.0d0*pc(4,ir,iw))
        endif

        return
        end

c subroutine to get P and S wave velocity, and derivatives
c at a given radius r
c
      subroutine viasp1(iw,ir,r,c,dc,d2c)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,d2c,dc,c
      INTEGER ir,iw
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn
c
      rn2=rn*rn
      if ((ir.gt.11).and.(r.ge.rsurf)) then
      c=apc(1,11,iw)+apc(2,11,iw)+apc(3,11,iw)+apc(4,11,iw)
      dc=rn*(apc(2,11,iw)+2.0d0*apc(3,11,iw)+3.0d0*apc(4,11,iw))
      d2c=rn2*(2.0d0*apc(3,11,iw)+6.0d0*apc(4,11,iw))
      else
      x=r*rn
      c=apc(1,ir,iw)+x*(apc(2,ir,iw)+x*(apc(3,ir,iw)+x*apc(4,ir,iw)))
      dc=rn*(apc(2,ir,iw)+x*(2.0d0*apc(3,ir,iw)+x*3.0d0*apc(4,ir,iw)))
      d2c=rn2*(2.0d0*apc(3,ir,iw)+x*6.0d0*apc(4,ir,iw))
      endif
      return
      end

c subroutine to get P and S wave velocity, and derivatives
c at a given radius r
c for ak135 model (136 layers)

      subroutine vak1351(iw,ir,r,c,dc,d2c)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 r,d2c,dc,c
      INTEGER ir,iw
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,rn

c isotropic model between 24.4 (6346.6) and 220 (6151) km
c
      if (r.ge.rsurf) then
         dc=akdc(nlak,iw)
         d2c=0.0d0
         c=akc(nlak,iw)
         return
      endif
      do i=1,nlak-1
         if (r.ge.akr(i).and.r.lt.akr(i+1)) then
            x=r-akr(i)
            dc=akdc(i,iw)
            d2c=0.0d0
            c=akc(i,iw)+x*akdc(i,iw)
         endif
      enddo
      return
      end

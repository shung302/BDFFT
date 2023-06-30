c     imesh_car.f: build regional cubic Cartesian mesh
c                  by transform the center of the region
c                  to lon,lat=(0,0) first
      implicit double precision(a-h,o-z)
      parameter(rE=6371.d0,pi=3.14159265d0,rad=180.d0/pi)
      parameter(levelxmax=7,levelymax=7,levelzmax=7)
      dimension trans(3,3)

      open(1,file='mesh.plot')
      open(2,file='mesh.config')
c      write(6,*)'enter the levels in x(long),y(lati),z(radi)directions'
c      read(5,*)levelx,levely,levelz
      write(6,*)'enter the dimensions in x(long),y(lati),z(radi)directions'
      read(5,*)nevelx,nevely,nevelz
      levelx=nint(log(real(nevelx))/log(2.)+0.5)
      levely=nint(log(real(nevely))/log(2.)+0.5)
      levelz=nint(log(real(nevelz))/log(2.)+0.5)
      write(6,*)'level in x,y,z directions=',levelx,levely,levelz 
      write(6,*)'enter the center position in lon, lat directions'
      read(5,*)onl00,atl00
      print *,'enter half-width lengths in lon, lat, and depths at base and top'
      read(5,*)onlL2,atlL2,dphL,dphL0
      onl0=onl00/rad
      atl0=atl00/rad
      if(atl00.lt.0.d0)then
       atlp=atl00+90.d0
       onlp=onl00
      else
       atlp=90.d0-atl00
       onlp=onl00+180.d0
       if(onlp.ge.360.d0)onlp=onlp-360.d0
      endif
      onlp=onlp/rad
      atlp=atlp/rad
      trans(1,1)=dcos(atl0)*dcos(onl0)
      trans(1,2)=dcos(atl0)*dsin(onl0)
      trans(1,3)=dsin(atl0)
      trans(3,1)=dcos(atlp)*dcos(onlp)
      trans(3,2)=dcos(atlp)*dsin(onlp)
      trans(3,3)=dsin(atlp)
      trans(2,1)=trans(3,2)*trans(1,3)-trans(3,3)*trans(1,2)
      trans(2,2)=trans(3,3)*trans(1,1)-trans(3,1)*trans(1,3)
      trans(2,3)=trans(3,1)*trans(1,2)-trans(3,2)*trans(1,1)
c      nx=2**(levelx-1)
c      ny=2**(levely-1)
c      nz=2**(levelz-1)
      nx=nevelx
      ny=nevely
      nz=nevelz
      dx=onlL2*2.d0/dble(nx-1)
      dy=atlL2*2.d0/dble(ny-1)
      dz=(dphL-dphL0)/dble(nz-1)
      write(2,'(6i5,3f15.5)')levelx,levely,levelz,nx,ny,nz
      write(2,'(6f15.5)')-onlL2,-atlL2,rE-dphL,dx,dy,dz
      write(2,'(3f15.5)')((trans(i,j),j=1,3),i=1,3)
c      write(2,'(i6)')nz
c      do i=1,nz
c       z=rE-dble(nz-i)*dz
c       write(2,'(f15.3)')z
c      enddo
c      write(2,'(2i6)')nx,ny
c      do j=1,ny
c       atl1=-atlL2+dble(j-1)*dy
c       atl2=atl1/rad
c       do i=1,nx
c        onl1=-onlL2+dble(i-1)*dx
c        onl2=onl1/rad
c        xp=dcos(atl2)*dcos(onl2)
c        yp=dcos(atl2)*dsin(onl2)
c        zp=dsin(atl2)
c        x=trans(1,1)*xp+trans(2,1)*yp+trans(3,1)*zp
c        y=trans(1,2)*xp+trans(2,2)*yp+trans(3,2)*zp
c        z=trans(1,3)*xp+trans(2,3)*yp+trans(3,3)*zp
c        r=dsqrt(x*x+y*y+z*z)
c        atdg=dasin(z/r)*rad
c        andg=datan2(y,x)*rad
c        if(andg.lt.0.d0)andg=andg+360.d0
c        write(2,'(4f15.5)')onl1,atl1,andg,atdg
c       end do
c      end do
      do j=1,ny
       atl1=-atlL2+dble(j-1)*dy
       atl2=atl1/rad
       write(1,'(5h> y= ,i3)')j
       onl1=-onlL2
       onl2=onl1/rad
       xp=dcos(atl2)*dcos(onl2)
       yp=dcos(atl2)*dsin(onl2)
       zp=dsin(atl2)
       x=trans(1,1)*xp+trans(2,1)*yp+trans(3,1)*zp
       y=trans(1,2)*xp+trans(2,2)*yp+trans(3,2)*zp
       z=trans(1,3)*xp+trans(2,3)*yp+trans(3,3)*zp
       r=dsqrt(x*x+y*y+z*z)
       atdg=dasin(z/r)*rad
       andg=datan2(y,x)*rad
       if(andg.lt.0.d0)andg=andg+360.d0
       write(1,'(2f15.5)')andg,atdg
c       write(*,*)onl1,atl1,andg,atdg
       onl1=-onlL2+dble(nx-1)*dx
       onl2=onl1/rad
       xp=dcos(atl2)*dcos(onl2)
       yp=dcos(atl2)*dsin(onl2)
       zp=dsin(atl2)
       x=trans(1,1)*xp+trans(2,1)*yp+trans(3,1)*zp
       y=trans(1,2)*xp+trans(2,2)*yp+trans(3,2)*zp
       z=trans(1,3)*xp+trans(2,3)*yp+trans(3,3)*zp
       r=dsqrt(x*x+y*y+z*z)
       atdg=dasin(z/r)*rad
       andg=datan2(y,x)*rad
       if(andg.lt.0.d0)andg=andg+360.d0
       write(1,'(2f15.5)')andg,atdg
c       write(*,*)onl1,atl1,andg,atdg
      end do
      do j=1,nx
       onl1=-onlL2+dble(j-1)*dx
       onl2=onl1/rad
       write(1,'(5h> x= ,i3)')j
       atl1=-atlL2
       atl2=atl1/rad
       xp=dcos(atl2)*dcos(onl2)
       yp=dcos(atl2)*dsin(onl2)
       zp=dsin(atl2)
       x=trans(1,1)*xp+trans(2,1)*yp+trans(3,1)*zp
       y=trans(1,2)*xp+trans(2,2)*yp+trans(3,2)*zp
       z=trans(1,3)*xp+trans(2,3)*yp+trans(3,3)*zp
       r=dsqrt(x*x+y*y+z*z)
       atdg=dasin(z/r)*rad
       andg=datan2(y,x)*rad
       if(andg.lt.0.d0)andg=andg+360.d0
c       write(*,*)onl1,atl1,andg,atdg
       write(1,'(2f15.5)')andg,atdg
       atl1=-atlL2+dble(ny-1)*dy
       atl2=atl1/rad
       xp=dcos(atl2)*dcos(onl2)
       yp=dcos(atl2)*dsin(onl2)
       zp=dsin(atl2)
       x=trans(1,1)*xp+trans(2,1)*yp+trans(3,1)*zp
       y=trans(1,2)*xp+trans(2,2)*yp+trans(3,2)*zp
       z=trans(1,3)*xp+trans(2,3)*yp+trans(3,3)*zp
       r=dsqrt(x*x+y*y+z*z)
       atdg=dasin(z/r)*rad
       andg=datan2(y,x)*rad
       if(andg.lt.0.d0)andg=andg+360.d0
       write(1,'(2f15.5)')andg,atdg
c       write(*,*)onl1,atl1,andg,atdg
      end do
      stop
      end




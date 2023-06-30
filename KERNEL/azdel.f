c subroutine packages to calculate delta, azimuth, and bounce points, etc...
c from Guy Masters (geographic coords convert to gecentric coords. when calculating
c epicentral distance and azimuth
c
      subroutine stadis(colat,colon,scolat,scolon,del,az)
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
c      real*8 co,si,caz,saz
c      data rad/57.29578/
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

c  first do eq coords.
c      t0=colat/rad
      t0=geocen(colat*rpd)    
      p0=colon*rpd
      c0=dcos(t0)
      s0=dsin(t0)
c  now do station coords.
c     t2=scolat*rpd
      t2=geocen(scolat*rpd)           
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
      if(az.lt.0.0d0) az=360.0d0 + az  
      return
      end
c
c subroutine packages to calculate delta, azimuth, and bounce points, etc...
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
c      real*8 co,si,caz,saz
c      data rad/57.29578/
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
c      write(6,*)'rpd and dpr at stadis=',rpd,dpr

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
      real*8 function geogrf(arg)
c input:
c   arg    = geocentric colatitude (radians)
c output:
c   geogrf = geographic colatitude (radians
c (n.b. fac=(1-f)**2)
c
      implicit real*8 (a-h,o-z)
      data fac/0.993305621334896d0/
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

c      geogrf=pi2-datan2(cos(arg),fac*dmax1(1.0d-30,dsin(arg))) 
      geogrf=pi2-datan(cos(arg)/(fac*dmax1(1.0d-30,dsin(arg)))) 
      return                       
      end
c
      subroutine givloc(colat,colon,del,az,t1,p1)
c input:
c   colat,colon = source colatitude and colongitude 
c   del         = distance of point in degrees
c   az          = azimuth of point from source
c output:
c   t1 = point latitude  ( + = N, - = S)
c   p1 = point longitude ( + = E, - = W)
      implicit real*8 (a-h,o-z)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
c      data rad/57.29578/
      delr=del*rpd
      azr=az*rpd
      t0=geocen(colat*rpd)    
      ctheta=dsin(delr)*dsin(t0)*dcos(azr) + dcos(t0)*dcos(delr)
      t1=dacos(ctheta)
      if (t0.eq.0.0d0) then
        p1=az
      elseif (t1.eq.0.0d0) then
        p1=0.0d0
      else
        sphi=dsin(delr)*dsin(azr)/dsin(t1)
        cphi=(dcos(delr) - dcos(t0)*ctheta)/(dsin(t0)*dsin(t1))
        p1=colon + datan2(sphi,cphi)*dpr
      endif
      t1=90.0-geogrf(t1)*dpr      
      if (p1.gt.360.0d0) p1 = p1 - 360.0d0 
      if (p1.gt.180.0d0) p1 = p1 - 360.0d0 
      return
      end

      real*8 function geocen(arg)
c input:
c   arg    = geographic colatitude (radians)
c output:
c   geocen = geocentric colatitude (radians)
c (n.b. fac=(1-f)**2)
      implicit real*8 (a-h,o-z)
      data fac/0.993305621334896d0/
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
c     data pi2,fac/1.570796326794895,0.993305621334896/
c      geocen=pi2-datan2(fac*dcos(arg),dmax1(1.0d-30,dsin(arg)))
      geocen=pi2-datan(fac*dcos(arg)/(dmax1(1.0d-30,dsin(arg))))
      return
      end

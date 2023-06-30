c calculate velocity and dv/dr at a given r
c see Dahlen and Tromp (1998)
c and derivatives of d(r,i,T)/d(phi) (equations 15.253, 15.255)
c and derivatives of d(dr/di0), (di/di0))/d(phi) (equations 15.257 and 258)

      subroutine derivs(ivin,ips,n,x,y,dydx,ir,iw,c,dc,d2c)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 y(n),x,dydx(n)
      integer ips
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /grad_model/ gc(2,3,2)
      COMMON /flat_model/ fc(2,3,2)

      iw=ips
      if (ir.eq.2) iw=1
      call velo1(ivin,iw,ir,y(1),c,dc,d2c)
      sii=dsin(y(2))
      coi=dcos(y(2))
      ri=1.0d0/y(1)
      dydx(1)=coi
      dydx(2)=sii*(dc/c-ri)
      dydx(3)=sii*ri

      return
      end

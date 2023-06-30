c calculate velocity and dv/dr at a given r
c and derivatives of d(r,i,T)/d(phi)

      subroutine derivs_pq(ivin,n,rayp,x,y,dydx,ir,iw)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      REAL*8 y(n),x,dydx(n),rayp,c,dc,d2c
      
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /grad_model/ gc(2,3,2)
      COMMON /flat_model/ fc(2,3,2)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak

      r=y(9)
      call velo1(ivin,iw,ir,r,c,dc,d2c)
      ri=1.0d0/r
      r2=r*r
      tai=dtan(y(10))
      cti=1.0d0/tai
      cti2=cti*cti
      c3=c*c*c
      dydx(1)=-y(3)*rayp*(d2c+dc*cti2*ri)/c
      dydx(2)=-y(4)*r*dc/(c3*rayp)
      dydx(3)=y(1)*r2/rayp
      dydx(4)=y(2)*r2/rayp 
      dydx(5)=-y(7)*rayp*(d2c+dc*cti2*ri)/c
      dydx(6)=-y(8)*r*dc/(c3*rayp)
      dydx(7)=y(5)*r2/rayp
      dydx(8)=y(6)*r2/rayp 
      dydx(9)=r*cti
      dydx(10)=r*dc/c-1.0d0

      return
      end

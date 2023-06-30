      SUBROUTINE driver_pq_rk4(ivin,ir,iw,n,rayp,y,dydx,x,h,yout,
     & derivs_pq)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INTEGER n,i
      REAL*8 h,x,dydx(n),y(n),yout(n)
      EXTERNAL derivs_pq
      REAL*8 h6,hh,xh,dym(NMX),dyt(NMX),yt(NMX)
      hh=h*0.5d0
      h6=h/6.0d0
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call derivs_pq(ivin,n,rayp,xh,yt,dyt,ir,iw)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call derivs_pq(ivin,n,rayp,xh,yt,dym,ir,iw)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call derivs_pq(ivin,n,rayp,x+h,yt,dyt,ir,iw)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
14    continue
      return
      END

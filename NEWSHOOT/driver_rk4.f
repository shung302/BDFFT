      SUBROUTINE driver_rk4(ivin,ips,n,y,dydx,x,h,yout,derivs)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INTEGER n,ips,i
      REAL*8 h,x,dydx(n),y(n),yout(n)
      EXTERNAL derivs
      REAL*8 h6,hh,xh,dym(NMX),dyt(NMX),yt(NMX)
      hh=h*0.5d0
      h6=h/6.0d0
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      call findr(ivin,yt(1),ir)
      call derivs(ivin,ips,n,xh,yt,dyt,ir,iw,vr,dvr,d2vr)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      call findr(ivin,yt(1),ir)
      call derivs(ivin,ips,n,xh,yt,dym,ir,iw,vr,dvr,d2vr)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      call findr(ivin,yt(1),ir)
      call derivs(ivin,ips,n,x+h,yt,dyt,ir,iw,vr,dvr,d2vr)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
14    continue
      return
      END

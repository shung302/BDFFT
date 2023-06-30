c subroutine to calculate integration over omega to get kernel for
c a given scatterer point
c
c
      subroutine dintw_kntbl_df(istf,a2,wb,we,wx,ww,nn)
      implicit doubleprecision (a-h,o-z)
      include 'parameters.inc'
      DIMENSION wx(*),ww(*),a2(*)
      complex*16 z105,z30,z60,z120,bigd2sh,bigd2sv
      COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),
     &   rkn(NKMX,MSMX,20)
      COMMON /ttkerndf/ nmu,dmu,rmutbl(NKMX),rkndf(NKMX,NKMX,20)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /bphi/ tdiff,sig
      COMMON /airyfcn/ dai(10),fdai(10),ai(10),fai(10),af2(10),daf2(10),
     & bigd2sh(NKMX,NTMX),bigd2sv(NKMX,NTMX)
      COMMON /zconsts/ z105,z30,z60,z120
      COMMON /rconsts/ oneth,twoth,c30,s30


c    real part of integration in the nominator
      do 31 im=1,nmu
         rmu=rmutbl(im)
      do 21 it=1,ntt
         tdiffo=ttdf(it)
         rkndf(it,im,istf)=0.0d0
      do 11 i=1,nn
         cwt=dcos(wx(i)*tdiffo)
         swt=dsin(wx(i)*tdiffo)
         fdf=dreal(bigd2sh(im,i))*cwt+dimag(bigd2sh(im,i))*swt
         rkndf(it,im,istf)=rkndf(it,im,istf)+ww(i)*fdf
11    continue
      rkndf(it,im,istf)=rkndf(it,im,istf)/(pi*rkd(istf))
      if (im.eq.1) then
      write(13,*)rmu,tdiffo,rkndf(it,im,istf)
      endif
21    continue
31    continue
      return
      end

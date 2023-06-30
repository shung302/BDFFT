c subroutine to calculate integration over omega to get kernel for
c a given scatterer point
c
c
      subroutine dintw_kntbl(istf,a2,wb,we,wx,ww,nn)
      implicit doubleprecision (a-h,o-z)
      include 'parameters.inc'
      DIMENSION wx(*),ww(*),a2(*)
      COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),
     &  rkn(NKMX,MSMX,20)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /bphi/ tdiff,sig

c      open(51,file='rkntbl.out')
c      write(51,'(''table of rkn/rkd'')')
c      close(51)
c      open(51,file='rkntbl.out',status='old',access='append')
      do 31 im=1,ms+1
         sig=dble(im-1)
      do 21 it=1,ntt
         rkn(it,im,istf)=0.0d0
         tdiff=ttdf(it)
      do 11 j=1,nn
      rkn(it,im,istf)=rkn(it,im,istf)+ww(j)*dfint(wx(j),a2(j))
11    continue
      rkn(it,im,istf)=rkn(it,im,istf)/(pi*rkd(istf))
c      write(51,111)im,tdiff,rkn(it,im,istf)
21    continue
31    continue
111   format(i4,2e15.6)
      close(51)
      return
      end

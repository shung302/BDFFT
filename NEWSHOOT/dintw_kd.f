      SUBROUTINE dintw_kd(istf,a2,a,b,wx,ww,nn)
      implicit doubleprecision (a-h,o-z)
      include 'parameters.inc'
      dimension ww(*),wx(*),a2(*)
      COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),
     &   rkn(NKMX,MSMX,20)

      INTEGER j
      ss=0.0d0
      do 11 j=1,nn
        ss=ss+ww(j)*a2(j)
11    continue
      rkd(istf)=ss
      return
      END

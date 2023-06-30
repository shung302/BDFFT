      subroutine getv_cmdl(frad,fcol,flong,amap)
      implicit real*4 (a-h,o-z)
      parameter (MXLH=36,MXLENY=(MXLH+1)**2,MXDEP=21,MXSPL=MXDEP+3,
     + MXMSZ=MXSPL*MXLENY,MAXP=1024)
      common /caltechblk/ wk1(MXLENY),wk2(MXLENY),wk3(MXLENY),
     + d0(MXLENY),spl(MXDEP),x(MXMSZ)
      common /caltechnum/ lmx,nsmx,ndmx,ndmn,natd

      character*80 mfl
      real*8 frad,fcol,flong,amap

      xlt=90.0d0-fcol
      xln=flong

c logitude between -180 and 180
      xlno=xln
      if (xln.gt.180.) xln=xln-360.

c    calculate Y_k(xlt,xln) at random depth level
      call ylm(xlt,xln,lmx,d0,wk1,wk2,wk3)

c    calculate spline coefficients
      nbeg=max(4,ndmn)
      nend=ndmx-3
      do i=nbeg,ndmx
       ind=(i-1)*natd+1 
       spl(i-3)=sdot(natd,d0,1,x(ind),1)
      enddo

      rcmb=3480.
      rmoho=6346.
      rearth=6371.
      r=frad
      xd=-1.+2.*(r-rcmb)/(rmoho-rcmb)
      dv=0.
      do ip=nbeg-3,nend
       dv=dv+spl(ip)*splh(ip-1,xd)
      enddo
      amap=dble(dv)
c      write(6,112) xlt,xlno,r,dv*100.
c112   format(f6.2,f9.2,f9.1,f9.2)

      return
      end

c -----------------------------------------------------

      subroutine ylm(xlat,xlon,lmax,y,wk1,wk2,wk3)
c
      complex temp,fac,dfac
      dimension wk1(1),wk2(1),wk3(1),y(1)
c
c     wk1,wk2,wk3 should be dimensioned at least (lmax+1)*4
c
      data radian/57.2957795/    ! 360./2pi
c
c     transform to spherical coordinates
      theta=(90.-xlat)/radian
      phi=xlon/radian
c
c    loop over l values
      ind=0
      lm1=lmax+1
      do 10 il1=1,lm1
      l=il1-1
      call legndr(theta,l,l,wk1,wk2,wk3)
c
      fac=(1.,0.)
      dfac=cexp(cmplx(0.,phi))
c
c    loop over m values
      do 20 im=1,il1
      temp=fac*cmplx(wk1(im),0.)
      ind=ind+1
      y(ind)=real(temp)
      if(im.eq.1) goto 20
      ind=ind+1
      y(ind)=aimag(temp)
   20 fac=fac*dfac   ! calculates exp(im phi)
c
   10 continue
      return
      end

c --------------------------------------------------------------------
      SUBROUTINE LEGNDR(THETA,L,M,X,XP,XCOSEC)
      DIMENSION X(*),XP(*),XCOSEC(*)
      DOUBLE PRECISION SMALL,SUM,COMPAR,CT,ST,FCT,COT,FPI,X1,X2,X3,
     1F1,F2,XM,TH,DFLOAT
      DATA FPI/12.56637062D0/
      DFLOAT(I)=FLOAT(I)
      SUM=0.D0
      LP1=L+1
      TH=THETA
      CT=DCOS(TH)
      ST=DSIN(TH)
      MP1=M+1
      FCT=DSQRT(DFLOAT(2*L+1)/FPI)
      SFL3=SQRT(FLOAT(L*(L+1)))
      COMPAR=DFLOAT(2*L+1)/FPI
      DSFL3=SFL3
      SMALL=1.D-16*COMPAR
      DO 1 I=1,MP1
      X(I)=0.
      XCOSEC(I)=0.
    1 XP(I)=0.
      IF(L.GT.1.AND.ABS(THETA).GT.1.E-5) GO TO 3
      X(1)=FCT
      IF(L.EQ.0) RETURN
      X(1)=CT*FCT
      X(2)=-ST*FCT/DSFL3
      XP(1)=-ST*FCT
      XP(2)=-.5D0*CT*FCT*DSFL3
      IF(ABS(THETA).LT.1.E-5) XCOSEC(2)=XP(2)
      IF(ABS(THETA).GE.1.E-5) XCOSEC(2)=X(2)/ST
      RETURN
    3 X1=1.D0
      X2=CT
      DO 4 I=2,L
      X3=(DFLOAT(2*I-1)*CT*X2-DFLOAT(I-1)*X1)/DFLOAT(I)
      X1=X2
    4 X2=X3
      COT=CT/ST
      COSEC=1./ST
      X3=X2*FCT
      X2=DFLOAT(L)*(X1-CT*X2)*FCT/ST
      X(1)=X3
      X(2)=X2
      SUM=X3*X3
      XP(1)=-X2
      XP(2)=DFLOAT(L*(L+1))*X3-COT*X2
      X(2)=-X(2)/SFL3
      XCOSEC(2)=X(2)*COSEC
      XP(2)=-XP(2)/SFL3
      SUM=SUM+2.D0*X(2)*X(2)
      IF(SUM-COMPAR.GT.SMALL) RETURN
      X1=X3
      X2=-X2/DSQRT(DFLOAT(L*(L+1)))
      DO 5 I=3,MP1
      K=I-1
      F1=DSQRT(DFLOAT((L+I-1)*(L-I+2)))
      F2=DSQRT(DFLOAT((L+I-2)*(L-I+3)))
      XM=K
      X3=-(2.D0*COT*(XM-1.D0)*X2+F2*X1)/F1
      SUM=SUM+2.D0*X3*X3
      IF(SUM-COMPAR.GT.SMALL.AND.I.NE.LP1) RETURN
      X(I)=X3
      XCOSEC(I)=X(I)*COSEC
      X1=X2
      XP(I)=-(F1*X2+XM*COT*X3)
    5 X2=X3
      RETURN
      END


c ------------------------------------------------------------------------------

      subroutine wint2ch(int,ch,lrch)

      character*2 ch

      if(int.lt.10) then
       lrch=1
       write(ch,'(i1)') int
      else if(int.ge.10.and.int.lt.100) then
       lrch=2
       write(ch,'(i2)') int
      else
       stop 'wint2ch; int.gt.99'
      endif

      end

c ------------------------------------------------
      function splh(ind,x)

      parameter (MXKNT=21)
      common/splhprm/spknt(MXKNT),qq0(MXKNT,MXKNT),qq(3,MXKNT,MXKNT)
      if(x.gt.1.or.x.lt.-1) then
       splh=0.
       goto 10
      endif

      splh=rsple(1,MXKNT,spknt(1),qq0(1,MXKNT-ind),qq(1,1,MXKNT-ind),x)

10    continue
      return
      end

c -------------------------------------------------

      FUNCTION RSPLE(I1,I2,X,Y,Q,S)
C
C C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
C
C   RSPLE RETURNS THE VALUE OF THE FUNCTION Y(X) EVALUATED AT POINT S
C   USING THE CUBIC SPLINE COEFFICIENTS COMPUTED BY RSPLN AND SAVED IN
C   Q.  IF S IS OUTSIDE THE INTERVAL (X(I1),X(I2)) RSPLE EXTRAPOLATES
C   USING THE FIRST OR LAST INTERPOLATION POLYNOMIAL.  THE ARRAYS MUST
C   BE DIMENSIONED AT LEAST - X(I2), Y(I2), AND Q(3,I2).
C
C                                                     -RPB
      DIMENSION X(*),Y(*),Q(3,*)
      DATA I/1/
      II=I2-1
C   GUARANTEE I WITHIN BOUNDS.
      I=MAX0(I,I1)
      I=MIN0(I,II)
C   SEE IF X IS INCREASING OR DECREASING.
      IF(X(I2)-X(I1))1,2,2
C   X IS DECREASING.  CHANGE I AS NECESSARY.
 1    IF(S-X(I))3,3,4
 4    I=I-1
      IF(I-I1)11,6,1
 3    IF(S-X(I+1))5,6,6
 5    I=I+1
      IF(I-II)3,6,7
C   X IS INCREASING.  CHANGE I AS NECESSARY.
 2    IF(S-X(I+1))8,8,9
 9    I=I+1
      IF(I-II)2,6,7
 8    IF(S-X(I))10,6,6
 10   I=I-1
      IF(I-I1)11,6,8
 7    I=II
      GO TO 6
 11   I=I1
C   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
 6    H=S-X(I)
      RSPLE=Y(I)+H*(Q(1,I)+H*(Q(2,I)+H*Q(3,I)))
      RETURN
      END

c ------------------------------------------------------------

      real function sdot(n,sx,incx,sy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end


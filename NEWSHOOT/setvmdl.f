c subroutine to determine which global models to be used for
c integration of kernels
c
      subroutine setvmdl(ips,mdl)
      implicit real (a-h,o-z)
      common/blocks/r1(30),r2(30),t1(180),t2(180),
     +    pmn(50000),pmx(50000),nshell,mshell,mlat(181),nlat
      common/blk$/block(50000),nblk
      parameter (MXLH=36,MXLENY=(MXLH+1)**2,MXDEP=21,MXSPL=MXDEP+3,
     + MXMSZ=MXSPL*MXLENY,MAXP=1024)
      common /caltechblk/ wk1(MXLENY),wk2(MXLENY),wk3(MXLENY),
     + d0(MXLENY),spl(MXDEP),x(MXMSZ)
      common /caltechnum/ lmx,nsmx,ndmx,ndmn,natd

      integer mdl
      character*40 mdlname*10,mdlblk

      goto (10,20,30,40,45,50) mdl
10    mdlname='SB4L18'
      mdlblk='SB4L18/block.model'
      call setup_smdl(mdlname,mdlblk)
      go to 50
20    mdlname='SB10L18'
      if (ips.eq.1) then
         mdlblk='SB10L18/block.model.p'
      elseif (ips.eq.2) then
         mdlblk='SB10L18/block.model.s'
      else
         mdlblk='SB10L18/block.model.c'
      endif
      call setup_smdl(mdlname,mdlblk)
      go to 50
30    if (ips.eq.1) then
         mdlname='P20RTS.sph'
      else
         mdlname='S20RTS.sph'
      endif
      call setup_cmdl(mdlname)
      go to 50
40    mdlname='SB4L18'
      mdlblk='SB4L18/rough.model'
      call setup_smdl(mdlname,mdlblk)
      go to 50
45    mdlname='SB4L18'
      mdlblk='SB4L18/erough.model'
      call setup_smdl(mdlname,mdlblk)
50    return
      end

cc ########################################################################cc
cc set up common blocks for scripps block models
cc ########################################################################cc

      subroutine setup_smdl(mdlname,mdlblk)
      character*128 filename,mdlname*10,mdlblk*40
      common/blocks/r1(30),r2(30),t1(180),t2(180),
     +    pmn(50000),pmx(50000),nshell,mshell,mlat(181),nlat
      common/blk$/block(50000),nblk

      lf=index(mdlname,' ')-1
      call blks(nblk,mdlname,lf)
      lf=index(mdlblk,' ')-1
      filename='/home/shung/GlobalModels/scripps/'//mdlblk(1:lf)
      open(9,file=filename, form='unformatted')
      read(9)n,phibar
c      print*,'n,nblk,phibar',n,nblk,phibar
      read(9) iter,(block(i),i=1,n),phibar,r
      close(9)
      do 20 i=1,n
        if(abs(block(i)).lt.1.e-15) block(i)=1.
   20   continue
      do 30 i=n+1,nblk
   30   block(i)=1.
      return
      end

      subroutine blks(nblk,fname,lf)
c*** routine to establish block parameters
      common/blocks/r1(30),r2(30),t1(180),t2(180),
     +    pmn(50000),pmx(50000),nshell,mshell,mlat(181),nlat
      character*256 filename,fname*10
      data rad/57.29578/
      filename='/home/shung/GlobalModels/scripps/'//fname(1:lf)//'/block.desc'
      open(12,file=filename)
c bsize=block size at equator
      read(12,*) bsize
      nlat=nint(180./bsize)
      bsize=180./nlat
c nshells = # shells
      read(12,*) nshell
      ib=0
      do 10 n=1,nshell
c rad1 and rad2 are lower and upper radii of shell
        read(12,*) rad1,rad2
        tmax=0.
        r1(n)=rad1
        r2(n)=rad2
        mlat(1)=0
        do 30 ii=1,nlat
        tmin=tmax
        tmax=tmin+bsize
        th=0.5*(tmin+tmax)/rad
        s1=sin(th)
        mlon=max(nint(360./bsize*s1),1)
        hsize=360./mlon
        pmax=0.
        do 35 jj=1,mlon
          pmin=pmax
          pmax=pmin+hsize
          ib=ib+1
          pmn(ib)=pmin
          pmx(ib)=pmax
   35     continue
        if(n.eq.1) then
           mlat(ii+1)=ib
           t1(ii)=tmin
           t2(ii)=tmax
        end if
   30   continue
      if(n.eq.1) mshell=ib
   10 continue
      nblk=ib
      close(12)
      return
      end

cc###############################################################cc
cc set up common block for Caltech spherical harmonics model
cc###############################################################cc

      subroutine setup_cmdl(mdlname)
      implicit real*4 (a-h,o-z)
      character*80 mfl,mdlname*10
      parameter (MXLH=36,MXLENY=(MXLH+1)**2,MXDEP=21,MXSPL=MXDEP+3,
     + MXMSZ=MXSPL*MXLENY,MAXP=1024)
      common /caltechblk/ wk1(MXLENY),wk2(MXLENY),wk3(MXLENY),
     + d0(MXLENY),spl(MXDEP),x(MXMSZ)
      common /caltechnum/ lmx,nsmx,ndmx,ndmn,natd

c    read sph model
      mfl='/home/shung/GlobalModels/caltech/'//mdlname(1:10)
      call rmod(mfl,x,lmx,nsmx,ndmx,ndmn)
      if(lmx.gt.MXLH) stop'lmx.gt.MXLH'
      if(ndmx.gt.MXSPL) stop'ndmx.gt.MXSPL'
      natd=(lmx+1)**2

c    Calculate the spline basis functions at a regular grid
      call splhsetup()

      return
      end

c------------------------------------------------
c read model coefficients (ASCII format)
      subroutine rmod(infl,x,lmx,nsmx,ndmx,ndmn)
      character*80 infl
      dimension x(*)

      open(10,file=infl,status='old')
      read(10,*) lmx,dum,nsmx

c--   JR:
c--   Hardwired are the following four parameters.
c--   See also Hendrik's original code.
c--

        lmx = 20
        nsmx = 24
        ndmx = 24
        ndmn = 4

      natd=(lmx+1)**2
      ind=(ndmn-1)*natd+1
      do i=ndmn,ndmx
       do j=0,lmx
        ind1=ind+2*j
        read(10,'(11e12.4)',end=100)(x(k),k=ind,ind1)
        ind=ind1+1
       Enddo
      Enddo

      goto 200

 100  stop 'incompatible sph header'
 200  continue
      close(10)

      end

c-------------------------------------
      subroutine splhsetup()
      parameter (MXKNT=21)
      common/splhprm/spknt(MXKNT),qq0(MXKNT,MXKNT),qq(3,MXKNT,MXKNT)
      data spknt/
     1   -1.00000,-0.78631,-0.59207,-0.41550,-0.25499
     1  ,-0.10909, 0.02353, 0.14409, 0.25367, 0.35329
     1  , 0.44384, 0.52615, 0.60097, 0.66899, 0.73081
     1  , 0.78701, 0.83810, 0.88454, 0.92675, 0.96512
     1  , 1.00000
     1    /
      dimension qqwk(3,MXKNT)
      do i=1,MXKNT
        do j=1,MXKNT
          if(i.eq.j) then
            qq0(j,i)=1.
          else
            qq0(j,i)=0.
          endif
        enddo
      enddo
      do i=1,MXKNT
        call rspln(1,MXKNT,spknt(1),qq0(1,i),qq(1,1,i),qqwk(1,1))
      enddo
      return
      end

c --------------------------------------------------------

      SUBROUTINE RSPLN(I1,I2,X,Y,Q,F)
C
C C$C$C$C$C$ CALLS ONLY LIBRARY ROUTINES C$C$C$C$C$
C
C   SUBROUTINE RSPLN COMPUTES CUBIC SPLINE INTERPOLATION COEFFICIENTS
C   FOR Y(X) BETWEEN GRID POINTS I1 AND I2 SAVING THEM IN Q.  THE
C   INTERPOLATION IS CONTINUOUS WITH CONTINUOUS FIRST AND SECOND 
C   DERIVITIVES.  IT AGREES EXACTLY WITH Y AT GRID POINTS AND WITH THE
C   THREE POINT FIRST DERIVITIVES AT BOTH END POINTS (I1 AND I2).
C   X MUST BE MONOTONIC BUT IF TWO SUCCESSIVE VALUES OF X ARE EQUAL
C   A DISCONTINUITY IS ASSUMED AND SEPERATE INTERPOLATION IS DONE ON
C   EACH STRICTLY MONOTONIC SEGMENT.  THE ARRAYS MUST BE DIMENSIONED AT
C   LEAST - X(I2), Y(I2), Q(3,I2), AND F(3,I2).  F IS WORKING STORAGE
C   FOR RSPLN.
C                                                     -RPB
C
      DIMENSION X(*),Y(*),Q(3,*),F(3,*),YY(3)
      EQUIVALENCE (YY(1),Y0)
      DATA SMALL/1.E-5/,YY/0.,0.,0./
      J1=I1+1
      Y0=0.
C   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL.
      IF(I2-I1)13,17,8
 8    A0=X(J1-1)
C   SEARCH FOR DISCONTINUITIES.
      DO 3 I=J1,I2
      B0=A0
      A0=X(I)
      IF(ABS((A0-B0)/AMAX1(A0,B0)).LT.SMALL) GO TO 4
 3    CONTINUE
 17   J1=J1-1
      J2=I2-2
      GO TO 5
 4    J1=J1-1
      J2=I-3
C   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
 5    IF(J2+1-J1)9,10,11
C   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
 10   J2=J2+2
      Y0=(Y(J2)-Y(J1))/(X(J2)-X(J1))
      DO 15 J=1,3
      Q(J,J1)=YY(J)
 15   Q(J,J2)=YY(J)
      GO TO 12
C   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
 11   A0=0.
      H=X(J1+1)-X(J1)
      H2=X(J1+2)-X(J1)
      Y0=H*H2*(H2-H)
      H=H*H
      H2=H2*H2
C   CALCULATE DERIVITIVE AT NEAR END.
      B0=(Y(J1)*(H-H2)+Y(J1+1)*H2-Y(J1+2)*H)/Y0
      B1=B0
C   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
      DO 1 I=J1,J2
      H=X(I+1)-X(I)
      Y0=Y(I+1)-Y(I)
      H2=H*H
      HA=H-A0
      H2A=H-2.*A0
      H3A=2.*H-3.*A0
      H2B=H2*B0
      Q(1,I)=H2/HA
      Q(2,I)=-HA/(H2A*H2)
      Q(3,I)=-H*H2A/H3A
      F(1,I)=(Y0-H*B0)/(H*HA)
      F(2,I)=(H2B-Y0*(2.*H-A0))/(H*H2*H2A)
      F(3,I)=-(H2B-3.*Y0*HA)/(H*H3A)
      A0=Q(3,I)
 1    B0=F(3,I)
C   TAKE CARE OF LAST TWO ROWS.
      I=J2+1
      H=X(I+1)-X(I)
      Y0=Y(I+1)-Y(I)
      H2=H*H
      HA=H-A0
      H2A=H*HA
      H2B=H2*B0-Y0*(2.*H-A0)
      Q(1,I)=H2/HA
      F(1,I)=(Y0-H*B0)/H2A
      HA=X(J2)-X(I+1)
      Y0=-H*HA*(HA+H)
      HA=HA*HA
C   CALCULATE DERIVITIVE AT FAR END.
      Y0=(Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
      Q(3,I)=(Y0*H2A+H2B)/(H*H2*(H-2.*A0))
      Q(2,I)=F(1,I)-Q(1,I)*Q(3,I)
C   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
      DO 2 J=J1,J2
      K=I-1
      Q(1,I)=F(3,K)-Q(3,K)*Q(2,I)
      Q(3,K)=F(2,K)-Q(2,K)*Q(1,I)
      Q(2,K)=F(1,K)-Q(1,K)*Q(3,K)
 2    I=K
      Q(1,I)=B1
C   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
 9    J2=J2+2
      DO 14 J=1,3
 14   Q(J,J2)=YY(J)
C   SEE IF THIS DISCONTINUITY IS THE LAST.
 12   IF(J2-I2)6,13,13
C   NO.  GO BACK FOR MORE.
 6    J1=J2+2
      IF(J1-I2)8,8,7
C   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
 7    DO 16 J=1,3
 16   Q(J,I2)=YY(J)
C   FINI.
 13   RETURN
      END


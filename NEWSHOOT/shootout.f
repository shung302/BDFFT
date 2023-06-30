c output the solutions 
c
      subroutine shootout(n,nid,fid,check,radistrs,gsp,e2g)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      DIMENSION e2g(3,3)
      INTEGER n,nid
      CHARACTER*80 fid,fout
      LOGICAL check
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      INCLUDE 'paths.inc'

      fout='shoot.'//fid(1:nid)
      write(6,*)radistrs,dsin(radistrs) 
c      write(6,*)fout
c      open(11,file=fout)
c      sii0=dsin(yp(2,1))
c      cti0=1.0d0/dtan(yp(2,1))
c      sidrs=dsin(radistrs)
c      te=90.
      do i=1,n
c         write(11,39)xp(i),yp(1,i),yp(2,i)*dpr,yp(3,i)*dpr,
c     &    vp(i),dvp(i),d2vp(i)
         re=yp(1,i)
         te=90.*rpd
         pe=yp(3,i)
         call rtp2xyz(re,te,pe,xe,ye,ze)
         call rotate(xe,ye,ze,e2g,xg,yg,zg)
         call xyz2rtp(xg,yg,zg,rg,tg,pg)
c         write(11,39)rg,90.-tg*dpr,pg*dpr,re,90.-te*dpr,pe*dpr
      enddo
c      close(11)
      check=.true.
c39    format(4f14.7,f8.2,2e14.6)
39    format(6f10.2)

      return
      end


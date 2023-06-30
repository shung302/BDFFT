c  find the nearest point on the central raypath for a given scatterer
c
      subroutine projq_df(iphase,nlags,lagb,lage,lag,ni,x,y,z,
     &   nq,ieta,xeta,yeta,zeta,q1,q2,rs,ts,ps)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
c      INTEGER nrbp,irbp(10)
      INTEGER nlags,lagb(10),lage(10)
      DIMENSION xeta(10),yeta(10),zeta(10),q(3,10),q1(10),q2(10),dq(10)
      DIMENSION qv(3),pk(3)
      DIMENSION ieta(10),ieta2(10),iprj(50),indx(10)
      COMMON /xyzpath/ xr(NSTPMX),yr(NSTPMX),zr(NSTPMX)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /diffr/ phip,phiq,tlp,tlq,tpq,cb,rb,rmu0c,crbt
      INCLUDE 'paths.inc'


      dmin=1.0d20

      iq=1
c find  the lags to be projeted
      dq(iq)=1.0d20
      do i=lagb(lag),lage(lag)
         d=(x-xr(i))**2+(y-yr(i))**2+(z-zr(i))**2
         if (d.lt.dq(iq)) then
            ieta(iq)=i
            dq(iq)=d
         endif
      enddo
c check if normal to wave vector
      ie=ieta(iq)
      xeta(iq)=xr(ie)
      yeta(iq)=yr(ie)
      zeta(iq)=zr(ie)
      q(1,iq)=x-xeta(iq)
      q(2,iq)=y-yeta(iq)
      q(3,iq)=z-zeta(iq)
      call pkv(yp(2,ie),yp(3,ie),pk)
      call productv(q(1,iq),pk,qvlg,pklg,prd,cs,acs)
      call decomq(iq,ieta,xeta,yeta,zeta,q,q1,q2)
c for diffracted segment, find traveltimes between projection point eta to points p and q
      if (lag.eq.2) then
      call xyz2rtp(xeta(iq),yeta(iq),zeta(iq),reta,teta,peta)
      tlp=(peta-phip)*rcmb/cb
      tlq=(phiq-peta)*rcmb/cb
c      write(6,'(''projq='',2i5,8e12.4,5f8.2)')
c     & iq,ieta(iq),xeta(iq),yeta(iq),zeta(iq),
c     & q1(iq),q2(iq),tlp,tlq,tpq,peta*dpr,phip*dpr,phiq*dpr,rcmb,cb
      endif
19    format(2i2,i5,e15.6,3f8.2,7f10.3)
29    format(2i2,i5,e15.6,3f8.2,7f10.3,'*')
39    format(2i2,i5,e15.6,3f8.2,7f10.3,'**')
      return
      end

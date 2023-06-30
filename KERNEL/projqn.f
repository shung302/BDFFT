c  find the nearest point on the central raypath for a given scatterer
c
      subroutine projqn(nlags,lagb,lage,lag,x,y,z,rs,
     &   iq,ieta,xeta,yeta,zeta,q1,q2)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
c      INTEGER nrbp,irbp(10)
      INTEGER nlags,lagb(10),lage(10)
      DIMENSION xeta(10),yeta(10),zeta(10),q(3,10),q1(10),q2(10),dq(10)
      DIMENSION qv(3),pk(3)
      DIMENSION ieta(10),ieta2(10),iprj(50),indx(10)
      COMMON /xyzpath/ xr(NSTPMX),yr(NSTPMX),zr(NSTPMX)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      INCLUDE 'paths.inc'


      dmin=1.0d20

      iq=1
c find  the lags to be projected
      lb=lag-1
      le=lag+1
      if (lag.eq.1) lb=1
      if (lag.eq.nlags) le=nlags
      do 200 lg=lb,le
         if (lg.eq.lag) go to 200
            iq=iq+1
            dq(iq)=1.0d20
            do 110 i=lagb(lg),lage(lg)
               d=(x-xr(i))**2+(y-yr(i))**2+(z-zr(i))**2
               if (d.lt.dq(iq)) then
                  ieta(iq)=i
                  dq(iq)=d
               endif
110         continue
c check if normal to wave vector
            ie=ieta(iq)
            xeta(iq)=xr(ie)
            yeta(iq)=yr(ie)
            zeta(iq)=zr(ie)
            q(1,iq)=x-xeta(iq)
            q(2,iq)=y-yeta(iq)
            q(3,iq)=z-zeta(iq)
c            call pkv(yp(2,ie),yp(3,ie),pk)
c            call productv(q(1,iq),pk,qvlg,pklg,prd,cs,acs)
c            if ( acs.gt.0.08d0 ) then
c  check if it's located at bounce points or near bounce points
             if ( ie.eq.lagb(lg).or.ie.eq.lage(lg) ) then
               ratio=(2.0d0*rcmb-rs)/rs
               xi=x*ratio
               yi=y*ratio
               zi=z*ratio
               dq(iq)=1.0d20
               do 120 i=lagb(lag),lage(lag)
                  d=(xi-xr(i))**2+(yi-yr(i))**2+(zi-zr(i))**2
                  if (d.lt.dq(iq)) then
                     ieta(iq)=i
                     dq(iq)=d
                  endif
120            continue
               ie=ieta(iq)
               xeta(iq)=xr(ie)
               yeta(iq)=yr(ie)
               zeta(iq)=zr(ie)
               q(1,iq)=xi-xeta(iq)
               q(2,iq)=yi-yeta(iq)
               q(3,iq)=zi-zeta(iq)
             endif
c          endif
200     continue
c        if (iq.eq.0) return
        call decomq(iq,ieta,xeta,yeta,zeta,q,q1,q2)
19	format(2i2,i5,e15.6,3f8.2,7f10.3)
29	format(2i2,i5,e15.6,3f8.2,7f10.3,'*')
39	format(2i2,i5,e15.6,3f8.2,7f10.3,'**')
       return
       end

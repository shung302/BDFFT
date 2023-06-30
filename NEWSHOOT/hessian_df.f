c calculate forward and backward hessian matrix (hmf,hmb)
      subroutine hessian_df(nid,fid,ivin,nvar,ystart,nlags,lagb,lage)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INCLUDE 'phases.inc'
      INTEGER nlags,lagb(10),lage(10)
      REAL*8 ystart(nvar)
      CHARACTER*80 fid,fout
c in path xp(1)=l (arc length in km)
c     yp(1,i)=radius
c     yp(2,i)=i
c     yp(3,i)=phi in radius
c         yp(1-11,nstpmx)=
c  1: r (radius in km)
c  2: ainc (incident angle in radius)
c  3: travel time (in sec)
c  4-7: (P1,P2,Q1,Q2) in forward direction starting in P1=P2=1 and Q1=Q2=0
c  8-11: (P1,P2,Q1,Q2) in foreward direction starting in P1=P2=0 and Q1=Q2=1 
c
      COMMON /rayint/ ips,idistrs,idepsrc,nsurf,ncmb,nocb,niob
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /diffr/ phip,phiq,tlp,tlq,tpq,cb,rb,rmu0c
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /grad_model/ gc(2,3,2)
      INCLUDE 'paths.inc'
      INCLUDE 'PQpaths.inc'
c
c initial conditions
c
      sirs=dsin(radistrs)/(rcmb*cb)
      fout='hess.'//fid(1:nid)
      open(12,file=fout)
c ystart(1)=P1=1, ystart(2)=Q1=0, ystart(3)=~P1=0, ystart(4)=~Q1=1
c pq(i=1-4,l)=(P1,Q1,~P1,~Q1)
c
      ystart(1)=1.0d0
      ystart(2)=0.0d0
      ystart(3)=0.0d0
      ystart(4)=1.0d0
      ystart(5)=radsrc
      ystart(6)=yp(2,1)
      call tracer_pqdf(ivin,nvar,ystart,lagb(1),lage(1))
      q1fp=pq(2,lage(1))
      q1tp=pq(4,lage(1))
      rbi=1.0d0/yp(1,lage(1))-dvp(lage(1))/vp(lage(1))
c q1ps= -cb*rs*cos(aincs)*b^-1*rb, where rb^-1=1/b-cb^dot/cb
      q1ps=-vp(lage(1))*radsrc*dcos(yp(2,lagb(1)))/(yp(1,lage(1))*rbi)

      ystart(1)=1.0d0
      ystart(2)=0.0d0
      ystart(3)=0.0d0
      ystart(4)=1.0d0
      ystart(5)=yp(1,lagb(3))
      ystart(6)=yp(2,lagb(3))
      call tracer_pqdf(ivin,nvar,ystart,lagb(3),lage(3))
      q1fd=pq(2,lage(3))
      q1td=pq(4,lage(3))

      rbi=1.0d0/yp(1,lagb(3))-dvp(lagb(3))/vp(lagb(3))
      q1rq=vp(lagb(3))*rsurf*dcos(yp(2,lage(3)))/(rcmb*rbi)

      write(6,*)'q1rq & q1ps=',q1rq,pq(2,lage(3)),q1ps,pq(2,lage(1))
      write(6,*)'q1=',q1fd,q1td,lage(3),q1fp,q1tp,lage(1),yp(2,1)*dpr,
     & yp(2,lagb(3))*dpr

      cb=vp(lage(1))
      dcb=dvp(lage(1))
      d2cb=d2vp(lage(1))

      sid=dsin(radistrs)
      sii0=dsin(yp(2,1))
      rconst=yp(1,1)*dcos(yp(2,1))
c      write(6,*)'sai0=',sai0,rconst
      do 200 k=1,lage(1)
c replace pq(5-8,k) as backward P1,P2,Q1,Q2 basd upon propagator matrix
         r=yp(1,k)
         r2=r*r
         sip=dsin(yp(3,k))
         sidp=dsin(radistrs-yp(3,k))
         pqnom=pq(1,k)*q1tp-pq(3,k)*q1fp
         pqden=pq(4,k)*q1fp-pq(2,k)*q1tp
         if (pq(2,k).eq.0.0d0.or.pqden.eq.0.0d0) then
            gn(k)=999.0d0
            write(6,*)yp(1,k),yp(2,k)*dpr,yp(3,k)*dpr
         else
            gn(k)=pq(1,k)/pq(2,k)+pqnom/pqden
         endif
         if (sip.eq.0.0d0.or.sidp.eq.0d0) then
            gh(k)=999.0d0
            write(6,*)yp(1,k),yp(2,k)*dpr,yp(3,k)*dpr
         else
            gh(k)=rayp*sid/(r2*sip*sidp)
         endif
         sii=dsin(yp(2,k))
         csi=dcos(yp(2,k))
         tai=sii/csi  
         cti=1.0d0/tai
c         q1test=-cb*r*csi/(rbi*rcmb)
         rcti0=r*csi*pq(1,k)+rayp2*dvp(k)*pq(2,k)/r-pq(2,k)/vp(k)
         write(12,49)yp(2,k)*dpr,yp(3,k)*dpr,yp(1,k),rcti0,rconst,
     &    gn(k),gh(k),pq(1,k),pq(2,k)
200   continue      

c for diffracted lag
      do 250 k=lagb(2),lage(2)
         r=yp(1,k)
c         r2=r*r
         sip=dsin(yp(3,k))
         sidp=dsin(radistrs-yp(3,k))
c         gh(k)=rayp*sid/(r2*sip*sidp)
         gh(k)=sirs/(sip*sidp)
         gn(k)=0.0d0
250   continue

      sii0=dsin(yp(2,lagb(3)))
      rconst=yp(1,lagb(3))*dcos(yp(2,lagb(3)))
c      sii0=dsin(yp(2,lage(3)))
c      rconst=yp(1,lage(3))*dcos(yp(2,lage(3)))
      do 300 k=lagb(3),lage(3)
         r=yp(1,k)
         r2=r*r
         sip=dsin(yp(3,k))
         sidp=dsin(radistrs-yp(3,k))
         pqnom=pq(1,k)*q1td-pq(3,k)*q1fd
         pqden=pq(4,k)*q1fd-pq(2,k)*q1td
         if (pq(2,k).eq.0.0d0.or.pqden.eq.0.0d0) then
            gn(k)=999.0d0
            write(6,*)yp(1,k),yp(2,k)*dpr,yp(3,k)*dpr
         else
            gn(k)=pq(1,k)/pq(2,k)+pqnom/pqden
         endif
         if (sip.eq.0.0d0.or.sidp.eq.0d0) then
            gh(k)=999.0d0
            write(6,*)yp(1,k),yp(2,k)*dpr,yp(3,k)*dpr
         else
            gh(k)=rayp*sid/(r2*sip*sidp)
         endif

c  check for accuracy in-plane spreading factor
c         sii=dsin(yp(2,k))
c         csi=dcos(yp(2,k))
c         tai=sii/csi
c         cti=1.0d0/tai
c         q1test=cb*r*csi/(rbi*rcmb)
c         rcti0=r*csi*pq(1,k)+rayp2*dvp(k)*pq(2,k)/r-pq(2,k)/vp(k)
c         write(12,49)yp(2,k)*dpr,yp(3,k)*dpr,yp(1,k),rcti0,rconst,
c     &    gn(k),gh(k),pq(1,k),pq(2,k)
c         write(12,49)yp(2,k)*dpr,yp(3,k)*dpr,yp(1,k),rcti0,rconst,vp(k),dvp(k),d2vp(k),
c     &  pq(1,k),pq(2,k)
300   continue
c replace the geometrical spreading factors at source and receivers 
c by the values at the neighboring points
c      gh(1)=gh(2)
c      gn(1)=gn(2)
c      gh(lage(3))=gh(lage(3)-1)
c      gn(lage(3))=gn(lage(3)-1)

      close(12)
49    format(3f10.4,8e14.6)
                                   
       return
       end 

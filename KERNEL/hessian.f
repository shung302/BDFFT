c calculate forward and backward hessian matrix (hmf,hmb)
      subroutine hessian(nid,fid,ivin,nvar,ystart,dh,ni,
     &   bigr,gs)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INCLUDE 'phases.inc'
      REAL*8 ystart(nvar)
      CHARACTER*80 fid,fout,fout1
c in path xp(1)=l (arc length in km)
c     yp(1,i)=radius
c     yp(2,i)=i
c     yp(3,i)=phi
c
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /rayint/ ips,idistrs,idepsrc,nsurf,ncmb,nocb,niob
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(12)
      COMMON /grad_model/ gc(2,3,2)
      COMMON /flat_model/ fc(2,3,2)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      INCLUDE 'paths.inc'
      INCLUDE 'PQpaths.inc'
c
c initial conditions
c
c      fout='pqfb.'//fid(1:nid)
c      open(12,file=fout)
c      fout='hess.'//fid(1:nid)
c      open(13,file=fout)
c      fout='sprd.'//fid(1:nid)
c      open(14,file=fout)
c ystart
c  1-4: (P1,P2,Q1,Q2) in forward direction starting in P1=P2=1 and Q1=Q2=0
c  5-8: (P1,P2,Q1,Q2) in backward direction starting in P1=P2=0 and Q1=Q2=1 
c  9: radius
c  10: incident angle

      ystart(1)=1.0d0
      ystart(2)=1.0d0
      ystart(3)=0.0d0
      ystart(4)=0.0d0
      ystart(5)=0.0d0
      ystart(6)=0.0d0
      ystart(7)=1.0d0
      ystart(8)=1.0d0
      ystart(9)=radsrc
      ystart(10)=yp(2,1)
      call tracer_pq(ivin,nvar,ystart,ni)

      q1fr=pq(3,ni)
      q2fr=pq(4,ni)
      q1tr=pq(7,ni)
      q2tr=pq(8,ni)


c      fout1='ckpq.'//fid(1:nid)
c      open(15,file=fout1)

      sii0=dsin(yp(2,1))
      rconst=radsrc*dcos(yp(2,1))
c      write(6,*)'sai0=',sai0,rconst
      do 200 k=1,ni
c replace pq(5-8,k) as backward P1,P2,Q1,Q2 basd upon propagator matrix
         pq(5,k)=pq(1,k)*q1tr-pq(5,k)*q1fr
         pq(6,k)=pq(2,k)*q2tr-pq(6,k)*q2fr
         pq(7,k)=pq(3,k)*q1tr-pq(7,k)*q1fr
         pq(8,k)=pq(4,k)*q2tr-pq(8,k)*q2fr

c forward hessian M11=P1/Q1 and M22=P3/Q3
         if (pq(3,k).eq.0.0d0) then
            hmf(1,k)=sign(999.d0,pq(1,k))
         else
            hmf(1,k)=pq(1,k)/pq(3,k)
         endif

         if (pq(4,k).eq.0.0d0) then
             hmf(2,k)=sign(999.d0,pq(2,k))
         else
             hmf(2,k)=pq(2,k)/pq(4,k)
         endif

c backward hessian M11=-P1/Q1 and M22=-P3/Q3
         if (pq(7,k).eq.0.0d0) then
            hmb(1,k)=sign(999.d0,pq(5,k))
         else
            hmb(1,k)=-pq(5,k)/pq(7,k)
         endif

         if (pq(8,k).eq.0.0d0) then
            hmb(2,k)=sign(999.d0,pq(6,k))
         else
            hmb(2,k)=-pq(6,k)/pq(8,k)
         endif                  

c  check for accuracy of geometrical spreading and P Q and hessian matrix M1, M2
         sip=dsin(yp(3,k))
         sipi=dsin(yp(3,k)+yp(2,k))
         r=yp(1,k)
         sii=dsin(yp(2,k))
         csi=dcos(yp(2,k))
         tai=sii/csi
         cti=csi/sii
         rp=(r*sii)/vp(k)*rpd
         gs2=r/sii0*sii*sip*yp(4,k)
         gs=dsqrt(dabs(gs2))
         bigr2=(pq(3,k)*pq(4,k))/(vp(1)*vp(1))
         bigr=dsqrt(dabs(bigr2))

         cti0=yp(5,k)*cti-(dvp(k)/vp(k)-1.0d0/r)*yp(4,k)
         rcti0=r*csi*pq(1,k)+rayp2*dvp(k)*pq(3,k)/r-pq(3,k)/vp(k)

         am2=sipi/(vp(k)*r*sip)
         ap2=radsrc*sipi/(rayp*vp(k))
         aq2=r*radsrc*sip/rayp

c         write(12,49)yp(3,k)*dpr,yp(1,k),(pq(j,k),j=1,8)      
c         write(13,49)yp(3,k)*dpr,yp(1,k),(hmf(j,k),j=1,2),(hmb(j,k),j=1,2)      
c         write(14,49)xp(k)*dpr,yp(1,k),bigr,bigr2,gs,gs2
c         write(15,49)yp(3,k)*dpr,yp(1,k),
c     &    (pq(j,k),j=1,4),ap2,aq2,am2,rcti0

200    continue
c       write(6,*)'geometrical spreading=',bigr,bigr2
c       close(12)
c       close(13)
c       close(14)
c       close(15)
49    format(2f12.6,8e14.6)
                                   
       return
       end 

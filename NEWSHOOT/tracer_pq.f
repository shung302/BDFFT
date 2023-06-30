      SUBROUTINE tracer_pq(ivin,nvar,ystart,ni)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      PARAMETER (NBISMX=100,EPS=1.0d-8)
      INCLUDE 'phases.inc'
      INTEGER nstep,nvar,NMAX
      REAL*8 x1,x2,ystart(NMX),y(NMX),dydx(NMX),yn(NMX)
      EXTERNAL derivs_pq
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /grad_model/ gc(2,3,2)
      COMMON /flat_model/ fc(2,3,2)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      INCLUDE 'paths.inc'
      INCLUDE 'PQpaths.inc'

c store initial conditoion in path xx(1) and y(1,1) and y(2,1)

c      write(6,*)(ystart(i),i=1,nvar),h,ncmb,nsurf,nvar 
      do 11 i=1,nvar-2
	 pq(i,1)=ystart(i)
	 y(i)=ystart(i)
 11   continue
      y(9)=ystart(9)
      y(10)=ystart(10)

      do 110 k=1,ni-1
         k1=k+1
         ir=irp(k)
         ir1=irp(k1)
         if (ir.ne.ir1.or.
     &      (yp(1,k1).eq.rcmb.or.yp(1,k1).eq.rsurf)) then
            rd=yp(1,k1)
            rdi=1.0d0/rd
            rd2=rd*rd
            c=vp(k)
            dc=dvp(k)
            d2c=d2vp(k)
            c1=vp(k1)
            dc1=dvp(k1)
            d2c1=d2vp(k1)
            csi=dcos(yp(2,k))
            if ((ir.eq.ir1).and.(rd.eq.rsurf.or.rd.eq.rcmb)) then
               csi1=-csi
               yn(1)=-y(1)
     &              +y(3)*((1.0d0/(csi*c1)-1.0d0/(csi1*c))/rd
     &              -rayp2*(dc1/csi-dc/csi1)/rd2)
               yn(2)=y(2)+y(4)*(csi1/c1-csi/c)/rd
               yn(3)=-y(3)
               yn(4)=y(4)
               yn(5)=-y(5)
     &              +y(7)*((1.0d0/(csi*c1)-1.0d0/(csi1*c))/rd
     &              -rayp2*(dc1/csi-dc/csi1)/rd2)
               yn(6)=y(6)+y(8)*(csi1/c1-csi/c)/rd
               yn(7)=-y(7)
               yn(8)=y(8)
            else
c hit internal discontinuities
               csi1=dcos(yp(2,k1))
               ratcs=csi1/csi
               yn(1)=y(1)/ratcs
     &              +y(3)*((1.0d0/(csi*c1)-1.0d0/(csi1*c))/rd
     &              -rayp2*(dc1/csi-dc/csi1)/rd2)
               yn(2)=y(2)+y(4)*(csi1/c1-csi/c)/rd
               yn(3)=y(3)*ratcs
               yn(4)=y(4)
               yn(5)=y(5)/ratcs
     &              +y(7)*((1.0d0/(csi*c1)-1.0d0/(csi1*c))/rd
     &              -rayp2*(dc1/csi-dc/csi1)/rd2)
               yn(6)=y(6)+y(8)*(csi1/c1-csi/c)/rd
               yn(7)=y(7)*ratcs
               yn(8)=y(8)
            endif
            yn(9)=rd
            yn(10)=yp(2,k1)
         else
            ir=irp(k)
            iw=iwp(k)
            x=yp(3,k)
            hi=yp(3,k1)-x
            call derivs_pq(ivin,nvar,rayp,x,y,dydx,ir,iw)
            call driver_pq_rk4(ivin,ir,iw,nvar,rayp,y,dydx,x,hi,yn,
     &           derivs_pq)
         endif         
c store the incremental steps for pq
c so the final point at nstep+1
         do i=1,nvar-2
            pq(i,k1)=yn(i)
            y(i)=yn(i)
         enddo
            y(9)=yn(9)
            y(10)=yn(10)
110   continue
      return
      end

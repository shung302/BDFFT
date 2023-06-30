c calculate kernel using paraxial approximation
c ci=c is the velocity at scatterer

      subroutine kernblk(ivin,nid,fid,ips,iphase,nlags,lagb,lage,ni,
     & nkpt,rbke,tbke,pbke,skblk,distrs,insflag,ntstf,tstar)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INCLUDE 'phases.inc'
      DIMENSION tstar(*),rbke(NBLKMX),tbke(NBLKMX),pbke(NBLKMX),skblk(20,NBLKMX)
      INTEGER nlags,lagb(10),lage(10)
      DIMENSION q1(10),q2(10),ieta(10),xeta(10),yeta(10),zeta(10)
      DIMENSION skq(10),tdiffq(10),sprq(10),rotm(3,3),sk(20)
      CHARACTER*80 fid,fout
      COMMON /xyzpath/ xr(NSTPMX),yr(NSTPMX),zr(NSTPMX)
      COMMON /bphi/ tdiff,sig
      COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),rkn(NKMX,MSMX,20)

      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      INCLUDE 'paths.inc'
      INCLUDE 'PQpaths.inc'
c      external dfint

c   
      ickq=0
      nl=ni
      ts0=pi2
      do 1000 k=1,nkpt
	 rs=rbke(k)
	 ts=tbke(k)
	 ps=pbke(k)
         call findr(ivin,rs,ir)
         call velo1(ivin,ips,ir,rs,c,dc,d2c)
         if (c.eq.0.0d0) call velo1(ivin,ips,ir+1,rs,c,dc,d2c)
         ci=0.5d0/c
         if (ts.eq.ts0) ickq=1
         call prjlag(radistrs,rs,ts,ps,nlags,lagb,lage,lag)
         call rtp2xyz(rs,ts,ps,xs,ys,zs)
         call projq(iphase,nlags,lagb,lage,lag,ni,xs,ys,zs,
     &       nq,ieta,xeta,yeta,zeta,q1,q2,rs,ts,ps)
c         call projq_r(iphase,nlags,lagb,lage,lag,ni,xs,ys,zs,
c     &                    nq,ieta,q1,q2,rs,ts,ps)

c        write(6,*)'r,t,p=',rs,ts*dpr,ps*dpr,xs,ys,zs,nq

         do istf=1,ntstf
            sk(istf)=0.0d0
         enddo

         do 500 iq=1,nq
            if (ieta(iq).eq.1.or.ieta(iq).eq.ni) go to 500
c           if (q1(iq)*q1(iq)+q2(iq)*q2(iq).gt.qmax) go to 500
c hmf1 = M11';  hmf2 = M22'
c hmb1 = M11''; hmb2 = M22''
c h11 = M11'+M11''
c h22 = M22'+M22''
            h11=hmf(1,ieta(iq))+hmb(1,ieta(iq))
            h22=hmf(2,ieta(iq))+hmb(2,ieta(iq))
            spr=dsqrt(dabs(h11*h22))
            tdiff=0.5d0*(h11*q1(iq)*q1(iq)+h22*q2(iq)*q2(iq))
	    if (tdiff.gt.ttmx) go to 500
            call signature(iphase,h11,h22,nsig,sig)
            call eval_kntstr(nsig,tdiff,spr,ntstf,sk)
500      continue

         do istf=1,ntstf
            sk(istf)=sk(istf)*ci
	    skblk(istf,k)=sk(istf)
         enddo
1000  continue 

19    format(f13.6,2f10.2,8e13.4)
39    format('*',i3,f8.2,9e13.5)
21    format(i3,10f10.3,f12.5,2e15.6)
      return
      end

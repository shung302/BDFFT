c calculate kernel using paraxial approximation
c
      subroutine kernela(ivin,nid,fid,ips,iphase,nlags,lagb,lage,ni,
     & dl,drtb,rtb,dtbd,tbd0,distrs,insflag,ntstf,tstar)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INCLUDE 'phases.inc'
      DIMENSION tstar(*)
c      INTEGER nrbp,irbp(10)
      INTEGER nlags,lagb(10),lage(10)
      DIMENSION q1(10),q2(10),ieta(10),xeta(10),yeta(10),zeta(10)
      DIMENSION skq(10),tdiffq(10),sprq(10),rotm(3,3),sk(20)
c      REAL*8 wx(1),ww(1),stfa2(1)
      CHARACTER*80 fid,fout
      COMMON /xyzpath/ xr(NSTPMX),yr(NSTPMX),zr(NSTPMX)
      COMMON /bphi/ tdiff,sig
      COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),
     &    rkn(NKMX,MSMX,20)

c      COMMON /stf/ alpha,tp, insflag, amp2(NTMX)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      INCLUDE 'paths.inc'
      INCLUDE 'PQpaths.inc'
c      external dfint
       data tdmax /100.0d0/

      fout='kern.'//fid(1:nid)
      nrecl=24
c      open(11,file=fout,status='unknown',access='direct',
c     &  recl=nrecl)
      fout='kern.'//fid(1:nid)
      open(11,file=fout)
      nrtb=int(rtb/drtb)+1
      write(6,*)nrtb,rtb,drtb
      ntb=int(360.0d0/dtbd)
      dtb=rpd*dtbd
      tb0=rpd*tbd0
      irec=0
      ickq=0
      nl=ni
      ni2=ni/2

      write(11,'(3i5,e15.8,2f10.2,f6.1,i4,28x,''(nl,nrtb,ntb,dl,drtb,dtbd,dist,depth)'')')
     & nl,nrtb,ntb,dl,drtb,dtbd,distrs,int(rsurf-radsrc)
      write(11,'(6f10.2,28x,''(lbeg,lend,rbeg,rend,tbeg,tend)'')')
     & 0.,xp(nl),0.0,rtb,tb0,tbd0+360.
      if (insflag.eq.1) then
      write(11,'(10f8.2,8x,''(tstart)'')')(tstar(i),i=1,ntstf)
      else
      write(11,'(10f8.2,8x,''(tpstf)'')')(tstar(i),i=1,ntstf)
      endif

      do 1100 lag=1,nlags
      do 1000 k=lagb(lag),lage(lag)
c 1/(twopi*c) where c uses reference velocity at central ray
         dl=xp(k+1)-xp(k)
         rl=yp(1,k)
         pl=yp(3,k)
         ail=yp(2,k)
         tl=pi2
c         cl=vp(k)
c         ci=0.50d0/cl
         call rtp2xyz(rl,tl,pl,xl,yl,zl)
         xeta(1)=xl
         yeta(1)=yl
         zeta(1)=zl
         ieta(1)=k
         call coordl(ail,pl,rotm)
         do 900 i=1,nrtb
            rx=dble(i)*drtb
            do 800 j=1,ntb
               tx=tb0+dble(j-1)*dtb
               q1(1)=rx*dcos(tx)
               q2(1)=rx*dsin(tx)
               nq=1
               call locates(xl,yl,zl,0.0d0,q1(1),q2(1),rotm,
     &           xs,ys,zs,rs,ts,ps)
c if outside of the mantle, skip
c               if (rs.gt.rsurf.or.rs.lt.rcmb) go to 800
               if (rs.gt.rsurf) go to 800

cc find background velocity
               call findr(ivin,rs,irs)
               call velo1(ivin,ips,irs,rs,cs,dcs,d2cs)
               ci=0.50d0/cs

               if (nlags.ne.1) 
     &            call projqn(nlags,lagb,lage,lag,xs,ys,zs,rs,
     &                 nq,ieta,xeta,yeta,zeta,q1,q2)
c	      write(6,*)'nq=',nq,(q1(iq),q2(iq),iq=1,nq)
          
               do istf=1,ntstf
               sk(istf)=0.0d0
               enddo

               do 500 iq=1,nq
                  if (ieta(iq).eq.1.or.ieta(iq).eq.ni) go to 500
c                  if (q1(iq)*q1(iq)+q2(iq)*q2(iq).gt.qmax) go to 500
c hmf1 = M11';  hmf2 = M22'
c hmb1 = M11''; hmb2 = M22''
c h11 = M11'+M11''
c h22 = M22'+M22''
                  h11=hmf(1,ieta(iq))+hmb(1,ieta(iq))
                  h22=hmf(2,ieta(iq))+hmb(2,ieta(iq))
                  spr=dsqrt(dabs(h11*h22))
                  tdiff=0.5d0*(h11*q1(iq)*q1(iq)+h22*q2(iq)*q2(iq))
                  if (tdiff.gt.ttmx) go to 500
c                  write(21,21)iq,rs,ts*dpr,ps*dpr,xs,ys,zs,
c     &               xeta(iq),yeta(iq),zeta(iq),q1(iq),tdiff,h11,h22
c                  write(6,*)'pq=',(pqf(ii,ieta(iq)),ii=1,4),
c     &            (pqb(ii,ieta(iq)),ii=1,4)
c                  write(6,*)'h11=',iq,h11,h22,tdiff,spr,rs,ts*dpr,ps*tpr

                  call signature(iphase,h11,h22,nsig,sig)
                  call eval_kntstr(nsig,tdiff,spr,ntstf,sk)
500            continue
               do istf=1,ntstf
               sk(istf)=sk(istf)*ci
               enddo
c            write(11,19)rs,ts*dpr,ps*dpr,(sk(istf),istf=1,ntstf),tdiff,spr
800         continue
900      continue
1000  continue 
1100  continue 
      close(11)
19    format(f13.6,2f10.2,8e13.4)
39    format('*',i3,f8.2,9e13.5)
21    format(i3,10f10.3,f12.5,2e15.6)
      return
      end

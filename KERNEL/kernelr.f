c calculate kernel using paraxial approximation
c ci=c is the velocity at scatterer

      subroutine kernelr(ivin,nid,fid,ips,iphase,nlags,lagb,lage,ni,
     & dra,dth,dph,doff,distrs,insflag,ntstf,tstar,g2e,e2g)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INCLUDE 'phases.inc'
      REAL*8 e2g(3,3),g2e(3,3)
      DIMENSION tstar(*)
c      INTEGER nrbp,irbp(10)
      INTEGER nlags,lagb(10),lage(10)
      DIMENSION q1(10),q2(10),ieta(10),xeta(10),yeta(10),zeta(10)
      DIMENSION skq(10),tdiffq(10),sprq(10),rotm(3,3),sk(20)
c      REAL*8 wx(1),ww(1),stfa2(1)
      CHARACTER*80 fid,fout
      COMMON /xyzpath/ xr(NSTPMX),yr(NSTPMX),zr(NSTPMX)
      COMMON /bphi/ tdiff,sig
      COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),rkn(NKMX,MSMX,20)

c      COMMON /stf/ alpha,tp, insflag, amp2(NTMX)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      INCLUDE 'paths.inc'
      INCLUDE 'PQpaths.inc'
c      external dfint

c      nrecl=24
c      open(11,file=fout,status='unknown',access='direct',
c     &  recl=nrecl)
      fout='kern.'//fid(1:nid)
      if (iphase.lt.15.or.(iphase.ge.21.and.iphase.lt.35)) then
         rend=rcmb
      else
         rend=0.0d0
      endif
      open(11,file=fout)
      nra=int((rsurf-rend)/dra)+2
      dra=(rsurf-rend)/dble(nra-1)
      nph=int((distrs+10.0d0+0.001d0*dph)/dph)+1
      if (dth.eq.0.0d0) then
          nth=1
          ts0=pi2
      else
          nth=int((doff+0.001d0*dth)/dth)+1
          ts0=pi2
      endif
      dphr=dph*rpd
      dthr=dth*rpd
      rs0=rsurf
      
      ps0=-5.0d0*rpd
      irec=0
      ickq=0

c   
      thend=90.0d0+dble(nth-1)*dth
      phend=-5.0d0+dble(nph-1)*dph
      write(11,'(3i5,e15.8,2f10.2,f6.1,i4,8x,''(nra,nth,nph,dra,dth,dph,dist,depth)'')')
     & nra,nth,nph,dra,dth,dph,distrs,int(rsurf-radsrc)
      write(11,'(6f10.2,8x,''(rbeg,rend,thbeg,thend,phbeg,phend)'')')
     & rs0,rend,90.0,thend,-5.0,phend
      if (insflag.eq.3) then
      write(11,'(10f8.2,8x,''(tstart)'')')(tstar(i),i=1,ntstf)
      elseif (insflag.eq.1.or.insflag.eq.2) then
      write(11,'(10f8.2,8x,''(tpstf)'')')(tstar(i),i=1,ntstf)
      elseif (insflag.eq.4) then
      write(11,'(2f8.2,8x,''(begin and end periods)'')')tstar(1),tstar(2)
      endif

      ickq=0
9     format(3i8,3f10.2,f13.6,2f10.2)
      nl=ni
      do 1000 k=1,nra
c         write(6,*)'k=',k
         rs=rs0-dble(k-1)*dra
c          rs=3971.0d0
c         if (rs.lt.rcmb) rs=rcmb
         call findr(ivin,rs,ir)
         call velo1(ivin,ips,ir,rs,c,dc,d2c)
         if (c.eq.0.0d0) call velo1(ivin,ips,ir+1,rs,c,dc,d2c)
c         if (rs.eq.rcmb) write(6,*)rs,ts*dpr,ps*dpr,c,dc,d2c
c         cst=0.5d0/(pi*c)
          ci=0.5d0/c
         do 900 j=1,nth
            ts=ts0+dble(j-1)*dthr
            if (ts.eq.ts0) ickq=1
            do 800 i=1,nph
               irec=irec+1
               ps=ps0+dble(i-1)*dphr
               call prjlag(radistrs,rs,ts,ps,nlags,lagb,lage,lag)
               call rtp2xyz(rs,ts,ps,xs,ys,zs)
               call projq(iphase,nlags,lagb,lage,lag,ni,xs,ys,zs,
     &                    nq,ieta,xeta,yeta,zeta,q1,q2,rs,ts,ps)
c              write(6,*)'r,t,p=',rs,ts*dpr,ps*dpr,xs,ys,zs,nq

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
c                 write(6,*)hmf(1,ieta(iq)),hmb(1,ieta(iq)),ieta(iq),hmf(2,ieta(iq)),hmb(2,ieta(iq))
c                 write(6,*)'h11=',h11,h22,tdiff,spr,rs,ts*dpr,ps*tpr

                  call signature(iphase,h11,h22,nsig,sig)
                  call eval_kntstr(nsig,tdiff,spr,ntstf,sk)
500            continue

               do istf=1,ntstf
                  sk(istf)=sk(istf)*ci
               enddo
c            write(11,19)rs,ts*dpr,ps*dpr,(sk(istf),istf=1,ntstf),tdiff,spr
            write(11,19)rs,ts*dpr,ps*dpr,(sk(istf),istf=1,ntstf)
c            write(6,19)rs,ts*dpr,ps*dpr,(sk(istf),istf=1,ntstf),tdiff,spr
800         continue
900      continue
1000  continue 

      close(11)
19    format(f13.6,2f10.2,8e13.4)
39    format('*',i3,f8.2,9e13.5)
21    format(i3,10f10.3,f12.5,2e15.6)
      return
      end

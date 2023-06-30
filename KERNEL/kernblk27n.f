ce calculate kernel using paraxial approximation
c ci=c is the velocity at scatterer

      subroutine kernblk27n(ivin,nid,fid,ips,iphase,nlags,lagb,lage,ni,
     & jsta,g2e,e2g,distrs,insflag,ntstf,tstar)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      INCLUDE 'phases.inc'
      INCLUDE 'paths.inc'
      INCLUDE 'PQpaths.inc'
      INCLUDE 'buildGx.inc'
      DIMENSION tstar(*)
      INTEGER nlags,lagb(10),lage(10)
      DIMENSION q1(10),q2(10),ieta(10),xeta(10),yeta(10),zeta(10)
      DIMENSION skq(10),tdiffq(10),sprq(10),rotm(3,3),sk(20)
      REAL*8 e2g(3,3),g2e(3,3)

      COMMON /xyzpath/ xr(NSTPMX),yr(NSTPMX),zr(NSTPMX)
c      COMMON /bphi/ tdiff,sig
      COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),rkn(NKMX,MSMX,20)

      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
c      external dfint

c   
c	write(6,*)'nkpt=',nkpt,distrs,nvolu,ntstf,nvertg,rC,delr
      ickq=0
c      nl=ni
      ts0=pi2
c	write(6,*)ntstf,nvolu,nlocat_g,jsta,nx,ny,nz,dx,dy,dz
c      do istf=1,ntstf
c         do i=1,nvolu
c            coef(i,istf)=0.0d0
c            write(91,*)'jsta=',jsta,i,coef(i,istf)
c         enddo
c      enddo

c      write(6,*)'iposit=',locat_g(nlocat_g),nlocat_g
      do 1000 ilocat=1,nlocat_g
c subroutine from chiao to calculate the position of gaussian-integration points
         iz=(locat_g(ilocat)-1)/((nx-1)*(ny-1))+1
         ixy=locat_g(ilocat)-(iz-1)*(nx-1)*(ny-1)
         iy=(ixy-1)/(nx-1)+1
         ix=ixy-(iy-1)*(nx-1)
         iposition(1)=(iz-1)*nx*ny+(iy-1)*nx+ix
         iposition(2)=(iz-1)*nx*ny+(iy-1)*nx+ix+1
         iposition(3)=(iz-1)*nx*ny+(iy)*nx+ix
         iposition(4)=(iz-1)*nx*ny+(iy)*nx+ix+1
         iposition(5)=(iz)*nx*ny+(iy-1)*nx+ix
         iposition(6)=(iz)*nx*ny+(iy-1)*nx+ix+1
         iposition(7)=(iz)*nx*ny+(iy)*nx+ix
         iposition(8)=(iz)*nx*ny+(iy)*nx+ix+1
c         write(6,*)'iposit=',(iposition(j),j=1,8),locat_g(ilocat),ix,iy,iz,ixy,nlocat_g
         lon0=x0+dble(ix-1)*dx
         lon1=lon0+dx
         lonc=(lon0+lon1)*0.5d0
         lat0=y0+dble(iy-1)*dy
         lat1=lat0+dy
         latc=(lat0+lat1)*0.5d0
         rad0=z0+dble(iz-1)*dz
         rad1=rad0+dz
         radc=(rad0+rad1)*0.5d0
         area0=dabs(dcos(pi/2.-lat1*rpd)-dcos(pi/2.-lat0*rpd))*0.5d0*
     *        dx/90.d0*pi
         do 900 iiz=1,3
          rad=radc+sample0(iiz)*dz*0.5d0
          area=area0*rad*rad
          wtz=weight0(iiz)*0.5d0
          do 800 iiy=1,3
           latp=latc+sample0(iiy)*dy*0.5d0
           wty=weight0(iiy)*0.5d0
           do 700 iix=1,3
            lonp=lonc+sample0(iix)*dx*0.5d0
            wtx=weight0(iix)*0.5d0
c get lat and long of geographic earth=(lon,lat)
            call tran(lon,lat,lonp,latp,xx,yy,zz,-1)
c	    write(6,*)'posit=',iix,iiy,iiz,lon,lat,lonp,latp
c start to calculate kernel
c first rotate to equatorial plane
c            write(6,*)lon,lat,rad
            rbg=rad
            tbg=(90.0d0-lat)*rpd
            pbg=lon*rpd
            call rtp2xyz(rbg,tbg,pbg,xbg,ybg,zbg)
            call rotate(xbg,ybg,zbg,g2e,xbe,ybe,zbe)
            call xyz2rtp(xbe,ybe,zbe,rs,ts,ps)
            call findr(ivin,rs,ir)
            call velo1(ivin,ips,ir,rs,c,dc,d2c)
c	    write(6,*)'r, theta, phi=',rbg,tbg,pbg,rs,ts,ps,c
            if (c.eq.0.0d0) call velo1(ivin,ips,ir+1,rs,c,dc,d2c)
                ci=0.5d0/c
            if (ts.eq.ts0) ickq=1
            call prjlag(radistrs,rs,ts,ps,nlags,lagb,lage,lag)
            call rtp2xyz(rs,ts,ps,xs,ys,zs)
            call projq(iphase,nlags,lagb,lage,lag,ni,xs,ys,zs,
     &           nq,ieta,xeta,yeta,zeta,q1,q2,rs,ts,ps)

c        write(6,*)'r,t,p=',rs,ts*dpr,ps*dpr,xs,ys,zs,nq

            do istf=1,ntstf
               sk(istf)=0.0d0
            enddo

            do 500 iq=1,nq
               if (ieta(iq).eq.1.or.ieta(iq).eq.ni) go to 500
c             if (q1(iq)*q1(iq)+q2(iq)*q2(iq).gt.qmax) go to 500
c hmf1 = M11';  hmf2 = M22'
c hmb1 = M11''; hmb2 = M22''
c h11 = M11'+M11''
c h22 = M22'+M22''
               h11=hmf(1,ieta(iq))+hmb(1,ieta(iq))
               h22=hmf(2,ieta(iq))+hmb(2,ieta(iq))
               spr=dsqrt(dabs(h11*h22))
               tdiff=0.5d0*(h11*q1(iq)*q1(iq)+h22*q2(iq)*q2(iq))
c	       write(*,*)'h11 & h22=',h11,h22,q1(iq),q2(iq),tdiff,spr
	       if (tdiff.gt.ttmx) then
c                  write(6,*)'sk=0',jsta,lon,lat,tdiff,spr
                  go to 500
               endif
               call signature(iphase,h11,h22,nsig,sig)
               call eval_kntstr(nsig,tdiff,spr,ntstf,sk)
500         continue
           
            do istf=1,ntstf
c	       write(6,*)'sk=',sk(istf)*ci,ci,area,dz,wtx,wty,wtz
               iii=(iiz-1)*9+(iiy-1)*3+iix
               wtw(iii,istf)=-0.5*sk(istf)*wtx*wty*wtz*area*dz
c               wtw(iii,istf)=-0.5*sk(istf)*wtx*wty*wtz*area*dz
c               write(6,*)'iii=',iii,wtw(iii,istf)
            end do
700        continue
800       continue
900      continue

c loop over all the source time functions
c$         do istf=1,ntstf
c$            do ii=1,8
c$            value(ii)=0.d0
c$          end do
c$            do iiz=1,3
c$            jjz=iiz+1
c$              if(jjz.gt.3)jjz=1
c$                kkz=iiz-1
c$              if(kkz.lt.1)kkz=3
c$              do iiy=1,3
c$                 jjy=iiy+1
c$                 if(jjy.gt.3)jjy=1
c$                 kky=iiy-1
c$                 if(kky.lt.1)kky=3
c$                 do iix=1,3
c$                  jjx=iix+1
c$                  if(jjx.gt.3)jjx=1
c$                  kkx=iix-1
c$                  if(kkx.lt.1)kkx=3
c$                  ipo=(iiz-1)*9+(iiy-1)*3+iix
c$                  value(1)=value(1)+wtw(ipo,istf)*wl(iiz)*wl(iiy)*wl(iix)
c$                  value(2)=value(2)+wtw(ipo,istf)*wl(kkz)*wl(iiy)*wl(iix)
c$                  value(3)=value(3)+wtw(ipo,istf)*wl(iiz)*wl(kky)*wl(iix)
c$                  value(4)=value(4)+wtw(ipo,istf)*wl(kkz)*wl(kky)*wl(iix)
c$                  value(5)=value(5)+wtw(ipo,istf)*wl(iiz)*wl(iiy)*wl(kkx)
c$                  value(6)=value(6)+wtw(ipo,istf)*wl(kkz)*wl(iiy)*wl(kkx)
c$                  value(7)=value(7)+wtw(ipo,istf)*wl(iiz)*wl(kky)*wl(kkx)
c$                  value(8)=value(8)+wtw(ipo,istf)*wl(kkz)*wl(kky)*wl(kkx)
c$                 end do
c$              end do
c$            end do
c$            do ii=1,8
cc$              write(92,*)'coef before=',jsta,iposition(ii), coef(iposition(ii),istf),value(ii)
c$                coef(iposition(ii),istf)=coef(iposition(ii),istf)+value(ii)
c$            end do
c$        end do
c      write(6,*)'ilocat',ilocat,jsta

c===================================================================
c bugs fixed by Peter Nelson and informed by Xiaofeng Liang, 04/27/2017
c see email from X. Liang 1/5/2017
c===================================================================

        do istf=1,ntstf
           do ii=1,8
            value(ii)=0.d0
           end do
           do iiz=1,3
             kkz=4-iiz
             do iiy=1,3
               kky=4-iiy
               do iix=1,3
                 kkx=4-iix
                 ipo=(iiz-1)*9+(iiy-1)*3+iix
                 value(1)=value(1)+wtw(ipo,istf)*wl(iiz)*wl(iiy)*wl(iix)
                 value(2)=value(2)+wtw(ipo,istf)*wl(iiz)*wl(iiy)*wl(kkx)
                 value(3)=value(3)+wtw(ipo,istf)*wl(iiz)*wl(kky)*wl(iix)
                 value(4)=value(4)+wtw(ipo,istf)*wl(iiz)*wl(kky)*wl(kkx)
                 value(5)=value(5)+wtw(ipo,istf)*wl(kkz)*wl(iiy)*wl(iix)
                 value(6)=value(6)+wtw(ipo,istf)*wl(kkz)*wl(iiy)*wl(kkx)
                 value(7)=value(7)+wtw(ipo,istf)*wl(kkz)*wl(kky)*wl(iix)
                 value(8)=value(8)+wtw(ipo,istf)*wl(kkz)*wl(kky)*wl(kkx)
               end do
             end do
           end do
           do ii=1,8
              coef(iposition(ii),istf)=coef(iposition(ii),istf)+value(ii)
           end do
       end do

1000  continue 

19    format(f13.6,2f10.2,8e13.4)
39    format('*',i3,f8.2,9e13.5)
21    format(i3,10f10.3,f12.5,2e15.6)
      return
      end

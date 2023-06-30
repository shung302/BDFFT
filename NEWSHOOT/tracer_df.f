      SUBROUTINE tracer_df(ivin,ips,nsurf,ncmb,noc,nvar,ystart,x1,x2,h,
     &  fvec,ni,isurf,icmb,iocb,iiob)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      PARAMETER (NBISMX=100,EPS=1.0d-8,TOL=1.0d-6)
      INCLUDE 'parameters.inc'
c      INTEGER nstep,nvar,NMAX,NSTPMX
      REAL*8 x1,x2,ystart(NMX),y(NMX),dydx(NMX),yn(NMX),yc(NMX)
      EXTERNAL derivs_df
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /grad_model/ gc(2,3,2)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn 
      INCLUDE 'paths.inc'

c store initial conditoion in path xp(1) and yp(1,1) and yp(2,1)

      isurf=0
      icmb=0
      iocb=0
      iiob=0
c      write(6,*)'nocb=',nocb
      do 11 i=1,nvar
	 y(i)=ystart(i)
	 ypdf(i,1)=y(i)
 11   continue
      hi=h
      xpdf(1)=x1
      x=x1
      call findr(ivin,y(1),ir1)
      call derivs_df(ivin,ips,nvar,x,y,dydx,ir1,iw,vr,dvr,d2vr)
      vpdf(1)=vr
      dvpdf(1)=dvr
      d2vpdf(1)=d2vr
      irpdf(1)=ir1
      iwpdf(1)=iw

      k=1
      do 110 while (k.le.nstpmx)
c determine the range of the radius "r=y(1)" to see if the continutiy 
c conditions need to be applied
c
         call driver_rk4(ivin,ips,nvar,y,dydx,x,hi,yn,derivs_df)
	 call findr(ivin,yn(1),ir2)
c	 write(6,*)'k=',k,ir2,x,yn(1),yn(2),dydx(1),dydx(2)
c check if the ray hits the discontinuity or reach the surface
         
         if ((ir2.eq.3) 
     &         .and.((yn(2)-pi2)*(y(2)-pi2).lt.0.0d0)) then
10          dr=rcmb-y(1)
            hi=dr/dcos(y(2)) 
20          call driver_rk4(ivin,ips,nvar,y,dydx,x,hi,yn,derivs_df) 
            if ((yn(2)-pi2)*(y(2)-pi2).lt.0.0d0) then
               hi=hi*0.5d0
               go to 20
            else
               k=k+1
               x=x+hi
               xpdf(k)=x
               do i=1,nvar
                  ypdf(i,k)=yn(i)
                  y(i)=yn(i)
               enddo
               vpdf(k)=vr
               dvpdf(k)=dvr
               d2vpdf(k)=d2vr
               irpdf(k)=ir1
               iwpdf(k)=iw
               if (rcmb-yn(1).lt.tol) then
                  return
               else
                  go to 10
               endif
            endif
         endif
         irdiff=abs(ir1-ir2)
         ird=min(ir1,ir2)
         if (irdiff.gt.1) then
            hi=hi*0.5d0
            go to 110
         endif
            
         if (irdiff.eq.1) then
            rd=rdis(ird)
            if (ir1.eq.ird) then
               if (ird.eq.2) then
                  vd1=vdis(1,ird,1)
	          vd2=vdis(2,ird,ips)
               else
                  vd1=vdis(1,ird,ips)
	          vd2=vdis(2,ird,ips)
               endif
            else
               if (ird.eq.2) then
                  vd1=vdis(2,ird,ips)
	          vd2=vdis(1,ird,1)
               else
                  vd1=vdis(2,ird,ips)
	          vd2=vdis(1,ird,ips)
               endif
            endif
c trace through discontinuity at ir2, must satisfy the continuity conditions
	    hh=h
	    do 50 kk=1,nbismx
               hh=0.50d0*hh
               call derivs_df(ivin,ips,nvar,x,y,dydx,ir1,iw,vr,dvr,d2vr)
               call driver_rk4(ivin,ips,nvar,y,dydx,x,hh,yn,derivs_df)
	       call findr(ivin,yn(1),ir)
               if (dabs(yn(1)-rd).lt.eps.and.ir.eq.ir1) then
                  k=k+1
                  ni=k
                  x=x+hh
                  xpdf(k)=x
	          do i=1,nvar
                     ypdf(i,k)=yn(i)
	          enddo
                  vpdf(k)=vr
                  dvpdf(k)=dvr
                  d2vpdf(k)=d2vr
                  irpdf(k)=ir1
                  iwpdf(k)=iw
c	 write(6,*)'disc=',k,ir2,x,yn(1),yn(2),dydx(1),dydx(2)
c hit the surface
                  if (rd.eq.rsurf) then
                     if (isurf.eq.nsurf) then
                        fvec=x2-ypdf(3,k)
c                        write(6,*)'fvec=',x2,x,fvec
                        return
                     else
                        isurf=isurf+1
                        call velo1(ivin,iw,ir1,yn(1),c,dc,d2c)
                        call velo1(ivin,iw,ir1,rd,c1,dc1,d2c1)
                        tai=dtan(yn(2))
                        tai1=-tai
                        rdi=1.0d0/rd
                        y(1)=rd
                        y(2)=pi-yn(2)
                        y(3)=yn(3)
                        call derivs_df(ivin,ips,nvar,x,y,dydx,ir1,
     &                       iw,vr,dvr,d2vr)
                        k=k+1
                        ni=k
                        xpdf(k)=x
                        do i=1,nvar
                           ypdf(i,k)=y(i)
                        enddo
                        vpdf(k)=vr
                        dvpdf(k)=dvr
                        d2vpdf(k)=d2vr
                        irpdf(k)=ir1
                        iwpdf(k)=iw
                        hi=h
                        go to 110
                     endif
c hit internal discontinuities
                  elseif (rd.ne.rsurf.and.rd.ne.rcmb) then
                     call velo1(ivin,iw,ir1,yn(1),c,dc,d2c)    
                     tai=dtan(yn(2))
                     k=k+1
                     ni=k
                     rdi=1.0d0/rd
                     y(1)=rd
                     si=vd2/vd1*dsin(yn(2))
                     if (dabs(si).ge.1.0d0) then
                        y(2)=pi2
                        y(3)=yn(3)
                        call velo1(ivin,iw,ir1,rd,c1,dc1,d2c1)
                        xpdf(k)=x
                        x2=x
                        irpdf(k)=ir1
                        do i=1,nvar
                           ypdf(i,k)=y(i)
                        enddo
                        vpdf(k)=vr
                        dvpdf(k)=dvr
                        d2vpdf(k)=d2vr
                        iwpdf(k)=iw
                        return
                     elseif (si.eq.0.0d0) then
                        y(2)=pi
                     else
                        ir1=ir2
                        if (yn(2).le.pi2) then
                           y(2)=dasin(si)
                        else
                           y(2)=pi-dasin(si)
                        endif
                     endif
                     y(3)=yn(3)
                     tai1=dtan(y(2))
                     call velo1(ivin,iw,ir1,rd,c1,dc1,d2c1) 
                     xpdf(k)=x
                     irpdf(k)=ir1
	             do i=1,nvar
                        ypdf(i,k)=y(i)
                     enddo
                     call derivs_df(ivin,ips,nvar,x,y,dydx,ir1,
     &                    iw,vr,dvr,d2vr)
                     vpdf(k)=vr
                     dvpdf(k)=dvr
                     d2vpdf(k)=d2vr
                     iwpdf(k)=iw
                     hi=h
                     go to 110
                   endif
               else 
c               if ((yn(2)-pi2)*(y(2)-pi2).le.0.0d0) then
c                  if ((ir.eq.3).and.((yn(2)-pi2)*(y(2)-pi2).le.0.0d0)) then
                   
                  if (ir.eq.ir1) then
c                     write(6,*)'hh=',hh,yn(2)*dpr
                     x=x+hh
                     do i=1,nvar
                        y(i)=yn(i)
                     enddo
                  endif
                  go to 50
               endif
50         continue
        endif
c end of continuation of ray at discontinuities
        k=k+1
        ni=k
        x=x+hi
        xpdf(k)=x
c store incremental integration of each step in xx(i+1),yp(1,i+1),yp(2,i+1)
c so the final point at nstep+1
        do i=1,nvar
          ypdf(i,k)=yn(i)
          y(i)=yn(i)
        enddo
        call derivs_df(ivin,ips,nvar,x,y,dydx,ir1,iw,vr,dvr,d2vr)
        vpdf(k)=vr
        dvpdf(k)=dvr
        d2vpdf(k)=d2vr
        irpdf(k)=ir1
        iwpdf(k)=iw
110	enddo

      return
      end


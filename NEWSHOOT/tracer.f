      SUBROUTINE tracer(ivin,ips,nsurf,ncmb,nocb,niob,nvar,ystart,x1,x2,h,fvec,
     &  ni,isurf,icmb,iocb,iiob)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      PARAMETER (NBISMX=100,EPS=1.0d-8)
      INCLUDE 'parameters.inc'
      INTEGER nstep,nvar,NMAX
      REAL*8 x1,x2,ystart(NMX),y(NMX),dydx(NMX),yn(NMX),yc(NMX)
      EXTERNAL derivs
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /grad_model/ gc(2,3,2)
      COMMON /flat_model/ fc(2,3,2)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn 
      INCLUDE 'paths.inc'

c store initial conditoion in path xp(1) and yp(1,1) and yp(2,1)

      isurf=0
      icmb=0
      iocb=0
      iiob=0
c      write(6,*)'noc=',nocb,niob,nvar,nsurf,ni,riob,rcmb
      do 11 i=1,nvar
	 y(i)=ystart(i)
	 yp(i,1)=y(i)
 11   continue
      hi=h
      xp(1)=x1
      x=x1
      call findr(ivin,y(1),ir1)
      call derivs(ivin,ips,nvar,x,y,dydx,ir1,iw,vr,dvr,d2vr)
      vp(1)=vr
      dvp(1)=dvr
      d2vp(1)=d2vr
      irp(1)=ir1
      iwp(1)=iw

      k=1
      do 110 while (k.le.nstpmx)
c determine the range of the radius "r=y(1)" to see if the continutiy 
c conditions need to be applied
c
         call driver_rk4(ivin,ips,nvar,y,dydx,x,hi,yn,derivs)
	 call findr(ivin,yn(1),ir2)
c	 write(6,*)'k=',k,ir2,x,yn(1),yn(2),dydx(1),dydx(2)
c check if the ray hits the discontinuity or reach the surface
            
         irdiff=abs(ir1-ir2)
         ird=min(ir1,ir2)
         if (irdiff.gt.1) then
            hi=hi*0.5d0
            go to 110
         endif
            
         if (irdiff.eq.1) then
            rd=rdis(ird)
c propagate upward 
            if (ir1.eq.ird) then
c outer core to mantle (2th layer-->3th layer)
               if (ird.eq.2) then
                  vd1=vdis(1,ird,1)
	          vd2=vdis(2,ird,ips)
c inner core to outer core
               elseif (ir2.eq.2) then
                  vd1=vdis(1,ird,ips)
                  vd2=vdis(2,ird,1)
               else
                  vd1=vdis(1,ird,ips)
	          vd2=vdis(2,ird,ips)
               endif
            else
c propagate downward
c mantle to outer core
               if (ird.eq.2) then
                  vd1=vdis(2,ird,ips)
	          vd2=vdis(1,ird,1)
c outer core to inner core
              elseif (ir1.eq.2) then
                  vd1=vdis(2,ird,1)
	          vd2=vdis(1,ird,ips)
              else
                  vd1=vdis(2,ird,ips)
	          vd2=vdis(1,ird,ips)
               endif
            endif
c            write(6,*)irdiff,ird,ir1,ir2,vd1,vd2
c
c trace through discontinuity at ir2, must satisfy the continuity conditions
	    hh=h
	    do 50 kk=1,nbismx
               hh=0.50d0*hh
               call derivs(ivin,ips,nvar,x,y,dydx,ir1,iw,vr,dvr,d2vr)
               call driver_rk4(ivin,ips,nvar,y,dydx,x,hh,yn,derivs)
	       call findr(ivin,yn(1),ir)
               if (dabs(yn(1)-rd).lt.eps.and.ir.eq.ir1) then
                  k=k+1
                  ni=k
                  x=x+hh
                  xp(k)=x
	          do i=1,nvar
                     yp(i,k)=yn(i)
	          enddo
                  vp(k)=vr
                  dvp(k)=dvr
                  d2vp(k)=d2vr
                  irp(k)=ir1
                  iwp(k)=iw
c	 write(6,*)'disc=',k,ir2,x,yn(1),yn(2),dydx(1),dydx(2)
c hit the surface
                  if (rd.eq.rsurf) then
                     if (isurf.eq.nsurf) then
                        fvec=x2-yp(3,k)
c                        write(6,*)'fvec=',x2,x,yp(3,k),fvec
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
                        call derivs(ivin,ips,nvar,x,y,dydx,ir1,
     &                       iw,vr,dvr,d2vr)
                        k=k+1
                        ni=k
                        xp(k)=x
                        do i=1,nvar
                           yp(i,k)=y(i)
                        enddo
                        vp(k)=vr
                        dvp(k)=dvr
                        d2vp(k)=d2vr
                        irp(k)=ir1
                        iwp(k)=iw
                        hi=h
                        go to 110
                     endif
c hit internal discontinuities
                  else
                     call velo1(ivin,iw,ir1,yn(1),c,dc,d2c)    
                     tai=dtan(yn(2))
                     k=k+1
                     ni=k
                     rdi=1.0d0/rd
                     y(1)=rd
                     if ((rd.eq.rcmb).and.(icmb.lt.ncmb).and.
     &                  yn(2).gt.pi2) then
                        y(2)=pi-yn(2)
                        icmb=icmb+1
                     elseif ((rd.eq.rcmb).and.(iocb.lt.nocb).
     &                  and.yn(2).lt.pi2) then
                        y(2)=pi-yn(2)
                        iocb=iocb+1
                     elseif ((rd.eq.riob).and.(iiob.lt.niob).
     &                  and.yn(2).gt.pi2) then
                        y(2)=pi-yn(2)
                        iiob=iiob+1
                     else
                        si=vd2/vd1*dsin(yn(2))
c                        write(6,*)'si=',si,vd2,vd1,yn(2)*dpr,rd
                        if (dabs(si).gt.1.0d0) then
                           y(2)=pi-yn(2)
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
                     endif
                     y(3)=yn(3)
                     tai1=dtan(y(2))
                     call velo1(ivin,iw,ir1,rd,c1,dc1,d2c1) 
                     xp(k)=x
                     irp(k)=ir1
	             do i=1,nvar
                        yp(i,k)=y(i)
                     enddo
                     call derivs(ivin,ips,nvar,x,y,dydx,ir1,
     &                    iw,vr,dvr,d2vr)
                     vp(k)=vr
                     dvp(k)=dvr
                     d2vp(k)=d2vr
                     iwp(k)=iw
                     hi=h
                     go to 110
                  endif
               else
                  if (ir.eq.ir1) then
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
        xp(k)=x
c store incremental integration of each step in xx(i+1),yp(1,i+1),yp(2,i+1)
c so the final point at nstep+1
        do i=1,nvar
          yp(i,k)=yn(i)
          y(i)=yn(i)
        enddo
        call derivs(ivin,ips,nvar,x,y,dydx,ir1,iw,vr,dvr,d2vr)
        vp(k)=vr
        dvp(k)=dvr
        d2vp(k)=d2vr
        irp(k)=ir1
        iwp(k)=iw
110	enddo

      return
      end


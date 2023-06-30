c determine velocity at discontinuities
c add ak135 model for ivelin=2
c 

      subroutine velodisc(ivelin,rayprm,grang)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      REAL*8 dvdr(2),vr(2),rayprm(2),grang(2)
      CHARACTER*5 mdn
      INCLUDE 'phases.inc'            
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /flat_model/ fc(2,3,2)
      COMMON /grad_model/ gc(2,3,2)
      COMMON /homo_model/ hc(3,2)
      COMMON /lyrs_model/ yc(11,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak

c determine velocity jumps at existing discontinuities
      if (ivelin.eq.1) then
         ndis=11
         do ir=1,ndis
            r=rdis(ir)
            call vprem(ir,r,dvdr,vr)
            vdis(1,ir,1)=vr(1)
            vdis(1,ir,2)=vr(2)
            call vprem(ir+1,r,dvdr,vr)
            vdis(2,ir,1)=vr(1)
            vdis(2,ir,2)=vr(2)
         enddo
      elseif (ivelin.eq.2) then
         open(22,file='ak135-c.model')
         read(22,*)mdn,nlak
         do ir=1,nlak
            read(22,*)akr(ir),rho,akc(ir,1),akc(ir,2)
         enddo
         jr=0
         do ir=1,nlak-1
            if (akr(ir).eq.akr(ir+1)) then
               jr=jr+1
               rdis(jr)=akr(ir)
               vdis(1,jr,1)=akc(ir,1)
               vdis(1,jr,2)=akc(ir,2)
               vdis(2,jr,1)=akc(ir+1,1)
               vdis(2,jr,2)=akc(ir+1,2)
               akdc(ir,1)=0.0d0
               akdc(ir,2)=0.0d0
            else
               akdc(ir,1)=(akc(ir+1,1)-akc(ir,1))/(akr(ir+1)-akr(ir))
               akdc(ir,2)=(akc(ir+1,2)-akc(ir,2))/(akr(ir+1)-akr(ir))
            endif
         enddo
         rdis(jr+1)=akr(nlak)
         vdis(1,jr+1,1)=akc(nlak,1)
         vdis(1,jr+1,2)=akc(nlak,2)
         vdis(2,jr+1,1)=akc(nlak,1)
         vdis(2,jr+1,2)=akc(nlak,2)
         akdc(nlak,1)=0.0d0
         akdc(nlak,2)=0.0d0
         ndis=jr+1
c         do i=1,ndis
c            write(6,*)rdis(i),(vdis(1,i,j),vdis(2,i,j),j=1,2)
c         enddo
         close(22)
       elseif (ivelin.eq.3) then
         open(22,file='ak135-o.model')
         read(22,*)mdn,nlak
         do ir=1,nlak
            read(22,*)akr(ir),rho,akc(ir,1),akc(ir,2)
         enddo
         jr=0
         do ir=1,nlak-1
            if (akr(ir).eq.akr(ir+1)) then
               jr=jr+1
               rdis(jr)=akr(ir)
               vdis(1,jr,1)=akc(ir,1)
               vdis(1,jr,2)=akc(ir,2)
               vdis(2,jr,1)=akc(ir+1,1)
               vdis(2,jr,2)=akc(ir+1,2)
               akdc(ir,1)=0.0d0
               akdc(ir,2)=0.0d0
            else
               akdc(ir,1)=(akc(ir+1,1)-akc(ir,1))/(akr(ir+1)-akr(ir))
               akdc(ir,2)=(akc(ir+1,2)-akc(ir,2))/(akr(ir+1)-akr(ir))
            endif
         enddo
         rdis(jr+1)=akr(nlak)
         vdis(1,jr+1,1)=akc(nlak,1)
         vdis(1,jr+1,2)=akc(nlak,2)
         vdis(2,jr+1,1)=akc(nlak,1)
         vdis(2,jr+1,2)=akc(nlak,2)
         akdc(nlak,1)=0.0d0
         akdc(nlak,2)=0.0d0
         ndis=jr+1
c         do i=1,ndis
c            write(6,*)rdis(i),(vdis(1,i,j),vdis(2,i,j),j=1,2)
c         enddo
         close(22)
      elseif (ivelin.eq.-4) then
         ndis=11
         do ir=1,ndis
            r=ardis(ir)
            call viasp(ir,r,dvdr,vr)
            vdis(1,ir,1)=vr(1)
            vdis(1,ir,2)=vr(2)
            call viasp(ir+1,r,dvdr,vr)
            vdis(2,ir,1)=vr(1)
            vdis(2,ir,2)=vr(2)
            rdis(ir)=r
         enddo
      elseif (ivelin.eq.0) then
         ndis=3
         rdis(3)=rsurf
         do ir=1,ndis
            if (ir.eq.3) then
            r=6371.0d0
            else
            r=rdis(ir)
            endif
            call vgrad(ir,r,dvdr,vr)
            vdis(1,ir,1)=vr(1)
            vdis(1,ir,2)=vr(2)
            call vgrad(ir+1,r,dvdr,vr)
            vdis(2,ir,1)=vr(1)
            vdis(2,ir,2)=vr(2)
          enddo
      elseif (ivelin.eq.-1) then
         ndis=3
         rdis(3)=rsurf
         do ir=1,ndis
            if (ir.eq.3) then             
            r=6371.0d0
            else
            r=rdis(ir)
            endif
            call vflat(ir,r,dvdr,vr)
            vdis(1,ir,1)=vr(1)
            vdis(1,ir,2)=vr(2)
            call vflat(ir+1,r,dvdr,vr)
            vdis(2,ir,1)=vr(1)
            vdis(2,ir,2)=vr(2)
          enddo
      elseif (ivelin.eq.-2) then
         ndis=3
         rdis(3)=rsurf
         do ir=1,ndis
            if (ir.eq.3) then
            r=6371.0d0
            else
            r=rdis(ir)
            endif
            call vhomo(ir,r,dvdr,vr)
            vdis(1,ir,1)=vr(1)
            vdis(1,ir,2)=vr(2)
            call vhomo(ir+1,r,dvdr,vr)
            vdis(2,ir,1)=vr(1)
            vdis(2,ir,2)=vr(2)
          enddo
      elseif (ivelin.eq.-3) then
         ndis=11
         do ir=1,ndis
            r=rdis(ir)
            call vlyrs(ir,r,dvdr,vr)
            vdis(1,ir,1)=vr(1)
            vdis(1,ir,2)=vr(2)
            call vlyrs(ir+1,r,dvdr,vr)
            vdis(2,ir,1)=vr(1)
            vdis(2,ir,2)=vr(2)
         enddo
      endif

c      write(6,*)'P wave'
c      write(6,*)(rdis(i),vdis(1,i,1),vdis(2,i,1),i=1,ndis)
c      write(6,*)'S wave'
c      write(6,*)(rdis(i),vdis(1,i,2),vdis(2,i,2),i=1,ndis)
                                                                         

c  determine the ray parameter for grazing ray at CMB
c for P wave
      call findr(ivelin,radsrc,irsrc)
      call velo(ivelin,irsrc,radsrc,dvdr,vr)
c      write(6,*)irsrc,radsrc,vr(1),vr(2),ivelin
      
      rayprm(1)=rcmb/vdis(2,2,1)
      rayprm(2)=rcmb/vdis(2,2,2)
      grang(1)=pi-dasin(rayprm(1)*vr(1)/radsrc)
      grang(2)=pi-dasin(rayprm(2)*vr(2)/radsrc)

c      write(6,*)'ray parameters & takeoff angles for grazing P & S rays'
c      write(6,19) (rayprm(i),grang(i)*dpr,i=1,2)
19    format('(ray parameter, takeoff angle) for grazing P & S rays:',/2x,
     & 2(e12.4,f10.1))

       if (ivelin.eq.1) then
       open(1,file='rayshoot.velout')
       do i=1,63710
          ra=i/10.
          call findr(ivelin,ra,ir)
          call vprem(ir,ra,dvdr,vr)
          write(1,'(3f12.2)')ra,vr(1),vr(2)
       enddo
       close(1)
       elseif (ivelin.eq.-4) then
       open(1,file='rayshoot.velout')
       do i=1,63710
          ra=i/10.
          call findr(ivelin,ra,ir)
          call viasp(ir,ra,dvdr,vr)
          write(1,'(3f12.2)')ra,vr(1),vr(2)
       enddo
       close(1)
       elseif (ivelin.ge.2) then
       open(1,file='rayshoot.velout')
       do i=1,63710
          ra=i/10.
          call vak135(ir,ra,dvdr,vr)
          write(1,'(3f12.2)')ra,vr(1),vr(2)
       enddo
       close(1)
       endif
       return
       end

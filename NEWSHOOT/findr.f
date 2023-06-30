c subroutine to find the range of r
c
      subroutine findr(ivin,r,ir)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      COMMON /prem_model/ pc(4,11,2),rdis(12),vdis(2,12,2)
      COMMON /iasp_model/ apc(4,11,2),ardis(11)
      COMMON /ak135_model/ akc(144,2),akdc(144,2),akr(144),nlak
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

c for flat and linear gradient model
      if (ivin.le.0.and.ivin.ge.-2) then
         if (r.le.rdis(1)) then
             ir=1
         elseif (r.gt.rdis(1).and.r.le.rdis(2)) then
             ir=2
         elseif (r.gt.rdis(2).and.r.le.rdis(11)) then
             ir=3
         else
             ir=4
         endif
c         if (r.eq.rcmb) ir=3
c PREM or 11 homogeneous layered models or IASP91
      elseif (ivin.eq.1.or.ivin.le.-3) then
c	 write(6,*)ivin,r,rdis(1),rdis(11)
         if (r.le.rdis(1)) then
         ir=1
         return
         endif

         if (r.gt.rdis(11)) then
         ir=12
         return
         endif
         do i=2,11
         if (r.gt.rdis(i-1).and.r.le.rdis(i)) then
            ir=i
            return
         endif
         enddo
c         if (r.eq.rcmb) ir=3
      elseif (ivin.eq.2) then
c for ak135 continent model
         if (r.le.rdis(1)) then
            ir=1
            return
         endif
         if (r.gt.rsurf) then
            ir=10
            return
         endif
         do i=2,9
         if (r.gt.rdis(i-1).and.r.le.rdis(i)) then
            ir=i
            return
         endif
         enddo
c         if (r.eq.rcmb) ir=3
      elseif (ivin.eq.3) then
c for ak135 spherical average model
         if (r.le.rdis(1)) then
            ir=1
            return
         endif
         if (r.gt.rsurf) then
            ir=13
            return
         endif
         do i=2,12
         if (r.gt.rdis(i-1).and.r.le.rdis(i)) then
            ir=i
            return
         endif
         enddo
c         if (r.eq.rcmb) ir=3
      endif
      return
      end

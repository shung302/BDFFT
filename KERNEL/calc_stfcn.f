c calculate source time function for given function
      subroutine calc_stfcn(tper,npts,wb,we,wx,ww,amp2)

      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      dimension wx(*),ww(*),amp2(*)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
 
c c generate tabulated abscissas and weights for gaussian-legendre
c integration between wb and we
      flow=twopi*wb
      fhigh=twopi*we
      call dgauleg(flow,fhigh,wx,ww,npts)
 
      do i=1,npts
         amp2(i)=wx(i)*wx(i)*1.0d0
      enddo        

      return
      end

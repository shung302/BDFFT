c calculate source time function for given function
      subroutine calc_stfsac(fsac,npts,wb,we,wx,ww,amp2)

      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      character*80 fsac
      parameter (ntmx=10000)
      real astf(ntmx),tstf(ntmx),dy2a(ntmx),dt,bt,dy1,dyn,xx,yy
      dimension wx(*),ww(*),amp2(*)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
 
c c generate tabulated abscissas and weights for gaussian-legendre
c integration between wb and we
      flow=twopi*wb
      fhigh=twopi*we
      call dgauleg(flow,fhigh,wx,ww,npts)
 
      call rsac1(fsac,astf,nt,bt,dt,ntmx,nerr)
      dt=twopi*dt
      bt=twopi*bt

	write(6,*)'nt=',bt,dt,nt,fsac,astf(10)
      do i=1,nt
         tstf(i)=bt+(i-1)*dt
      enddo
c use cubic spline to get stf at a given frequency wx
      dy1=(astf(2)-astf(1))/dt
      dyn=(astf(nt)-astf(nt-1))/dt
      call spline(tstf,astf,nt,dy1,dyn,dy2a)

      do i=1,npts
         xx=real(wx(i))
         call splint(tstf,astf,dy2a,nt,xx,yy)
         amp2(i)=dble(yy)
c	 write(6,*)twopi/xx,amp2(i)
      enddo        

      return
      end

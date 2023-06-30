c use bisection method to do 2-point ray tracing
c
      subroutine bisec(ivin,iraytype,nvar,dh,x1,x2,ystart,nbrk,dbrk,
     &  ni,nid,fid,gsp,e2g)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INCLUDE 'parameters.inc'
      PARAMETER (EPS=1.0d-8,TOL=1.0d-10)
      REAL*8 ystart(NMX)
      INTEGER iraytype,ips,idistrs,idepsrc,nsurf,ncmb
      DIMENSION e2g(3,3)
      CHARACTER*80 fid
      LOGICAL hit
      EXTERNAL derivs
      COMMON /rayint/ ips,idistrs,idepsrc,nsurf,ncmb,nocb,niob
      COMMON /raysrc/ radsrc,radistrs,rayp,rayp2
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      INCLUDE 'paths.inc'

        hit=.false.
        h=dh
c start bisection
        nbis=nbrk-1
        do 51 while (hit.eqv..false.)
           dbis=dbrk
c        write(6,*)'nbis=',nbis,iraytype      
           call tracer(ivin,ips,nsurf,ncmb,nocb,niob,nvar,ystart,x1,x2,h,fvec1,
     &     ni,isurf,icmb,iocb,iiob)
c           write(6,*)'fvec1 in tracer=',ni,fvec1,xp(ni),yp(1,ni),yp(2,ni),yp(3,ni)*dpr,
c     &      ystart(2)*dpr
           if (dabs(fvec1).lt.eps.and.(icmb.eq.ncmb))
     &        call shootout(ni,nid,fid,hit,radistrs,gsp,e2g)
           do 21 ib=1,nbis
              ystart(2)=ystart(2)+dbis
              call tracer(ivin,ips,nsurf,ncmb,nocb,niob,nvar,ystart,x1,x2,h,fvec2,
     &        ni,isurf,icmb,iocb,iiob)
c              write(6,*)'fvec2 in tracer=',ni,fvec2,xp(ni),yp(1,ni),yp(2,ni),yp(3,ni)*dpr,
c     &           ystart(2)*dpr
              if (dabs(fvec2).lt.eps) then
                 call shootout(ni,nid,fid,hit,radistrs,gsp,e2g)
                 go to 51
              elseif (fvec1*fvec2.lt.0.0d0) then
                  ystart(2)=ystart(2)-dbis
                  if (dabs(dbis).lt.tol) then
                     call shootout(ni,nid,fid,hit,radistrs,gsp,e2g)
                  else
                    dbrk=dbrk/dble(nbrk)
c                    write(6,*)'brake',ystart(2)*dpr,(ystart(2)+dbis)*dpr
                  endif
                  go to 51
              elseif (icmb.ne.ncmb) then
                 ystart(2)=ystart(2)-dbis
                 dbrk=dbrk/dble(nbrk)
                 go to 51
              else
                  fvec1=fvec2
              endif
21         continue
31      dbrk=dbrk/dble(nbrk)
51      enddo

        return
        end

c evaluate nominator of kernel from linear interpolation of tables
c
        subroutine eval_kntstrcl(nsig,tdiff,spr,nstf,ci,sk)
        implicit real*8 (a-h,o-z)
        dimension sk(1)
        integer jsign,msign
        include 'parameters.inc'
        COMMON /ttkern/ ms,ntt,dtt,ttmx,ttdf(NKMX),rkd(20),
     &   rkn(NKMX,MSMX,20)

        t=dabs(tdiff)
c if difference of maslov indices = 0
        if (nsig.eq.0) then
           ifac=dsign(1.0d0,tdiff)
           im=1
           go to 1
        endif
c otherwise, determine sign of kernel from dT and dM >0
        ns=abs(nsig)
        msign=sign(1,nsig)
        jsign=dsign(1.0d0,tdiff)
        im=ns+1
        ifac=(msign)**(ns)*(jsign)**(im)
c        write(61,11)tdiff,nsig,ns,msign,jsign,im,ifac
11      format(f10.4,6i5)
c       if(t.ge.ttdf(ntt)) then
c          do istf=1,nstf
c             sk(istf)=sk(istf)
c          enddo
c          return
c        endif
1       do i=2,ntt
           if (t.lt.ttdf(i)) then
              t0=(t-ttdf(i-1))/dtt
              do istf=1,nstf
                 sk(istf)=sk(istf)
     &           -ifac*(t0*(rkn(i,im,istf)-rkn(i-1,im,istf))
     &                   +rkn(i-1,im,istf))*spr*ci
              enddo
              return
           endif
        enddo

        return
        end


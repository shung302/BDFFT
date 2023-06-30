c________________________________________________________________
c      LSQR package
c________________________________________________________________
c          call linit
c          call lldrow (coef,jdx,ncoef,rhs)
c          call slsqr(nit)
c          call lgtsol(dd,errx,nnn)
c
c
c
c
c------------------------------
      subroutine linit

c ... initialize for conjugate gradient least squares

      include 'lsqrcom'
      n = 0
      m = 0
      nel=0
      return
      end

      subroutine lldrow(coef,jdx,ncoef,rhs)

c ... add one equation of data into row of matrix
c     for conjugate gradient method (slsqr) of solution

      integer ncoef,i,jdx(ncoef)
      double precision  coef(ncoef),rhs
      include 'lsqrcom'
c      write(*,*)ncoef,rhs,n,m,nel
      if (ncoef.gt.0) then

c      modified to handle the case while jdx is nor sorted in ascending order
c       i.e. the follwing loop save the work of invoking SORT2 in the main
c       program
c
c       n = max0(n,jdx(ncoef))

c        write(*,*)ncoef,rhs,n,m,nel
c      pause
       nc_temp=0
       do i=1,ncoef
        if(coef(i).ne.0.d0)then
         n=max0(n,jdx(i))
        else
         nc_temp=nc_temp+1
        endif
       enddo
        m = m+1
        b(m) = rhs
        na(m) = ncoef-nc_temp
        do i=1,ncoef
       if(coef(i).ne.0.d0)then
         nel = nel+1
         if (nel.gt.big) then
          write (*,*) 'sorry, out of room in coeff. array for lldrow'
          write(*,*)big,nel
          stop
         end if
         ra(nel) = coef(i)
         ja(nel) = jdx(i)
       endif
      enddo
      endif
c      write(*,*)ncoef,rhs,n,m,nel
      return
      end

      subroutine slsqr(itct)

c
c ... conjugate gradient routine for non-square least squares
c     formal parameters:
c       atol = ATOL from Paige and Saunders
c       btol = BTOL           "
c       conlim = CONLIM       "
c       itct = iteration count (returned)

      include 'lsqrcom'
      double precision atol,btol,conlim
c      parameter(atol=1.d-6,btol=1.d-6,conlim=0.d0)
c Y.Shen, 08/23/04
      parameter(atol=1.d-7,btol=1.d-7,conlim=0.d0)

      integer i,itct,maxit
      double precision alpha,beta,phibar,rhobar,phi,rho,c,s,theta
      double precision temp,anorm,dotp,test1,test2,rnorm
      external dotp,lnrliz

c     maxit should be set approximatly 4*n

      parameter(maxit=30000)

c ... initial housekeeping

      itct = 0
      anorm = dotp(ra,ra,nel)
      anorm = dsqrt(anorm)

c ... set initial values for start of iterations
      write(*,*)n,m,nel
      do i=1,m
        u(i) = b(i)
      enddo
      call lnrliz(beta,u,m)
c ... compute A transpose u and put results into v

      call laty(v,u)
      call lnrliz(alpha,v,n)

c ... move v into w

      do i=1,n
        w(i) = v(i)
      enddo

c ... zero out x and sig vectors

      do i=1,n
        x(i) = 0.d0
        sig(i) = 0.d0
      enddo

c ... assign phibar and rhobar

      phibar = beta
      rhobar = alpha

c****************************************
c     top of main iteration loop        *
c****************************************

150   continue

c ... compute Av, subtract alpha*u, and assign to u

      call lax(q,v)
      call lupdat(u,q,'-',alpha,u,m)
      call lnrliz(beta,u,m)

c ... compute A transpose u, subtract beta*v, and assign to v

      call laty(q,u)
      call lupdat(v,q,'-',beta,v,n)
      call lnrliz(alpha,v,n)
c

c ... do steps (a) through (g) in Paige and Saunders

      rho = dsqrt(rhobar**2+beta**2)
      if(rho.eq.0.d0)then
       return
      else
       c = rhobar/rho
       s = beta/rho
      endif
      theta = s*alpha
      rhobar = -c*alpha
      phi = c*phibar
      phibar = s*phibar

c ... update solution vector and w vector

      do i=1,n
        temp = w(i)/rho
        x(i) = x(i)+phi*temp
        w(i) = v(i)-theta*temp
        sig(i) = sig(i)+temp**2
      enddo

c ... loop back if convergence not attained
c ... test for convergence
c
c       calculate r (tmp1)
c
c      call lax(tmp1,x)
c      call lupdat(tmp1,b,'-',1.0,tmp1,m)
c      rnorm=dsqrt(dotp(tmp1,tmp1,n))
      rnorm=phibar

c     criteria 1 in P&S p.54

      test1=rnorm/anorm/dsqrt(dotp(x,x,n))


c ... A(transpose)r (tmp2)
c
c      call laty(tmp2,tmp1)
c
c     criteria 2 in P&S p.54
c
cc
c      test2=dsqrt(dotp(tmp2,tmp2,n))/anorm/rnorm
      test2=alpha*dabs(c)/anorm

      itct = itct+1
c      write(*,*)itct,test1,atol,test2
      if(mod(itct,500).eq.0)then
      write(*,'(19hiter., test1,test2:,4i11,3d13.4)')
     ,   itct,m,n,nel,test1,test2,phibar
      end if
      if((test1.le.atol).or.(test2.le.btol).or.(itct.eq.maxit))then
        write (*,'(19hiter., test1,test2:,4i11,3d13.4)')
     ,   itct,m,n,nel,test1,test2,phibar
        res2 = phibar
        return
      else
        go to 150
      endif
      end


      subroutine lnrliz(scale,vect,n)

c ... normalize the vector 'vect' to unit magnitude and
c     return the scale factor in 'scale'

      integer n,i
      double precision scale,vect(n),dotp
      external dotp
      scale = dotp(vect,vect,n)
      scale = dsqrt(scale)
      if(scale.eq.0.d0)return
      do i=1,n
        vect(i) = vect(i)/scale
      enddo
      return
      end



      subroutine lupdat(a,b,flag,scal,c,n)

c ... update vector A by adding scaled version of vector B
c     character flag indicates addition or subtraction

      integer n,i
      double precision a(n),b(n),c(n),scal
      character*2 flag
      if (flag(1:1).eq.'+') then
        do i=1,n
          a(i) = b(i)+c(i)*scal
      enddo
      else
        do i=1,n
          a(i) = b(i)-c(i)*scal
      enddo
      end if
      return
      end


      double precision function dotp(a,b,n)
      double precision a(n),b(n)
      dotp=0.d0
      do i=1,n
      dotp=dotp+a(i)*b(i)
      enddo
      return
      end


      subroutine lgtsol(solvec,serr,len)

c ... extract solution vector from conjugate gradient least squares
c ... formal parameters:
c       solvec = solution vector
c       serr = standard error estimates

      include 'lsqrcom'
      double precision solvec(*),serr(*)
      double precision temp
      integer i,len

      temp = res2/dble(max0(m-len,1))
      do i=1,len
        solvec(i) = x(i)
        serr(i) = dsqrt(temp*sig(i))
      enddo
      return
      end

      subroutine lax(yy,xx)

c ... Form matrix product AX and put result in vector Y
c     X is input in the dummy array "xx" which is assumed
cc       to  be of length "n";  the output Y is returned in
c       the dummy array "yy" of length "m"

      double precision xx(*),yy(*),sum
      integer i,j,l,l1,l2
      include 'lsqrcom'
      l2 = 0
      do i=1,m
        sum = 0.d0
        l1 = l2+1
        l2 = l2+na(i)
        do l=l1,l2
          j = ja(l)
          sum = sum+ra(l)*xx(j)
      enddo
        yy(i) = sum
      enddo
      return
      end

      subroutine laty(xx,yy)

c ... form matrix product A(transpose)Y and put result in vector X
c     The input vector Y of length "m" is passed in the dummy
c       array "yy" and the results are returned in the dummy array "xx"
c       which is of length "n"

      double precision xx(*),yy(*),yi
      integer i,j,l,l1,l2
      include 'lsqrcom'
      l2 = 0
      do i=1,n
        xx(i) = 0.0d0
      enddo
      do i=1,m
        yi = yy(i)
        l1 = l2+1
        l2 = l2+na(i)
        do l=l1,l2
          j = ja(l)
          xx(j) = xx(j)+ra(l)*yi
      enddo
      enddo
      return
      end





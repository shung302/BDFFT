c_________________________________________________________________________
c     solver_noSE.f:  
c
c        output files:
c         "try.sum"-- summary of fitting results including
c          ndata,nx,ny,nz,VARIANCE REDUCTION, MISFIT, MODEL VARIANCE and MODEL NORM
c         "try.xyz"-- the solved solution of m which is organized in arry
c          (((1,nx),1,ny),1,nz)
c
c        input files and parameters:
c         input_matrix-- G matrix formulated in 3D Cartecian. The nx,ny,nz are
c         # of regular grid intervals in X-,Y-,Z-direction.  They should all
c         equal or greater than 1 and are of the form of nx=2**(nlevel(1)-1),
c         ny=2**(nlevel(2)-1), nz=2**(nlevel(3)-1). i.e. choose nlevel(3)=1
c         means a 2D XY inversion...The total d.o.f of pursued vector m is then
c         nvolu=nx*ny*nz.
c         the parameter <mode (imultis)> controls several inversion options that
c         mode=1--> Simple damped LS inversion of Gm=d according to the original
c                   pixel parameterization
c             =2--> MultiScaled inversion; 3D linear interpolatory wavelet
c                   check solver_LSQR_pixel.f for 3D Haar or cubic splines
c                   in subroutine <face_lift> (bior1.1 or bior3.1)
c                   check solver_LSQR_node.f for linear interpolatory wavelet
c             =3--> Convulutional quelling that behaves similar as MINIMUN CURVATURE
c                   but is more flexible by convoling with  GAUSSIAN function.
c                   The parameter <sigma0> controls the one sigma # of grids of the
c                   imposed a priori smoothness preference, bigger <sigma0> enforces
c                   smoother solution!
c
c         Be sure to compare with <solver_SVD.f> for other insights
c
c        !!!!!!! make sure that declaration in <lsqrcom> is consistent !!!!!!!!!!
c
c   08/13 Y. Shen: change mdata and station location files
c		   change the defination of variance (fit in the program).
c
      character*150 input_matrix
      character*10 mode0(4),mode
      character*5  dumstn
      parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,          mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,          mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,          mbase=1,mvolu=mivolu*mbase,mdata=100000,mstation=1000,
     ,          mvar=mvolu)
      integer nbase,nlevel(3),isign,nlevel_syn(3),niedge(3),
     ,        nivolu,iivolu0,iivolu,ivolu,nvolu,niface,
     ,        i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_tmp(3),i3,ii0,
     ,        nlevel2D(2),niedge2D(2),nlevel_syn2D(2)
      integer ncoef,jdx(mvar),jdx0(mvar)
      real coef0(mvar),rhs,image0(mvar)
      double precision image(mvolu),coef(mvar),thres,errx(mvar),
     ,       cov,fit,fit_tmp,lamb,sigma0,vmax,vmin,ttmax,norm,tmp,
     ,       mean,vv1(mvar),weig_sta,sigmaz,vmin_sta,vmax_sta,
     ,       image2D(mx*my),chi2
      real noise_level, rhs_noise(mdata)
      real gasdev, err, tstd
      data idum /0/

c      double precision gtgdiag(mvar)
c
      parameter(thres=1.d-7,data_thres=100.0, gtgthres=1.d0)
      common/basics/nlevel,nivolu,nvolu,niface,niedge
      common/basics_C_2D/nlevel2D,niedge2D
      data mode0/'Simple    ','Multiscale','Quelling  ','HybridVQLM'/

      write(*,*)' mode:'
      write(*,*)'     0 = Resynthesizing to different level using'
      write(*,*)'          previous results!'
      write(*,*)'     1 = Simple damping;'
      write(*,*)'     2 = MultiScaled transformation:'
      write(*,*)'     3 = Convolution quelling.'
      write(*,*)'     4 = vertically quelling and laterally multiscale.'
      read(*,*)imultis
      if(imultis.lt.0.or.imultis.gt.4)then
       stop' irrelevant mode !'
      elseif(imultis.gt.0)then
       write(*,*)' Dampping factor='
       read(*,*)lamb
      write(*,*) lamb
      endif
 
      print *,' file storing GRAM-MATRIX(eg.G_64_64)='
      read(*,'(a)')input_matrix
      write(6,*)input_matrix
      open(1,file=input_matrix,form='unformatted',err=1001)
c      read(1,err=1001,end=1001)ndat,nlevel(1),nlevel(2),nlevel(3)
      read(1,err=1001,end=1001)ndat,nx,ny,nz
      write(6,*)ndat,nx,ny,nz

      nlevel(1)=nint(log(real(nx-1))/log(2.0))
      nlevel(2)=nint(log(real(ny-1))/log(2.0))
      nlevel(3)=nint(log(real(nz-1))/log(2.0))
      nx0=2**(nlevel(1)-1)+1
      ny0=2**(nlevel(2)-1)+1
      nz0=2**(nlevel(3)-1)+1
      if (nx0.ne.nx.or.ny0.ne.ny.or.nz0.ne.nz) then
         write(6,*)'grid dimensions is not power of 2 +1'
         write(6,*)nx0,ny0,nz0,nx,ny,nz,(nlevel(i),i=1,3)
      endif

      if(ndat.gt.mdata.or.nlevel(1).gt.mlevelx.or.
     .                    nlevel(2).gt.mlevely.or.
     .                    nlevel(3).gt.mlevelz) then
       write(*,*)ndat,nlevel(1),nlevel(2),nlevel(3)
       write(*,*)mdata,mlevelx,mlevely,mlevelz
       stop' inconsistent dim declaration!'
      end if

      do i=1,3
       nlevel_syn(i)=nlevel(i)
       nlevel_tmp(i)=nlevel(i)
      end do
c      nx=2**(nlevel(1)-1)+1
c      ny=2**(nlevel(2)-1)+1
c      nz=2**(nlevel(3)-1)+1
      niface=nx*ny
      niedge(1)=nx
      niedge(2)=ny
      niedge(3)=nz
      do i=1,2
       nlevel2D(i)=nlevel(i)
       niedge2D(i)=niedge(i)
       nlevel_syn2D(i)=nlevel_syn(i)
      enddo
      nivolu=niface*nz
      nvolu=nivolu
      nvar=nvolu
      write(*,*)ndat,(nlevel(i),i=1,3),nx,ny,nz,nvolu
      if(imultis.eq.0)then
       open(79,file='results',status='old',form='unformatted',
     ,                        err=1002)
       read(79,err=1002)ii,nsta,lamb
       mode=mode0(ii)
       write(*,*)' Your levels are:',(nlevel(i),i=1,3)
       write(*,*)' Now what level do you want to synthesize uo to?'
       read(*,*)(nlevel_syn(i),i=1,3)
      elseif(imultis.ge.3)then
       if(imultis.eq.3)then
        write(*,*)' The 2 sigmas for Gaussian quelling:'
        write(*,*)'  # of grids in both lateral and vertical direction='
        read(*,*)sigma0,sigmaz
       elseif(imultis.eq.4)then
        write(*,*)'  # of grids for quelling in the vertical direction:'
        read(*,*)sigmaz
        sigma0=0.
       endif
       sigma0=dabs(sigma0)
       sigmaz=dabs(sigmaz)
       call Gaussian(sigma0,sigmaz)
       mode=mode0(imultis)
       open(79,file='results',form='unformatted',err=1002)
       write(79,err=1002)imultis,nsta,lamb
      else
       mode=mode0(imultis)
       open(79,file='results',form='unformatted',err=1002)
       write(79,err=1002)imultis,nsta,lamb
      endif
c reading noise level
      write(*,*)'add noise level (gaussian random errors)'
      read(*,*)noise_level
      noise_level=noise_level/0.05
      if(imultis.ne.0)then
       rhs_min=1.e7
       rhs_max=-1.e7
       ndat_act=0
       do i=1,ndat
        read(1,err=1001)ncoef,(jdx0(j),coef0(j),j=1,ncoef),rhs
c        write(6,*)i,coef0(1),rhs
c add noise
        rhs_noise(i)=rhs*noise_level*gasdev(idum)*0.25
        rhs=rhs+rhs_noise(i)
c        write(6,*)i,rhs,rhs_noise(i)
        if(abs(rhs).le.data_thres)then
         do j=1,nvar
           image(j)=0.d0
         end do
         do j=1,ncoef
          k=jdx0(j)
           image(k)=dble(coef0(j))
         end do
         if(imultis.eq.2)then
          call face_lift_Hybrid_Nbt(image,2,nlevel_syn)
c          call face_lift_Pt3(image,2,nlevel_syn)
         elseif(imultis.eq.3)then
          call quell(image)
         elseif(imultis.eq.4)then
          call quell(image)
          do iiz=1,nz
           iiz0=(iiz-1)*nx*ny
           do iiy=1,ny
            iiy0=(iiy-1)*nx
            do iix=1,nx
             image2D(iiy0+iix)=image(iiz0+iiy0+iix)
            enddo
           enddo
           call face_liftLC(image2D,2,nlevel_syn2D)
           do iiy=1,ny
            iiy0=(iiy-1)*nx
            do iix=1,nx
             image(iiz0+iiy0+iix)=image2D(iiy0+iix)
            enddo
           enddo
          enddo
         end if
         do j=1,nvar
           ttmax=dmax1(ttmax,dabs(image(j)))
         end do
         ncoef=0
         do j=1,nvar
           tmp=image(j)
c         if(dabs(tmp).gt.thres*ttmax)then
          if(dabs(tmp).gt.thres)then
           ncoef=ncoef+1
           jdx(ncoef)=j
           coef(ncoef)=tmp
          end if
         end do
         call lldrow (coef,jdx,ncoef,dble(rhs))
         rhs_min=amin1(rhs_min,rhs)
         rhs_max=amax1(rhs_max,rhs)
         ndat_act=ndat_act+1
        endif
c        write(17,'(2i8,3f10.3)')i,ndat_act,rhs,rhs_min,rhs_max
       end do
       write(*,*)ndat_act,rhs_min,rhs_max
c      Dampping.....
       do i=1,nvar
        coef(1)=lamb
        jdx(1)=i
        call lldrow (coef,jdx,1,0.d0)
       end do
       close(1)
       call slsqr(nit)
       call lgtsol(vv1,errx,nvar)

       if(imultis.eq.2)then
        do ii=1,nvolu
         image(ii)=vv1(ii)
        enddo
        call face_lift_Hybrid_Nbt(image,-1,nlevel_syn)
c        call face_lift_Pt3(image,-1,nlevel_syn)
        do ii=1,nvolu
         vv1(ii)=image(ii)
        enddo
       elseif(imultis.eq.3)then
        do ii=1,nvolu
         image(ii)=vv1(ii)
        enddo
        call quell(image)
        do ii=1,nvolu
         vv1(ii)=image(ii)
        enddo
       elseif(imultis.eq.4)then
        do ii=1,nvolu
         image(ii)=vv1(ii)
        enddo
        do iiz=1,nz
         iiz0=(iiz-1)*nx*ny
         do iiy=1,ny
          iiy0=(iiy-1)*nx
          do iix=1,nx
           image2D(iiy0+iix)=image(iiz0+iiy0+iix)
          enddo
         enddo
         call face_liftLC(image2D,-1,nlevel_syn2D)
         do iiy=1,ny
          iiy0=(iiy-1)*nx
          do iix=1,nx
           image(iiz0+iiy0+iix)=image2D(iiy0+iix)
          enddo
         enddo
        enddo
        call quell(image)
        do ii=1,nvolu
         vv1(ii)=image(ii)
        enddo
       endif
       write(79)(vv1(i),i=1,nvar)
       write(79)(errx(i),i=1,nvar)
      else
       close(1)
       read(79,err=1002)(vv1(i),i=1,nvar)
       read(79,err=1002)(errx(i),i=1,nvar)
       do ii=1,nvolu
        image(ii)=vv1(ii)
       enddo
       call face_lift_Hybrid_Nbt(image,1,nlevel_tmp)
       call face_lift_Hybrid_Nbt(image,-1,nlevel_syn)
c       call face_lift_Pt3(image,1,nlevel_tmp)
c       call face_lift_Pt3(image,-1,nlevel_syn)
       do ii=1,nvolu
        vv1(ii)=image(ii)
       enddo
      endif
      close(79)

c     Quelled or Multiscaled, the fitting should be examined
c      in the original space

      open(1,file=input_matrix,form='unformatted')
      read(1)ndat,nlevel(1),nlevel(2),nlevel(3)

	write(*,*) ndat, nlevel(1),nlevel(2),nlevel(3)

      fit=0.d0
      fit_tmp=0.d0
      do i=1,ndat
       read(1)ncoef,(jdx0(j),coef0(j),j=1,ncoef),rhs
       rhs=rhs+rhs_noise(i)
       if(abs(rhs).le.data_thres)then
        do j=1,nvar
         image0(j)=0.
        end do
        do j=1,ncoef
         k=jdx0(j)
         image0(k)=coef0(j)
        end do 
        fit_tmp=fit_tmp+dble(rhs)*dble(rhs)
        tmp=0.d0
        do j=1,nvar
         tmp=tmp+dble(image0(j))*vv1(j)
        end do
c	write(*,*) '> ',i,rhs,tmp
        fit=fit+(tmp-dble(rhs))*(tmp-dble(rhs))
       endif

      end do
c      fit=(1.-dsqrt(fit/fit_tmp))*100.
c  Y. Shen, 8/13/2004
c  I define variance as deltT*deltT, thus variance reduction is
c  (1 - {sum[deltT(final)]**2/sum[deltT(original)]**2})*100.  
c
	fit = (1. - (fit/fit_tmp))*100.0
        chi2=fit/dble(ndat)

      close(1)
c      write(*,*)fit

c  Y.Shen, 08/16/04
c  much of the model space is unsampled or poorly samples.  It is
c  more meaningful to evaluate model norm only at the grids that are
c  sampled as measured by Gt*G (see Hung, Shen, Chiao, 2004).
c
c  This should be done separately (BUILDG/gtgdiag.f) and the result to be read in
c       open(81,file='/home/yang/FFT/africa/BUILDG/diag.out')
c       do k=1,nvar
c       read(81,*) ndum,gtgdiag(k)
c       end do
c       close(81)
c  calculate model norm etc.
c
      mean=0.d0
      vmin=1.d7
      vmax=-1.d7
      cov=0.d0
      norm=0.d0
      vmin_sta=1.d7
      vmax_sta=-1.d7
c      do j=1,nvar
c Y.Shen, 8/14/04
c       nk=0
       do j=1,nvolu
c       if (gtgdiag(j).gt.gtgthres) then
c       nk=nk+1
       cov=cov+errx(j)*errx(j)
       norm=norm+vv1(j)*vv1(j)
       mean=mean+vv1(j)
        vmin=dmin1(vv1(j),vmin)
        vmax=dmax1(vv1(j),vmax)
c       endif
        end do
c      norm=dsqrt(norm)/dble(nvar)
c      cov=dsqrt(cov)/dble(nvar)
c      mean=mean/dble(nvar)
       norm = dsqrt(norm/dble(nvolu))
       cov  = dsqrt(cov/dble(nvolu))
       mean=mean/dble(nvolu)
c       norm = dsqrt(norm/dble(nk))
c       cov  = dsqrt(cov/dble(nk))
c       mean=mean/dble(nk)
c       write(*,*) 'number of points with gtgdiag > ',gtgthres, nk

      open(78,file='try.xyz')
      write(78,'(8e12.3)')(vv1(j),j=1,nvolu)
      close(78)

      open(77,file='try.sum')
c      write(77,*)' Your last run was under mode:'
      write(77,*)'*************************************'
      write(77,'(a)')mode
      write(77,*)'*************************************'
c      write(77,*)' Remember that model variance can only be compared'
c      write(77,*)'  under the same mode !!'
      write(77,*)' '
      write(77,100)
100   format(10h var_redu.,10h model var,10h modelnorm,
     ,       6x,4hmean,7x,3hmin,7x,3hmax,6x,4hSmin,6x,4hSmax,6x,4hChi2)
      write(77,'(8e10.2,e12.4)')fit,cov,norm,mean,vmin,vmax,vmin_sta,vmax_sta,chi2
      write(77,*)' '
      write(77,'(18h Dampping factor: ,e10.3)')lamb
      if(imultis.eq.3)write(77,'(16h Sigma0,sigmaz= ,2f6.1)')sigma0,sigmaz
      close(77)
c      write(*,*)' Contents in file try.sum:'
c      write(*,*)' Your last run was under mode:'
      write(*,*)'*************************************'
      write(*,'(a)')mode
      write(*,*)'*************************************'
c      write(*,*)' Remember that model variance can only be compared'
c      write(*,*)'  under the same mode !!'
      write(*,*)' '
      write(*,100)
      write(*,'(8e10.2,e12.4)')fit,cov,norm,mean,vmin,vmax,vmin_sta,vmax_sta,chi2
      write(*,*)' '
      write(*,*)' Dampping factor: ',lamb
      if(imultis.ge.3)write(*,*)' Sigma0, sigmaz= ',sigma0,sigmaz

      stop
1001  stop' ERROR reading or opening G-file!'
1002  stop' ERROR reading or opening old results file!'
400   stop' ERROR reading or opening station_position !'
      end


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





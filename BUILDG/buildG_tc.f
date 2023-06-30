c     buildG_ray.f: build Gram matrix in LSQR fashion 
c               for regional Cartesian tomography
c
      character*150 ib,fray,fdat,fscr,fsta
      character*1 ia
      parameter(mlx=7,mly=7,mlz=7,
     ,          mnx=2**(mlx-1)+1,mny=2**(mly-1)+1,
     ,          mnz=2**(mlz-1)+1,
     ,          mvolu=mnx*mny*mnz,mdata=200000,mpath=2000,
     ,          mstation=1000,mevent=3000,mvar=mvolu+mstation+mevent)
      dimension npaths(mpath),locat(mpath),
     ,          locat_g(mvolu),judge(mvolu),jdx(mvar)
      double precision path(mpath)
      double precision lat0,lon0,lat1,lon1,dp0,dp1,thres
      double precision coef_ref(mvar),coef(mvar),coef_LSQR(mvar)
      parameter(thres=1.d-8)
      double precision x0,y0,z0,dx,dy,dz
      double precision trans(3,3)
      dimension ist(mdata),ist0(mdata),iev(mevent)
      real res(mdata),ttres,et,errtt(mdata)
      real coef_real(mvar)
      real tmp
      character*10 xtmp
      common/meshb1/x0,y0,z0,dx,dy,dz,nx,ny,nz,nface,nvolu
      common/meshb2/trans
      data coef/mvar*0.0d0/
      include 'lsqrcom1'

      write(6,*)'enter shootrays file'
      read(5,'(a)') fray
      write(6,*)'enter data and stations # file (P.xxx.datsta)'
      read(5,'(a)') fdat
      open(72,file=fdat,err=200)
      read(72,*,err=200)ndata,nstation,nevent
      close(72)
      write(6,*)'station terms included (1=yes)?'
      read(5,*)istacorflag
      if (istacorflag.eq.1) then
         write(6,*)'station and event correction file (P.xxx.stacor)'
         read(5,'(a)') fscr
         open(72,file=fscr,err=400)
         do i=1,ndata
            read(72,*,err=400,end=400)
     ,      idata,iev(i),jsta,ist(i),t1,t2
         enddo
         close(72)
         if(idata.ne.ndata)then
           write(*,*)idata,ndata,nevent
           stop' inconsistent setup'
           endif
           write(*,*)ndata,nevent,nstation
      endif
      write(6,*)'enter the file with stations used'
      read(5,'(a)') fsta
      open(72,file=fsta)
      hmean=0.
      read(72,*)nsta
      do i=1,nsta
         read(72,*)xtmp,xtmp,tmp,tmp,selv
         hmean=hmean+selv
      enddo
        hmean=hmean/nsta*0.001
        write(6,*)'average bathymetery of stations in km='
        write(6,*)hmean
      close(72)
      open(72,file='mesh.config',err=71)
      read(72,'(6i5)',err=71,end=71)lx,ly,lz,nx,ny,nz
      read(72,'(6f15.5)',err=71,end=71)x0,y0,z0,dx,dy,dz
      read(72,'(3f15.5)',err=71,end=71)((trans(i,j),j=1,3),i=1,3)
      close(72)
c      call linit
      nface=nx*ny
      nvolu=nface*nz
      if (istacorflag.eq.1) then
      nvar=nvolu+nstation+nevent
      else
      nvar=nvolu
      endif

      open (1,file=fray,err=100)

      open (3,file='Ray.summary')
      open (21,file='Kernel.pos')
      open(11,file='G_ray_d',form='unformatted')
      write(11)ndata,nx,ny,nz,nstation,nevent
      open(12,file='Gw_ray_d',form='unformatted')
      write(12)ndata,nx,ny,nz,nstation,nevent
      ndata0=ndata
      ndata=0
      do idat=1,mdata
       read(1,'(a80)',err=100,end=100)ib
       if(ib(1:1).ne.'>')goto 100
       write(21,'(a)')ib
       read(1,'(i3)',err=100,end=100)nsta
       write(21,'(i3)')nsta
       do jsta=1,nsta
        do i=1,nvolu
          coef(i)=0.d0
        end do
        nseg=0
        read(1,'(a)',err=100,end=100)ib
c        read(ib,'(29x,2f10.4,f7.1)')atl,onl,dph
c read event, station positions and residuals
        read(ib,'(29x,f9.4,2f10.4,f7.1,38x,f9.4,f9.5,f9.3,f9.4)')
     &   hs,atl,onl,dph,ttres,sigmatt,ttprd,tc
        tc=(hs-hmean)/hs*tc
c        write(6,*)'tc=',hmean,hs,ttres,tc
c reset error of tt to be 0.05 if error<0.05
        if (sigmatt.lt.0.05) sigmatt=0.05
c write the reset error in the Kernel.pos file,
c then don't need to do weighting while calculating kernel values.
        write(21,'(a112,f9.5,f9.3,f9.4)')ib(1:112),sigmatt,ttprd,tc
        dph=6371.0-dph
        nlocat_g=0
        do ivolu=1,nvolu
         judge(ivolu)=0
        end do
1       continue
        read(1,'(f10.3,2f10.4)',err=300,end=300)dp,at,on
        dp=6371.0-dp
        if(nseg.eq.0)then
         lat0=dble(at)
         lon0=dble(on)
         dp0=dble(dp)
         ats1=at
         ons1=on
         dps1=dp
         read(1,'(f10.3,2f10.4)',err=100,end=100)dp,at,on
         dp=6371.0-dp
         lat1=dble(at)
         lon1=dble(on)
         dp1=dble(dp)
         id=0
        else
         lat0=lat1
         lon0=lon1
         dp0=dp1
         lat1=dble(at)
         lon1=dble(on)
         dp1=dble(dp)
         ats2=at
         ons2=on
         dps2=dp
        endif
        nseg=nseg+1
        call inside(lon0,lat0,dp0,lon1,lat1,dp1,id,npath,npaths,path,
     ,              nlocat,locat)
c       accumulating ray path for G_ray
        if(npath.gt.0)then
         do ipath=1,npath
           coef(npaths(ipath))=coef(npaths(ipath))+path(ipath)
         end do
        endif       
c       find gaussian points for G_kernel
        do ilocat=1,nlocat
         iz=(locat(ilocat)-1)/((nx-1)*(ny-1))+1
         ixy=locat(ilocat)-(iz-1)*(nx-1)*(ny-1)
         iy=(ixy-1)/(nx-1)+1
         ix=ixy-(iy-1)*(nx-1)
         ix1=ix-4
         if(ix1.lt.1)ix1=1
         ix2=ix+4
         if(ix2.gt.nx)ix2=nx
         iy1=iy-4
         if(iy1.lt.1)iy1=1
         iy2=iy+4
         if(iy2.gt.ny)iy2=ny
         iz1=iz-4
         if(iz1.lt.1)iz1=1
         iz2=iz+4
         if(iz2.gt.nz)iz2=nz
         do iiz=iz1,iz2
          i1=(iiz-1)*(nx-1)*(ny-1)
          do iiy=iy1,iy2
           i2=(iiy-1)*(nx-1)
           do iix=ix1,ix2
            i3=i1+i2+iix
            if(judge(i3).eq.0)then
             nlocat_g=nlocat_g+1
             locat_g(nlocat_g)=i3
             judge(i3)=1
            endif
           end do
          end do
         end do
        end do
        goto 1
300     continue
        write(21,'(i8)')nlocat_g
        write(21,'(12i8)')(locat_g(ilocat),ilocat=1,nlocat_g)
        ndata=ndata+1
         if(ndata.gt.ndata0) then
           write(6,*)ndata,ndata0
           stop' inconsistent ray counts!'
         endif
         ncoef=0
         res(ndata)=ttres-tc
         errtt(ndata)=sigmatt
         do i=1,nvolu
          if(dabs(coef(i)).gt.thres)then
           ncoef=ncoef+1
           jdx(ncoef)=i
           coef_LSQR(ncoef)=coef(i)
          endif
         end do

c         call lldrow(coef_LSQR,jdx,ncoef)
         if(ndata/100*100.eq.ndata)write(*,*)ndata,idat,jsta,nseg,ncoef
c          write(*,*)ndata,idat,jsta,nseg,ncoef
         write(11)
     ,    ncoef,(jdx(k),sngl(coef_LSQR(k)),k=1,ncoef),
     ,    res(ndata)
         write(3,'(i6,2f15.5,f9.3,2f15.5,f9.3)')
     ,    ndata,ons1,ats1,dps1,ons2,ats2,dps2
c         write(3,'(a)')ib
        backspace(1)
       end do
      end do
100   continue
c      write(*,*)ndata,idat,jsta,nseg
c read G_raw and build G_ray_d
c read G_raw and build G_ray_d
        if (ndata.ne.ndata0) print *,'inconsistent # of data',ndata,ndata0
        rewind(11)
        read(11) ndata,nx,ny,nz,nstation,nevent
        do l=1,ndata
           et=errtt(l)
           read(11)ncoef,(jdx(k),coef_real(k),k=1,ncoef),restt
           write(12)ncoef,(jdx(k),coef_real(k)/et,k=1,ncoef),res(l)/et
        enddo

        close(11)
	close(12)
        close(3)

c      write(72,'(i8)')ndata
      stop
71    stop' Error in opening or reading mesh.config file !'
200   stop' Error in opening or reading P.l.ice.datsta !'
400   stop' Error in opening or reading P.l.ice.stacor !'
      end


      subroutine lldrow(coef,jdx,ncoef)
c ... add one equation of data into row of matrix
c     for conjugate gradient method (slsqr) of solution
      dimension jdx(ncoef)
      double precision  coef(ncoef)
      parameter(mlx=7,mly=7,mlz=7,
     ,          mnx=2**(mlx-1)+1,mny=2**(mly-1)+1,
     ,          mnz=2**(mlz-1)+1,
     ,          mvolu=mnx*mny*mnz,mdata=200000,mpath=2000,mst=1000,mevt=3000,
     ,          mvar=mvolu+mst+mevt)
      include 'lsqrcom1'

      if (ncoef.gt.0) then
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
       end do
      endif
      return
      end


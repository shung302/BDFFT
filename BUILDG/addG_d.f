c     buildG_ray.f: build Gram matrix in LSQR fashion 
c               for regional Cartesian tomography
c
      character*125 fgk(20),fref(20),fdatsta
      character*1 ia
      parameter(mlevelx=7,mlevely=7,mlevelz=7,
     ,          mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,          mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,          mbase=1,mvolu=mivolu*mbase,mdata=100000,mstation=500,
     ,          mvar=mvolu+mstation)
      dimension npaths(mdata),locat(mdata),
     ,          locat_g(mvolu),judge(mvolu),jdx(mvar)
      real coef_real(mvar)
      parameter(thres=1.d-7)
      dimension res_tt(mdata)
      data nstations /500/

      open(72,file='mesh.config',err=71)
      read(72,'(6i5)',err=71,end=71)lx,ly,lz,nx,ny,nz
      read(72,'(6f15.5)',err=71,end=71)x0,y0,z0,dx,dy,dz
c      write(6,*)nx,ny,nz
      nface=nx*ny
      nvolu=nface*nz
      close(72)

      ngk=0
      write(6,*)'enter P.xxx.datsta file'
      read(5,'(a)') fdatsta
      write(6,*)'add G matrix for several freq bands'
1     write(6,*)'enter G_ray_d or G_kernel_dfile or none'
      read(5,'(a)') fgk(ngk+1)
      if ( fgk(ngk+1)(1:4) .ne. 'none') then 
          ngk=ngk+1
          go to 1
      endif
     
      open(1,file=fdatsta)
      read(1,*)ndata,nstation,nevent
      ndata_count=0
      do i=1,ngk
         open(12,file=fgk(i),form='unformatted')
         read(12)ndata0,nx0,ny0,nz0,nstation0,nevent0
	 write(6,*)fgk(i),ndata0,nx0,ny0,nz0
         do l=1,ndata0
            read(12)ncoef,(jdx(k),coef_real(k),k=1,ncoef),res_tt(l)
         enddo
c           write(6,*)(coef_real(k),k=1,10)
         close(12)
         ndata_count=ndata_count+ndata0
      enddo
      write(6,*)'total data=',ndata_count,nx,ny,nz
      if (ndata_count.ne.ndata) then
         write(6,*)'# of data is not consitent'
         go to 71
      endif

      open(22,file='G_d_add',form='unformatted')
      write(22)ndata,nx,ny,nz,nstation,nevent
      nd=0
      do i=1,ngk
         open(12,file=fgk(i),form='unformatted')
         read(12)ndata0,nx0,ny0,nz0,nstation,nevent
         write(6,*)ndata0,nx0,ny0,nz0,nstation,nevent
         do l=1,ndata0
            read(12)ncoef,(jdx(k),coef_real(k),k=1,ncoef),res_tt(l)
            nd=nd+1
            write(22)ncoef,(jdx(k),coef_real(k),k=1,ncoef),res_tt(l)
         enddo
         close(12)
      enddo
      if (nd.ne.ndata) write(6,*)'warning nd and ndata not equal=',nd,ndata
      close(22)

71    stop
      end

c     addG_sta_evtcor: add station and event corrections to the built-Gram matrix in LSQR fashion 
c               for regional Cartesian tomography
c
      character*125 ib,fscr,fray,fdat,fg,outfile
      character*1 ia
      parameter(mlx=7,mly=7,mlz=7,
     ,          mnx=2**(mlx-1)+1,mny=2**(mly-1)+1,mnz=2**(mlz-1)+1,
     ,          mvolu=mnx*mny*mnz,mdata=200000,mpath=2000,mevent=10000,
     ,          mstation=1000,mvar=mvolu+mevent+mstation)
      dimension npaths(mpath),locat(mpath),
     ,          locat_g(mvolu),judge(mvolu),jdx(mvar)
      real coef_real(mvar)
      dimension ist(mdata),ist0(mdata),iev(mdata),res(mdata),res0(mdata)
      parameter(thres=1.d-8)
      dimension res_tt(mdata)
     
      write(6,*)'enter station correction file (P.xxx.stacor)'
      read(5,'(a)') fscr
      write(6,*)'enter data and stations # file (P.xxx.datsta)'
      read(5,'(a)') fdat
      write(6,*)'add event and station correction in G matrix'
      write(6,*)'enter G file (G_ray_d_Pxxxx or G_kernel_d_Pxxxx)'
      read(5,'(a)') fg
      lb=index(fg,'G')
      le=index(fg,'_d_')+3
      ll=index(fg,' ')-1

      open(72,file=fdat,err=200)
      read(72,*,err=200)ndata,nstation,nevent
      write(6,*)'ndata,nstation,nevent=',ndata,nstation,nevent
      close(72)
      open(72,file=fscr,err=400)
      do i=1,ndata
       read(72,*,err=400,end=400)
     ,      idata,iev(i),j,ist(i),t1,t2
       res(i)=t1
      enddo
      close(72)

      open(72,file='mesh.config',err=71)
      read(72,'(6i5,3f15.5)',err=71,end=71)lx,ly,lz,nx,ny,nz
      read(72,'(6f15.5)',err=71,end=71)x0,y0,z0,dx,dy,dz
      close(72)

      nface=nx*ny
      nvolu=nface*nz
      nvar=nvolu+nstation+nevent

      open(12,file=fg,form='unformatted')
      read(12)ndata0,nx,ny,nz,nstation,nevent
      write(6,*)fg
      write(6,*)ndata0,nx,ny,nz,nstation,nevent
      if (ndata.ne.ndata0) then
         write(6,*)'data # not equal',ndata,ndata0
         goto 400
      endif
      outfile=fg(lb:le-1)//'se_'//fg(le:ll)
c	write(6,*)lb,le,ll,outfile
      open(22,file=outfile,form='unformatted')
      write(22)ndata,nx,ny,nz,nstation,nevent
      do l=1,ndata
         read(12)ncoef,(jdx(k),coef_real(k),k=1,ncoef),res_tt(l)
c         write(6,*)ncoef,res_tt(l)
         ncoef=ncoef+1
         jdx(ncoef)=nvolu+ist(l)
         coef_real(ncoef)=1.
         ncoef=ncoef+1
         jdx(ncoef)=nvolu+nstation+iev(l)
         coef_real(ncoef)=1.
         write(22)ncoef,(jdx(k),coef_real(k),k=1,ncoef),res_tt(l)
      enddo
      close(12)
      close(22)
      stop
71    stop'error in reading mesh.config'
200   stop'error in opening file'
400   stop'Error in reading G_ray_d or G_kernel file'
      stop
      end

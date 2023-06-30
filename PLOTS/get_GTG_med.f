c     buildG_ray.f: build Gram matrix in LSQR fashion 
c               for regional Cartesian tomography
c
      character*250 ib,input_matrix,outfile,outfile1
      character*1 ia
      parameter(mlevelx=7,mlevely=7,mlevelz=7,
     ,          mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,          mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,          mbase=1,mvolu=mivolu*mbase,mdata=100000,mstation=200,
     ,          mvar=mvolu+mstation)
      dimension npaths(mdata),locat(mdata),
     ,          locat_g(mvolu),judge(mvolu),jdx(mvar),ncounts(mvolu)
      real coef_real(mvar),plg2(mvolu),res_tt(mdata),thres,plglog(mvolu)
     ,     ,plg(mvolu),plg2log(mvolu)
      dimension indx(mvolu)
      parameter(thres=1.e-7)

      open(72,file='mesh.config',err=71)
      read(72,'(6i5)',err=71,end=71)lx,ly,lz,nx,ny,nz
      read(72,'(6f15.5)',err=71,end=71)x0,y0,z0,dx,dy,dz
      write(6,*)nx,ny,nz
      nface=nx*ny
      nvolu=nface*nz
      close(72)

      write(6,*)'give G matrix name='
      read(5,'(a)') input_matrix
      write(6,*)'output GTG name='
      read(5,'(a)') outfile
      write(6,*)'output ray count name='
      read(5,'(a)') outfile1

      ndata=0
      open(12,file=input_matrix,form='unformatted')
      read(12)ndat,nx0,ny0,nz0
      write(*,*) ndat, nx0, ny0, nz0
        write(6,*)nx,ny,nz
c      if (nx0.ne.nx.or.ny0.ne.ny.or.nz0.ne.nz) then
c       stop' inconsistent dim declaration!'
c      endif

      do i=1,nvolu
         plg2(i)=0.
         plg(i)=0.
         plg2log(i)=0.
         plglog(i)=0.
         ncounts(i)=0
      enddo
      do l=1,ndat
         read(12)ncoef,(jdx(k),coef_real(k),k=1,ncoef),res_tt(l)
c         write(*,*) ncoef,res_tt(l)
c	 write(71,*)ncoef,res_tt(l),(coef_real(k),k=1,10)
c110	 continue
         do k=1,ncoef
            plg2(jdx(k))=plg2(jdx(k))+coef_real(k)*coef_real(k)
            if (abs(coef_real(k)).ge.thres) then
               ncounts(jdx(k))=ncounts(jdx(k))+1 
            endif
         enddo
      enddo
      close(12)
      open(22,file=outfile)
c      do i=1,nvolu
c           write(6,*)i,plg2(i)
c      enddo

      call indexx(mvolu,nvolu,plg2,indx)
      nmed=int((nvolu+1)/2)
      pmed=plg2(indx(nmed))

      write(6,*)'median=',nmed, pmed
      open(22,file=outfile)
      write(22,'(8e12.3)')(plg2(i)/pmed,i=1,nvolu)
      close(22)
      open(22,file=outfile1)
      write(22,'(8i12)')(ncounts(i),i=1,nvolu)
      close(22)
71    stop
      end

c sorting
      SUBROUTINE indexx(nmax,n,arr,indx)
      INTEGER n,indx(nmax),M,NSTACK
      REAL arr(nmax)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) print *, 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END


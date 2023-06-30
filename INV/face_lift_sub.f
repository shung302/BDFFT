c___________________________________________________________________________
c
c  CONTENTS:
c
c    face_lift_Nhelix: node based, lifted linear, helix transformed
c                     (instead of tensor products)
c       HELIX version of multidimensional WT using helix 1D WT
c        note that the present image is organized as (((1,nx),ny),nz)
c        when passed in.  It is then reorganized from the longest strip
c        that means the basic strip consists in the direction with the
c        highest nlevel and correspondingly, the periodic B.C. will be
c        set. For example, nlevel_x=5, nlevel_y=4, nllevel_z=3,then
c        the implicitly assumed periodic B.C. on is on YZ-edges.
c        (compare with the Tensor product version in <face_lift_Nt1>)
c        (N in _Nhelix means node based <instead of pixel-based>)

c    face_lift_Pt1: pixel based Harr,   tensor products.
c        (common/basics/nlevel,nivolu,nvolu,niface,niedge)

c    face_lift_Pt3: pixel based Spline, tensor products.
c        (common/basics/nlevel,nivolu,nvolu,niface,niedge)

c    face_lift_Pb1: pixel based Harr, 3D cube (not tensor products),
c                    need to run TRANSPLANT first.
c        (common/basics/nlevel,nivolu,nvolu,niface,niedge)

c      TRANSPLANT : configurating 3D cubes
c        (common/basics/nlevel,nivolu,nvolu,niface,niedge)

c    face_lift_Nt1: node based, lifted linear, tensor products
c                    (in fact, rotate X,Y,Z succesively transforming)
c                    notice that nx=2**(lx-1)+1 not 2**(lx-1) in all
c                     the pixel based scheme.
c        (common/basics/nlevel,nivolu,nvolu,niface,niedge)

c    face_lift_Hybrid_Nbt: Hybrid,  laterally  Nb1 (triangular meshes) +
c                           vertically Nt1 (All lifted Linear node based).
c         common/basics/nlevel,nivolu,nvolu,niface,niedge
c         common/basics_C_2D/nlevel2D,niedge2D

c      face_liftLC:  2D,laterally  Nb1 (triangular meshes)
c         common/basics_C_2D/nlevel2D,niedge2D

c    face_lift_Hybrid_Nbt_spherical: same as above but laterally on
c                                    spherical surface (geodesic triangles
c                                    seeded from icosahedron).
c         common/basics_spherical/nlevel,nvolu,niedge
c         common/base0/nbase,nlevelS,niface,nivert,nedge_head,
c     ,                nface,nvert,nvert_global,icc,icq
c         common/base1/face,face_global,indx,face_conf,
c     ,                face_index,node_conf

c      face_liftL: 2D spherical geodesic triangles seeded from icosahedron
c         common/base0/nbase,nlevelS,niface,nivert,nedge_head,
c     ,                nface,nvert,nvert_global,icc,icq
c         common/base1/face,face_global,indx,face_conf,
c     ,                face_index,node_conf
c
c     Notice that both of above (face_lift_Hybrid_Nbt & 
c      face_lift_Hybrid_Nbt_spherical) rotate seperately through each level
c      within lateral and vertical direction first than switch to the other
c      direction;
c      Another two in the following, on the other hand, go through each
c       direction first, then cycle through each level
c

c     face_lift_Nb_tera: similar to face_lift_Pt3 but is node-based,
c      instead of pixel-based, 3D tessellation built upon tetrahedron
c      elements

c     quell(image): 3D Cartesian convolutional quelling based on
c      Gaussian func. preprocessor needed=<Gaussian>
c          common/basics/nlevel,nivolu,nvolu,niface,niedge

c       Gaussian(sigma,sigmaz): preprocessor, building 3D weighting function
c        for convolutional quelling based on Gaussian func.
c          common/basics/nlevel,nivolu,nvolu,niface,niedge

c     quellS(nn,y): 3D spherical convolutional quelling
c       (nn=nvolu,y=image), need to run <prequell> first.
c          common/basics_spherical/nlevel,nvolu,niedge

c       prequell(vertg_x,sigma0)
c          common/basics_spherical/nlevel,nvolu,niedge
c          common/base0/nbase,nlevelS,niface,nivert,nedge_head,
c     ,                 nface,nvert,nvert_global,icc,icq
c          common/base1/face,face_global,indx,face_conf,
c     ,                 face_index,node_conf

c___________________________________________________________________________

c    call face_lift(image,isign,nlevel_syn)

c    isign:  1,  IMAGE analysis (decomposition)
c           -1,  IMAGE synthesis (reconstruction)
c            0,  doing nothing
c            2,  PATH-integral analysis
c           -2,  PATH-integral synthesis

c    notice that in the present tomographic study, each row of
c     the G-matrix (each path scoring within highest level faces)
c     is first decomposed (2); the resulted model image after
c     inversion is then synthesized (-1)
c
c    notice also that, face_conf contains the reorgnized block
c     numbering in scale hirachi, when the original numbering is
c     the usual:
c
c         (((1,nx),ny,nz)
c
c    compare with the pyrimal form in spherical wavelet (or other
c       plannar triangulation grid)
c
c    New version of decomposition and synthesizing is performed successively
c     in the order of (((1,nx),1,ny),1,nz).  In other words, preprocessing
c     by transplant is dropped, the reorganized index FACE_CONF is no more
c     needed.  Also, it is now possible to have different refinement level
c     in X-,Y-,and Z- direction  <2000/11/10>.

c    biorthogonal bior3.1 wavelet <2001/4/3>

c    change the order of processing (rotate X,Y,Z first within each level)
c    <2001/4/13> (lift_b31.f vs. lift_b31_level.f)

c    New version of LINEAR VERTEX-BASED INTERPOLATORY wavelets             
c       <2002/1/24>                                                        
c       !notice that for reconstruction (synthesizing),                    
c       !the second (lifting) phase has to be undone first                 
c       !(the 2nd phase is the lower procedure in Fig.2 of                 
c       !Sweldens' 5 minute tour)                                          
c       ! i.e.,                                                            
c       !the order of operation in Fig.2 of Sweldens's 5                   
c       !minute tour is switched for DECOMPOSITION and                     
c       !SYNTHESIZING

c    Due to the undesired characteristics of tensor products, it is decided
c     that 2D triangular mesh and 3D tetrahedral tellesation is preferable
c     when invoked in inverse formulation
c     this version is the lateral triangular 2D version
c     that might be invoked with hybrid vertical queeling + lateral multiscale
c     <2002/8/25>

c___________________________________________________________________________


c___________________________________________________________________________
       subroutine face_lift_Pt1(image,isign,nlevel_syn)
       implicit none
c
c      pixel based Harr, tensor products
c
 
       integer mlevelx, mlevely, mlevelz, mlevelxyz
       integer mivolu, mbase, mvolu, nbase
       parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,           mlevelxyz=10,
     ,           mivolu=2**(mlevelx-1)*2**(mlevely-1)*2**(mlevelz-1),
     ,           mbase=1,mvolu=mivolu*mbase,nbase=mbase)
       integer nlevel(3),isign,nlevel_syn(3),niedge(3),
     ,         nivolu,iivolu0,iivolu,ivolu,nvolu,niface,
     ,         i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_syntmp(3),i3,ii0
       double precision image(mvolu)
       double precision y(mivolu),tmp(2**(mlevelxyz-1))
       common/basics/nlevel,nivolu,nvolu,niface,niedge

c      bior_1.1
       integer leng_filt
       parameter(leng_filt=2)
       integer iab(leng_filt),iab1(leng_filt),iab2(leng_filt)
      double precision hi_D1(leng_filt),lo_D1(leng_filt),ab(leng_filt),
     ,       hi_R1(leng_filt),lo_R1(leng_filt),filter(leng_filt,2),
     ,       hi_D2(leng_filt),lo_D2(leng_filt),ab1(leng_filt),
     ,       hi_R2(leng_filt),lo_R2(leng_filt),ab2(leng_filt)
       data lo_D1/  0.5d0,  0.5d0/
       data hi_D1/ -0.5d0,  0.5d0/
       data lo_R1/  1.0d0,  1.0d0/
       data hi_R1/  1.0d0, -1.0d0/
       data lo_D2/  1.0d0,  1.0d0/
       data hi_D2/ -1.0d0,  1.0d0/
       data lo_R2/  0.5d0,  0.5d0/
       data hi_R2/  0.5d0, -0.5d0/

c      bior_3.1
c       parameter(leng_filt=4)
c       integer iab(leng_filt),iab1(leng_filt),iab2(leng_filt)
c      double precision hi_D1(leng_filt),lo_D1(leng_filt),ab(leng_filt),
c     ,       hi_R1(leng_filt),lo_R1(leng_filt),filter(leng_filt,2),
c     ,       hi_D2(leng_filt),lo_D2(leng_filt),ab1(leng_filt),
c     ,       hi_R2(leng_filt),lo_R2(leng_filt),ab2(leng_filt)
c       data lo_D1/  -0.25d0,   0.75d0,   0.75d0,  -0.25d0/
c       data hi_D1/ -0.125d0,  0.375d0, -0.375d0,  0.125d0/
c       data lo_R1/   0.25d0,   0.75d0,   0.75d0,   0.25d0/
c       data hi_R1/   -0.5d0,   -1.5d0,    1.5d0,    0.5d0/
c       data lo_D2/   0.25d0,   0.75d0,   0.75d0,   0.25d0/
c       data hi_D2/    0.5d0,    1.5d0 ,  -1.5d0,   -0.5d0/
c       data lo_R2/  -0.25d0,   0.75d0,   0.75d0,  -0.25d0/
c       data hi_R2/  0.125d0, -0.375d0,  0.375d0, -0.125d0/

c      local variables 
       integer k, n
       integer i0, itmp
       integer ix0, ixs, ixs1
       integer iy0, iys, iys1
       integer iz0, izs, izs1
       integer ibase
       integer leng0, leng_filt0
       integer niedge12, nlevelT

       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       elseif(isign.gt.0)then
        do i=1,leng_filt
         if(isign.eq.1)then
          filter(i,1)=lo_D1(i)
          filter(i,2)=hi_D1(i)
         else
          filter(i,1)=lo_D2(i)
          filter(i,2)=hi_D2(i)
         endif
        enddo
       elseif(isign.lt.0)then
        do i=1,leng_filt
         if(isign.eq.-1)then
          filter(i,1)=lo_R1(i)
          filter(i,2)=hi_R1(i)
         else
          filter(i,1)=lo_R2(i)
          filter(i,2)=hi_R2(i)
         endif
        enddo
       end if
       niedge12=niedge(1)*niedge(2)
       leng0=leng_filt/2
       if(leng0*2.ne.leng_filt)leng0=leng0+1
       leng_filt0=leng_filt
       nlevelT=max0(nlevel(1),nlevel(2),nlevel(3))
 
       do i0=1,nbase
        iivolu0=(i0-1)*nivolu
 

        if(isign.gt.0)then
 
c        Decomposition or Analyzing....
c        S(j,k)=SUM_n[lo_D(2k+leng0-n)*S(j-1,n)]
c        d(j,k)=SUM_n[hi_D(2k+leng0-n)*S(j-1,n)]

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          y(iivolu)=image(ivolu)
         end do
 
         do j=2,nlevelT
          i1=2**(j-1)
          ic=i1/2

c         Starting from X-decomposition.... 
          if(j.le.nlevel(1))then
           i2=niedge(1)/i1
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do iys=1,niedge(2)
             iy0=(iys-1)*niedge(1)
             ixs1=iz0+iy0
             do itmp=1,niedge(1)
              tmp(itmp)=y(ixs1+itmp)
             end do
             do k=1,i2
              call forward(k,ic,niedge(1),leng0,leng_filt0,
     ,                     ab1,filter,tmp)
              ibase=(k-1)*i1+1
              itmp=ixs1+ibase
              y(itmp)=ab1(1)
              y(itmp+ic)=ab1(2)
             end do
            end do
           end do
          end if

c         Followed by Y-decomposition....
          if(j.le.nlevel(2))then
           i2=niedge(2)/i1
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do ixs=1,niedge(1)
             ix0=ixs
             iys1=iz0+ix0
             do itmp=1,niedge(2)
              tmp(itmp)=y(iys1+(itmp-1)*niedge(1))
             end do
             do k=1,i2
              call forward(k,ic,niedge(2),leng0,leng_filt0,
     ,                     ab1,filter,tmp)
              ibase=(k-1)*i1+1
              itmp=iys1+(ibase-1)*niedge(1)
              y(itmp)=ab1(1)
              y(itmp+ic*niedge(1))=ab1(2)
             end do
            end do
           end do
          end if 

c         And followed by Z-decomposition....
          if(j.le.nlevel(3))then
           i2=niedge(3)/i1
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            do ixs=1,niedge(1)
             ix0=ixs
             izs1=iy0+ix0
             do itmp=1,niedge(3)
              tmp(itmp)=y(izs1+(itmp-1)*niedge12)
             end do
             do k=1,i2
              call forward(k,ic,niedge(3),leng0,leng_filt0,
     ,                     ab1,filter,tmp)
              ibase=(k-1)*i1+1
              itmp=izs1+(ibase-1)*niedge12
              y(itmp)=ab1(1)
              y(itmp+ic*niedge12)=ab1(2)
             end do
            end do
           end do
          end if

         end do
 
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do
 
        elseif(isign.lt.0)then
 
c        Synthesizing.....
c        S(j-1,n)=SUM_k[lo_R(n-2k+leng0+1)*S(j,k)+
c                       hi_R(n-2k+leng0+1)*d(j,k)]

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          y(iivolu)=image(ivolu)
         end do

c        first take out all contribution above the specfiled synthesizing level
         nlevel_syntmp(3)=nlevel_syn(3)-1
         i1=nlevel(3)-nlevel_syntmp(3)
         i2=2**(i1-1)
         i3=niedge(3)/i2
         if(nlevel_syntmp(3).lt.0)then
          write(*,*)' Strange Z-synthesizing level !',nlevel_syn(3)
          stop
         else
          do iys=1,niedge(2)
           iy0=(iys-1)*niedge(1)
           do ixs=1,niedge(1)
            ix0=ixs
            izs1=iy0+ix0
            do i=1,i3
             ii0=(i-1)*i2
             do j=2,i2
              y(izs1+(ii0+j-1)*niedge12)=0.d0
             end do
            end do
           end do
          end do
         endif
         nlevel_syntmp(2)=nlevel_syn(2)-1
         i1=nlevel(2)-nlevel_syntmp(2)
         i2=2**(i1-1)
         i3=niedge(2)/i2
         if(nlevel_syntmp(2).lt.0)then
          write(*,*)' Strange Y-synthesizing level !',nlevel_syn(2)
          stop
         elseif(i2.ge.2)then
          do izs=1,niedge(3)
           iz0=(izs-1)*niedge12
           do ixs=1,niedge(1)
            ix0=ixs
            iys1=iz0+ix0
            do i=1,i3
             ii0=(i-1)*i2
             do j=2,i2
              y(iys1+(ii0+j-1)*niedge(1))=0.d0
             end do
            end do
           end do
          end do
         endif
         nlevel_syntmp(1)=nlevel_syn(1)-1
         i1=nlevel(1)-nlevel_syntmp(1)
         i2=2**(i1-1)
         i3=niedge(1)/i2
         if(nlevel_syntmp(1).lt.0)then
          write(*,*)' Strange X-synthesizing level !',nlevel_syn(1)
          stop
         elseif(i2.ge.2)then
          do izs=1,niedge(3)
           iz0=(izs-1)*niedge12
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            ixs1=iz0+iy0
            do i=1,i3
             ii0=(i-1)*i2
             do j=2,i2
              y(ixs1+ii0+j)=0.d0
             end do
            end do
           end do
          end do
         endif

         do j=nlevelT,2,-1
          i1=2**(j-1)
          ic=i1/2

c         The order of synthesizing reverses the order of analyzing
c         Starting from Z_synthesizing.....
          if(j.le.nlevel(3))then
           i2=niedge(3)/ic
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            do ixs=1,niedge(1)
             ix0=ixs
             izs1=iy0+ix0
             do itmp=1,niedge(3)
              tmp(itmp)=y(izs1+(itmp-1)*niedge12)
             enddo 
             do n=1,i2
              call reverse(n,i1,ic,niedge(3),leng0,leng_filt0,
     ,                    ab1,filter,tmp)
              ibase=(n-1)*ic+1
              itmp=izs1+(ibase-1)*niedge12
              y(itmp)=ab1(1)+ab1(2)
             end do
            end do
           end do
          endif

c         Followed by Y-synthesizing...
          if(j.le.nlevel(2))then
           i2=niedge(2)/ic
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do ixs=1,niedge(1)
             ix0=ixs
             iys1=iz0+ix0
             do itmp=1,niedge(2)
              tmp(itmp)=y(iys1+(itmp-1)*niedge(1))
             enddo
             do n=1,i2
              call reverse(n,i1,ic,niedge(2),leng0,leng_filt0,
     ,                    ab1,filter,tmp)
              ibase=(n-1)*ic+1
              itmp=iys1+(ibase-1)*niedge(1)
              y(itmp)=ab1(1)+ab1(2)
             end do
            end do
           end do
          endif

c         and then followed by X-synthesizing...
          if(j.le.nlevel(1))then
           i2=niedge(1)/ic
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do iys=1,niedge(2)
             iy0=(iys-1)*niedge(1)
             ixs1=iz0+iy0
             do itmp=1,niedge(1)
              tmp(itmp)=y(ixs1+itmp)
             enddo
             do n=1,i2
              call reverse(n,i1,ic,niedge(1),leng0,leng_filt0,
     ,                    ab1,filter,tmp)
              ibase=(n-1)*ic+1
              itmp=ixs1+ibase
              y(itmp)=ab1(1)+ab1(2)
             end do
            end do
           end do
          endif

         end do
 
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do
 
        endif
 
       end do
 
       return
       end

c___________________________________________________________________________
       subroutine face_lift_Pt3(image,isign,nlevel_syn)
       implicit none
c
c      pixel based Spline, tensor products
c
 
       integer mlevelx, mlevely, mlevelz, mlevelxyz
       integer mivolu, mbase, mvolu, nbase
       parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,           mlevelxyz=10,
     ,           mivolu=2**(mlevelx-1)*2**(mlevely-1)*2**(mlevelz-1),
     ,           mbase=1,mvolu=mivolu*mbase,nbase=mbase)
       integer nlevel(3),isign,nlevel_syn(3),niedge(3),
     ,         nivolu,iivolu0,iivolu,ivolu,nvolu,niface,
     ,         i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_syntmp(3),i3,ii0
       double precision image(mvolu)
       double precision y(mivolu),tmp(2**(mlevelxyz-1))
       common/basics/nlevel,nivolu,nvolu,niface,niedge

c      bior_1.1
c       parameter(leng_filt=2)
c       integer iab(leng_filt),iab1(leng_filt),iab2(leng_filt)
c      double precision hi_D1(leng_filt),lo_D1(leng_filt),ab(leng_filt),
c     ,       hi_R1(leng_filt),lo_R1(leng_filt),filter(leng_filt,2),
c     ,       hi_D2(leng_filt),lo_D2(leng_filt),ab1(leng_filt),
c     ,       hi_R2(leng_filt),lo_R2(leng_filt),ab2(leng_filt)
c       data lo_D1/  0.5d0,  0.5d0/
c       data hi_D1/ -0.5d0,  0.5d0/
c       data lo_R1/  1.0d0,  1.0d0/
c       data hi_R1/  1.0d0, -1.0d0/
c       data lo_D2/  1.0d0,  1.0d0/
c       data hi_D2/ -1.0d0,  1.0d0/
c       data lo_R2/  0.5d0,  0.5d0/
c       data hi_R2/  0.5d0, -0.5d0/

c      bior_3.1
       integer leng_filt
       parameter(leng_filt=4)
       integer iab(leng_filt),iab1(leng_filt),iab2(leng_filt)
      double precision hi_D1(leng_filt),lo_D1(leng_filt),ab(leng_filt),
     ,       hi_R1(leng_filt),lo_R1(leng_filt),filter(leng_filt,2),
     ,       hi_D2(leng_filt),lo_D2(leng_filt),ab1(leng_filt),
     ,       hi_R2(leng_filt),lo_R2(leng_filt),ab2(leng_filt)
       data lo_D1/  -0.25d0,   0.75d0,   0.75d0,  -0.25d0/
       data hi_D1/ -0.125d0,  0.375d0, -0.375d0,  0.125d0/
       data lo_R1/   0.25d0,   0.75d0,   0.75d0,   0.25d0/
       data hi_R1/   -0.5d0,   -1.5d0,    1.5d0,    0.5d0/
       data lo_D2/   0.25d0,   0.75d0,   0.75d0,   0.25d0/
       data hi_D2/    0.5d0,    1.5d0 ,  -1.5d0,   -0.5d0/
       data lo_R2/  -0.25d0,   0.75d0,   0.75d0,  -0.25d0/
       data hi_R2/  0.125d0, -0.375d0,  0.375d0, -0.125d0/

c local variable
       integer i0, k, n
       integer ibase, itmp
       integer ix0, ixs, ixs1
       integer iy0, iys, iys1
       integer iz0, izs, izs1
       integer leng0, leng_filt0
       integer niedge12
       integer nlevelT

       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       elseif(isign.gt.0)then
        do i=1,leng_filt
         if(isign.eq.1)then
          filter(i,1)=lo_D1(i)
          filter(i,2)=hi_D1(i)
         else
          filter(i,1)=lo_D2(i)
          filter(i,2)=hi_D2(i)
         endif
        enddo
       elseif(isign.lt.0)then
        do i=1,leng_filt
         if(isign.eq.-1)then
          filter(i,1)=lo_R1(i)
          filter(i,2)=hi_R1(i)
         else
          filter(i,1)=lo_R2(i)
          filter(i,2)=hi_R2(i)
         endif
        enddo
       end if
       niedge12=niedge(1)*niedge(2)
       leng0=leng_filt/2
       if(leng0*2.ne.leng_filt)leng0=leng0+1
       leng_filt0=leng_filt
       nlevelT=max0(nlevel(1),nlevel(2),nlevel(3))
 
       do i0=1,nbase
        iivolu0=(i0-1)*nivolu
 

        if(isign.gt.0)then
 
c        Decomposition or Analyzing....
c        S(j,k)=SUM_n[lo_D(2k+leng0-n)*S(j-1,n)]
c        d(j,k)=SUM_n[hi_D(2k+leng0-n)*S(j-1,n)]

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          y(iivolu)=image(ivolu)
         end do
 
         do j=2,nlevelT
          i1=2**(j-1)
          ic=i1/2

c         Starting from X-decomposition.... 
          if(j.le.nlevel(1))then
           i2=niedge(1)/i1
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do iys=1,niedge(2)
             iy0=(iys-1)*niedge(1)
             ixs1=iz0+iy0
             do itmp=1,niedge(1)
              tmp(itmp)=y(ixs1+itmp)
             end do
             do k=1,i2
              call forward(k,ic,niedge(1),leng0,leng_filt0,
     ,                     ab1,filter,tmp)
              ibase=(k-1)*i1+1
              itmp=ixs1+ibase
              y(itmp)=ab1(1)
              y(itmp+ic)=ab1(2)
             end do
            end do
           end do
          end if

c         Followed by Y-decomposition....
          if(j.le.nlevel(2))then
           i2=niedge(2)/i1
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do ixs=1,niedge(1)
             ix0=ixs
             iys1=iz0+ix0
             do itmp=1,niedge(2)
              tmp(itmp)=y(iys1+(itmp-1)*niedge(1))
             end do
             do k=1,i2
              call forward(k,ic,niedge(2),leng0,leng_filt0,
     ,                     ab1,filter,tmp)
              ibase=(k-1)*i1+1
              itmp=iys1+(ibase-1)*niedge(1)
              y(itmp)=ab1(1)
              y(itmp+ic*niedge(1))=ab1(2)
             end do
            end do
           end do
          end if 

c         And followed by Z-decomposition....
          if(j.le.nlevel(3))then
           i2=niedge(3)/i1
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            do ixs=1,niedge(1)
             ix0=ixs
             izs1=iy0+ix0
             do itmp=1,niedge(3)
              tmp(itmp)=y(izs1+(itmp-1)*niedge12)
             end do
             do k=1,i2
              call forward(k,ic,niedge(3),leng0,leng_filt0,
     ,                     ab1,filter,tmp)
              ibase=(k-1)*i1+1
              itmp=izs1+(ibase-1)*niedge12
              y(itmp)=ab1(1)
              y(itmp+ic*niedge12)=ab1(2)
             end do
            end do
           end do
          end if

         end do
 
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do
 
        elseif(isign.lt.0)then
 
c        Synthesizing.....
c        S(j-1,n)=SUM_k[lo_R(n-2k+leng0+1)*S(j,k)+
c                       hi_R(n-2k+leng0+1)*d(j,k)]

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          y(iivolu)=image(ivolu)
         end do

c        first take out all contribution above the specfiled synthesizing level
         nlevel_syntmp(3)=nlevel_syn(3)-1
         i1=nlevel(3)-nlevel_syntmp(3)
         i2=2**(i1-1)
         i3=niedge(3)/i2
         if(nlevel_syntmp(3).lt.0)then
          write(*,*)' Strange Z-synthesizing level !',nlevel_syn(3)
          stop
         else
          do iys=1,niedge(2)
           iy0=(iys-1)*niedge(1)
           do ixs=1,niedge(1)
            ix0=ixs
            izs1=iy0+ix0
            do i=1,i3
             ii0=(i-1)*i2
             do j=2,i2
              y(izs1+(ii0+j-1)*niedge12)=0.d0
             end do
            end do
           end do
          end do
         endif
         nlevel_syntmp(2)=nlevel_syn(2)-1
         i1=nlevel(2)-nlevel_syntmp(2)
         i2=2**(i1-1)
         i3=niedge(2)/i2
         if(nlevel_syntmp(2).lt.0)then
          write(*,*)' Strange Y-synthesizing level !',nlevel_syn(2)
          stop
         elseif(i2.ge.2)then
          do izs=1,niedge(3)
           iz0=(izs-1)*niedge12
           do ixs=1,niedge(1)
            ix0=ixs
            iys1=iz0+ix0
            do i=1,i3
             ii0=(i-1)*i2
             do j=2,i2
              y(iys1+(ii0+j-1)*niedge(1))=0.d0
             end do
            end do
           end do
          end do
         endif
         nlevel_syntmp(1)=nlevel_syn(1)-1
         i1=nlevel(1)-nlevel_syntmp(1)
         i2=2**(i1-1)
         i3=niedge(1)/i2
         if(nlevel_syntmp(1).lt.0)then
          write(*,*)' Strange X-synthesizing level !',nlevel_syn(1)
          stop
         elseif(i2.ge.2)then
          do izs=1,niedge(3)
           iz0=(izs-1)*niedge12
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            ixs1=iz0+iy0
            do i=1,i3
             ii0=(i-1)*i2
             do j=2,i2
              y(ixs1+ii0+j)=0.d0
             end do
            end do
           end do
          end do
         endif

         do j=nlevelT,2,-1
          i1=2**(j-1)
          ic=i1/2

c         The order of synthesizing reverses the order of analyzing
c         Starting from Z_synthesizing.....
          if(j.le.nlevel(3))then
           i2=niedge(3)/ic
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            do ixs=1,niedge(1)
             ix0=ixs
             izs1=iy0+ix0
             do itmp=1,niedge(3)
              tmp(itmp)=y(izs1+(itmp-1)*niedge12)
             enddo 
             do n=1,i2
              call reverse(n,i1,ic,niedge(3),leng0,leng_filt0,
     ,                    ab1,filter,tmp)
              ibase=(n-1)*ic+1
              itmp=izs1+(ibase-1)*niedge12
              y(itmp)=ab1(1)+ab1(2)
             end do
            end do
           end do
          endif

c         Followed by Y-synthesizing...
          if(j.le.nlevel(2))then
           i2=niedge(2)/ic
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do ixs=1,niedge(1)
             ix0=ixs
             iys1=iz0+ix0
             do itmp=1,niedge(2)
              tmp(itmp)=y(iys1+(itmp-1)*niedge(1))
             enddo
             do n=1,i2
              call reverse(n,i1,ic,niedge(2),leng0,leng_filt0,
     ,                    ab1,filter,tmp)
              ibase=(n-1)*ic+1
              itmp=iys1+(ibase-1)*niedge(1)
              y(itmp)=ab1(1)+ab1(2)
             end do
            end do
           end do
          endif

c         and then followed by X-synthesizing...
          if(j.le.nlevel(1))then
           i2=niedge(1)/ic
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do iys=1,niedge(2)
             iy0=(iys-1)*niedge(1)
             ixs1=iz0+iy0
             do itmp=1,niedge(1)
              tmp(itmp)=y(ixs1+itmp)
             enddo
             do n=1,i2
              call reverse(n,i1,ic,niedge(1),leng0,leng_filt0,
     ,                    ab1,filter,tmp)
              ibase=(n-1)*ic+1
              itmp=ixs1+ibase
              y(itmp)=ab1(1)+ab1(2)
             end do
            end do
           end do
          endif

         end do
 
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do
 
        endif
 
       end do
 
       return
       end



       subroutine forward(k,ic,niedge,leng0,leng_filt,
     ,                    ab1,filter,tmp)
       implicit none

       integer n, k, ic, niedge, leng0, leng_filt
       integer mlevelxyz
       integer ifb, ifilt, ift1, itmp
       parameter(mlevelxyz=10)
       double precision ab1(2),filter(leng_filt,2),
     ,                  tmp(2**(mlevelxyz-1))

       ifb=2*k+leng0-leng_filt
       ab1(1)=0.d0
       ab1(2)=0.d0
       do ifilt=1,leng_filt
        n=ifb+ifilt-1
        ift1=(n-1)*ic+1
        if(ift1.ge.1.and.ift1.le.niedge)then
         itmp=2*k+leng0-n
         ab1(1)=ab1(1)+filter(itmp,1)*tmp(ift1)
         ab1(2)=ab1(2)+filter(itmp,2)*tmp(ift1)
        end if
       end do
       return
       end



       subroutine reverse(n,i1,ic,niedge,leng0,leng_filt,
     ,                    ab1,filter,tmp)
       implicit none
 
       integer k, n, i1, ic, niedge, leng0, leng_filt
       integer mlevelxyz
       integer ifilt, ift1, ift2, itmp, kb
       parameter(mlevelxyz=10)
       double precision ab1(2),filter(leng_filt,2),
     ,                  tmp(2**(mlevelxyz-1))

       kb=(n+leng0-leng_filt)/2
       ab1(1)=0.d0
       ab1(2)=0.d0
       do ifilt=1,leng_filt/2+1
        k=kb+ifilt
        ift1=(k-1)*i1+1
        ift2=ift1+ic
        if(ift1.ge.1.and.ift1.le.niedge)then
         itmp=n-2*k+leng0+1
         if(itmp.ge.1.and.itmp.le.leng_filt)then
          ab1(1)=ab1(1)+filter(itmp,1)*tmp(ift1)
          ab1(2)=ab1(2)+filter(itmp,2)*tmp(ift2)
         end if
        end if
       end do
       return
       end



c________________________________________________________________________
       subroutine face_lift_Pb1(image,isign,nlevel_syn0)
       implicit none
c
c      pixel based Harr, 3D cube (not tensor products),need to run
c        transplant first
c
       integer mlevel, mivolu, mbase, mvolu, nbase
       parameter(mlevel=8,mivolu=8**(mlevel-1),
     ,           mbase=1,mvolu=mivolu*mbase,nbase=mbase)
       integer face_conf(mivolu)
       integer nlevel0(3),isign,nlevel_syn0(3),niedge0(3),
     ,         nivolu,iivolu0,iivolu,ivolu,nvolu,niface,
     ,         i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_syntmp,i3,ii0
       real*8 image(mvolu)
       real*8 y(mivolu)
       integer ll(8)
       real*8 yy(8),ff(8)
       real*4 T1(8,8),T2(8,8)
       common/lift/face_conf
       common/basics/nlevel0,nivolu,nvolu,niface,niedge0
 
       data T1/
c     /  0.125,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,
c     ,  0.125, 0.4375,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,
c     ,  0.125,-0.0625, 0.4375,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,
c     ,  0.125,-0.0625,-0.0625, 0.4375,-0.0625,-0.0625,-0.0625,-0.0625,
c     ,  0.125,-0.0625,-0.0625,-0.0625, 0.4375,-0.0625,-0.0625,-0.0625,
c     ,  0.125,-0.0625,-0.0625,-0.0625,-0.0625, 0.4375,-0.0625,-0.0625,
c     ,  0.125,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625, 0.4375,-0.0625,
c     ,  0.125,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625,-0.0625, 0.4375/
     /  0.125,-0.125,-0.125, 0.125,-0.125, 0.125, 0.125,-0.125,
     ,  0.125, 0.125,-0.125,-0.125,-0.125,-0.125, 0.125, 0.125,
     ,  0.125,-0.125, 0.125,-0.125,-0.125, 0.125,-0.125, 0.125,
     ,  0.125, 0.125, 0.125, 0.125,-0.125,-0.125,-0.125,-0.125,
     ,  0.125,-0.125,-0.125, 0.125, 0.125,-0.125,-0.125, 0.125,
     ,  0.125, 0.125,-0.125,-0.125, 0.125, 0.125,-0.125,-0.125,
     ,  0.125,-0.125, 0.125,-0.125, 0.125,-0.125, 0.125,-0.125,
     ,  0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125/
       data T2/
c     /           1.,-2.,-2.,-2.,-2.,-2.,-2.,-2.,
c     ,           1., 2., 0., 0., 0., 0., 0., 0.,
c     ,           1., 0., 2., 0., 0., 0., 0., 0.,
c     ,           1., 0., 0., 2., 0., 0., 0., 0.,
c     ,           1., 0., 0., 0., 2., 0., 0., 0.,
c     ,           1., 0., 0., 0., 0., 2., 0., 0.,
c     ,           1., 0., 0., 0., 0., 0., 2., 0.,
c     ,           1., 0., 0., 0., 0., 0., 0., 2./
     /  1.0  ,-1.0  ,-1.0  , 1.0  ,-1.0  , 1.0  , 1.0  ,-1.0  ,
     ,  1.0  , 1.0  ,-1.0  ,-1.0  ,-1.0  ,-1.0  , 1.0  , 1.0  ,
     ,  1.0  ,-1.0  , 1.0  ,-1.0  ,-1.0  , 1.0  ,-1.0  , 1.0  ,
     ,  1.0  , 1.0  , 1.0  , 1.0  ,-1.0  ,-1.0  ,-1.0  ,-1.0  ,
     ,  1.0  ,-1.0  ,-1.0  , 1.0  , 1.0  ,-1.0  ,-1.0  , 1.0  ,
     ,  1.0  , 1.0  ,-1.0  ,-1.0  , 1.0  , 1.0  ,-1.0  ,-1.0  ,
     ,  1.0  ,-1.0  , 1.0  ,-1.0  , 1.0  ,-1.0  , 1.0  ,-1.0  ,
     ,  1.0  , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  , 1.0  /
 
c local variables
       integer i0, niedge, nlevel, nlevel_syn
       
c   notice that transpose(T2)=inverse(T1)
 
       nlevel=nlevel0(1)
       nlevel_syn=nlevel_syn0(1)
       niedge=niedge0(1)
c       write(*,'(8e12.5)')(image(j),j=15361,16385)
c       write(*,*)nlevel,nlevel_syn,niedge,nivolu,niface


       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       endif
 
       do i0=1,nbase
        iivolu0=(i0-1)*nivolu
 
        if(isign.gt.0)then
c        Decomposition or Analyzing....
         do iivolu=1,nivolu
          ivolu=iivolu0+face_conf(iivolu)
          y(iivolu)=image(ivolu)
         end do
 
         do i=2,nlevel
          i1=8**(i-1)
          ic=i1/8
          i2=nivolu/i1
          do l=1,i2
           ll(1)=(l-1)*i1+1
           ff(1)=y(ll(1))
           do j=2,8
            ll(j)=ll(1)+ic*(j-1)
            ff(j)=y(ll(j))
           enddo
           if(isign.eq.1)then
c           operation T1
            do ltmp=1,8
             yy(ltmp)=0.d0
             do ktmp=1,8
              yy(ltmp)=yy(ltmp)+dble(T1(ltmp,ktmp))*ff(ktmp)
             end do
            end do
           elseif(isign.eq.2)then
c           operation T2
            do ltmp=1,8
             yy(ltmp)=0.d0
             do ktmp=1,8
              yy(ltmp)=yy(ltmp)+dble(T2(ltmp,ktmp))*ff(ktmp)
             end do
            end do
           end if
           do ltmp=1,8
            y(ll(ltmp))=yy(ltmp)
           enddo
          end do
         end do
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do
 
        elseif(isign.lt.0)then
 
c        Synthesizing.....
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          y(iivolu)=image(ivolu)
         end do
c         do iqz=1,niedge
c          write(*,*)iqz
c          do iqy=1,niedge
c           write(*,'(8f10.4)')
c     ,      (image(face_conf((iqz-1)*niface+(iqy-1)*niedge+iqx)),
c     ,      iqx=1,niedge)
c          enddo
c         enddo
 
         nlevel_syntmp=nlevel_syn-1
         if(nlevel_syntmp.lt.0)then
          do iivolu=1,nivolu
           y(iivolu)=0.d0
          end do
         else
          i1=nlevel-nlevel_syntmp
          i2=8**(i1-1)
          i3=nivolu/i2
          do i=1,i3
           ii0=(i-1)*i2
           do j=2,i2
            y(ii0+j)=0.d0
           end do
          end do
          do i=nlevel,2,-1
           i1=8**(i-1)
           ic=i1/8
           i2=nivolu/i1
           do l=1,i2
            ll(1)=(l-1)*i1+1
            ff(1)=y(ll(1))
            do j=2,8
             ll(j)=ll(1)+ic*(j-1)
             ff(j)=y(ll(j))
            enddo
            if(isign.eq.-1)then
c            operation inv(T1) [or rather, transpose (T2)]
             do ltmp=1,8
              yy(ltmp)=0.d0
              do ktmp=1,8
               yy(ltmp)=yy(ltmp)+dble(T2(ktmp,ltmp))*ff(ktmp)
              end do
             end do
            elseif(isign.eq.-2)then
c            operation inv(T2) [or rather, transpose (T1)]
             do ltmp=1,8
              yy(ltmp)=0.d0
              do ktmp=1,8
               yy(ltmp)=yy(ltmp)+dble(T1(ktmp,ltmp))*ff(ktmp)
              end do
             end do
            end if
            do ltmp=1,8
             y(ll(ltmp))=yy(ltmp)
            enddo
           end do
          end do
         end if
 
         do iivolu=1,nivolu
          ivolu=iivolu0+face_conf(iivolu)
          image(ivolu)=y(iivolu)
         end do
 
        endif
 
       end do
 
c       write(*,'(8e12.5)')(image(j),j=15361,16385)
c       write(*,*)nlevel,nlevel_syn,niedge,nivolu,nifac
       return
       end



c________________________________________________________________
       subroutine transplant
       implicit none
       integer mlevel, mivolu
       parameter(mlevel=8,mivolu=8**(mlevel-1))
       integer nlevel0(3),nivolu,niface,niedge0(3),nvolu
       integer il,ic,ic1,iv,ix,ix1,ix2,iy,iy1,iy2,iz,iz1,iz2
       integer ii(8),face_conf(mivolu),face_tmp(mivolu)
c local variables
       integer i, j, i8, iiedge, iiface, imv, in
       integer ncc, niedge, nlevel

       common/lift/face_conf
       common/basics/nlevel0,nivolu,nvolu,niface,niedge0
 
       nlevel=nlevel0(1)
       niedge=niedge0(1)

       do i=1,nivolu
        face_conf(i)=i
       end do
       ic=2
       ic1=1
 
       if(nlevel.gt.2)then
        do i=nlevel,3,-1
         iiface=4**(i-1)
         iiedge=2**(i-1)
         il=nlevel-i+1
         iv=8**(il-1)
         in=niedge/(2**il)
         ncc=0
 
         do iz=1,in
          iz1=(iz-1)*ic+1
          iz2=iz1+ic1
          do iy=1,in
           iy1=(iy-1)*ic+1
           iy2=iy1+ic1
           do ix=1,in
            ix1=(ix-1)*ic+1
            ix2=ix1+ic1
            ii(1)=(iz1-1)*iiface+(iy1-1)*iiedge+ix1
            ii(2)=(iz1-1)*iiface+(iy1-1)*iiedge+ix2
            ii(3)=(iz1-1)*iiface+(iy2-1)*iiedge+ix1
            ii(4)=(iz1-1)*iiface+(iy2-1)*iiedge+ix2
            ii(5)=(iz2-1)*iiface+(iy1-1)*iiedge+ix1
            ii(6)=(iz2-1)*iiface+(iy1-1)*iiedge+ix2
            ii(7)=(iz2-1)*iiface+(iy2-1)*iiedge+ix1
            ii(8)=(iz2-1)*iiface+(iy2-1)*iiedge+ix2
            do i8=1,8
             do imv=1,iv
              ncc=ncc+1
              face_tmp(ncc)=face_conf((ii(i8)-1)*iv+imv)
             end do
            end do
           end do
          end do
         end do
 
         do j=1,nivolu
          face_conf(j)=face_tmp(j)
c          write(*,*)j,face_conf(j)
         end do
 
        end do
       endif
 
       return
       end


c___________________________________________________________________________
       subroutine face_lift_Nt1(image,isign,nlevel_syn)
       implicit none
c
c      node based, lifted linear, tensor products (in fact, rotate X,Y,Z
c                                  succesively transforming)
c       notice that nx=2**(lx-1)+1 not 2**(lx-1) in all the pixel based scheme
c
       integer mlevelx, mlevely, mlevelz, mlevelxyz, mx, my, mz
       integer mbase, mvolu, nbase, mivolu
       integer i0, ieffect, itmp
       integer ix0, ixs, ixs1
       integer iy0, iys, iys1
       integer iz0, izs, izs1
       integer niedge12, nlevelT

       parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,           mlevelxyz=10,
     ,           mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,           mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,           mbase=1,mvolu=mivolu*mbase,nbase=mbase)
       integer nlevel(3),isign,nlevel_syn(3),niedge(3),
     ,         nivolu,iivolu0,iivolu,ivolu,nvolu,niface,
     ,         i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_syntmp(3),i3,ii0
       double precision image(mvolu)
       double precision y(mivolu),tmp(mivolu)


       common/basics/nlevel,nivolu,nvolu,niface,niedge

       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       end if
       niedge12=niedge(1)*niedge(2)
       nlevelT=max0(nlevel(1),nlevel(2),nlevel(3))

       do i0=1,nbase
        iivolu0=(i0-1)*nivolu

        if(isign.gt.0)then

c        Decomposition or Analyzing....

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          y(iivolu)=image(ivolu)
         end do

         do j=2,nlevelT
          i1=2**(j-1)
          ic=i1/2

c         Starting from X-decomposition....
          if(j.le.nlevel(1))then
           i2=(niedge(1)-1)/i1
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do iys=1,niedge(2)
             iy0=(iys-1)*niedge(1)
             ixs1=iz0+iy0
             do itmp=1,niedge(1)
              tmp(itmp)=y(ixs1+itmp)
             end do
             call forward_node(i1,i2,ic,tmp,isign,1)
             do itmp=1,niedge(1)
              y(ixs1+itmp)=tmp(itmp)
             end do
            end do
           end do
          end if

c         Followed by Y-decomposition....
          if(j.le.nlevel(2))then
           i2=(niedge(2)-1)/i1
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do ixs=1,niedge(1)
             ix0=ixs
             iys1=iz0+ix0
             do itmp=1,niedge(2)
              tmp(itmp)=y(iys1+(itmp-1)*niedge(1))
             end do
             call forward_node(i1,i2,ic,tmp,isign,1)
             do itmp=1,niedge(2)
              y(iys1+(itmp-1)*niedge(1))=tmp(itmp)
             end do
            end do
           end do
          end if

c         And followed by Z-decomposition....
          if(j.le.nlevel(3))then
           i2=(niedge(3)-1)/i1
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            do ixs=1,niedge(1)
             ix0=ixs
             izs1=iy0+ix0
             do itmp=1,niedge(3)
              tmp(itmp)=y(izs1+(itmp-1)*niedge12)
             end do
             call forward_node(i1,i2,ic,tmp,isign,1)
             do itmp=1,niedge(3)
              y(izs1+(itmp-1)*niedge12)=tmp(itmp)
             end do
            end do
           end do
          end if
         end do

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do

        elseif(isign.lt.0)then

c        Synthesizing.....

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          y(iivolu)=image(ivolu)
         end do

         do j=nlevelT,2,-1
          i1=2**(j-1)
          ic=i1/2

c         The order of synthesizing reverses the order of analyzing
c         Starting from Z_synthesizing.....
          if(j.le.nlevel(3))then
           i2=(niedge(3)-1)/i1
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            do ixs=1,niedge(1)
             ix0=ixs
             izs1=iy0+ix0
             do itmp=1,niedge(3)
              tmp(itmp)=y(izs1+(itmp-1)*niedge12)
             enddo
             if(j.le.nlevel_syn(3))then
              ieffect=1
             else
              ieffect=0
             endif
             call forward_node(i1,i2,ic,tmp,isign,ieffect)
             do itmp=1,niedge(3)
              y(izs1+(itmp-1)*niedge12)=tmp(itmp)
             enddo
            end do
           end do
          endif

c         Followed by Y-synthesizing...
          if(j.le.nlevel(2))then
           i2=(niedge(2)-1)/i1
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do ixs=1,niedge(1)
             ix0=ixs
             iys1=iz0+ix0
             do itmp=1,niedge(2)
              tmp(itmp)=y(iys1+(itmp-1)*niedge(1))
             enddo
             if(j.le.nlevel_syn(2))then
              ieffect=1
             else
              ieffect=0
             endif
             call forward_node(i1,i2,ic,tmp,isign,ieffect)
             do itmp=1,niedge(2)
              y(iys1+(itmp-1)*niedge(1))=tmp(itmp)
             end do
            end do
           end do
          endif

c         and then followed by X-synthesizing...
          if(j.le.nlevel(1))then
           i2=(niedge(1)-1)/i1
           do izs=1,niedge(3)
            iz0=(izs-1)*niedge12
            do iys=1,niedge(2)
             iy0=(iys-1)*niedge(1)
             ixs1=iz0+iy0
             do itmp=1,niedge(1)
              tmp(itmp)=y(ixs1+itmp)
             enddo
             if(j.le.nlevel_syn(1))then
              ieffect=1
             else
              ieffect=0
             endif
             call forward_node(i1,i2,ic,tmp,isign,ieffect)
             do itmp=1,niedge(1)
              y(ixs1+itmp)=tmp(itmp)
             end do
            end do
           end do
          endif

         end do

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do

        endif

       end do

       return
       end



       subroutine forward_node(i1,i2,ic,tmp,isign,ieffect)
       implicit none

       integer i1, i2, ic, isign,ieffect
       integer mlevelx, mlevely, mlevelz, mx, my, mz
       integer mivolu, mbase, mvolu, nbase

       parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,           mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,           mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,           mbase=1,mvolu=mivolu*mbase,nbase=mbase)
       double precision tmp(mivolu)

       integer k, n1, n2, n3

       do k=1,i2
        n1=(k-1)*i1+1
        n2=n1+ic
        n3=n2+ic
c       1st phase of lifting for DECOMPOSING
        if(isign.eq.1)then
         tmp(n2)=tmp(n2)-0.5d0*tmp(n1)-0.5d0*tmp(n3)
        elseif(isign.eq.2)then
         tmp(n1)=tmp(n1)+0.5d0*tmp(n2)
         tmp(n3)=tmp(n3)+0.5d0*tmp(n2)
c       undone 2nd phase of lifting first during SYNTHESIZING
        elseif(isign.eq.-1)then
         if(ieffect.eq.1)then
          tmp(n1)=tmp(n1)-0.25d0*tmp(n2)
          tmp(n3)=tmp(n3)-0.25d0*tmp(n2)
         else
          tmp(n2)=0.d0
         endif
        elseif(isign.eq.-2)then
         if(ieffect.eq.1)then
          tmp(n2)=tmp(n2)+0.25d0*tmp(n1)+0.25d0*tmp(n3)
         else
          tmp(n2)=   0.d0+0.25d0*tmp(n1)+0.25d0*tmp(n3)
         endif
        endif
       end do

       do k=1,i2
        n1=(k-1)*i1+1
        n2=n1+ic
        n3=n2+ic
c       2nd lifting phase for DECOMPOSING
        if(isign.eq.1)then
         tmp(n1)=tmp(n1)+0.25d0*tmp(n2)
         tmp(n3)=tmp(n3)+0.25d0*tmp(n2)
        elseif(isign.eq.2)then
         tmp(n2)=tmp(n2)-0.25d0*tmp(n1)-0.25d0*tmp(n3)
c       undone 1st lifting phase after 2nd effects are removed
        elseif(isign.eq.-1)then
         if(ieffect.eq.1)then
          tmp(n2)=tmp(n2)+0.5d0*tmp(n1)+0.5d0*tmp(n3)
         else
          tmp(n2)=   0.d0+0.5d0*tmp(n1)+0.5d0*tmp(n3)
         endif
        elseif(isign.eq.-2)then
         if(ieffect.eq.1)then
          tmp(n1)=tmp(n1)-0.5d0*tmp(n2)
          tmp(n3)=tmp(n3)-0.5d0*tmp(n2)
         else
          tmp(n2)=0.d0
         endif
        endif
       end do

       return
       end


c___________________________________________________________________________
       subroutine face_lift_Hybrid_Nbt_spherical(image,isign,nlevel_syn)
       implicit none
c***************************************************************************
c   isign:  1,  IMAGE analysis (decomposition)                             *
c          -1,  IMAGE synthesis (reconstruction)                           *
c           0,  doing nothing                                              *
c           2,  PATH-integral analysis                                     *
c          -2,  PATH-integral synthesis                                    *
c                                                                          *
c   notice that in the present tomographic study, each row of              *
c    the G-matrix (each path scoring within highest level faces)           *
c    is first decomposed (2); the resulted model image after               *
c    inversion is then synthesized (-1)                                    * 
c                                                                          *
c    compare with the pyrimal form in spherical wavelet (or other          *
c       plannar triangulation grid)                                        *
c                                                                          *
c   New version of decomposition and synthesizing is performed successively*
c    in the order of (((1,nx),1,ny),1,nz).  In other words, preprocessing  *
c    by transplant is dropped, the reorganized index FACE_CONF is no more  *
c    needed.  Also, it is now possible to have different refinement level  *
c    in X-,Y-,and Z- direction  <2000/11/10>.                              *
c                                                                          *
c   New version of LINEAR VERTEX-BASED INTERPOLATORY wavelets              *
c       <2002/1/24>                                                        *
c       !notice that for reconstruction (synthesizing),                    *
c       !the second (lifting) phase has to be undone first                 *
c       !(the 2nd phase is the lower procedure in Fig.2 of                 *
c       !Sweldens' 5 minute tour)                                          *
c       ! i.e.,                                                            *
c       !the order of operation in Fig.2 of Sweldens's 5                   *
c       !minute tour is switched for DECOMPOSITION and                     *
c       !SYNTHESIZING                                                      *
c***************************************************************************
 
       integer mbase, mlevel, miface, mface, medge_head, mivert, mvert
       integer mlevelV, mvolu 
       parameter(mbase=20,mlevel=5,miface=4**(mlevel-1),
     ,           mface=miface*mbase,
     ,           medge_head=2**(mlevel-1)+1,
     ,           mivert=((medge_head+1)*medge_head)/2,
     ,           mvert=20*mivert)
       parameter(mlevelV=6,
     ,           mvolu=(2**(mlevelV-1)+1)*mvert)
       integer nbase,nlevel(2),isign,nlevel_syn(2),niedge(2),
     ,         nivolu,iivolu0,iivolu,ivolu,nvolu,niface,nlevel1,
     ,         i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_syntmp(2),i3,ii0
       integer nfase,nlevelS,nedge_head,nvert_global,nivert,nvert,
     ,         icc,icq(mlevel-1,2),face(3,miface)
       integer face_global(3,mface),indx(mvert),face_conf(miface),
     ,         face_index(mlevel-1),node_conf(miface*(mlevel-1),6)
       real*8 image(mvolu)
       real*8 y(mvolu),tmpV(2**(mlevelV-1)+1),tmpL(mvert)
       integer iab(2)
       common/basics_spherical/nlevel,nvolu,niedge
       common/base0/nbase,nlevelS,niface,nivert,nedge_head,
     ,              nface,nvert,nvert_global,icc,icq
       common/base1/face,face_global,indx,face_conf,
     ,              face_index,node_conf
 
        integer itmp
        integer ix
        integer iys
        integer izs, iz0
        integer lx
        integer n1, n2, n3, nface

c       write(*,*)nbase,nlevel(1),nlevel_syn(1),nlevel(2),nlevel_syn(2)
c       write(*,*)nvolu,niedge(1),niedge(2),isign
c       return

       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       endif
       nvolu=niedge(1)*niedge(2)

c       niedge(1)=nvert_global
c       niedge(2)=2**(nlevel(2)-1)+1
c      compare with the face-based bio-Harr where
c       niedge(1)=4**(nlevel(1)-1)*nbase
c       niedge(2)=2**(nlevel(2)-1)
 
       do iivolu=1,nvolu
        y(iivolu)=image(iivolu)
       end do

       if(isign.gt.0)then
c       Decomposition or Analyzing....

c       Starting from lateral spherical surface decomposition....
        do izs=1,niedge(2)
         iz0=(izs-1)*niedge(1)
         do itmp=1,niedge(1)
          tmpL(itmp)=y(iz0+itmp)
         enddo
         nlevel1=nlevel(1)
         call face_liftL(tmpL,isign,nlevel_syn(1))
         do itmp=1,niedge(1)
          y(iz0+itmp)=tmpL(itmp)
         enddo
        enddo

c       Followed by vertical radial decomposition....
        do iys=1,niedge(1)

         do itmp=1,niedge(2)
          tmpV(itmp)=y(iys+(itmp-1)*niedge(1))
         enddo
         do ix=2,nlevel(2)
          i1=2**(ix-1)
          ic=i1/2
          i2=(niedge(2)-1)/i1

c         first phase of lifting
          do lx=1,i2
           n1=(lx-1)*i1+1
           n2=n1+ic
           n3=n2+ic
           if(isign.eq.1)then
            tmpV(n2)=tmpV(n2)-0.5d0*tmpV(n1)-0.5d0*tmpV(n3)
           elseif(isign.eq.2)then
            tmpV(n1)=tmpV(n1)+0.5d0*tmpV(n2)
            tmpV(n3)=tmpV(n3)+0.5d0*tmpV(n2)
           endif
          enddo

c         2nd phase
          do lx=1,i2
           n1=(lx-1)*i1+1
           n2=n1+ic
           n3=n2+ic
           if(isign.eq.1)then
            tmpV(n1)=tmpV(n1)+0.25d0*tmpV(n2)
            tmpV(n3)=tmpV(n3)+0.25d0*tmpV(n2)
           elseif(isign.eq.2)then
            tmpV(n2)=tmpV(n2)-0.25d0*tmpV(n1)-0.25d0*tmpV(n3)
           endif
          enddo
         enddo
         do itmp=1,niedge(2)
          y(iys+(itmp-1)*niedge(1))=tmpV(itmp)
         enddo

        enddo

       elseif(isign.lt.0)then
c       Synthesizing.....

c       The order of synthesizing reverses the order of analyzing
c       Starting from vertical radial synthesizing.....
        nlevel_syntmp(2)=nlevel_syn(2)-1
        if(nlevel_syntmp(2).lt.0)then
         write(*,*)' Strange radial-synthesizing level !',nlevel_syn(2)
         stop
        else
         do iys=1,niedge(1)

          do itmp=1,niedge(2)
           tmpV(itmp)=y(iys+(itmp-1)*niedge(1))
          enddo
          i1=nlevel(2)-nlevel_syntmp(2)
          i2=2**(i1-1)
          i3=niedge(2)/i2
          do i=1,i3
           ii0=(i-1)*i2
           if(i2.ge.2)then
            do j=2,i2
             tmpV(ii0+j)=0.d0
            end do
           endif
          end do
          do i=nlevel(2),2,-1
           i1=2**(i-1)
           ic=i1/2
           i2=(niedge(2)-1)/i1

c          first phase
           do l=1,i2
            n1=(l-1)*i1+1
            n2=n1+ic
            n3=n2+ic
            if(isign.eq.-1)then
             tmpV(n1)=tmpV(n1)-0.25d0*tmpV(n2)
             tmpV(n3)=tmpV(n3)-0.25d0*tmpV(n2)
            elseif(isign.eq.-2)then
             tmpV(n2)=tmpV(n2)+0.25d0*tmpV(n1)+0.25d0*tmpV(n3)
            endif
           enddo

c          2nd phase
           do l=1,i2
            n1=(l-1)*i1+1
            n2=n1+ic
            n3=n2+ic
            if(isign.eq.-1)then
             tmpV(n2)=tmpV(n2)+0.5d0*tmpV(n1)+0.5d0*tmpV(n3)
            elseif(isign.eq.-2)then
             tmpV(n1)=tmpV(n1)-0.5d0*tmpV(n2)
             tmpV(n3)=tmpV(n3)-0.5d0*tmpV(n2)
            endif
           enddo
          enddo
          do itmp=1,niedge(2)
           y(iys+(itmp-1)*niedge(1))=tmpV(itmp)
          enddo
         end do
        endif

c       Followed by  lateral spherical surface synthesizing...
        nlevel_syntmp(1)=nlevel_syn(1)-1
        if(nlevel_syntmp(1).lt.0)then
         write(*,*)' Strange X-synthesizing level !',nlevel_syn(1)
         stop
        else
         do izs=1,niedge(2)
          iz0=(izs-1)*niedge(1)
          do itmp=1,niedge(1)
           tmpL(itmp)=y(iz0+itmp)
          enddo
          nlevel1=nlevel(1)
          call face_liftL(tmpL,isign,nlevel_syn(1))
          do itmp=1,niedge(1)
           y(iz0+itmp)=tmpL(itmp)
          enddo
         enddo
        endif

       endif

       do iivolu=1,nvolu
        image(iivolu)=y(iivolu)
       end do
 
       return
       end

c______________________________________________________________________
       subroutine face_liftL(image,isign,nlevel_syn)
       implicit none
c**********************************************************************
c    isign:  1,  IMAGE analysis  (decomposition)                      *
c           -1,  IMAGE synthesis (reconstruction)                     *
c            0,  doing nothing                                        *
c            2,  PATH-integral analysis                               *
c           -2,  PATH-integral synthesis                              *
c     <2002/1/17>                                                     *
c**********************************************************************

       integer isign, nlevel_syn

       integer mlevel, miface, mface, medge_head, mivert, mvert
       parameter(mlevel=5,miface=4**(mlevel-1),mface=20*miface,
     ,           medge_head=2**(mlevel-1)+1,
     ,           mivert=((medge_head+1)*medge_head)/2,
     ,           mvert=20*mivert)
       logical ifill(mvert,2)
       integer i,j,k,l,iface,nedge_head,nlevel,nvert_global,
     ,         niface,nivert,nface,nbase,nvert_base,icc,
     ,         indx(mvert),face_global(3,mface),neigh,npoin,
     ,         face_index(mlevel-1),face_conf(miface),nelem,
     ,         node_conf(miface*(mlevel-1),6),icq(mlevel-1,2)
       integer face0(3,20),face(3,miface),llocal(3,3)
       real*8 image(mvert),data(mface*10)
 
       integer ibase, ibase0, ifaceb, ifacef, ilevel
       integer itmp, l1,l2, l3, n1, n2, n3, number_face, nvert

       common/base0/nbase,nlevel,niface,nivert,nedge_head,
     ,              nface,nvert,nvert_global,icc,icq
       common/base1/face,face_global,indx,face_conf,
     ,              face_index,node_conf

       do i=1,nvert_global
        do j=1,2
         ifill(i,j)=.TRUE.
        end do
       end do

       do ilevel=1,nlevel-1
        if(isign.eq.1.or.isign.eq.2)then
         number_face=face_index(ilevel)
         ifaceb=icq(ilevel,1)
         ifacef=icq(ilevel,2)
        elseif(isign.eq.-1.or.isign.eq.-2)then
         number_face=face_index(nlevel-ilevel)
         ifaceb=icq(nlevel-ilevel,1)
         ifacef=icq(nlevel-ilevel,2)
        else
         return
        endif

c       !first phase for decomposition and synthesizing
c       !that invokes the unlifted HAT
        do ibase=1,nbase
         ibase0=(ibase-1)*nivert
         do iface=ifaceb,ifacef
          llocal(1,1)=node_conf(iface,1)
          llocal(1,2)=node_conf(iface,6)
          llocal(1,3)=node_conf(iface,2)
          llocal(2,1)=node_conf(iface,2)
          llocal(2,2)=node_conf(iface,4)
          llocal(2,3)=node_conf(iface,3)
          llocal(3,1)=node_conf(iface,3)
          llocal(3,2)=node_conf(iface,5)
          llocal(3,3)=node_conf(iface,1)
          do itmp=1,3
           l1=llocal(itmp,1)
           l2=llocal(itmp,2)
           l3=llocal(itmp,3)
           n1=indx(ibase0+l1)
           n2=indx(ibase0+l2)
           n3=indx(ibase0+l3)
           if(ifill(n2,1))then
            if(isign.eq.1)then
             image(n2)=image(n2)
     -                -0.5d0*image(n1)-0.5d0*image(n3)
            elseif(isign.eq.2)then
             image(n1)=image(n1)+0.5d0*image(n2)
             image(n3)=image(n3)+0.5d0*image(n2)
c            !notice that for reconstruction (synthesizing),
c            !the second (lifting) phase has to be undone first
c            !(the 2nd phase is the lower procedure in Fig.2 of
c            !Sweldens' 5 minute tour)
            elseif(isign.eq.-1)then
             if(ilevel.le.nlevel_syn-1)then
              image(n1)=image(n1)-0.25d0*image(n2)
              image(n3)=image(n3)-0.25d0*image(n2)
             else
              image(n2)=0.d0
             endif
            elseif(isign.eq.-2)then
             if(ilevel.le.nlevel_syn-1)then
              image(n2)=image(n2)
     +                 +0.25d0*image(n1)+0.25d0*image(n3)
             else
              image(n2)=0.d0
     +                 +0.25d0*image(n1)+0.25d0*image(n3)
             endif
            endif
            ifill(n2,1)=.FALSE.
           endif
          end do
         end do
        end do

c       !second phase of the LIFTING process, notice that
c       !the order of operation in Fig.2 of Sweldens's 5
c       !minute tour is switched for DECOMPOSITION and
c       !SYNTHESIZING
        do ibase=1,nbase
         ibase0=(ibase-1)*nivert
         do iface=ifaceb,ifacef
          llocal(1,1)=node_conf(iface,1)
          llocal(1,2)=node_conf(iface,6)
          llocal(1,3)=node_conf(iface,2)
          llocal(2,1)=node_conf(iface,2)
          llocal(2,2)=node_conf(iface,4)
          llocal(2,3)=node_conf(iface,3)
          llocal(3,1)=node_conf(iface,3)
          llocal(3,2)=node_conf(iface,5)
          llocal(3,3)=node_conf(iface,1)
          do itmp=1,3
           l1=llocal(itmp,1)
           l2=llocal(itmp,2)
           l3=llocal(itmp,3)
           n1=indx(ibase0+l1)
           n2=indx(ibase0+l2)
           n3=indx(ibase0+l3)
           if(ifill(n2,2))then
            if(isign.eq.1)then
             image(n1)=image(n1)+0.25d0*image(n2)
             image(n3)=image(n3)+0.25d0*image(n2)
            elseif(isign.eq.2)then
             image(n2)=image(n2)
     -                -0.25d0*image(n1)-0.25d0*image(n3)
            elseif(isign.eq.-1)then
             if(ilevel.le.nlevel_syn-1)then
              image(n2)=image(n2)
     +                 +0.5d0*image(n1)+0.5d0*image(n3)
             else
              image(n2)=0.d0
     +                 +0.5d0*image(n1)+0.5d0*image(n3)
             endif
            elseif(isign.eq.-2)then
             if(ilevel.le.nlevel_syn-1)then
              image(n1)=image(n1)-0.5d0*image(n2)
              image(n3)=image(n3)-0.5d0*image(n2)
             else
              image(n2)=0.d0
             endif
            endif
            ifill(n2,2)=.FALSE.
           endif
          end do
         end do
        end do

       end do

       return
       end


c___________________________________________________________________________
      subroutine face_lift_Hybrid_Nbt (image, isign, nlevel_syn)
       implicit none
c
c     Hybrid: laterally  Nb1 (triangular meshes) +
c             vertically Nt1 (All lifted Linear node based)
c
       integer mlevelx, mlevely, mlevelz, mx, my, mz, mivolu
       integer mbase, mvolu, nbase
       parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,           mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,           mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,           mbase=1,mvolu=mivolu*mbase,nbase=mbase)
       integer nlevel(3),isign,nlevel_syn(3),niedge(3),
     ,         nivolu,iivolu0,iivolu,ivolu,nvolu,niface,
     ,         nnx(3),nny(3),nnz(3),local(27),
     ,         nlevel2D(2),niedge2D(2),nlevel_syn2D(2),
     ,         i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_syntmp(3),i3,ii0

       integer i0, itmp
       integer ix, ixs, ixys, iy0, iys, iz0, izs
       integer lx, n1, n2, n3
       integer niedge12, nlevel_syntmpL, nlevelT

       double precision image(mvolu)
       double precision y(mivolu),tmpV(mz),tmpL(mx*my)

       common/basics/nlevel,nivolu,nvolu,niface,niedge
       common/basics_C_2D/nlevel2D,niedge2D

       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       end if
       niedge12=niedge(1)*niedge(2)
       nlevelT=min0(nlevel(1),nlevel(2))
       do i=1,2
        nlevel2D(i)=nlevel(i)
        niedge2D(i)=niedge(i)
        nlevel_syn2D(i)=nlevel_syn(i)
       enddo

       do i0=1,nbase
        iivolu0=(i0-1)*nivolu
        do iivolu=1,nivolu
         ivolu=iivolu0+iivolu
         y(iivolu)=image(ivolu)
        end do

        if(isign.gt.0)then
c        Decomposition or Analyzing....

c        Starting from lateral XY decomposition
         do izs=1,niedge(3)
          iz0=(izs-1)*niedge12
          do iys=1,niedge(2)
           iy0=(iys-1)*niedge(1)
           do ixs=1,niedge(1)
            tmpL(iy0+ixs)=y(iz0+iy0+ixs)
           enddo
          enddo
          call face_liftLC(tmpL,isign,nlevel_syn2D)
          do iys=1,niedge(2)
           iy0=(iys-1)*niedge(1)
           do ixs=1,niedge(1)
            y(iz0+iy0+ixs)=tmpL(iy0+ixs)
           enddo
          enddo
         enddo

c        Followed by vertical radial decomposition....
         do ixys=1,niedge12

          do itmp=1,niedge(3)
           tmpV(itmp)=y(ixys+(itmp-1)*niedge12)
          enddo
          do ix=2,nlevel(3)
           i1=2**(ix-1)
           ic=i1/2
           i2=(niedge(3)-1)/i1

c          first phase of lifting
           do lx=1,i2
            n1=(lx-1)*i1+1
            n2=n1+ic
            n3=n2+ic
            if(isign.eq.1)then
             tmpV(n2)=tmpV(n2)-0.5d0*tmpV(n1)-0.5d0*tmpV(n3)
            elseif(isign.eq.2)then
             tmpV(n1)=tmpV(n1)+0.5d0*tmpV(n2)
             tmpV(n3)=tmpV(n3)+0.5d0*tmpV(n2)
            endif
           enddo

c          2nd phase
           do lx=1,i2
            n1=(lx-1)*i1+1
            n2=n1+ic
            n3=n2+ic
            if(isign.eq.1)then
             tmpV(n1)=tmpV(n1)+0.25d0*tmpV(n2)
             tmpV(n3)=tmpV(n3)+0.25d0*tmpV(n2)
            elseif(isign.eq.2)then
             tmpV(n2)=tmpV(n2)-0.25d0*tmpV(n1)-0.25d0*tmpV(n3)
            endif
           enddo
          enddo
          do itmp=1,niedge(3)
           y(ixys+(itmp-1)*niedge12)=tmpV(itmp)
          enddo

         enddo

        elseif(isign.lt.0)then
c        Synthesizing.....

c        The order of synthesizing reverses the order of analyzing
c        Starting from vertical radial synthesizing.....
         nlevel_syntmp(3)=nlevel_syn(3)-1
         if(nlevel_syntmp(3).lt.0)then
          write(*,*)' Strange vertical synthesizing level !'
          write(*,*)nlevel_syn(3)
          stop
         else
          do ixys=1,niedge12

           do itmp=1,niedge(3)
            tmpV(itmp)=y(ixys+(itmp-1)*niedge12)
           enddo
           i1=nlevel(3)-nlevel_syntmp(3)
           i2=2**(i1-1)
           i3=niedge(3)/i2
           do i=1,i3
            ii0=(i-1)*i2
            if(i2.ge.2)then
             do j=2,i2
              tmpV(ii0+j)=0.d0
             end do
            endif
           end do
           do ix=nlevel(3),2,-1
            i1=2**(ix-1)
            ic=i1/2
            i2=(niedge(3)-1)/i1
 
c           first phase of lifting
            do lx=1,i2
             n1=(lx-1)*i1+1
             n2=n1+ic
             n3=n2+ic
             if(isign.eq.-1)then
              tmpV(n1)=tmpV(n1)-0.25d0*tmpV(n2)
              tmpV(n3)=tmpV(n3)-0.25d0*tmpV(n2)
             elseif(isign.eq.-2)then
              tmpV(n2)=tmpV(n2)+0.25d0*tmpV(n1)+0.25d0*tmpV(n3)
             endif
            enddo

c           2nd phase
            do lx=1,i2
             n1=(lx-1)*i1+1
             n2=n1+ic
             n3=n2+ic
             if(isign.eq.-1)then
              tmpV(n2)=tmpV(n2)+0.5d0*tmpV(n1)+0.5d0*tmpV(n3)
             elseif(isign.eq.-2)then
              tmpV(n1)=tmpV(n1)-0.5d0*tmpV(n2)
              tmpV(n3)=tmpV(n3)-0.5d0*tmpV(n2)
             endif
            enddo
           enddo
           do itmp=1,niedge(3)
            y(ixys+(itmp-1)*niedge12)=tmpV(itmp)
           enddo
          enddo
         endif

c        Followed by  lateral spherical surface synthesizing...
         nlevel_syntmpL=min0(nlevel_syn(1),nlevel_syn(2))
         if(nlevel_syntmpL.lt.1)then
          write(*,*)' Strange laterally synthesizing level !',nlevel_syntmpL
          stop
         else
          do izs=1,niedge(3)
           iz0=(izs-1)*niedge12
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            do ixs=1,niedge(1)
             tmpL(iy0+ixs)=y(iz0+iy0+ixs)
            enddo
           enddo
           call face_liftLC(tmpL,isign,nlevel_syn2D)
           do iys=1,niedge(2)
            iy0=(iys-1)*niedge(1)
            do ixs=1,niedge(1)
             y(iz0+iy0+ixs)=tmpL(iy0+ixs)
            enddo
           enddo
          enddo
         endif

        endif

        do iivolu=1,nivolu
         ivolu=iivolu0+iivolu
         image(ivolu)=y(iivolu)
        end do

       end do
       
       return
       end


c______________________________________________________________________
       subroutine face_liftLC(image,isign,nlevel_syn)
       implicit none
c
c       2D node-based transformation based on implicit triangulation
c        formed for ((1,nx),ny)
c
c               7-8-9
c            Y  | |/|
c            ^  4-5-6
c            |  |/| |
c            |  1-2-3
c            |
c            ---->X
c
       integer  mlevelx, mlevely, mx, my, miface, mbase, mface, nbase
       parameter(mlevelx=8,mlevely=8,
     ,           mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,           miface=mx*my,
     ,           mbase=1,mface=miface*mbase,nbase=mbase)
       integer nlevel(2),isign,nlevel_syn(2),niedge(2),
     ,         nnx(3),nny(3),local(9)
       integer i0, i1, i2x, i2y
       integer ic, ip, iface, iiface, iiface0
       integer ix, ix0
       integer iy, iy0
       integer j, kx, ky
       integer niface, nlevel_synT, nlevelT

       double precision image(mface)
       double precision tmp(miface)
       common/basics_C_2D/nlevel,niedge

       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       end if
       niface=niedge(1)*niedge(2)
       nlevelT=min0(nlevel(1),nlevel(2))
       nlevel_synT=min0(nlevel_syn(1),nlevel_syn(2))

       do i0=1,nbase
        iiface0=(i0-1)*niface
        do iiface=1,niface
         iface=iiface0+iiface
         tmp(iiface)=image(iface)
        enddo

        if(isign.gt.0)then
c        Decomposition or Analyzing....

         do j=2,nlevelT
          i1=2**(j-1)
          ic=i1/2
          i2x=(niedge(1)-1)/i1
          i2y=(niedge(2)-1)/i1

c         1st phase of lifting for DECOMPOSING
          do ky=1,i2y
           nny(1)=(ky-1)*i1+1
           nny(2)=nny(1)+ic
           nny(3)=nny(2)+ic
           do kx=1,i2x
            nnx(1)=(kx-1)*i1+1
            nnx(2)=nnx(1)+ic
            nnx(3)=nnx(2)+ic
            do iy=1,3
             iy0=(nny(iy)-1)*niedge(1)
             do ix=1,3
              ix0=nnx(ix)
              ip=(iy-1)*3+ix
              local(ip)=iy0+ix0
             end do
            end do
            if(isign.eq.1)then
             tmp(local(5))=tmp(local(5))
     -                    -0.5d0*tmp(local(1))-0.5d0*tmp(local(9))
             tmp(local(6))=tmp(local(6))
     -                    -0.5d0*tmp(local(3))-0.5d0*tmp(local(9))
             tmp(local(8))=tmp(local(8))
     -                    -0.5d0*tmp(local(7))-0.5d0*tmp(local(9))
             if(kx.eq.1)then
              tmp(local(4))=tmp(local(4))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(7))
             elseif(ky.eq.1)then
              tmp(local(2))=tmp(local(2))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(3))
             endif
            elseif(isign.eq.2)then
             tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(5))
             tmp(local(9))=tmp(local(9))+0.5d0*tmp(local(5))
             tmp(local(3))=tmp(local(3))+0.5d0*tmp(local(6))
             tmp(local(9))=tmp(local(9))+0.5d0*tmp(local(6))
             tmp(local(7))=tmp(local(7))+0.5d0*tmp(local(8))
             tmp(local(9))=tmp(local(9))+0.5d0*tmp(local(8))
             if(kx.eq.1)then
              tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(4))
              tmp(local(7))=tmp(local(7))+0.5d0*tmp(local(4))
             elseif(ky.eq.1)then
              tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(2))
              tmp(local(3))=tmp(local(3))+0.5d0*tmp(local(2))
             endif
            endif
           enddo
          enddo

c         2nd phase of lifting for DECOMPOSING
          do ky=1,i2y
           nny(1)=(ky-1)*i1+1
           nny(2)=nny(1)+ic
           nny(3)=nny(2)+ic
           do kx=1,i2x
            nnx(1)=(kx-1)*i1+1
            nnx(2)=nnx(1)+ic
            nnx(3)=nnx(2)+ic
            do iy=1,3
             iy0=(nny(iy)-1)*niedge(1)
             do ix=1,3
              ix0=nnx(ix)
              ip=(iy-1)*3+ix
              local(ip)=iy0+ix0
             end do
            end do
            if(isign.eq.2)then
             tmp(local(5))=tmp(local(5))
     -                    -0.25d0*tmp(local(1))-0.25d0*tmp(local(9))
             tmp(local(6))=tmp(local(6))
     -                    -0.25d0*tmp(local(3))-0.25d0*tmp(local(9))
             tmp(local(8))=tmp(local(8))
     -                    -0.25d0*tmp(local(7))-0.25d0*tmp(local(9))
             if(kx.eq.1)then
              tmp(local(4))=tmp(local(4))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(7))
             elseif(ky.eq.1)then
              tmp(local(2))=tmp(local(2))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(3))
             endif
            elseif(isign.eq.1)then
             tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(5))
             tmp(local(9))=tmp(local(9))+0.25d0*tmp(local(5))
             tmp(local(3))=tmp(local(3))+0.25d0*tmp(local(6))
             tmp(local(9))=tmp(local(9))+0.25d0*tmp(local(6))
             tmp(local(7))=tmp(local(7))+0.25d0*tmp(local(8))
             tmp(local(9))=tmp(local(9))+0.25d0*tmp(local(8))
             if(kx.eq.1)then
              tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(4))
              tmp(local(7))=tmp(local(7))+0.25d0*tmp(local(4))
             elseif(ky.eq.1)then
              tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(2))
              tmp(local(3))=tmp(local(3))+0.25d0*tmp(local(2))
             endif
            endif
           enddo
          enddo
         enddo     

        elseif(isign.lt.0)then
c        Synthesizing.....

         do j=nlevelT,2,-1
          i1=2**(j-1)
          ic=i1/2
          i2x=(niedge(1)-1)/i1
          i2y=(niedge(2)-1)/i1

c         undo 2nd phase of DECOMPOSING when synthesizing
          do ky=1,i2y
           nny(1)=(ky-1)*i1+1
           nny(2)=nny(1)+ic
           nny(3)=nny(2)+ic
           do kx=1,i2x
            nnx(1)=(kx-1)*i1+1
            nnx(2)=nnx(1)+ic
            nnx(3)=nnx(2)+ic
            do iy=1,3
             iy0=(nny(iy)-1)*niedge(1)
             do ix=1,3
              ix0=nnx(ix)
              ip=(iy-1)*3+ix
              local(ip)=iy0+ix0
             end do
            end do
            if(isign.eq.-2)then
             if(j.le.nlevel_synT)then
              tmp(local(5))=tmp(local(5))
     -                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(9))
              tmp(local(6))=tmp(local(6))
     -                     +0.25d0*tmp(local(3))+0.25d0*tmp(local(9))
              tmp(local(8))=tmp(local(8))
     -                     +0.25d0*tmp(local(7))+0.25d0*tmp(local(9))
              if(kx.eq.1)then
               tmp(local(4))=tmp(local(4))
     -                      +0.25d0*tmp(local(1))+0.25d0*tmp(local(7))
              elseif(ky.eq.1)then
               tmp(local(2))=tmp(local(2))
     -                      +0.25d0*tmp(local(1))+0.25d0*tmp(local(3))
              endif
             else
              tmp(local(5))=0.d0
     -                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(9))
              tmp(local(6))=0.d0
     -                     +0.25d0*tmp(local(3))+0.25d0*tmp(local(9))
              tmp(local(8))=0.d0
     -                     +0.25d0*tmp(local(7))+0.25d0*tmp(local(9))
              if(kx.eq.1)then
               tmp(local(4))=0.d0
     -                      +0.25d0*tmp(local(1))+0.25d0*tmp(local(7))
              elseif(ky.eq.1)then
               tmp(local(2))=0.d0
     -                      +0.25d0*tmp(local(1))+0.25d0*tmp(local(3))
              endif
             endif
            elseif(isign.eq.-1)then
             if(j.le.nlevel_synT)then
              tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(5))
              tmp(local(9))=tmp(local(9))-0.25d0*tmp(local(5))
              tmp(local(3))=tmp(local(3))-0.25d0*tmp(local(6))
              tmp(local(9))=tmp(local(9))-0.25d0*tmp(local(6))
              tmp(local(7))=tmp(local(7))-0.25d0*tmp(local(8))
              tmp(local(9))=tmp(local(9))-0.25d0*tmp(local(8))
              if(kx.eq.1)then
               tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(4))
               tmp(local(7))=tmp(local(7))-0.25d0*tmp(local(4))
              elseif(ky.eq.1)then
               tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(2))
               tmp(local(3))=tmp(local(3))-0.25d0*tmp(local(2))
              endif
             else
              tmp(local(5))=0.d0
              tmp(local(6))=0.d0
              tmp(local(8))=0.d0
              if(kx.eq.1)then
               tmp(local(4))=0.d0
              elseif(ky.eq.1)then
               tmp(local(2))=0.d0
              endif
             endif
            endif
           enddo
          enddo

c         undo 1st phase of DECOMPOSING when synthesizing
          do ky=1,i2y
           nny(1)=(ky-1)*i1+1
           nny(2)=nny(1)+ic
           nny(3)=nny(2)+ic
           do kx=1,i2x
            nnx(1)=(kx-1)*i1+1
            nnx(2)=nnx(1)+ic
            nnx(3)=nnx(2)+ic
            do iy=1,3
             iy0=(nny(iy)-1)*niedge(1)
             do ix=1,3
              ix0=nnx(ix)
              ip=(iy-1)*3+ix
              local(ip)=iy0+ix0
             end do
            end do
            if(isign.eq.-1)then
             if(j.le.nlevel_synT)then
              tmp(local(5))=tmp(local(5))
     -                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(9))
              tmp(local(6))=tmp(local(6))
     -                     +0.5d0*tmp(local(3))+0.5d0*tmp(local(9))
              tmp(local(8))=tmp(local(8))
     -                     +0.5d0*tmp(local(7))+0.5d0*tmp(local(9))
              if(kx.eq.1)then
               tmp(local(4))=tmp(local(4))
     -                      +0.5d0*tmp(local(1))+0.5d0*tmp(local(7))
              elseif(ky.eq.1)then
               tmp(local(2))=tmp(local(2))
     -                      +0.5d0*tmp(local(1))+0.5d0*tmp(local(3))
              endif
             else
              tmp(local(5))=0.d0
     -                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(9))
              tmp(local(6))=0.d0
     -                     +0.5d0*tmp(local(3))+0.5d0*tmp(local(9))
              tmp(local(8))=0.d0
     -                     +0.5d0*tmp(local(7))+0.5d0*tmp(local(9))
              if(kx.eq.1)then
               tmp(local(4))=0.d0
     -                      +0.5d0*tmp(local(1))+0.5d0*tmp(local(7))
              elseif(ky.eq.1)then
               tmp(local(2))=0.d0
     -                      +0.5d0*tmp(local(1))+0.5d0*tmp(local(3))
              endif
             endif
            elseif(isign.eq.-2)then
             if(j.le.nlevel_synT)then
              tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(5))
              tmp(local(9))=tmp(local(9))-0.5d0*tmp(local(5))
              tmp(local(3))=tmp(local(3))-0.5d0*tmp(local(6))
              tmp(local(9))=tmp(local(9))-0.5d0*tmp(local(6))
              tmp(local(7))=tmp(local(7))-0.5d0*tmp(local(8))
              tmp(local(9))=tmp(local(9))-0.5d0*tmp(local(8))
              if(kx.eq.1)then
               tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(4))
               tmp(local(7))=tmp(local(7))-0.5d0*tmp(local(4))
              elseif(ky.eq.1)then
               tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(2))
               tmp(local(3))=tmp(local(3))-0.5d0*tmp(local(2))
              endif
             else
              tmp(local(5))=0.d0
              tmp(local(6))=0.d0
              tmp(local(8))=0.d0
              if(kx.eq.1)then
               tmp(local(4))=0.d0
              elseif(ky.eq.1)then
               tmp(local(2))=0.d0
              endif
             endif
            endif
           enddo
          enddo
         enddo

        endif

        do iiface=1,niface
         iface=iiface0+iiface
         image(iface)=tmp(iiface)
        enddo

       enddo
       return
       end


c___________________________________________________________________________
      subroutine face_lift_Nb_tetra (image, isign, nlevel_syn)
       implicit none
c
c   Due to the undesired characteristics of tensor products, it is decided
c    that 2D triangular mesh and 3D tetrahedral tellesation is preferable
c    when invoked in inverse formulation
c    this version is the lateral triangular 2D version
c    that might be invoked with hybrid vertical queeling + lateral multiscale
c    <2002/8/25>
c
c   change the order of processing (rotate X,Y,Z first within each level)
c   <2001/4/13> (lift_b31.f vs. lift_b31_level.f)
c
c   biorthogonal bior3.1 wavelet <2001/4/3>
c
c   isign:  1,  IMAGE analysis (decomposition)
c          -1,  IMAGE synthesis (reconstruction)
c           0,  doing nothing
c           2,  PATH-integral analysis
c          -2,  PATH-integral synthesis
c
c      notice also that, face_conf contains the reorgnized block
c      numbering in scale hirachi, when the original numbering is
c      the usual:
c
c         (((1,nx),ny,nz)
c
c    compare with the pyrimal form in spherical wavelet (or other
c     plannar triangulation grid)
c
       integer mlevelx, mlevely, mlevelz, mx, my, mz
       integer mivolu, mbase, mvolu, nbase
       parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,           mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,           mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,           mbase=1,mvolu=mivolu*mbase,nbase=mbase)
       integer nlevel(3),isign,nlevel_syn(3),niedge(3),
     ,         nivolu,iivolu0,iivolu,ivolu,nvolu,niface,
     ,         nnx(3),nny(3),nnz(3),local(27),
     ,         i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_syntmp(3),i3,ii0
       integer ip, i0, ix, ix0, i2x, iy, iy0, i2y, iz, iz0, i2z
       integer kx, ky, kz
       integer niedge12, nlevelT

       double precision image(mvolu)
       double precision y(mivolu),tmp(mivolu)
       common/basics/nlevel,nivolu,nvolu,niface,niedge

       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       end if
       niedge12=niedge(1)*niedge(2)
       nlevelT=min0(nlevel(1),nlevel(2),nlevel(3))

       do i0=1,nbase
        iivolu0=(i0-1)*nivolu

        if(isign.gt.0)then

c        Decomposition or Analyzing....

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          tmp(iivolu)=image(ivolu)
         end do

         do j=2,nlevelT
          i1=2**(j-1)
          ic=i1/2

          i2x=(niedge(1)-1)/i1
          i2y=(niedge(2)-1)/i1
          i2z=(niedge(3)-1)/i1
c         1st phase of lifting for DECOMPOSING
          do kz=1,i2z
           nnz(1)=(kz-1)*i1+1
           nnz(2)=nnz(1)+ic
           nnz(3)=nnz(2)+ic
           do ky=1,i2y
            nny(1)=(ky-1)*i1+1
            nny(2)=nny(1)+ic
            nny(3)=nny(2)+ic
            do kx=1,i2x
             nnx(1)=(kx-1)*i1+1
             nnx(2)=nnx(1)+ic
             nnx(3)=nnx(2)+ic
             do iz=1,3
              iz0=(nnz(iz)-1)*niedge12
              do iy=1,3
               iy0=(nny(iy)-1)*niedge(1)
               do ix=1,3
                ix0=nnx(ix)
                ip=(iz-1)*9+(iy-1)*3+ix
                local(ip)=iz0+iy0+ix0
               end do
              end do
             end do
             if(isign.eq.1)then
              tmp(local(2))=tmp(local(2))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(3))
              tmp(local(4))=tmp(local(4))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(7))
              tmp(local(10))=tmp(local(10))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(19))
              if(kx.eq.1)then
               tmp(local(13))=tmp(local(13))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(25))
              elseif(kx.eq.i2x)then
               tmp(local(15))=tmp(local(15))
     -                     -0.5d0*tmp(local(9))-0.5d0*tmp(local(21))
              else
               tmp(local(13))=tmp(local(13))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(25))
     -                     -0.5d0*tmp(local(7))-0.5d0*tmp(local(19))
              endif
              if(ky.eq.1)then
               tmp(local(11))=tmp(local(11))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(21))
              elseif(ky.eq.i2y)then
               tmp(local(17))=tmp(local(17))
     -                     -0.5d0*tmp(local(9))-0.5d0*tmp(local(25))
              else
               tmp(local(11))=tmp(local(11))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(21))
     -                     -0.5d0*tmp(local(3))-0.5d0*tmp(local(19))
              endif
              if(kz.eq.1)then
               tmp(local(5))=tmp(local(5))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(9))
              elseif(kz.eq.i2z)then
               tmp(local(23))=tmp(local(23))
     -                     -0.5d0*tmp(local(21))-0.5d0*tmp(local(25))
              else
               tmp(local(5))=tmp(local(5))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(9))
     -                     -0.5d0*tmp(local(3))-0.5d0*tmp(local(7))
              endif
              tmp(local(14))=tmp(local(14))
     -                     -0.5d0*tmp(local(1))-0.5d0*tmp(local(9))
     -                     -0.5d0*tmp(local(3))-0.5d0*tmp(local(7))
     -                     -0.5d0*tmp(local(19))-0.5d0*tmp(local(27))
     -                     -0.5d0*tmp(local(21))-0.5d0*tmp(local(25))
             elseif(isign.eq.2)then
              tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(2))
              tmp(local(3))=tmp(local(3))+0.5d0*tmp(local(2))
              tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(4))
              tmp(local(7))=tmp(local(7))+0.5d0*tmp(local(4))
              tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(10))
              tmp(local(19))=tmp(local(19))+0.5d0*tmp(local(10))
              if(kx.eq.1)then
               tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(13))
               tmp(local(25))=tmp(local(25))+0.5d0*tmp(local(13))
              elseif(kx.eq.i2x)then
               tmp(local(9))=tmp(local(9))+0.5d0*tmp(local(15))
               tmp(local(21))=tmp(local(21))+0.5d0*tmp(local(15))
              else
               tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(13))
               tmp(local(7))=tmp(local(7))+0.5d0*tmp(local(13))
               tmp(local(19))=tmp(local(19))+0.5d0*tmp(local(13))
               tmp(local(25))=tmp(local(25))+0.5d0*tmp(local(13))
              endif
              if(ky.eq.1)then
               tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(11))
               tmp(local(21))=tmp(local(21))+0.5d0*tmp(local(11))
              elseif(ky.eq.i2y)then
               tmp(local(9))=tmp(local(9))+0.5d0*tmp(local(17))
               tmp(local(25))=tmp(local(25))+0.5d0*tmp(local(17))
              else
               tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(11))
               tmp(local(3))=tmp(local(3))+0.5d0*tmp(local(11))
               tmp(local(19))=tmp(local(19))+0.5d0*tmp(local(11))
               tmp(local(21))=tmp(local(21))+0.5d0*tmp(local(11))
              endif
              if(kz.eq.1)then
               tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(5))
               tmp(local(9))=tmp(local(9))+0.5d0*tmp(local(5))
              elseif(kz.eq.i2z)then
               tmp(local(21))=tmp(local(21))+0.5d0*tmp(local(23))
               tmp(local(25))=tmp(local(25))+0.5d0*tmp(local(23))
              else
               tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(5))
               tmp(local(3))=tmp(local(3))+0.5d0*tmp(local(5))
               tmp(local(7))=tmp(local(7))+0.5d0*tmp(local(5))
               tmp(local(9))=tmp(local(9))+0.5d0*tmp(local(5))
              endif
              tmp(local(1))=tmp(local(1))+0.5d0*tmp(local(14))
              tmp(local(3))=tmp(local(3))+0.5d0*tmp(local(14))
              tmp(local(7))=tmp(local(7))+0.5d0*tmp(local(14))
              tmp(local(9))=tmp(local(9))+0.5d0*tmp(local(14))
              tmp(local(19))=tmp(local(19))+0.5d0*tmp(local(14))
              tmp(local(21))=tmp(local(21))+0.5d0*tmp(local(14))
              tmp(local(25))=tmp(local(25))+0.5d0*tmp(local(14))
              tmp(local(27))=tmp(local(27))+0.5d0*tmp(local(14))
             endif
            end do
           end do
          end do

c         2nd phase of lifting for DECOMPOSING
          do kz=1,i2z
           nnz(1)=(kz-1)*i1+1
           nnz(2)=nnz(1)+ic
           nnz(3)=nnz(2)+ic
           do ky=1,i2y
            nny(1)=(ky-1)*i1+1
            nny(2)=nny(1)+ic
            nny(3)=nny(2)+ic
            do kx=1,i2x
             nnx(1)=(kx-1)*i1+1
             nnx(2)=nnx(1)+ic
             nnx(3)=nnx(2)+ic
             do iz=1,3
              iz0=(nnz(iz)-1)*niedge12
              do iy=1,3
               iy0=(nny(iy)-1)*niedge(1)
               do ix=1,3
                ix0=nnx(ix)
                ip=(iz-1)*9+(iy-1)*3+ix
                local(ip)=iz0+iy0+ix0
               end do
              end do
             end do
c            undone 2nd phase of lifting first during SYNTHESIZING
             if(isign.eq.1)then
              tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(2))
              tmp(local(3))=tmp(local(3))+0.25d0*tmp(local(2))
              tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(4))
              tmp(local(7))=tmp(local(7))+0.25d0*tmp(local(4))
              tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(10))
              tmp(local(19))=tmp(local(19))+0.25d0*tmp(local(10))
              if(kx.eq.1)then
               tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(13))
               tmp(local(25))=tmp(local(25))+0.25d0*tmp(local(13))
              elseif(kx.eq.i2x)then
               tmp(local(9))=tmp(local(9))+0.25d0*tmp(local(15))
               tmp(local(21))=tmp(local(21))+0.25d0*tmp(local(15))
              else
               tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(13))
               tmp(local(7))=tmp(local(7))+0.25d0*tmp(local(13))
               tmp(local(19))=tmp(local(19))+0.25d0*tmp(local(13))
               tmp(local(25))=tmp(local(25))+0.25d0*tmp(local(13))
              endif
              if(ky.eq.1)then
               tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(11))
               tmp(local(21))=tmp(local(21))+0.25d0*tmp(local(11))
              elseif(ky.eq.i2y)then
               tmp(local(9))=tmp(local(9))+0.25d0*tmp(local(17))
               tmp(local(25))=tmp(local(25))+0.25d0*tmp(local(17))
              else
               tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(11))
               tmp(local(3))=tmp(local(3))+0.25d0*tmp(local(11))
               tmp(local(19))=tmp(local(19))+0.25d0*tmp(local(11))
               tmp(local(21))=tmp(local(21))+0.25d0*tmp(local(11))
              endif
              if(kz.eq.1)then
               tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(5))
               tmp(local(9))=tmp(local(9))+0.25d0*tmp(local(5))
              elseif(kz.eq.i2z)then
               tmp(local(21))=tmp(local(21))+0.25d0*tmp(local(23))
               tmp(local(25))=tmp(local(25))+0.25d0*tmp(local(23))
              else
               tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(5))
               tmp(local(3))=tmp(local(3))+0.25d0*tmp(local(5))
               tmp(local(7))=tmp(local(7))+0.25d0*tmp(local(5))
               tmp(local(9))=tmp(local(9))+0.25d0*tmp(local(5))
              endif
              tmp(local(1))=tmp(local(1))+0.25d0*tmp(local(14))
              tmp(local(3))=tmp(local(3))+0.25d0*tmp(local(14))
              tmp(local(7))=tmp(local(7))+0.25d0*tmp(local(14))
              tmp(local(9))=tmp(local(9))+0.25d0*tmp(local(14))
              tmp(local(19))=tmp(local(19))+0.25d0*tmp(local(14))
              tmp(local(21))=tmp(local(21))+0.25d0*tmp(local(14))
              tmp(local(25))=tmp(local(25))+0.25d0*tmp(local(14))
              tmp(local(27))=tmp(local(27))+0.25d0*tmp(local(14))
             elseif(isign.eq.2)then
              tmp(local(2))=tmp(local(2))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(3))
              tmp(local(4))=tmp(local(4))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(7))
              tmp(local(10))=tmp(local(10))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(19))
              if(kx.eq.1)then
               tmp(local(13))=tmp(local(13))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(25))
              elseif(kx.eq.i2x)then
               tmp(local(15))=tmp(local(15))
     -                     -0.25d0*tmp(local(9))-0.25d0*tmp(local(21))
              else
               tmp(local(13))=tmp(local(13))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(25))
     -                     -0.25d0*tmp(local(7))-0.25d0*tmp(local(19))
              endif
              if(ky.eq.1)then
               tmp(local(11))=tmp(local(11))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(21))
              elseif(ky.eq.i2y)then
               tmp(local(17))=tmp(local(17))
     -                     -0.25d0*tmp(local(9))-0.25d0*tmp(local(25))
              else
               tmp(local(11))=tmp(local(11))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(21))
     -                     -0.25d0*tmp(local(3))-0.25d0*tmp(local(19))
              endif
              if(kz.eq.1)then
               tmp(local(5))=tmp(local(5))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(9))
              elseif(kz.eq.i2z)then
               tmp(local(23))=tmp(local(23))
     -                     -0.25d0*tmp(local(21))-0.25d0*tmp(local(25))
              else
               tmp(local(5))=tmp(local(5))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(9))
     -                     -0.25d0*tmp(local(3))-0.25d0*tmp(local(7))
              endif
              tmp(local(14))=tmp(local(14))
     -                     -0.25d0*tmp(local(1))-0.25d0*tmp(local(9))
     -                     -0.25d0*tmp(local(3))-0.25d0*tmp(local(7))
     -                     -0.25d0*tmp(local(19))-0.25d0*tmp(local(27))
     -                     -0.25d0*tmp(local(21))-0.25d0*tmp(local(25))
             endif
            end do
           end do
          end do

         end do

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do

        elseif(isign.lt.0)then

c        Synthesizing.....

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          y(iivolu)=image(ivolu)
         end do

         do j=nlevelT,2,-1
          i1=2**(j-1)
          ic=i1/2

c         The order of synthesizing reverses the order of analyzing
          do kz=1,i2z
           nnz(1)=(kz-1)*i1+1
           nnz(2)=nnz(1)+ic
           nnz(3)=nnz(2)+ic
           do ky=1,i2y
            nny(1)=(ky-1)*i1+1
            nny(2)=nny(1)+ic
            nny(3)=nny(2)+ic
            do kx=1,i2x
             nnx(1)=(kx-1)*i1+1
             nnx(2)=nnx(1)+ic
             nnx(3)=nnx(2)+ic
             do iz=1,3
              iz0=(nnz(iz)-1)*niedge12
              do iy=1,3
               iy0=(nny(iy)-1)*niedge(1)
               do ix=1,3
                ix0=nnx(ix)
                ip=(iz-1)*9+(iy-1)*3+ix
                local(ip)=iz0+iy0+ix0
               end do
              end do
             end do
             if(isign.eq.-1)then
              tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(2))
              tmp(local(3))=tmp(local(3))-0.25d0*tmp(local(2))
              tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(4))
              tmp(local(7))=tmp(local(7))-0.25d0*tmp(local(4))
              tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(10))
              tmp(local(19))=tmp(local(19))-0.25d0*tmp(local(10))
              if(kx.eq.1)then
               tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(13))
               tmp(local(25))=tmp(local(25))-0.25d0*tmp(local(13))
              elseif(kx.eq.i2x)then
               tmp(local(9))=tmp(local(9))-0.25d0*tmp(local(15))
               tmp(local(21))=tmp(local(21))-0.25d0*tmp(local(15))
              else
               tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(13))
               tmp(local(7))=tmp(local(7))-0.25d0*tmp(local(13))
               tmp(local(19))=tmp(local(19))-0.25d0*tmp(local(13))
               tmp(local(25))=tmp(local(25))-0.25d0*tmp(local(13))
              endif
              if(ky.eq.1)then
               tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(11))
               tmp(local(21))=tmp(local(21))-0.25d0*tmp(local(11))
              elseif(ky.eq.i2y)then
               tmp(local(9))=tmp(local(9))-0.25d0*tmp(local(17))
               tmp(local(25))=tmp(local(25))-0.25d0*tmp(local(17))
              else
               tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(11))
               tmp(local(3))=tmp(local(3))-0.25d0*tmp(local(11))
               tmp(local(19))=tmp(local(19))-0.25d0*tmp(local(11))
               tmp(local(21))=tmp(local(21))-0.25d0*tmp(local(11))
              endif
              if(kz.eq.1)then
               tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(5))
               tmp(local(9))=tmp(local(9))-0.25d0*tmp(local(5))
              elseif(kz.eq.i2z)then
               tmp(local(21))=tmp(local(21))-0.25d0*tmp(local(23))
               tmp(local(25))=tmp(local(25))-0.25d0*tmp(local(23))
              else
               tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(5))
               tmp(local(3))=tmp(local(3))-0.25d0*tmp(local(5))
               tmp(local(7))=tmp(local(7))-0.25d0*tmp(local(5))
               tmp(local(9))=tmp(local(9))-0.25d0*tmp(local(5))
              endif
              tmp(local(1))=tmp(local(1))-0.25d0*tmp(local(14))
              tmp(local(3))=tmp(local(3))-0.25d0*tmp(local(14))
              tmp(local(7))=tmp(local(7))-0.25d0*tmp(local(14))
              tmp(local(9))=tmp(local(9))-0.25d0*tmp(local(14))
              tmp(local(19))=tmp(local(19))-0.25d0*tmp(local(14))
              tmp(local(21))=tmp(local(21))-0.25d0*tmp(local(14))
              tmp(local(25))=tmp(local(25))-0.25d0*tmp(local(14))
              tmp(local(27))=tmp(local(27))-0.25d0*tmp(local(14))
             elseif(isign.eq.-2)then
              tmp(local(2))=tmp(local(2))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(3))
              tmp(local(4))=tmp(local(4))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(7))
              tmp(local(10))=tmp(local(10))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(19))
              if(kx.eq.1)then
               tmp(local(13))=tmp(local(13))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(25))
              elseif(kx.eq.i2x)then
               tmp(local(15))=tmp(local(15))
     +                     +0.25d0*tmp(local(9))+0.25d0*tmp(local(21))
              else
               tmp(local(13))=tmp(local(13))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(25))
     +                     +0.25d0*tmp(local(7))+0.25d0*tmp(local(19))
              endif
              if(ky.eq.1)then
               tmp(local(11))=tmp(local(11))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(21))
              elseif(ky.eq.i2y)then
               tmp(local(17))=tmp(local(17))
     +                     +0.25d0*tmp(local(9))+0.25d0*tmp(local(25))
              else
               tmp(local(11))=tmp(local(11))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(21))
     +                     +0.25d0*tmp(local(3))+0.25d0*tmp(local(19))
              endif
              if(kz.eq.1)then
               tmp(local(5))=tmp(local(5))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(9))
              elseif(kz.eq.i2z)then
               tmp(local(23))=tmp(local(23))
     +                     +0.25d0*tmp(local(21))+0.25d0*tmp(local(25))
              else
               tmp(local(5))=tmp(local(5))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(9))
     +                     +0.25d0*tmp(local(3))+0.25d0*tmp(local(7))
              endif
              tmp(local(14))=tmp(local(14))
     +                     +0.25d0*tmp(local(1))+0.25d0*tmp(local(9))
     +                     +0.25d0*tmp(local(3))+0.25d0*tmp(local(7))
     +                     +0.25d0*tmp(local(19))+0.25d0*tmp(local(27))
     +                     +0.25d0*tmp(local(21))+0.25d0*tmp(local(25))
             endif
            end do
           end do
          end do

          do kz=1,i2z
           nnz(1)=(kz-1)*i1+1
           nnz(2)=nnz(1)+ic
           nnz(3)=nnz(2)+ic
           do ky=1,i2y
            nny(1)=(ky-1)*i1+1
            nny(2)=nny(1)+ic
            nny(3)=nny(2)+ic
            do kx=1,i2x
             nnx(1)=(kx-1)*i1+1
             nnx(2)=nnx(1)+ic
             nnx(3)=nnx(2)+ic
             do iz=1,3
              iz0=(nnz(iz)-1)*niedge12
              do iy=1,3
               iy0=(nny(iy)-1)*niedge(1)
               do ix=1,3
                ix0=nnx(ix)
                ip=(iz-1)*9+(iy-1)*3+ix
                local(ip)=iz0+iy0+ix0
               end do
              end do
             end do
             if(isign.eq.-1)then
              tmp(local(2))=tmp(local(2))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(3))
              tmp(local(4))=tmp(local(4))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(7))
              tmp(local(10))=tmp(local(10))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(19))
              if(kx.eq.1)then
               tmp(local(13))=tmp(local(13))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(25))
              elseif(kx.eq.i2x)then
               tmp(local(15))=tmp(local(15))
     +                     +0.5d0*tmp(local(9))+0.5d0*tmp(local(21))
              else
               tmp(local(13))=tmp(local(13))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(25))
     +                     +0.5d0*tmp(local(7))+0.5d0*tmp(local(19))
              endif
              if(ky.eq.1)then
               tmp(local(11))=tmp(local(11))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(21))
              elseif(ky.eq.i2y)then
               tmp(local(17))=tmp(local(17))
     +                     +0.5d0*tmp(local(9))+0.5d0*tmp(local(25))
              else
               tmp(local(11))=tmp(local(11))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(21))
     +                     +0.5d0*tmp(local(3))+0.5d0*tmp(local(19))
              endif
              if(kz.eq.1)then
               tmp(local(5))=tmp(local(5))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(9))
              elseif(kz.eq.i2z)then
               tmp(local(23))=tmp(local(23))
     +                     +0.5d0*tmp(local(21))+0.5d0*tmp(local(25))
              else
               tmp(local(5))=tmp(local(5))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(9))
     +                     +0.5d0*tmp(local(3))+0.5d0*tmp(local(7))
              endif
              tmp(local(14))=tmp(local(14))
     +                     +0.5d0*tmp(local(1))+0.5d0*tmp(local(9))
     +                     +0.5d0*tmp(local(3))+0.5d0*tmp(local(7))
     +                     +0.5d0*tmp(local(19))+0.5d0*tmp(local(27))
     +                     +0.5d0*tmp(local(21))+0.5d0*tmp(local(25))
             elseif(isign.eq.-2)then
              tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(2))
              tmp(local(3))=tmp(local(3))-0.5d0*tmp(local(2))
              tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(4))
              tmp(local(7))=tmp(local(7))-0.5d0*tmp(local(4))
              tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(10))
              tmp(local(19))=tmp(local(19))-0.5d0*tmp(local(10))
              if(kx.eq.1)then
               tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(13))
               tmp(local(25))=tmp(local(25))-0.5d0*tmp(local(13))
              elseif(kx.eq.i2x)then
               tmp(local(9))=tmp(local(9))-0.5d0*tmp(local(15))
               tmp(local(21))=tmp(local(21))-0.5d0*tmp(local(15))
              else
               tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(13))
               tmp(local(7))=tmp(local(7))-0.5d0*tmp(local(13))
               tmp(local(19))=tmp(local(19))-0.5d0*tmp(local(13))
               tmp(local(25))=tmp(local(25))-0.5d0*tmp(local(13))
              endif
              if(ky.eq.1)then
               tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(11))
               tmp(local(21))=tmp(local(21))-0.5d0*tmp(local(11))
              elseif(ky.eq.i2y)then
               tmp(local(9))=tmp(local(9))-0.5d0*tmp(local(17))
               tmp(local(25))=tmp(local(25))-0.5d0*tmp(local(17))
              else
               tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(11))
               tmp(local(3))=tmp(local(3))-0.5d0*tmp(local(11))
               tmp(local(19))=tmp(local(19))-0.5d0*tmp(local(11))
               tmp(local(21))=tmp(local(21))-0.5d0*tmp(local(11))
              endif
              if(kz.eq.1)then
               tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(5))
               tmp(local(9))=tmp(local(9))-0.5d0*tmp(local(5))
              elseif(kz.eq.i2z)then
               tmp(local(21))=tmp(local(21))-0.5d0*tmp(local(23))
               tmp(local(25))=tmp(local(25))-0.5d0*tmp(local(23))
              else
               tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(5))
               tmp(local(3))=tmp(local(3))-0.5d0*tmp(local(5))
               tmp(local(7))=tmp(local(7))-0.5d0*tmp(local(5))
               tmp(local(9))=tmp(local(9))-0.5d0*tmp(local(5))
              endif
              tmp(local(1))=tmp(local(1))-0.5d0*tmp(local(14))
              tmp(local(3))=tmp(local(3))-0.5d0*tmp(local(14))
              tmp(local(7))=tmp(local(7))-0.5d0*tmp(local(14))
              tmp(local(9))=tmp(local(9))-0.5d0*tmp(local(14))
              tmp(local(19))=tmp(local(19))-0.5d0*tmp(local(14))
              tmp(local(21))=tmp(local(21))-0.5d0*tmp(local(14))
              tmp(local(25))=tmp(local(25))-0.5d0*tmp(local(14))
              tmp(local(27))=tmp(local(27))-0.5d0*tmp(local(14))
             endif
            end do
           end do
          end do

         end do

         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=y(iivolu)
         end do

        endif

       end do

       return
       end


c____________________________________________________________________
      subroutine Gaussian(sigma,sigmaz)
      implicit none
      integer isigma, isigmaz, itt
      integer i, ii, ix, iy, iz

      integer nlevel(3),nivolu,nvolu,niface,niedge(3)
      integer iw(5000)
      double precision w(5000),xx,yy,zz,wsum,sigma,sigmaz
      common/basics/nlevel,nivolu,nvolu,niface,niedge
      common/gauss/w,iw,isigma,isigmaz,itt

      wsum=0.d0
      isigma=idint(sigma)+1
      isigmaz=idint(sigmaz)+1
      if(sigma.eq.0.d0)isigma=0
      if(sigmaz.eq.0.d0)isigmaz=0
      itt=2*isigma+1
      do iz=-isigmaz,isigmaz
       if(sigmaz.eq.0.d0)then
        zz=0.d0
       else
        zz=(dabs(dble(iz))/sigmaz)**2
       endif
       do iy=-isigma,isigma
        if(sigma.eq.0.d0)then
         yy=0.d0
        else
         yy=(dabs(dble(iy))/sigma)**2
        endif
        do ix=-isigma,isigma
         if(sigma.eq.0.d0)then
          xx=0.d0
         else
          xx=(dabs(dble(ix))/sigma)**2
         endif
         ii=(iz+isigmaz)*itt*itt+(iy+isigma)*itt+(ix+isigma)+1
         if(ii.gt.5000)stop' Too many quelling coefficients!'
         iw(ii)=iz*niface+iy*niedge(1)+ix
         w(ii)=dexp(-(xx+yy+zz))
         wsum=wsum+w(ii)
        end do
       end do
      end do
      ii=(2*isigmaz+1)*itt*itt
      do i=1,ii
       w(i)=w(i)/wsum
      end do
      write(*,*)isigma,itt,ii
      return
      end



c________________________________________________________________________
       subroutine quell(image)
       implicit none
       integer mlevelx, mlevely, mlevelz, mx, my, mz
       integer mivolu, mbase, mvolu, nbase
       integer isigma,isigmaz,itt
       parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,           mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,           mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,           mbase=1,mvolu=mivolu*mbase,nbase=1)
       integer nlevel(3),nivolu,nvolu,niface,niedge(3)
       integer i0, ii, ivolu, iivolu, iivolu0
       integer ix, ixb, ixe
       integer iy, iyb, iye
       integer iz, izb, ize
       integer jj, jx, jy, jz
       integer iw(5000)
       double precision w(5000)
       double precision image(mvolu)
       double precision y0(mivolu),y(mivolu)
       common/basics/nlevel,nivolu,nvolu,niface,niedge
       common/gauss/w,iw,isigma,isigmaz,itt

       do i0=1,nbase
        iivolu0=(i0-1)*nivolu

        do iivolu=1,nivolu
         ivolu=iivolu0+iivolu
         y0(iivolu)=image(ivolu)
         y(iivolu)=0.d0
        end do

        do iz=1,niedge(3)
         izb=iz-isigmaz
         if(izb.lt.1)izb=1
         ize=iz+isigmaz
         if(ize.gt.niedge(3))ize=niedge(3)
         do iy=1,niedge(2)
          iyb=iy-isigma
          if(iyb.lt.1)iyb=1
          iye=iy+isigma
          if(iye.gt.niedge(2))iye=niedge(2)
          do ix=1,niedge(1)
           ixb=ix-isigma
           if(ixb.lt.1)ixb=1
           ixe=ix+isigma
           if(ixe.gt.niedge(1))ixe=niedge(1)
           ii=(iz-1)*niface+(iy-1)*niedge(1)+ix
           do jz=izb-iz,ize-iz
            do jy=iyb-iy,iye-iy
             do jx=ixb-ix,ixe-ix
              jj=(jz+isigmaz)*itt*itt+(jy+isigma)*itt+(jx+isigma)+1
              y(ii)=y(ii)+y0(ii+iw(jj))*w(jj)
             end do
            end do
           end do
          end do
         end do
        end do

        do iivolu=1,nivolu
         ivolu=iivolu0+iivolu
         image(ivolu)=y(iivolu)
        end do

       end do
       return
       end


c_______________________________________________________________________
      subroutine quellS(nn,y)
      implicit none
      integer nn
      integer mbase, mlevel, miface, mface
      integer mneigh, medge_head, mivert,mvert
      integer mlevelV, mvolu

      parameter(mbase=20,mlevel=5,miface=4**(mlevel-1),
     ,          mface=miface*mbase,mneigh=2000,
     ,          medge_head=2**(mlevel-1)+1,
     ,          mivert=((medge_head+1)*medge_head)/2,
     ,          mvert=20*mivert)
      parameter(mlevelV=6,
     ,          mvolu=(2**(mlevelV-1)+1)*mvert)
      integer windex(mvert,mneigh+1),wi(8*mvert),
     ,        niedge(2),nlevel(2)
      integer i, j, k0, k1, nface, iface, iface0, nlayer, ilayer, ir, iw
      integer maxneigh, minneigh
      integer nrr, nsigma, nsigmar, nvolu

      real*8 weigh(mvert,mneigh),ww(8*mvert),wr,weight_sum
      real*8 y(nn),ytmp(nn),sigma0
      common/GAUSS2/sigma0,weigh,windex
      common/basics_spherical/nlevel,nvolu,niedge

      nface=niedge(1)
      nrr=niedge(2)
      if(nface*nrr.ne.nn)stop' Inconsistent variable in quell!!'

      do i=1,nn
       ytmp(i)=0.d0
      enddo
      minneigh=mneigh
      maxneigh=0

      nsigma=int(sigma0)+1
      nsigmar=1
      do ir=1,nrr
       iface0=(ir-1)*nface
       k0=ir-nsigmar
       if(k0.lt.1)k0=1
       k1=ir+nsigmar
       if(k1.gt.nrr)k1=nrr
       nlayer=k1-k0+1
       do iface=1,nface
        i=iface0+iface
        weight_sum=0.d0
        iw=0
        do ilayer=k0,k1
         wr=dble(iabs(ilayer-ir))/sigma0
         wr=1.d0/dexp(wr*wr)
         do j=1,windex(iface,1)
          weight_sum=weight_sum+weigh(iface,j)*wr
          iw=iw+1
          wi(iw)=(ilayer-1)*nface+windex(iface,j+1)
          ww(iw)=weigh(iface,j)*wr
         end do
        end do
        minneigh=min0(minneigh,iw)
        maxneigh=max0(maxneigh,iw)
        weight_sum=1.d0/weight_sum
        do j=1,iw
         ytmp(wi(j))=ytmp(wi(j))+y(i)*ww(j)*weight_sum
        end do
       end do
      end do

      do i=1,nn
       y(i)=ytmp(i)
      enddo
c      write(*,*)minneigh,maxneigh
      return
      end



c________________________________________________________________________
      subroutine prequell(vertg_x,sigma0)
      implicit none
      integer mbase, mlevel, miface, mface, mneigh
      integer medge_head, mivert, mvert
      parameter(mbase=20,mlevel=5,miface=4**(mlevel-1),
     ,          mface=mbase*miface,mneigh=2000,
     ,          medge_head=2**(mlevel-1)+1,
     ,          mivert=((medge_head+1)*medge_head)/2,
     ,          mvert=20*mivert)
      integer nbase,nlevelS,niface,nivert,nedge_head,
     ,        nface,nvert,nvert_global,icc,icq(mlevel-1,2),
     ,        nlevel(2),niedge(2),neigh,npoin,nelem,face(3,miface)
      integer windex(mvert,mneigh+1),ntbl(mvert,10),windex_tmp
      integer face_global(3,mface),indx(mvert),face_conf(miface),
     ,        face_index(mlevel-1),node_conf(miface*(mlevel-1),6)
      integer i,j, k, l, m, n, ii, isigma, nvolu
      integer lpoin, neighm, neighn
      real*4 vertg_x(mvert,3)
      real*8 weigh(mvert,mneigh),weigh_sum
      real*8 sigma0,sigma,sigma1,sigma00
      real*8 x0,y0,z0,x1,y1,z1,d_tmp,ddd,tmp
      common/basics_spherical/nlevel,nvolu,niedge
      common/base0/nbase,nlevelS,niface,nivert,nedge_head,
     ,             nface,nvert,nvert_global,icc,icq
      common/base1/face,face_global,indx,face_conf,
     ,             face_index,node_conf
      common/GAUSS2/sigma00,weigh,windex

      sigma00=dabs(sigma0)
      if(sigma0.eq.0.d0)then
       do i=1,nvert_global
        windex(i,1)=1
        windex(i,2)=i
        weigh(i,1)=1.d0
       enddo
       write(*,*)' No quelling !!'
       return
      endif

      neigh=0
      npoin=nvert_global
      nelem=nface
      do i=1,npoin
       ntbl(i,1)=0
       do j=1,nelem
        do k=1,3
         if(i.eq.face_global(k,j))then
          ntbl(i,1)=ntbl(i,1)+1
          ntbl(i,ntbl(i,1)+1)=j
         endif
        end do
       end do
       neigh=max0(neigh,ntbl(i,1))
      end do

      isigma=int(dabs(sigma00))+1
      neighn=mneigh
      neighm=1
      do i=1,npoin
       windex(i,1)=1
       windex(i,2)=i
       do j=1,isigma
        windex_tmp=windex(i,1)
        do k=1,windex(i,1)
         do l=2,ntbl(windex(i,k+1),1)+1
          do m=1,3
           lpoin=ntbl(l,m)
           do n=1,windex_tmp
            if(windex(i,n+1).eq.lpoin)goto 1
           enddo
           windex_tmp=windex_tmp+1
           if(windex_tmp.ge.mneigh)then
            stop' neighbooring too large within prequell!'
           endif
           windex(i,windex_tmp+1)=lpoin
1          continue
          end do
         end do
         windex(i,1)=windex_tmp
        end do
       end do
       neighn=amin0(neighn,windex(i,1))
       neighm=amax0(neighm,windex(i,1))
      end do

      sigma0=sigma00
      x0=vertg_x(1,1)
      y0=vertg_x(1,2)
      z0=vertg_x(1,3)
      tmp=1.d10
      do i=2,windex(1,1)+1
       ii=windex(1,i)
       x1=vertg_x(ii,1)
       y1=vertg_x(ii,2)
       z1=vertg_x(ii,3)
       d_tmp=x0*x1+y0*y1+z0*z1
       ddd=dacos(d_tmp)*6371.d0
       tmp=dmin1(tmp,ddd)
      end do
      sigma=sigma0*tmp
      write(*,*)' sigma in KM: ',sigma,neighn,neighm

      do i=1,npoin
       x0=vertg_x(i,1)
       y0=vertg_x(i,2)
       z0=vertg_x(i,3)
       weigh_sum=0.d0
       do j=2,windex(i,1)+1
        k=windex(i,j)
        x1=vertg_x(k,1)
        y1=vertg_x(k,2)
        z1=vertg_x(k,3)
        d_tmp=x0*x1+y0*y1+z0*z1
        ddd=dacos(d_tmp)*6371.d0
        d_tmp=ddd/sigma
        weigh(i,j-1)=1.d0/dexp(d_tmp*d_tmp)
        weigh_sum=weigh_sum+weigh(i,j-1)
       end do
       do j=1,windex(i,1)
        weigh(i,j)=weigh(i,j)/weigh_sum
       end do
      end do

      return
      end



c___________________________________________________________________________
       subroutine face_lift_Nhelix(image,isign,nlevel_syn)
       implicit none
c
c       HELIX version of multidimensional WT using helix 1D WT
c        note that the present image is organized as (((1,nx),ny),nz)
c        when passed in.  It is then reorganized from the longest strip
c        that means the basic strip consists in the direction with the
c        highest nlevel and correspondingly, the periodic B.C. will be
c        set. For example, nlevel_x=5, nlevel_y=4, nllevel_z=3,then
c        the implicitly assumed periodic B.C. on is on YZ-edges.
c        (compare with the Tensor product version in <face_lift_Nt1>)
c        (N in _Nhelix means node based <instead of pixel-based>)
c
c       notice that nx=2**(lx-1)+1 not 2**(lx-1) in all the pixel based scheme
c
       integer mlevelx, mlevely, mlevelz, mlevelxyz
       integer mx, my, mz, mivolu, mbase, mvolu, nbase
       parameter(mlevelx=8,mlevely=8,mlevelz=7,
     ,           mlevelxyz=10,
     ,           mx=2**(mlevelx-1)+1,my=2**(mlevely-1)+1,
     ,           mz=2**(mlevelz-1)+1,mivolu=mx*my*mz,
     ,           mbase=1,mvolu=mivolu*mbase,nbase=mbase)
       integer nlevel(3),isign,nlevel_syn(3),niedge(3),
     ,         nivolu,iivolu0,iivolu,ivolu,nvolu,niface,ipos_r(6),
     ,         i,i1,ic,i2,l,j,ltmp,ktmp,nlevel_syntmp(3),i3,ii0
       integer i0, id_r, ieffect, ipos, ix, iy, iz, niedge12, niedge_r
       integer nlevelT, nlevelT_syn, nx, nxy, ny, nyz, nz, nzx
       double precision image(mvolu)
       double precision tmp(mivolu),tmp0(mivolu)
       common/basics/nlevel,nivolu,nvolu,niface,niedge

       if(isign.eq.0)then
        return
       elseif(iabs(isign).gt.2)then
        write(*,*)' wrong mode entering face_lift!'
        stop
       end if
       niedge12=niedge(1)*niedge(2)
       nlevelT=max0(nlevel(1),nlevel(2),nlevel(3))
       nlevelT_syn=max0(nlevel_syn(1),nlevel_syn(2),nlevel_syn(3))
       niedge_r=2**(nlevelT-1)+1
       nx=niedge(1)
       ny=niedge(2)
       nz=niedge(3)
       nxy=nx*ny
       nyz=ny*nz
       nzx=nz*nx
           if(nlevel(1).ge.nlevel(2).and.nlevel(2).ge.nlevel(3))then
        id_r=1
       elseif(nlevel(1).ge.nlevel(3).and.nlevel(3).ge.nlevel(2))then
        id_r=2
       elseif(nlevel(2).ge.nlevel(1).and.nlevel(1).ge.nlevel(3))then
        id_r=3
       elseif(nlevel(2).ge.nlevel(3).and.nlevel(3).ge.nlevel(1))then
        id_r=4
       elseif(nlevel(3).ge.nlevel(1).and.nlevel(1).ge.nlevel(2))then
        id_r=5
       elseif(nlevel(3).ge.nlevel(2).and.nlevel(2).ge.nlevel(1))then
        id_r=6
       endif

       do i0=1,nbase
        iivolu0=(i0-1)*nivolu

        if(isign.gt.0)then
c        Decomposition or Analyzing....
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          tmp0(iivolu)=image(ivolu)
         end do
         do iz=1,nz
          do iy=1,ny
           do ix=1,nx
                 ipos=(iz-1)*nxy+(iy-1)*nx+ix
            ipos_r(1)=ipos
            ipos_r(2)=(iy-1)*nzx+(iz-1)*nx+ix
            ipos_r(3)=(iz-1)*nxy+(ix-1)*ny+iy
            ipos_r(4)=(ix-1)*nyz+(iz-1)*ny+iy
            ipos_r(5)=(ix-1)*nyz+(iy-1)*nz+iz
            ipos_r(6)=(iy-1)*nzx+(ix-1)*nz+iz
            tmp(ipos_r(id_r))=tmp0(ipos)
           end do
          end do
         end do
         do j=2,nlevelT
          i1=2**(j-1)
          ic=i1/2
          i2=(niedge_r-1)/i1

          call forward_node(i1,i2,ic,tmp,isign,1)

          do iz=1,nz
           do iy=1,ny
            do ix=1,nx
                  ipos=(iz-1)*nxy+(iy-1)*nx+ix
             ipos_r(1)=ipos
             ipos_r(2)=(iy-1)*nzx+(iz-1)*nx+ix
             ipos_r(3)=(iz-1)*nxy+(ix-1)*ny+iy
             ipos_r(4)=(ix-1)*nyz+(iz-1)*ny+iy
             ipos_r(5)=(ix-1)*nyz+(iy-1)*nz+iz
             ipos_r(6)=(iy-1)*nzx+(ix-1)*nz+iz
             tmp0(ipos)=tmp(ipos_r(id_r))
            end do
           end do
          end do
         end do
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=tmp0(iivolu)
         end do

        elseif(isign.lt.0)then

c        Synthesizing.....
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          tmp0(iivolu)=image(ivolu)
         end do
         do iz=1,nz
          do iy=1,ny
           do ix=1,nx
                 ipos=(iz-1)*nxy+(iy-1)*nx+ix
            ipos_r(1)=ipos
            ipos_r(2)=(iy-1)*nzx+(iz-1)*nx+ix
            ipos_r(3)=(iz-1)*nxy+(ix-1)*ny+iy
            ipos_r(4)=(ix-1)*nyz+(iz-1)*ny+iy
            ipos_r(5)=(ix-1)*nyz+(iy-1)*nz+iz
            ipos_r(6)=(iy-1)*nzx+(ix-1)*nz+iz
            tmp(ipos_r(id_r))=tmp0(ipos)
           end do
          end do
         end do
         do j=nlevelT,2,-1
          i1=2**(j-1)
          ic=i1/2
          i2=(niedge_r-1)/i1

          if(j.le.nlevelT_syn)then
           ieffect=1
          else
           ieffect=0
          endif
          call forward_node(i1,i2,ic,tmp,isign,ieffect)

          do iz=1,nz
           do iy=1,ny
            do ix=1,nx
                  ipos=(iz-1)*nxy+(iy-1)*nx+ix
             ipos_r(1)=ipos
             ipos_r(2)=(iy-1)*nzx+(iz-1)*nx+ix
             ipos_r(3)=(iz-1)*nxy+(ix-1)*ny+iy
             ipos_r(4)=(ix-1)*nyz+(iz-1)*ny+iy
             ipos_r(5)=(ix-1)*nyz+(iy-1)*nz+iz
             ipos_r(6)=(iy-1)*nzx+(ix-1)*nz+iz
             tmp0(ipos)=tmp(ipos_r(id_r))
            end do
           end do
          end do
         end do
         do iivolu=1,nivolu
          ivolu=iivolu0+iivolu
          image(ivolu)=tmp0(iivolu)
         end do

        endif

       end do

       return
       end



C gaussian distribution
      FUNCTION gasdev(idum)
       implicit none
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

C random generator
      FUNCTION ran1(idum)
       implicit none
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END


      FUNCTION MOD1(I,N)
       implicit none
       integer i, n, MOD1
       integer ii
C
C     Yields MOD(I,N) with the exception that if I=K*N, result is N.
C
      MOD1=0.0
      IF (I.EQ.0) RETURN
C
      IF (N.EQ.0) THEN
         WRITE(6,*)
     $  'WARNING:  Attempt to take MOD(I,0) in FUNCTION MOD1.'
         RETURN
      ENDIF
      II = I+N-1
      MOD1 = MOD(II,N)+1
      RETURN
      END


      FUNCTION NDIV1(I,N)
       implicit none
       integer i, n, ndiv1 
       integer ii
C
C     Yields I/N with the exception that if I=K*N, result is I/N-1.
C
      NDIV1=0
      IF (I.EQ.0) RETURN
C
      IF (N.EQ.0) THEN
         WRITE(6,*)
     $  'WARNING:  Attempt to take I/0 in FUNCTION NDIV1.'
         RETURN
      ENDIF
      II = I+N-1
      NDIV1 = II/N - 1
      RETURN
      END

                                   

c are returned by real*4 function fval
c e.g. amap=fval(srad,scol,slong)
c where fval, srad, scol, slong have to be declared real*4
c the returned value is dvs/vs.
c
c there are two file names hardwired into the subroutines:
c 1) block.desc (to read the model description file)
c 2) block.model (the binary block model)
c for linux machine, the block.model file needs to be swab the bytes first
c run: swab -4 < block.model > block.model.swab
c
      subroutine getv_smdl(frad,fcol,flong,amap)
      implicit real*4(a-h,o-z)
      common/blocks/r1(30),r2(30),t1(180),t2(180),
     +    pmn(50000),pmx(50000),nshell,mshell,mlat(181),nlat
      common/blk$/block(50000),nblk
c single precision stuff
      real*8 frad,fcol,flong,amap
      real*4 srad,scol,slong

c loop over angout, the angle to be output
c choose radius between 3480 (CMB) and Moho (6348)
c
c choose longitude going from 0 to 360 deg
      srad=frad
c colatitude
      scol=fcol
      if (flong.lt.0.0d0) then
         slong=flong+360.0d0
      else
         slong=flong
      endif
      amap=fval(srad,scol,slong)
c this makes percentage dvs/vs
c      write(6,112)scol,slong,srad,float(100.*amap)
112   format(f6.2,f9.2,f9.1,f9.2)

      return
      end

c################################################################################


      real*4 function fval(r0,t0,p0)
      common/blk$/block(50000),nblk
      call fblk(r0,t0,p0,ib)
      fval=block(ib)
      return
      end

      subroutine fblk(r0,t0,p0,ib)
c*** routine to find block number given radial index and colat and long
      common/blocks/r1(30),r2(30),t1(180),t2(180),
     +    pmn(50000),pmx(50000),nshell,mshell,mlat(181),nlat
      do 10 i=1,nshell
      if(r0.ge.r1(i).and.r0.le.r2(i)) then
         n=i
         goto 20
      end if
   10 continue
      ib=-1
      return
   20 do 30 i=1,nlat
      if(t0.ge.t1(i).and.t0.le.t2(i)) then
         j1=(n-1)*mshell+mlat(i)+1
         j2=(n-1)*mshell+mlat(i+1)
         goto 40
      end if
   30 continue
      ib=-1
      return
   40 do 50 i=j1,j2
      if(p0.ge.pmn(i).and.p0.le.pmx(i)) then
         ib=i
         return
      end if
   50 continue
      ib=-1
      return
      end

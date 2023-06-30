c find scatterer's position around the tube
c
      subroutine locates(xl,yl,zl,xt,yt,zt,rotm,xs,ys,zs,rs,ts,ps)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      DIMENSION rotm(3,3)

      xs=xl+rotm(1,1)*xt+rotm(1,2)*yt+rotm(1,3)*zt
      ys=yl+rotm(2,1)*xt+rotm(2,2)*yt+rotm(2,3)*zt
      zs=zl+rotm(3,1)*xt+rotm(3,2)*yt+rotm(3,3)*zt
      call xyz2rtp(xs,ys,zs,rs,ts,ps)
      return
      end
   

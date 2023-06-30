      parameter(mlx=7,mly=7,mlz=7,
     ,          mnx=2**(mlx-1)+1,mny=2**(mly-1)+1,mnz=2**(mlz-1)+1,
     ,          mvolu=mnx*mny*mnz,mdata=14000,mpath=240,
     ,          thres=1.d-7)
      dimension jdx(mvolu),locat_g(mvolu),iposition(8)
      double precision lat,lon,latp,lonp,
     ,                 lat0,lon0,lat1,lon1,latc,lonc,lonl,latl,
     ,                 rad0,rad1,radc,radl
      double precision coef_ref(mvolu,20),coef(mvolu,20),coef_LSQR(mvolu)
      double precision x0,y0,z0,dx,dy,dz
      double precision trans(3,3)
      real coef_real(mvolu)
      common/meshb1/x0,y0,z0,dx,dy,dz,nx,ny,nz,nface,nvolu,nlocat_g,locat_g
      common/coeffs/ coef_ref,coef,coef_LSQR
      common/meshb2/trans
      include 'lsqrcom1'
      double precision sample0(3),weight0(3),wtx,wty,wtz,
     ,                 wtw(27,1),value(8),area0,area,wl(3)
      data sample0 /-0.77459667d0,0.d0,0.77459667d0/
c      data wl /(1.d0+sample0(1))*0.5d0,0.5d0,
c     ,        (1.d0+sample0(3))*0.5d0/
      data wl /.112701665d0,0.5d0,
     ,        .887298335d0/
      data weight0/0.55555556d0,0.88888889d0,0.55555556d0/


      parameter(mlevel=5,mbase=20,medge=2**(mlevel-1)+1,
     ,          mface=4**(mlevel-1)*mbase,
     ,          mvert=(medge*(medge+1))/2*mbase,
     ,          mlevelR=6,medgeR=2**(mlevelR-1)+1,
     ,          mvoluf=mface*(medgeR-1),mvolu=mvert*medgeR)
      character*120 title
      integer faceg(3,mface),jdx(mvolu),
     ,        local(6),locatG(mvoluf)
      real*8 gauss(mface,4,2),kernel_value(12,1),weight(3),rthedaph(12,3)
      real*8 image(mvolu,1),coef(mvolu),coef_value(6),
     ,       ttmp,thres,thres1,tmp(3,3)
      parameter(thres=1.d-5)

      common /buildg12/ locatG,gauss,image,coef,faceg
     &       ,nvertg,nface,nedgeR,nvolu,rC,delr

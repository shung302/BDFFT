c determine the numbers of bounce points at surface and core-mantle boundary
c
      subroutine bpt(iraytype,ips,iupdown,nsurf,ncmb,nocb,niob)
      integer iraytype,ips,iupdown,nsurf,ncmb,nocb,niob

      ips=1
      iupdown=1
      nocb=0
      niob=0
      if (iraytype.gt.20) ips=2
      if (iraytype.eq.10.or.iraytype.eq.11.or.iraytype.eq.30
     &   .or.iraytype.eq.31) iupdown=2
      if (iraytype.eq.1.or.iraytype.eq.21) then
c for P or S
         nsurf=0
         ncmb=0
      elseif (iraytype.eq.2.or.iraytype.eq.22) then
c  for PP or SS
         nsurf=1
         ncmb=0
      elseif (iraytype.eq.3.or.iraytype.eq.23) then
c  for PPP or SSS
         nsurf=2
	 ncmb=0
      elseif (iraytype.eq.4.or.iraytype.eq.24) then
c  for PcP or ScS
         nsurf=0
	 ncmb=1
      elseif (iraytype.eq.5.or.iraytype.eq.25) then
c  for PcP2 or ScS2
         nsurf=1
	 ncmb=2
      elseif (iraytype.eq.6.or.iraytype.eq.26) then
c  for PcP3 or ScS3
         nsurf=2
	 ncmb=3
      elseif (iraytype.eq.7.or.iraytype.eq.27) then
c  for PcP4 or ScS4
         nsurf=3
	 ncmb=4
      elseif (iraytype.eq.8.or.iraytype.eq.28) then
c  for PcP5 or ScS5
         nsurf=4
         ncmb=5
      elseif (iraytype.eq.9.or.iraytype.eq.29) then
c  for PcP6 or ScS6
         nsurf=5
         ncmb=6
      elseif (iraytype.eq.10.or.iraytype.eq.30) then
c  for pP or sS (only for depth phases)
         nsurf=1
         ncmb=0
      elseif (iraytype.eq.11.or.iraytype.eq.31) then
c  for pPP or sSS (only for depth phases)
         nsurf=2
         ncmb=0
      elseif (iraytype.eq.12.or.iraytype.eq.32) then
c  for Pd or Sd (only for diffracted phases)
         nsurf=0
	 ncmb=999
      elseif (iraytype.eq.13.or.iraytype.eq.33) then
c  for PS or SP
         nsurf=1
         ncmb=0
      elseif (iraytype.eq.14.or.iraytype.eq.34) then
c  for PcS or ScP
         nsurf=0
         ncmb=1
      elseif (iraytype.eq.15.or.iraytype.eq.35) then
c  for PKS or SKP
         nsurf=0
         ncmb=0
      elseif (iraytype.eq.16.or.iraytype.eq.17.or.iraytype.eq.36) then
c  for PKPab, PKPbc, or SKSac
         nsurf=0
         ncmb=0
      elseif (iraytype.eq.18.or.iraytype.eq.38) then
c  for PKPdf or SKSdf
         nsurf=0
         ncmb=0
      elseif (iraytype.eq.19.or.iraytype.eq.37.or.iraytype.eq.39) then
c  for PKKPdf, SKKSac, or SKKSdf
         nsurf=0
         ncmb=0
         nocb=1
      elseif (iraytype.eq.20.or.iraytype.eq.40) then
c  for PKPcd (PKiKP) and SKiKP
         nsurf=0
         ncmb=0
         nicb=1
      endif

      return
      end

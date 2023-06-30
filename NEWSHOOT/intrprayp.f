c look up the table of ray parameter for starting incident angle
c
      subroutine intrprayp(iraytype,distrs,deps,radsrc,c,
     &  ndi,nds,ditbl,dstbl,rptbl,gai,ai0,ai)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      real*8 rptbl(ndi,nds),ditbl(ndi),dstbl(nds),ai(*)
      COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn
      deg5=5.0d0*rpd
      deg10=10.0d0*rpd
      do j=2,nds
         if (deps.le.dstbl(j)) then
            deps1=dstbl(j-1)
            deps2=dstbl(j)
            jds=j
            go to 10
         endif
      enddo

10    do i=2,ndi
         if (distrs.le.ditbl(i)) then
            distrs1=ditbl(i-1)
            distrs2=ditbl(i)
            idi=i
            go to 20
	 endif
      enddo
20    continue 

c      write(6,*)'idi,jds',idi,jds,deps,distrs,rptbl(idi,jds),distrs1,distrs2,rptbl(idi-1,jds-1),ndi
c      do i=1,ndi
c      write(6,*)ditbl(i),rptbl(i,jds)
c      enddo
c      write(6,*)ditbl(idi),dstbl(jds),rptbl(idi,jds)

      sai1=rptbl(idi-1,jds-1)*c*dpr/radsrc
      sai2=rptbl(idi,jds-1)*c*dpr/radsrc
      sai3=rptbl(idi,jds)*c*dpr/radsrc
      sai4=rptbl(idi-1,jds)*c*dpr/radsrc
      t=(distrs-distrs1)/(distrs2-distrs1)
      u=(deps-deps1)/(deps2-deps1)
      sai=(1-t)*(1-u)*sai1+t*(1-u)*sai2+t*u*sai3+(1-t)*u*sai4
c      write(6,*)idi,jds,rptbl(idi-1,jds-1),rptbl(idi,jds-1),rptbl(idi-1,jds),rptbl(idi,jds)
      if (iraytype.eq.10.or.iraytype.eq.11.or.
     $    iraytype.eq.30.or.iraytype.eq.31) then
      ai0=dasin(sai)
      ai(1)=dasin(sai1)
      ai(3)=dasin(sai3)
      ai(2)=dasin(sai2)
      ai(4)=dasin(sai4)
c      write(6,*)'ai=',(ai(i)*dpr,i=1,4)
      else
      ai0=pi-dasin(sai)
      ai(1)=pi-dasin(sai1)
      ai(3)=pi-dasin(sai3)
      ai(2)=pi-dasin(sai2)
      ai(4)=pi-dasin(sai4)
      endif

      if (iraytype.eq.1.or.iraytype.eq.21) then
         ai(1)=ai(1)-deg5
         ai(2)=ai(2)+deg5
         if (ai(2).gt.gai) ai(2)=gai
      elseif (iraytype.eq.2.or.iraytype.eq.22) then
         ai(1)=ai(1)-deg10
         ai(2)=ai(2)+deg10
         if (ai(2).gt.gai) ai(2)=gai
      elseif (iraytype.eq.4.or.iraytype.eq.24) then
         ai(1)=ai(1)+deg5
         ai(2)=ai(2)-deg5
         if (ai(2).lt.gai) ai(2)=gai
      else
         ai(1)=ai(1)-2.*rpd
         ai(2)=ai(2)+2.*rpd
         if (ai(2).ge.pi) ai(2)=pi-0.00001
      endif
      
      return
      end

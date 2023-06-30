c subroutine to read in rayparameter tables
c
	subroutine rayptbl(iraytype,ndi,nds,ddi,dds,di0,ds0,
     &       irpflag,ditbl,dstbl,rptbl)
        IMPLICIT DOUBLEPRECISION (a-h,o-z)
        data imax /100000/, jmax /1000/
        dimension rptbl(ndi,nds),ditbl(ndi),dstbl(nds)
        INCLUDE 'phases.inc'
        character*80 fnm

        lp=index(pha(iraytype),' ')-1
        fnm='iasp/rayp.'//pha(iraytype)(1:lp)//'.out'
c        write(6,*)fnm
        open(11,file=fnm,status='old',err=400)

        irpflag=1
        do i=1,nds
           dstbl(i)=ds0+dble(i-1)*dds
        enddo
        do i=1,ndi
           ditbl(i)=di0+dble(i-1)*ddi
        enddo

c        write(6,*)d1,d2,ndi,nds
        do 200 ii=1,imax
           read(11,*,end=300)d3,d4
           if (d3.eq.d4) then
              i=int((d3-ds0)/dds)+1
           else
              j=int((d3-di0)/ddi)+1
              rptbl(j,i)=d4
           endif
c           write(6,'(e12.5,2f10.2,2i8)')rptbl(j,i),d3,d4,j,i
200	continue
300     close(11)
        goto 500

400     write(*,*)'ERROR:rayptbl: no iasp/rayp* file found - stop'
        stop

500	return
	end

c make color table from file color.rgb
c
	write(6,*)'input the range xmin to xmax'
	read(5,*) x,y

	open(1,file='/home/shung/bin/svel12.rgb')
	open(2,file='svel12.cpt')
	read(1,*)
	read(1,*)nc
	dc=(y-x)/nc
	do i=1,nc
	   read(1,*)d1,ir1,ig1,ib1,d2,ir2,ig2,ib2
	   if (i.eq.1) then
	   ir0=ir1
	   ig0=ig1
	   ib0=ib1
           endif
	   y1=x+(i-1)*dc
	   y2=y1+dc
	   write(2,99)y1,ir1,ig1,ib1,y2,ir2,ig2,ib2
	enddo
        write(2,199)ir0,ig0,ib0
	write(2,299)ir2,ig2,ib2
	write(2,399)128,128,128
99      format(f12.4,3i8,f12.4,3i8)
199     format('B',8x,3i8)
299     format('F',8x,3i8)
399     format('N',8x,3i8)

	stop
	end


c rotate coordinate in one to another given a rotation matrix
	subroutine rotate(xo,yo,zo,rotm,xn,yn,zn)
	implicit REAL*8 (a-h,o-z)
        dimension rotm(3,3)

	xn=xo*rotm(1,1)+yo*rotm(1,2)+zo*rotm(1,3)
	yn=xo*rotm(2,1)+yo*rotm(2,2)+zo*rotm(2,3)
	zn=xo*rotm(3,1)+yo*rotm(3,2)+zo*rotm(3,3)
c	vlen=sqrt(xn**2+yn**2+zn**2)
c	xn=xn/vlen
c	yn=yn/vlen
c	zn=zn/vlen

	return
	end
	

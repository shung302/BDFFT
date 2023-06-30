c functions to call
c  function to calculate function value at nominator for a given omega
c
        function dfint1(w)
        implicit doubleprecision (a-h,o-z)
        common /bphi/ tdiff,sig  
        common /stf/ alpha,tp
        COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn


        w2=w*w
        w5=w2*w2*w
        dfint1=w5*dexp(-w2/alpha*0.50d0)
     &        *dsin(w*tdiff-sig*pi2)

c	 write(6,*)'dfint1=',w,dd,T1,T2,T,M1,M2,M,dexp(-w2/alpha*0.50d0),dsin(w*(T1+T2-T)-(M1+M2-M)*pi2)
        return
        end

c function to calculate function value at denominator for a given omega
c
        function dfint2(w)

        implicit doubleprecision (a-h,o-z)
        common /stf/ alpha,tp

        w2=w*w
        w4=w2*w2
        dfint2=w4*dexp(-w2/alpha*0.50d0)

        return
        end

c function to calculate function value at denominator for a given omega
c for mixing two closely-arriving phases
c
        function dfint3(w,dit,rsp)

        implicit doubleprecision (a-h,o-z)
        common /stf/ alpha,tp

        w2=w*w
        w4=w2*w2
	fac=1.0d0+rsp*rsp+2.0d0*rsp*dcos(w*dit)
        dfint3=w4*dexp(-w2/alpha*0.50d0)*fac

        return
        end


      FUNCTION MOD1(I,N)
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

c find root for PcPp
	function dfaincs1(x)

	implicit doubleprecision (a-h,o-z)
	common /facoeffs/ c1,c2,c3,b2,bc2,a

	xs=x*x
c	write(6,*)'b2, bc=',c1,c2,c3,b2,bc2,a,x
c	write(6,*)'sqrt=',1.0d0-bc2*xs,1.0d0-b2*xs
	dfaincs1=c1+a*x-c2*dsqrt(1.0d0-bc2*xs)-dsqrt(1.0d0-b2*xs)

	return
	end

c find root for PcPn
	function dfaincs2(x)

	implicit doubleprecision (a-h,o-z)
	common /facoeffs/ c1,c2,c3,b2,bc2,a

	xs=x*x
	dfaincs2=c3+a*x-c2*dsqrt(1.0d0-bc2*xs)+dsqrt(1.0d0-b2*xs)
	return
	end

c find root for pP when zd != 0
c faincs3=sin(x)
cc        function dfaincs3(x)

c        implicit doubleprecision (a-h,o-z)
c        common /facoeffs/ c1,c2,c3,b2,bc2,a
c        common /facoeffsd/ b12,hhd2,bn
c
c	x2=x*x
cc       write(6,*)'b12 and hhd2 and a in faincs3=',b12,hhd2,a
c        dfaincs3=1.0d0-(a*x+dsqrt(1.0d0-x2)-2.0d0*bn*dsqrt(1.0d0-hhd2*x2))**2
c     &          -b12*x2
c        return
c        end


c find root for PP when zd !=0
        function dfaincs4(x)

        implicit doubleprecision (a-h,o-z)
        common /facoeffs/ c1,c2,c3,b2,bc2,a
        common /facoeffsd/ b12,hhd2,bn

        x2=x*x
c       write(6,*)'b12 and hhd2 and a in faincs3=',b12,hhd2,a
        dfaincs4=1.0d0-(a*x-dsqrt(1.0d0-x2)
     &          -2.0d0*bn*dsqrt(1.0d0-hhd2*x2))**2
     &          -b12*x2
        return
        end


c find root for pP when zd != 0
        function dfaincs3(x)

        implicit doubleprecision (a-h,o-z)
        common /facoeffs/ c1,c2,c3,b2,bc2,a
        common /facoeffsd/ b12,hhd2,bn

	sa=dsin(x)
	ca=dcos(x)
	sa2=sa*sa
c	write(6,*)'b12 and hhd2 and a in faincs3=',b12,hhd2,a
        dfaincs3=1.0d0-(a*sa+ca-2.0d0*bn*dsqrt(1.0d0-hhd2*sa2))**2
     &          -b12*sa2
        return
        end
c
c
cc find root for PP when zd !=0
c        function dfaincs4(x)
c
c        implicit doubleprecision (a-h,o-z)
c        common /facoeffs/ c1,c2,c3,b2,bc2,a
c        common /facoeffsd/ b12,hhd2,bn
c
c        sa=dsin(x)
c        ca=dcos(x)
c        sa2=sa*sa
c        dfaincs4=1.0d0-(a*sa-ca-2.0d0*bn*dsqrt(1.0d0-hhd2*sa2))**2
c     &          -b12*sa2
c        return
c        end

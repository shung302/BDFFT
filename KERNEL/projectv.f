c project a vector onto b vector and return with
c the projection vector's length as bb
c and residual vector as a
      subroutine projectv(a,b,bb)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      dimension a(3),b(3)

      ab=0.0d0
      bb=0.0d0
      do 23192 i=1,3
      ab=ab+a(i)*b(i)
      bb=bb+b(i)*b(i)
23192 continue
      if(bb.ne.0.0d0) go to 23194
      bb=1.0d0
23194 continue
      bs=dsqrt(bb)
      aa=0.0d0
      bb=ab/bs
c      ab=ab/bb
c      bb=dabs(ab/bs)
c      do 23196 i=1,3
c      a(i)=a(i)-ab*b(i)
c      aa=aa+a(i)*a(i)
c23196 continue
      return
      end
c      if(aa.ne.0.0d0) go to 23198
c      a(1)=b(2)
c      a(2)=-b(1)
c      aa=a(1)*a(1)+a(2)*a(2)
c23198 continue
c      aa=dsqrt(aa)

c      if(aa.ne.0.0d0) go to 23200
c      aa=1.0d0
c23200 continue
c      do 23202 i=1,3
c      a(i)=a(i)/aa
c23202 continue
c      return
c      end


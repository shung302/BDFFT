c calculat product of two vectors
c return with their lengths, product and cosine of angle 
c betwen these two vectors
      subroutine productv(qv,pk,qvlg,pklg,prd,cs,acs)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      DIMENSIOn pk(3),qv(3)

      prd=0.0d0
      qvlg=0.0d0
      pklg=0.0d0
      do i=1,3
         prd=prd+qv(i)*pk(i)
         qvlg=qvlg+qv(i)*qv(i)
         pklg=pklg+pk(i)*pk(i)
      enddo
      if (prd.ne.0.0d0) go to 10
      cs=0.0d0
      acs=0.0d0
      prd=0.0d0
      return
10    if (qvlg.ne.0.0d0) go to 20
      cs=0.0d0
      acs=0.0d0
      prd=0.0d0
      return
20    qvlg=dsqrt(qvlg)
      cs=prd/(qvlg*pklg)
      acs=dabs(cs)

      return
      end

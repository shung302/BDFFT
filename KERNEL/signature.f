c calculate sig(hf+hb)

      subroutine signature(iphase,h11,h22,nsig,sig)
      implicit doubleprecision (a-h,o-z)
      integer nsig

c      if (iphase.ne.2.and.iphase.ne.3.and.
c     &   iphase.ne.22.and.iphase.ne.23) then
c       sig=0.0d0
c       return
c      endif
      if (iphase.eq.1.or.iphase.eq.21) then
         nsig=0
         sig=0.0d0
         return
      endif

      npos=0
      nneg=0

      if (h11.lt.0.0d0) nneg=nneg+1
      if (h11.gt.0.0d0) npos=npos+1
      if (h22.lt.0.0d0) nneg=nneg+1
      if (h22.gt.0.0d0) npos=npos+1
      nsig=(npos-nneg-2)/2
      sig=dble(nsig)
      return
      end 

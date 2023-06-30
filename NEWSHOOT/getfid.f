c
      subroutine getfid(ivin,ipha,irs,idep,rs,dep,nid,fid)
      IMPLICIT DOUBLEPRECISION (a-h,o-z)
      INTEGER ipha,irs,idep,nid
      CHARACTER*80 fid
      INCLUDE 'models.inc'
      INCLUDE 'phases.inc'   

      im=ivin+5
      drs=rs-dble(irs)
c	write(6,*)'here',irs,ipha,mdl(im),idep,rs,deps,drs
      if (idep.lt.10) then

         if (drs.ne.0.) then
            if (rs.lt.10.) then
               write(fid,'(a4,''.r'',f3.1,''.d'',i1,''.'',a6)')
     &         mdl(im),rs,idep,pha(ipha)
            elseif (rs.ge.10..and.rs.lt.100.) then
               write(fid,'(a4,''.r'',f4.1,''.d'',i1,''.'',a6)')
     &         mdl(im),rs,idep,pha(ipha)
            elseif (rs.ge.100..and.rs.lt.1000.) then
               write(fid,'(a4,''.r'',f5.1,''.d'',i1,''.'',a6)')
     &         mdl(im),rs,idep,pha(ipha)
            endif
         else
            if (irs.lt.10) then
               write(fid,'(a4,''.r'',i1,''.d'',i1,''.'',a6)')
     &         mdl(im),irs,idep,pha(ipha)
            elseif (irs.ge.10.and.irs.lt.100) then
               write(fid,'(a4,''.r'',i2,''.d'',i1,''.'',a6)')
     &         mdl(im),irs,idep,pha(ipha)
            elseif (irs.ge.100.and.irs.lt.1000) then
               write(fid,'(a4,''.r'',i3,''.d'',i1,''.'',a6)')
     &         mdl(im),irs,idep,pha(ipha)
            endif
         endif

      elseif (idep.ge.10.and.idep.lt.100) then

         if (drs.ne.0.) then
            if (rs.lt.10.) then
               write(fid,'(a4,''.r'',f3.1,''.d'',i2,''.'',a6)')
     &         mdl(im),rs,idep,pha(ipha)
            elseif (rs.ge.10..and.rs.lt.100.) then
               write(fid,'(a4,''.r'',f4.1,''.d'',i2,''.'',a6)')
     &          mdl(im),rs,idep,pha(ipha)
            elseif (rs.ge.100..and.rs.lt.1000.) then
               write(fid,'(a4,''.r'',f5.1,''.d'',i2,''.'',a6)')
     &          mdl(im),rs,idep,pha(ipha)
            endif
         else
            if (irs.lt.10) then
               write(fid,'(a4,''.r'',i1,''.d'',i2,''.'',a6)')
     &          mdl(im),irs,idep,pha(ipha)
            elseif (irs.ge.10.and.irs.lt.100) then
               write(fid,'(a4,''.r'',i2,''.d'',i2,''.'',a6)')
     &          mdl(im),irs,idep,pha(ipha)
            elseif (irs.ge.100.and.irs.lt.1000) then
               write(fid,'(a4,''.r'',i3,''.d'',i2,''.'',a6)')
     &          mdl(im),irs,idep,pha(ipha)
            endif
         endif

      elseif (idep.ge.100.and.idep.lt.1000) then

         if (drs.ne.0.) then
            if (rs.lt.10.) then
               write(fid,'(a4,''.r'',f3.1,''.d'',i3,''.'',a6)')
     &          mdl(im),rs,idep,pha(ipha)
            elseif (rs.ge.10..and.rs.lt.100.) then
               write(fid,'(a4,''.r'',f4.1,''.d'',i3,''.'',a6)')
     &          mdl(im),rs,idep,pha(ipha)
            elseif (rs.ge.100..and.rs.lt.1000.) then
               write(fid,'(a4,''.r'',f5.1,''.d'',i3,''.'',a6)')
     &          mdl(im),rs,idep,pha(ipha)
            endif
         else
            if (irs.lt.10) then
               write(fid,'(a4,''.r'',i1,''.d'',i3,''.'',a6)')
     &          mdl(im),irs,idep,pha(ipha)
            elseif (irs.ge.10.and.irs.lt.100) then
               write(fid,'(a4,''.r'',i2,''.d'',i3,''.'',a6)')
     &          mdl(im),irs,idep,pha(ipha)
            elseif (irs.ge.100.and.irs.lt.1000) then
               write(fid,'(a4,''.r'',i3,''.d'',i3,''.'',a6)')
     &          mdl(im),irs,idep,pha(ipha)
            endif
         endif

      endif                              

      nid=index(fid,' ')-1

      return
      end

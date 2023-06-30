c-----------------------------------------------------------------------
c
      subroutine eulerf(sth,sph,rth,rph,alpha,beta,gamma,
     & del,pth,pph,g2e,e2g)
c   Finds the Euler angles alpha,beta,gamma that rotate the coordinate axes
c   so that the z-axis is at the pole of the source-receiver great circle
c   (s x r), and the x-axis is at the source. See Edmonds' Angular Momentum
c   in Quantum Mechanics, page 7 for the angle conventions.
c     input: sth,sph = source coordinates in radians
c            rth,rph = receiver coordinates in radians
c     output: alpha,beta,gamma = euler angles in radians which rotate the
c             original coordinate system to the one with the source-receiver
c             great circle on the equator, the source at (PI/2,0). The minor
c             arc to the receiver is in the positive phi direction.
c            del = source-receiver separation in radians.
c            pth,pph = source-receiver great circle pole location.
c g2e: rotation matrix from 
      implicit real*8(a-h,o-z)
      dimension e2g(3,3),g2e(3,3)
c      data pi/3.14159265358979d0/
       pi=dasin(1.0d0)*2.0d0
c  Get cartesian coordinates for source and receiver
      call cart(sth,sph,sx,sy,sz)
      call cart(rth,rph,rx,ry,rz)
      del = dacos(sx*rx + sy*ry + sz*rz)
c cross --> (px,py,pz)=(sx,sy,sz)x(rx,ry,rz)
      call cross(sx,sy,sz,rx,ry,rz,px,py,pz)
      pth = datan2(dsqrt(px*px+py*py),pz)
      if(px.eq.0. .and. py.eq.0.) then
c        special case of pole at z or -z
         pph = 0.
      else
         pph = datan2(py,px)
      endif
      alpha = pph
      beta = pth
c  the x'' axis (call it t) is at pth+pi/2,pph
      ttheta = pth + pi/2.
      call cart(ttheta,pph,tx,ty,tz)
c  the third Euler angle, gamma, rotates x'' to the source s.
      gamma = dacos(sx*tx + sy*ty + sz*tz)
c  form q = x'' x s to check the sign of gamma (q/|q| = +-p/|p|)
      call cross(tx,ty,tz,sx,sy,sz,qx,qy,qz)
      sgn = px*qx + py*qy + pz*qz
      if(sgn .lt. 0.) gamma = -gamma
      ca=dcos(alpha)
      sa=dsin(alpha)
      cb=dcos(beta)
      sb=dsin(beta)
      cg=dcos(gamma)
      sg=dsin(gamma)

c g2e: rotation matrix from x' to x (great circle plane to equator plane)
      g2e(1,1)=ca*cb*cg-sa*sg
      g2e(2,1)=-ca*cb*sg-sa*cg
      g2e(3,1)=ca*sb
      g2e(1,2)=sa*cb*cg+ca*sg
      g2e(2,2)=-sa*cb*sg+ca*cg
      g2e(3,2)=sa*sb
      g2e(1,3)=-sb*cg
      g2e(2,3)=sb*sg
      g2e(3,3)=cb

c g2e': rotation matrix from x to x' (equator plane to great circle plane)
      call transpose(3,3,g2e,e2g)
c      e2g(1,1)=g2e(1,1)
c      e2g(2,1)=sa*cb*cg+ca*sg
c      e2g(3,1)=-sb*cg
c      e2g(1,2)=-ca*cb*sg-sa*cg
c      e2g(2,2)=g2e(2,2)
c      e2g(3,2)=sa*sb
c      e2g(1,3)=g2e(1,3)
c      e2g(2,3)=g2e(2,3)
c      e2g(3,3)=g2e(3,3)
      
      return
      end
c get transpose of the matrix
      subroutine transpose(n,m,a,b)
      implicit real*8 (a-h,o-z)
      dimension a(n,m),b(m,n)

      do i=1,n
      do j=1,m
         b(i,j)=a(j,i)
      enddo
      enddo
      return
      end
     

c
c--------------------------------------------------------------------------
c
      subroutine cart(thet,phi,x,y,z)
      implicit real*8(a-h,o-z)
      s=dsin(thet)
      x=s*dcos(phi)
      y=s*dsin(phi)
      z=dcos(thet)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine cross(sx,sy,sz,rx,ry,rz,px,py,pz)
      implicit real*8(a-h,o-z)
      px = sy*rz - sz*ry
      py = sz*rx - sx*rz
      pz = sx*ry - sy*rx
      return
      end
c
c----------------------------------------------------------------------



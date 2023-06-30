c functions to call
c  function to calculate function value at nominator for a given omega
c
        function dfint(w,a2)
        implicit doubleprecision (a-h,o-z)
        real*8 dfint
        common /bphi/ tdiff,sig  
        COMMON /consts/ pi,pi2,twopi,rpd,dpr,rsurf,rcmb,riob,rn

        dfint=w*a2*dsin(w*tdiff-sig*pi2)

        return
        end


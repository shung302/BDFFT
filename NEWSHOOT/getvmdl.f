c subroutine to determine which global models to be used for
c integration of kernels
c
      subroutine getvmdl(mdl,frad,fcol,flong,dv)
      implicit real*8(a-h,o-z)

      goto (10,20,30,40,45,50) mdl
10    call getv_smdl(frad,fcol,flong,dv) 
      go to 50
20    call getv_smdl(frad,fcol,flong,dv)
      go to 50
30    call getv_cmdl(frad,fcol,flong,dv)
      go to 50
40    call getv_smdl(frad,fcol,flong,dv)
      go to 50
45    call getv_smdl(frad,fcol,flong,dv)
50    return
      end

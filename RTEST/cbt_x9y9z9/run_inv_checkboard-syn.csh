#!/bin/csh
set TBTDIR = ../..
set BDIR = $TBTDIR/BUILDG
set PD = $TBTDIR/INV
cp $PD/solver_noSE_wgtnoise_perc .
set PD = `pwd`
cp $BDIR/mesh.config .
# build G and synthetic travel time data
# make sure mesh.config in the running directory
set damping = ( 100 200 500 )
# create checkerboard anomaly
# IASP91 model
set ivelin = 2
set mdl = ak1c.tc
# P (1) or S (2) wave
set ips = 2
# velocity perturbations
set dv = -0.04
# 1-sigma error for travel time measurement
set ifdemean = 0
# checkboard width
set log = log.kern_syn
date>$log
foreach nxc ( 9 )
set nyc = $nxc
set nzc = $nxc
#set nzc = `echo $nxc | awk ' { print $1*2 } '`
set pattern = cbt_x"$nxc"y"$nyc"z"$nzc"

#goto ray
# inversion kernel G matrix
foreach PH ( S.hm.hic )
set FG = ./G_kernel_d_syn_$PH.$mdl.$pattern
foreach noise ( 0.05 0 0.1 )
goto multis
simple:
# for simple DLS
echo DLS solutions!
foreach damp ( `echo $damping` )
$PD/solver_noSE_wgtnoise_perc>>$log<<!
1
$damp
$FG
$noise
!
mv try.xyz try.simple_kern_syn_$damp.$PH.$mdl.$pattern.n${noise}
mv try.sum sum.simple_kern_syn_$damp.$PH.$mdl.$pattern.n${noise}
mv results results.simple_kern_syn_$damp.$PH.$mdl.$pattern.n${noise}
date >> $log
end #foreach damp

quelling:
echo Quelling solutions!
foreach sigmaz ( 2 1 )
set sigmah = $sigmaz
set sigma = "$sigmaz"-"$sigmah"
foreach damp ( `echo $damping` )
$PD/solver_noSE_wgtnoise_perc<<!>>$log
3
$damp
$FG
$sigmah $sigmaz
$noise
!
mv try.xyz try.quell_kern_syn_"$damp"-"$sigma".$PH.$mdl.$pattern.n${noise}
mv try.sum sum.quell_kern_syn_"$damp"-"$sigma".$PH.$mdl.$pattern.n${noise}
mv results results.quell_syn_"$damp"-"$sigma".$PH.$mdl.$pattern.n${noise}
date >>$log
end #foreach sigmaz
end #foreach damp

multis:
echo Multiscale solutions!
foreach damp ( `echo $damping` )
$PD/solver_noSE_wgtnoise_perc<<!>>& $log
2
$damp
$FG
$noise
!
mv try.xyz try.multis_kern_syn_"$damp".$PH.$mdl.$pattern.n${noise}
mv try.sum sum.multis_kern_syn_"$damp".$PH.$mdl.$pattern.n${noise}
mv results results.multis_kern_syn_"$damp".$PH.$mdl.$pattern.n${noise}
end #foreach damp
goto skip1

hybrid:
echo Hybrid solutions!
foreach sigmaz ( 2 1 )
foreach damp ( `echo $damping` )
$PD/solver_noSE_wgtnoise_perc<<!>>& $log
4
$damp
$FG
$sigmaz
$noise
!
mv try.xyz try.hybrid_kern_syn_$damp-$sigmaz.$PH.$mdl.$pattern.n${noise}
mv try.sum sum.hybrid_kern_syn_$damp-$sigmaz.$PH.$mdl.$pattern.n${noise}
mv results results.hybrid_kern_syn_$damp-$sigmaz.$PH.$mdl.$pattern.n${noise}

end #foreach sigmaz
end #foreach damp

skip1:

end #foreach noise
end #PH
exit
#end #nxc
ray:
set log = log.ray_syn
date>$log
foreach PH ( S.m.hic )
set FG = ./G_ray_d_syn_$PH.$mdl.$pattern
foreach noise ( 0.05 0 0.1 )
#goto rmultis

rsimple:
echo simple DLS solutions!
foreach damp ( `echo $damping` )
$PD/solver_noSE_wgtnoise_perc>>$log<<!
1
$damp
$FG
$noise
!
mv try.xyz try.simple_ray_syn_$damp.$PH.$mdl.$pattern.n${noise}
mv try.sum sum.simple_ray_syn_$damp.$PH.$mdl.$pattern.n${noise}
mv results results.simple_ray_syn_$damp.$PH.$mdl.$pattern.n${noise}
date >>$log
end #foreach damp

rquelling:
echo Quelling solution!
foreach sigmaz ( 2 1 )
set sigmah = $sigmaz
set sigma = "$sigmaz"-"$sigmah"
foreach damp ( `echo $damping` )
$PD/solver_noSE_wgtnoise_perc<<!>>& $log
3
$damp
$FG
$sigmah $sigmaz
$noise
!
set sigma = "$sigmaz"-"$sigmah"
mv try.xyz try.quell_ray_syn_"$damp"-"$sigma".$PH.$mdl.$pattern.n${noise}
mv try.sum sum.quell_ray_syn_"$damp"-"$sigma".$PH.$mdl.$pattern.n${noise}
mv results results.quell_ray_syn_"$damp"-"$sigma".$PH.$mdl.$pattern.n${noise}
date >>$log
end #foreach damp
end #foreach sigma

rmultis:
echo Multiscale solutions!
foreach damp ( `echo $damping` )
$PD/solver_noSE_wgtnoise_perc<<!>>& log.solver
2
$damp
$FG
$noise
!
mv try.xyz try.multis_ray_syn_"$damp".$PH.$mdl.$pattern.n${noise}
mv try.sum sum.multis_ray_syn_"$damp".$PH.$mdl.$pattern.n${noise}
mv results results.multis_ray_syn_"$damp".$PH.$mdl.$pattern.n${noise}
end #foreach damp

#goto skip2
rhybrid:
echo Hybrid solution!
foreach damp ( `echo $damping` )
foreach sigmaz ( 2 1 )
$PD/solver_noSE_wgtnoise_perc<<!>>& $log
4
$damp
$FG
$sigmaz
$noise
!
mv try.xyz try.hybrid_ray_syn_$damp-$sigmaz.$PH.$mdl.$pattern.n${noise}
mv try.sum sum.hybrid_ray_syn_$damp-$sigmaz.$PH.$mdl.$pattern.n${noise}
mv results results.hybrid_ray_syn_$damp-$sigmaz.$PH.$mdl.$pattern.n${noise}
date >>$log
end #foreach sigmaz
end #foreach damp
skip2:

end #foreach noise
end #PH

end # of foreach nxc

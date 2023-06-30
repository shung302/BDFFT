#!/bin/csh
# inversion using ray G matrix with station correction terms
#
# parameters:
# mode: 1 = simple damping; 2 = multiscale transformation;
#	3 = convolution quelling; 4 = vertically quelling & laterally multiscale
# damping factor
# G matrix
# number of stations and weight
set TBTDIR = ../../..
set PD = .
set GDIR = ../../../KERNEL

set mdl = ak1c.tc
set area = hic
if ( ! ( -d WGT ) ) then
mkdir WGT
endif
\cp ../solver_WSE WGT
cd WGT
date > log.solver

foreach phase ( S )
foreach fq ( m )
set F = $phase.$fq.$area
set STFILE = $TBTDIR/$F.stations_used
cp $STFILE stations_used
set EVFILE = $TBTDIR/$F.events_used
head -1 $EVFILE >  events_used
more +2 $EVFILE | awk '{print $8, $9 } ' >> events_used

foreach sigmaz ( 1 2 )
foreach ws ( 1 10 )
set we = $ws
#foreach damp ( 20000 15000 10000 30000 50000 8000 5000 3000 1000 10000 100000 )
foreach damp ( 50 100 200 300 500 1000 2000 3000 4000 5000 7000 10000 20000 50000 )
#echo "damp = $damp , output = try.hybrid_kern_se"$ws"-"$we"_$damp.$F.${mdl}"
echo "Inversioning....... Damping:$damp"
set weight = "$ws"-"$we"
set sigma = "$sigmaz"
if ( -f try.hybrid_kern_se"$weight"_"$damp"-"$sigma".$F.$mdl ) continue
$PD/solver_WSE<<!>>& log.solver
4
$damp
$GDIR/Gw_kernel_d_se_$F.$mdl
$ws
$we
$sigmaz
!
mv try.xyz try.hybrid_kern_se"$weight"_"$damp"-"$sigma".$F.$mdl
mv try.sum sum.hybrid_kern_se"$weight"_"$damp"-"$sigma".$F.$mdl
mv results results.hybrid_kern_se"$weight"_"$damp"-"$sigma".$F.$mdl
mv station_correction stacor.hybrid_kern_se"$weight"_"$damp"-"$sigma".$F.$mdl
mv event_correction evtcor.hybrid_kern_se"$weight"_"$damp"-"$sigma".$F.$mdl
end #foreach damp
end #foreach ws
end #foreach sigmaz

end #foreach fq
end #foreach phase

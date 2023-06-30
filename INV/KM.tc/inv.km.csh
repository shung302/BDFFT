#!/bin/csh
# inversion using kern G matrix with station correction terms
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
foreach fq ( hm h )
set F = $phase.$fq.$area
set STFILE = $TBTDIR/$F.stations_used
cp $STFILE stations_used
set EVFILE = $TBTDIR/$F.events_used
head -1 $EVFILE >  events_used
more +2 $EVFILE | awk '{print $8, $9 } ' >> events_used

foreach ws ( 1 10 )
set we = $ws
#foreach damp ( 8000 10000 12000 13000 14000 11000 7000 5000 2000 )
foreach damp ( 50 100 200 300 500 1000 2000 3000 4000 5000 7000 10000 20000 50000 )
#echo "damp = $damp , output = try.simple_kern_se"$ws-$we"_$damp.$F.$mdl"
echo "Inversioning....... Damping:$damp"
set weight = "$ws-$we"
if ( -f try.multis_kern_se"$weight"_"$damp".$F.$mdl ) continue
$PD/solver_WSE<<!>>& log.solver
2
$damp
$GDIR/Gw_kernel_d_se_$F.$mdl
$ws
$we
!
mv try.xyz try.multis_kern_se"$weight"_"$damp".$F.$mdl
mv try.sum sum.multis_kern_se"$weight"_"$damp".$F.$mdl
mv results results.multis_kern_se"$weight"_"$damp".$F.$mdl
mv station_correction stacor.multis_kern_se"$weight"_"$damp".$F.$mdl
mv event_correction evtcor.multis_kern_se"$weight"_"$damp".$F.$mdl
end #foreach damp
end #foreach ws

end #foreach fq
end #foreach phase

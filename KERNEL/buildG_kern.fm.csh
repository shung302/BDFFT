#!/bin/csh
# input file: mesh.config 
set BDIR = ../BUILDG
set area = hic
# add phase and frequencies together
foreach mdl ( ak1c )
# add h and m and l frequencies together
goto addsta
foreach phase ( S )
echo "add all the frequencies together for $phase"
# weighted
#$BDIR/addG_d <<!
#../$phase.hm.$area.datsta
#Gw_kernel_d_$phase.h.$area.$mdl
#Gw_kernel_d_$phase.m.$area.$mdl
#none
#!
#mv G_d_add Gw_kernel_d_$phase.hm.$area.$mdl
## unweighted
#$BDIR/addG_d <<!
#../$phase.hm.$area.datsta
#G_kernel_d_$phase.h.$area.$mdl
#G_kernel_d_$phase.m.$area.$mdl
#none
#!
#mv G_d_add G_kernel_d_$phase.hm.$area.$mdl
# weighted G with topo-corrected d
$BDIR/addG_d <<!
../$phase.hm.$area.datsta
Gw_kernel_d_$phase.h.$area.$mdl.tc
Gw_kernel_d_$phase.m.$area.$mdl.tc
none
!
mv G_d_add Gw_kernel_d_$phase.hm.$area.$mdl.tc
# unweighted G with topo-corrected d
$BDIR/addG_d <<!
../$phase.hm.$area.datsta
G_kernel_d_$phase.h.$area.$mdl.tc
G_kernel_d_$phase.m.$area.$mdl.tc
none
!
mv G_d_add G_kernel_d_$phase.hm.$area.$mdl.tc
end

######################################################
# buildG with station and event corrections terms
# add station and event terms in G_kernel_d or Gw_kernel_d
######################################################
addsta:
foreach phase ( S )
foreach fq ( h m )
set PH = $phase.$fq.$area
# weighted
#$BDIR/addG_sta_evtcor<<!
#../$phase.$fq.$area.stacor
#../$phase.$fq.$area.datsta
#Gw_kernel_d_$phase.$fq.$area.$mdl
#!
##unweighted
#$BDIR/addG_sta_evtcor<<!
#../$phase.$fq.$area.stacor
#../$phase.$fq.$area.datsta
#G_kernel_d_$phase.$fq.$area.$mdl
#!
# weighted with tc
$BDIR/addG_sta_evtcor<<!
../$phase.$fq.$area.stacor
../$phase.$fq.$area.datsta
Gw_kernel_d_$phase.$fq.$area.$mdl.tc
!
# unweighted
$BDIR/addG_sta_evtcor<<!
../$phase.$fq.$area.stacor
../$phase.$fq.$area.datsta
G_kernel_d_$phase.$fq.$area.$mdl.tc
!
end
end

end # end of mdl

#!/bin/csh
# set your path if necessary
set BDIR = .
set DDIR = ..
set area = hic
set ccmdl = tc
#goto addsta
foreach mdl ( ak1c )
# input file: mesh.config 
# buildG 
foreach phase ( S )
foreach fq ( m )
###########################################
# build G matrix with no weight
###########################################
set PH = $phase.$fq.$area
echo ../NEWSHOOT/shootrays.$PH.$mdl
# without correction
$BDIR/buildG<<!
../NEWSHOOT/shootrays.$PH.$mdl
$DDIR/$PH.datsta
0
!
mv G_ray_d G_ray_d_$PH.$mdl
mv Gw_ray_d Gw_ray_d_$PH.$mdl
mv Kernel.pos Kw_pos_$PH.$mdl
mv Ray.summary Rayw.summary.$PH.$mdl

corr_tc:
# build ray-based G matrix and apply topography correction on measured tt
set PH = $phase.$fq.$area
echo ../NEWSHOOT/shootrays.$PH.$mdl
# with correction
$BDIR/buildG_tc<<!
../NEWSHOOT/shootrays.$PH.$mdl
$DDIR/$PH.datsta
0
$DDIR/stations_used
!
mv G_ray_d G_ray_d_$PH.$mdl.tc
mv Gw_ray_d Gw_ray_d_$PH.$mdl.tc
mv Kernel.pos Kw_pos_$PH.$mdl.tc
mv Ray.summary Rayw.summary.$PH.$mdl.tc
end # end of foreach fq
end # end of foreach phase
end # end of foreach mdl

######################################################
# buildG with station and event corrections terms
# add station and event terms in G_ray_d or Gw_ray_d
######################################################
addsta:
foreach mdl ( ak1c )
foreach phase ( S )
foreach fq ( m )
set PH = $phase.$fq.$area
$BDIR/addG_sta_evtcor<<!
$DDIR/$PH.stacor
$DDIR/$PH.datsta
Gw_ray_d_$PH.$mdl
!
$BDIR/addG_sta_evtcor<<!
$DDIR/$PH.stacor
$DDIR/$PH.datsta
G_ray_d_$PH.$mdl
!
$BDIR/addG_sta_evtcor<<!
$DDIR/$PH.stacor
$DDIR/$PH.datsta
Gw_ray_d_$PH.$mdl.$ccmdl
!
$BDIR/addG_sta_evtcor<<!
$DDIR/$PH.stacor
$DDIR/$PH.datsta
G_ray_d_$PH.$mdl.$ccmdl
!
end
end
end  # end of foreach mdl

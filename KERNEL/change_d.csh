#!/bin/csh
set BPD = ../BUILDG
set OUTDIR = .
set INDIR = ../BUILDG
set area = hic
if ( ! ( -f mesh.config ) ) then
cp $INDIR/mesh.config .
endif
foreach refmdl ( ak1c )
foreach cc ( tc )
foreach phase ( S )
foreach fq ( h m )
foreach gm ( Gw G )
# change d in a given G matrix, the one to be changed
set GIN1 = ${gm}_kernel_d_${phase}.$fq.$area.$refmdl
# G matrix that provides new d
set GIN2 = ${gm}_ray_d_${phase}.$fq.$area.$refmdl.$cc
# output G
set GOUT = ${gm}_kernel_d_${phase}.$fq.$area.$refmdl.$cc
$BPD/changeG_d<<!
$OUTDIR/$GIN1
$INDIR/$GIN2
$OUTDIR/$GOUT
!
end
end
end
end
end

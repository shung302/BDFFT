set TBTDIR = ..
set KDIR = $TBTDIR/KERNEL
set BDIR = $TBTDIR/BUILDG
set RDIR = $TBTDIR/BUILDG
set BPD =  ../BUILDG
set area = hic
set ccmdl = ak1c.tc
cp $BDIR/mesh.config .
#cp $KDIR/ak135-c.model .
# build G and synthetic travel time data
# make sure mesh.config in the running directory
# create checkerboard anomaly
# ak135-continental model (ivelin = 2 )
# iasp model
set ivelin = 2
# P (1) or S (2) wave
set ips = 2
if ( $ips == 1 ) then
set phase = P
else
set phase = S
endif
# velocity perturbations
set dv = -0.04
# 1-sigma error for travel time measurement
set tstd = 0.1
set ifdemean = 0
# checkboard width
#set nxc = 2
#set nyc = 2
#set nzc = 2

foreach nxc ( 9 )
set nyc = $nxc
set nzc = $nxc
#set nzc = `echo $nxc | awk ' { print $1*2 } '`
# checkerboard starting from the top
set pattern = cgs_x"$nxc"y"$nyc"z"$nzc"
if ( ! ( -d $pattern ) ) mkdir $pattern
foreach fq ( hm ) 
#goto rayg
# create synthetic kernl-based tt
set PH = $phase.$fq.$area.$ccmdl
set FG = $KDIR/G_kernel_d_$PH
$BPD/synthetic_ray_checkz_gauss <<!
$ivelin
$ips
$FG
$tstd
$dv
$ifdemean
$nxc $nyc $nzc
!
mv try.slowness.syn  try.slowness.syn.$pattern
mv try.velo.syn  try.velo.syn.$pattern
mv $KDIR/G_kernel_d_syn_$PH $pattern/G_kernel_d_syn_$PH.$pattern
end	# end foreach fq

rayg:
foreach fq ( m )
set PH = $phase.$fq.$area.$ccmdl
set FG = $RDIR/G_ray_d_$PH
echo $FG
$BPD/synthetic_ray_checkz_gauss <<!
$ivelin
$ips
$FG
$tstd
$dv
$ifdemean
$nxc $nyc $nzc
!
mv try.slowness.syn  try.slowness.syn.$pattern
rm -f try.velo.syn
#mv try.velo.syn  try.velo.syn.$pattern
mv $RDIR/G_ray_d_syn_$PH $pattern/G_ray_d_syn_$PH.$pattern
end # end of foreach fq

end # of foreach nxc

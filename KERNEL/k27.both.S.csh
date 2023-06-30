#!/bin/csh
set area = hic
set band = ( h m )
set f1 = ( 0.1 0.03 )
set f2 = ( 0.5 0.1 )

# ivelin =  1  for PREM
# ivelin = -4 for IASP91
# ivelin = 2 for ak135-continent model
# ivelin = 3 for ak135 global average model
foreach ivelin ( 2 )
if ( $ivelin == 1 ) then
set mdl = prem
else if ( $ivelin == 3 ) then
set mdl = ak1o
else if ( $ivelin == 2 ) then
set mdl = ak1c
else
set mdl = iasp
endif
echo $mdl

# the directory with built meshes
set TBTDIR = ..
set PD = .

# the directory include *datsta files and waveform spectra files
set TBTEVT = ..

cp $TBTDIR/BUILDG/mesh.config .

date > run.out
@ i = 0

foreach phase ( 21 )
if ( $phase == 21 ) then
set pha = S
endif
if ( $phase == 24 ) then
set pha = ScS
endif
if ( $phase == 36 ) then
set pha = SKSac
endif

foreach fq ( h m )
@ i ++
set flow = $f1[$i]
set fhigh = $f2[$i]

echo $fq $flow - $fhigh Hz $phase $pha 
echo $TBTEVT/$pha.$fq.$area.datsta

# direct P wave
# with topo correction
echo "Now generating G matrix of finite-frequency kernel for $pha phase, $fq band"
if ( -f G_kernel_d_$pha.$fq.$area.$mdl.tc ) continue

$PD/kernbg27_both << ! > run.out
1
$ivelin
$phase
10. 10.
10., 5., .5
1000
4 1
$flow $fhigh
$TBTDIR/waveform_stf/$pha.stack.$fq.t.am.ascii
1
$TBTDIR/BUILDG/Kw_pos_$pha.$fq.$area.$mdl.tc
$TBTDIR/$pha.$fq.$area.datsta
!
mv kobs.$mdl.G_kernel_d		G_kernel_d_"$pha".$fq.$area.$mdl.tc
mv kobs.$mdl.Gw_kernel_d	Gw_kernel_d_"$pha".$fq.$area.$mdl.tc

end  # foreach fq
end  # foreach of phase
end # foreach of ivelin
# change topo-corrected tt in d vector
#change_d.csh

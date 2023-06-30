#!/bin/csh
#gmtset BASEMAP_TYPE plain
gmtset ANOT_FONT_SIZE 12p LABEL_FONT_SIZE 12p
gmtset LABEL_FONT Helvetica ANNOT_FONT_PRIMARY Helvetica ANNOT_FONT_SECONDARY Helvetica HEADER_FONT Helvetica
set SUTURE = mtibet.lineation
# csh plt.cmd 2 100 0.01
echo "input the max value of velocity perturbations: "
set cmax = $1
echo "input damping factor of the model: "
set damp0 = $2
echo "input the threhold of sensitivity matrix to be plotted: "
set perc = $3
echo "input type of model parameterization: e.g., KM.tc, ...: "
set PARA = $4
echo "input type of reference model: e.g., prem, iasp, iasp.tc, prem.tc, ...: "
set refmdl = $5
echo "input weight of station and event corrections: "
set ws = $6
echo "input frequencies of tt measurements, e.g., h, hml, ...:  "
set fq = $7
set idemean = 0
set phase = S
set RK = `echo $PARA | awk '{print substr($1,1,1)}'`
#if ( $RK == "K" ) then
#set fq = hml
#else
#set fq = h
#endif
set weight = ${ws}-$ws
set mdl = `echo $refmdl | awk -F. '{print $1}'`
# set grid interval in degree
set COLOR = purple

set stsym = ( t0.10 s0.10 d0.10 t0.10 i0.10 a0.10 )
set stcol = ( yellow brown 245/150/0 0/200/0 cyan purple )

set secor = se
set area = hic

set TBTDIR = ..
set TBTEVT = ..
set PDIR = .
set BDIR = $TBTDIR/BUILDG

\cp $BDIR/mesh.config .
set PROG = preplot_ray_map100
set PROG1 = preplot_ray_map100_nan

# create color table for dvs
set CPT = dvs.cpt
set ddd = `echo $cmax | awk '{ if ( $1 >= 2. ) {printf("%6.1f\n",$1/2) } \
else if ( $1 < 2 && $1 > 1 )  {printf("%6.1f\n",$1/3)} else {printf("%6.1f\n",$1/2)}} '`
#goto makectbl
svcpt12_table_cont <<!>> junk
-$cmax $cmax
!
mv svel12.cpt $CPT

set width = 2.4
set swathwidth = 1.5

#set im = 2              #AK135-Continent model
#set im = 1              #PREM model
if ( $mdl == "prem" ) then
set im = 1
else if ( $mdl == "iasp" ) then
set im = -4
else if ( $mdl == "ak1c" ) then
set im = 2
cp ../NEWSHOOT/ak135* .
else if ( $mdl == "ak1o" ) then
set im = 3
endif
set ips = 1             #P(=1)

set RANGE = 81/91.5/25.5/36
set lon1 = `echo $RANGE | awk -F/ ' { print $1 } '`
set lon2 = `echo $RANGE | awk -F/ ' { print $2 } '`
set lonc = `echo $lon1 $lon2 | awk ' { print ($1+$2)*0.5 } '`
set lat1 = `echo $RANGE | awk -F/ ' { print $3 } '`
set lat2 = `echo $RANGE | awk -F/ ' { print $4 } '`
set latc = `echo $lat1 $lat2 | awk ' { print ($1+$2)*0.5 } '`
set latf = `echo $lat1 $lat2 | awk ' { print $1-1 } '`
set PROJC = ( 90 90 )
set PROJ = `echo $PROJC | awk '{printf("S%s/%s/%si",$1,$2,w) }' w=$width`
set PROJ = `echo $RANGE | awk -F/ '{printf("S%s/%s/%si",($1+$2)/2,90,w) }' w=$width`

set disc = ( 5711 5961 )
set im = 2              #AK135-Continent model
set ips = 1             #P(=1)
set idemean = 0         #0=no demean

set DLON = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $4 } '`
set DLON2 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print 0.5*$4*'$swathwidth' }'`
echo $DLON2
#set DLON2 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $4 } '`
set DLAT = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $5 } '`
set DDEP = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set NLON = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $4 } '`
set NLAT = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $5 } '`
set NDEP = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set DDEP0 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3 } '`
set DDEP1 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3+('$NDEP'-1)*'$DDEP' } '`
echo DDEP = "$DDEP"  NDEP = "$NDEP"

# for P wave model
set K = `echo $PARA | awk '{print substr($1,1,1) }'`
set P = `echo $PARA | awk '{print substr($1,2,1)}'`
set WDIR = WGT
set dtype =  $area.$refmdl
if ( $K == "R" ) then
set GK = ray
set gtgfile = $TBTDIR/BUILDG/Gw_ray_d_$phase.$fq.$area.$refmdl.gtg
else 
set GK = kern
set gtgfile = $TBTDIR/KERNEL/Gw_kernel_d_$phase.$fq.$area.$refmdl.gtg
endif
if ( ! -f $gtgfile ) then
set FG = `echo $gtgfile | awk -F.gtg '{print $1}'`
$PDIR/get_GTG_med <<!> log
$FG
G.gtg
G.counts
!
# outout velocity images
mv G.gtg $gtgfile
endif
if ( $P == "M" ) then
set GP = multis
set damp = $damp0
else if ( $P == "H" ) then
set GP = hybrid
set damp = ${damp0}-1
else if ( $P == "Q" ) then
set GP = quell
set damp = ${damp0}-1-1
else
set GP = simple
set damp = ${damp0}
endif
set GM = ${GP}_$GK
set TRYDIR = $TBTDIR/INV/$PARA/$WDIR
set wg = `echo $WDIR | tr '[A-Z]' '[a-z]' | awk ' { print substr($1,1,3) }'`
set pstype = $GM.$refmdl.$phase.$fq.${wg}
set stlst = $TBTEVT/$phase.$fq.$area.stations_used
if ( $idemean == 1 ) then
set PS = dvs.$pstype.${weight}_$damp.${perc}.dm.4m2x.ps
else
set PS = dvs.$pstype.${weight}_$damp.${perc}.4m2x.ps
endif
echo $PS

set FF = try."$GM"_"$secor""$weight"_$damp.$phase.$fq.$dtype
set SF = $TRYDIR/stacor."$GM"_"$secor""$weight"_$damp.$phase.$fq.$dtype

set F = $TRYDIR/$FF
echo $F

# outout velocity images

set FV3D =  image.3D."$GM"_"$secor""$weight"_$damp.$phase.$fq.$dtype.nan.$perc
if ( ! ( -f $FV3D ) ) then
\cp /dev/null image.3D
$PDIR/findgrid_sqrtgtg_NaN<<!
$gtgfile
$perc
!
$PDIR/setmodel_NaN<<!
$F
!
cp try.xyz $F.nan.$perc
echo $F.nan.$perc
@ layer = 1
while ( $layer <= $NDEP )
# rad : depth
set rad = `echo $NDEP $layer $DDEP $DDEP1 | awk ' { printf("%10.5f\n",$4-($1-$2)*$3) } '`
set depth = `echo $rad | awk ' { print 6371-$1 } '`
$PDIR/$PROG1 <<!>> junk
$idemean
$im
$ips
$F.nan.$perc
$layer
!
mv slow.slice image.lyr$layer
awk ' { if ( $3 != "NaN" ) { print $1, $2, '$rad', $3 } } ' image.lyr$layer >> image.3D
@ layer ++
end # of while
mv image.3D $FV3D
echo $FV3D
endif

# ----- plot mapview -----
set RMAP =  81/91.5/25.5/35.9

set tx = `echo $RMAP | awk -F/ ' { print ($1+$2)/2 } '`
set ty = `echo $RMAP | awk -F/ ' { print $4+1.1 } '`
set PROJM = `echo $RMAP | awk -F/ '{printf("S%s/%s/%si",($1+$2)/2,90,w) }' w=$width`

@ i = 0
foreach layer ( 31 29 27 22 )
@ i ++
#set depth = `echo $NDEP $layer $DDEP | awk ' { printf("%8.0f\n", ($1-$2)*$3) } '`
set rad = `echo $NDEP $layer $DDEP $DDEP1 | awk ' { printf("%10.1f\n",$4-($1-$2)*$3) } '`  
set depth = `echo $rad | awk ' { print 6371.-$1 } ' `
$PDIR/$PROG1 <<!>> junk
$idemean
$im
$ips
$F.nan.$perc
$layer
!
triangulate slow.slice `minmax -I1 slow.slice` -I0.1/0.1 -Gslice.grd -E > tria.out
#blockmean slow.slice `minmax -I0.1 slow.slice` -I0.1/0.1 > slow.tmp
#surface slow.tmp `minmax -I0.1 slow.slice` -I0.1/0.1 -Gslice.grd
if ( $i == 1 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::Wens -X0.6 -Y5 -K > $PS
else if ( $i == 4 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wEns -X2.7 -O -K >> $PS
else
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wens -X2.7 -O -K >> $PS
endif

psclip ../BUILDG/hiclimb_clip.xy -J$PROJ -R$RMAP -O -K >> $PS
grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -O -K -Ba2f1/a2f1::wens >> $PS
psclip -C -O -K >> $PS
if ( $i == 1 ) then
more +2 $stlst | awk '{print $4, $3} ' | psxy -R -J -St0.1 -W0.5p,0 -O -K >> $PS
cat $SF | awk ' { if ( $3 >= 0 ) { print $2, $1, $3*xsc } }' xsc=1000 | \
psxy  -R -J -O -K -Sx -W2/red >> $PS
cat $SF | awk ' { if ( $3 < 0 ) { print $2, $1, -$3*xsc } }' xsc=1000 | \
psxy  -R -J -O -K -Sc -W2/blue >> $PS
endif
psxy $SUTURE -J -R -M -H1 -A -W.5/0/0/0p -O -K >> $PS
pstext -R -J -O -K -N <<EOF>> $PS
81.5 26 12 0 1 ML $depth km
EOF

end # of foreach layer

# plot cross section
#
#psscale -D-3/-0.7/2.4/0.17h -C$CPT -B"$cmax"::/:"@~d@~lnV@-S@- (%)": -V -E -O -N >> $PS

set yscale = 1.5
set slicewidth = 2.6
set blayer = 18

set COLOR = purple

set stsym = ( t0.10 s0.10 d0.10 t0.10 i0.10 a0.10 )
set stcol = ( yellow brown 245/150/0 0/200/0 cyan purple )

set secor = se

set GCP1 = ( 83.6/30.3 86.5/27.95 )
set GCP2 = ( 90.2/28.6 83.5/35.1 )

set plane = ( xzp yzp )

set RANGE = 81/91.5/25.5/36
set width = 2.2

set disc = ( 5711 5961 )

set t1 = ( "A"  "B"  "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" )
set t2 = ( "A'" "B'" "C'" "D'" "E'" "F'" "G'" "H'" "I'" "J'" "K'" "L'" )
set dtc = 0.02

@ i = 0
foreach SL ( $GCP1 )
@ i ++
set GC1 = $SL
set GC2 = $GCP2[$i]
@ layer = $blayer
set rad = `echo $NDEP $layer $DDEP $DDEP1 | awk ' { printf("%10.5f\n",$4-($1-$2)*$3) } '`
echo $rad
cat $FV3D | awk '{ if ( $3 >= '$rad' ) { print $0 } } ' > image.3D
set FV = image.proj
set kmdeg = `echo $DDEP1 | awk '{print 2*3.1416*$1/360.}'`
set DLON2km = `echo $DLON2 | awk '{print $1*'$kmdeg' }'`
echo $DLON2km $kmdeg
project image.3D -C$GC1 -E$GC2 -W-$DLON2km/$DLON2km -F$plane[$i] -Q -Lw > $FV
echo $FV
cat $FV | awk '{ if ( $2 == 6371 ) {print $0} } ' > surf.proj
cat $FV | awk '{print $1, 6371-$2, $3 }' > $FV.tmp
mv $FV.tmp $FV
sort -k 4 -g surf.proj > surf.sort
set xmin = `head -1 surf.sort | awk '{print $4 } '`
set xmax = `tail -1 surf.sort | awk '{print $4 } '`
set range = `minmax -C $FV | awk ' { print $1 "/" $2 "/" $3 "/" $4 }'`
set rangens = $range
if ( $i == 1 ) then
set rangens = $range
echo "NS: $rangens"
else
set rangeew = $range
echo "EW : $rangeew"
endif

if ( $slicewidth != 0 ) then
if ( $i == 1 ) then
set scaleunit = `echo $xmin $xmax | awk '{print '$slicewidth'/($2-$1) } '`
set xwidth = $slicewidth
set ywidth = `echo $range | awk -F/ '{print ($4-$3)*'$scaleunit' } '`
else
set xwidth = `echo $xmin $xmax | awk '{print ($2-$1)*'$scaleunit' } '`
set ywidth = `echo $range | awk -F/ '{print ($4-$3)*'$scaleunit' } '`
endif
else 
if ( $i == 1 ) then
set xwidth = 1.8037
set ywidth = 1.15
else
set xwidth =  2.24767
set ywidth = 1.15
endif
endif

set nx = `echo $range | awk -F/ '{print int(($2-$1)/dx+1) } ' dx=$DLON`
set DD = `echo $range | awk -F/ '{ print ($2-$1)/n } ' n=$nx`
triangulate $FV -R$range -I$DD/$DDEP -Gslice.grd -E > tria.out
set dist2 = `echo $range | awk -F/ ' { print ($1+$2)*0.5 } '`
set sizearrow = 0.2i
set scalarrow = 0.03i/0.08i/0.04i
set SCALE = X${xwidth}i/-${ywidth}i

switch ($i)

case 1:
# E-W profile
psbasemap -R$range -J$SCALE -Ba2f1/a100f50:"Depth (km)":Wens -K -X-6.5 -Y-2.8 -O >> $PS
grdimage slice.grd -R -J -Ba2f1/a100f50::Wens -C$CPT -O -K -Q >> $PS
#grdcontour slice.grd -R -J -Ccontmdl1.d -O -K -W1.5/255 -B2/200we >> $PS
#grdcontour slice.grd -R -J -Ccontmdl.d -O -K -W1.5/255t12_7_12_7:7 -B2/200we >> $PS
#grdcontour slice.grd -R -J -Ccontmdl2.d -O -K -W1.5/255to -B2/200we >> $PS
set ty = `echo $range | awk -F/ ' { print $4+40 } '`
foreach deg ( 84 86 88 90 )
pstext -R -J -O -K -N <<EOF>>$PS
$deg $ty 12 0 0 MC $deg\312
EOF
end

pstext -R -J -O -K -N <<EOF>>$PS
89.6 350 12 0 2 MC V@-S@-
EOF
breaksw

case 2:
psbasemap -R$range -J$SCALE -Ba2f1/a100f50:"Depth (km)":Wens -O -K -X3.7 -V >> $PS
grdimage slice.grd -R$range -J$SCALE -Ba2f1/a100f50::Wens -C$CPT -O -K -Q -V >> $PS
#grdcontour slice.grd -R -J -Ccontmdl1.d -O -K -W1.5/255 -B2/200we >> $PS
#grdcontour slice.grd -R -J -Ccontmdl.d -O -K -W1.5/255t12_7_12_7:7 -B2/200we >> $PS
#grdcontour slice.grd -R -J -Ccontmdl2.d -O -K -W1.5/255to -B2/200we >> $PS
echo $range
set ty = `echo $range | awk -F/ ' { print $4+40 } '`
foreach deg ( 28 30 32 34 )
pstext -R -J -O -K -N <<EOF>>$PS
$deg $ty 12 0 0 MC $deg\312
EOF
end
pstext -R -J -O -K -N -G0 <<EOF>>$PS
33.85 350 12 0 2 MC V@-S@-
EOF

breaksw

default:
breaksw
echo "none"
breaksw

endsw

set tx1 = `echo $range | awk -F/ ' { print $1-0.1 } '`
set tx2 = `echo $range | awk -F/ ' { print $2+0.1 } '`
pstext -R -J -O -K -N <<EOF>>$PS
$tx1 -50 12 0 6 MC $t1[$i]
$tx2 -50 12 0 6 MC $t2[$i]
EOF
end  # of foreach SL
psscale -D-0.5/-0.6/3/0.17h -C$CPT -B"$cmax"::/:"@~d@~lnV@-S@- (%)": -V -E -O -N >> $PS

/bin/rm image.lyr* slice.grd slow.slice slow.tmp log junk try.xyz tria.out surf.proj surf.sort G.counts

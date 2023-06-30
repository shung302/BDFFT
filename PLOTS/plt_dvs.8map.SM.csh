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
set PROJ = M${width}i

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

set DDEP = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set NLON = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $4 } '`
set NLAT = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $5 } '`
set NDEP = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set DDEP0 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3 } '`
set DDEP1 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3+('$NDEP'-1)*'$DDEP' } '`
#echo DDEP = "$DDEP"  NDEP = "$NDEP"

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
#if ( ! -f $gtgfile ) then
set FG = `echo $gtgfile | awk -F.gtg '{print $1}'`
$PDIR/get_GTG_med <<!> log
$FG
G.gtg
G.counts
!
# outout velocity images
mv G.gtg $gtgfile
#endif
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
set PS = dvs.$pstype.${weight}_$damp.${perc}.dm.8m.ps
else
set PS = dvs.$pstype.${weight}_$damp.${perc}.8m.ps
endif
echo $PS

set FF = try."$GM"_"$secor""$weight"_$damp.$phase.$fq.$dtype
set SF = $TRYDIR/stacor."$GM"_"$secor""$weight"_$damp.$phase.$fq.$dtype

set F = $TRYDIR/$FF
echo $F

# outout velocity images

set FV3D =  image.3D."$GM"_"$secor""$weight"_$damp.$phase.$fq.$dtype.nan.$perc

#if ( ! ( -f $FV3D ) ) then
#if ( $perc == 0 ) then
#cp $F $F.nan.$perc
#else
\cp /dev/null image.3D
$PDIR/findgrid_sqrtgtg_NaN<<!
$gtgfile
$perc
!
$PDIR/setmodel_NaN<<!
$F
!
cp try.xyz $F.nan.$perc
#endif

# ----- plot mapview -----
set RMAP =  81/91.5/25.5/35.9

set tx = `echo $RMAP | awk -F/ ' { print ($1+$2)/2 } '`
set ty = `echo $RMAP | awk -F/ ' { print $4+1.1 } '`

@ i = 0
foreach layer ( 32 30 29 27 26 25 22 18 )
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
#triangulate slow.slice `minmax -I1 slow.slice` -I0.2/0.2 -Gslice.grd -E > tria.out
blockmean slow.slice `minmax -I0.1 slow.slice` -I0.1/0.1 > slow.tmp
surface slow.tmp `minmax -I0.1 slow.slice` -I0.1/0.1 -Gslice.grd
if ( $i <= 4 ) then
if ( $i == 1 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::Wens -X0.6 -Y5 -K > $PS
else if ( $i == 4 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wEns -X2.7 -O -K >> $PS
else
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wens -X2.7 -O -K >> $PS
endif
else
if ( $i == 5 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::WenS -X-8.1 -Y-2.8 -O -K >> $PS
else if ( $i == 8 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wEnS -X2.7 -O -K >> $PS
else
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wenS -X2.7 -O -K >> $PS
endif
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
psscale -D-3/-0.7/2.4/0.17h -C$CPT -B"$cmax"::/:"@~d@~lnV@-S@- (%)": -V -E -O -N >> $PS
ps2pdf $PS dvs.pdf
/bin/rm slice.grd slow.slice slow.tmp log junk try.xyz tria.out

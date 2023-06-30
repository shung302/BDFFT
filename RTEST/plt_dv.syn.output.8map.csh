#!/bin/csh
#gmtset BASEMAP_TYPE plain
gmtset ANOT_FONT_SIZE 12p LABEL_FONT_SIZE 12p
gmtset LABEL_FONT Helvetica ANNOT_FONT_PRIMARY Helvetica ANNOT_FONT_SECONDARY Helvetica HEADER_FONT Helvetica
# csh plt.cmd 2 100 0.01
echo "input the max value of velocity perturbations: "
set cmax = $1
echo "input damping factor of the model: "
set damp0 = $2
echo "input the noise level added in data: "
set noise = $3
echo "input type of model parameterization: e.g., KM.tc, ...: "
set PARA = $4
echo "input type of reference model: e.g., prem, iasp, iasp.tc, prem.tc, ...: "
set ccmdl = $5
set mdl = `echo $ccmdl | awk -F. '{print $1}'`
echo "input frequencies of tt measurements, e.g., h, hml, ...:  "
set fq = $6
echo "input sigma in horizontal and vertical direction, e.g., 1, 2, 3, ...:  "
set sigma = $7
echo "synthetic pattern, e.g., cbt_x4y4z4 "
set pattern = $8
set idemean = 0
set pha = $9
set RK = `echo $PARA | awk '{print substr($1,1,1)}'`
# set grid interval in degree
set COLOR = purple

set stsym = ( t0.10 s0.10 d0.10 t0.10 i0.10 a0.10 )
set stcol = ( yellow brown 245/150/0 0/200/0 cyan purple )

set area = hic

set RMAP =  81/91.5/25.5/35.9
set TBTDIR = ..
set TBTEVT = ..
set PDIR = $TBTDIR/PLOTS
set BDIR = $TBTDIR/BUILDG

\cp $BDIR/mesh.config .
set PROG = preplot_ray_map100
set PROG1 = preplot_ray_map100_nan

# create color table for dvp
set CPT = dv.cpt
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
else if ( $mdl == "ak1o" ) then
set im = 3
endif
if ( $pha == p ) then
set ips = 1             #P(=1)
set phase = P
else
set ips = 2
set phase = S
endif

set DDEP = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set NLON = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $4 } '`
set NLAT = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $5 } '`
set NDEP = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set DDEP0 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3 } '`
set DDEP1 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3+('$NDEP'-1)*'$DDEP' } '`
#echo DDEP = "$DDEP"  NDEP = "$NDEP"
set K = `echo $PARA | awk '{print substr($1,1,1) }'`
set P = `echo $PARA | awk '{print substr($1,2,1)}'`
if ( $K == "R" ) then
set GK = ray
else
set GK = kern
endif


# for P wave model
set dtype =  $area.$ccmdl

if ( $P == "M" ) then
set GP = multis
set damp = $damp0
else if ( $P == "H" ) then
set GP = hybrid
set damp = ${damp0}-${sigma}
else if ( $P == "Q" ) then
set GP = quell
set damp = ${damp0}-${sigma}-${sigma}
else
set GP = simple
set damp = ${damp0}
endif
set GM = ${GP}_$GK
set TRYDIR = $TBTDIR/RTEST/$pattern
set pstype = $GM.$ccmdl.$phase.$fq
set stlst = $TBTEVT/$phase.$fq.$area.stations_used
set PS = dv$pha.$pstype.$damp.${pattern}.n${noise}.8m.ps
echo $PS

set FF = try."$GM"_syn_$damp.$phase.$fq.$dtype.$pattern.n${noise}

set F = $TRYDIR/$FF

# outout velocity images


# ----- plot mapview -----
set tx = `echo $RMAP | awk -F/ ' { print $1+0.5 } '`
set ty = `echo $RMAP | awk -F/ ' { print $3+0.8 } '`

@ i = 0
foreach layer ( 32 30 28 27 26 23 20 16 )
@ i ++
set rad = `echo $NDEP $layer $DDEP $DDEP1 | awk ' { printf("%10.1f\n",$4-($1-$2)*$3) } '`                                           
set depth = `echo $rad | awk ' { print 6371.-$1 } ' `
$PDIR/$PROG1 <<!>> junk
$idemean 
$im
$ips
$F
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

psclip $BDIR/hiclimb_clip.xy -J$PROJ -R$RMAP -O -K >> $PS
grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -O -K -Ba2f1/a2f1::wens >> $PS
psclip -C -O -K >> $PS
#grdcontour slice.grd -R -J -C$CONDIR/contmdl1.d -O -K -W1.5/255 -B1/100we >> $PS
#grdcontour slice.grd -R -J -C$CONDIR/contmdl1.5.d -O -K -W1.5/255t12_7_12_7:7 -B1/100we >> $PS
#grdcontour slice.grd -R -J -C$CONDIR/contmdl2.d -O -K -W1.5/255 -B1/100we >> $PS

if ( $i == 1 ) then
more +2 $stlst | awk '{print $4, $3} ' | psxy -R -J -St0.15 -W1p,0 -O -K >> $PS
endif
#grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -O -K -Ba2f1/a2f1::wens >> $PS
#grdcontour slice.grd -J -R -O -K -A1 -C1 -W1p/255 >> $PS 
#psclip -J$PROJ -R$RMAP -C -O -K >> $PS
#psxy $BDIR/mesh.xy1 -R -J -W1p,0,. -M -O -K >> $PS

pstext -R -J -O -K -N <<EOF>> $PS
$tx $ty 14 0 1 ML $depth km
EOF

end # of foreach layer
psscale -D-3/-0.7/2.4/0.17h -C$CPT -B"$cmax"::/:"@~d@~lnV${pha} (%)": -V -E -O -N >> $PS
/bin/rm junk slow.* slice.grd

#!/bin/csh
#gmtset BASEMAP_TYPE plain
gmtset ANOT_FONT_SIZE 12p LABEL_FONT_SIZE 12p
gmtset LABEL_FONT Helvetica ANNOT_FONT_PRIMARY Helvetica ANNOT_FONT_SECONDARY Helvetica HEADER_FONT Helvetica
#echo "input the max value of GTG: "
set SUTURE = mtibet.lineation
set cmax = 1.5
set cmin = -1
set phase = S
set area = hic
echo "input kernel or G-matrix type: ray or kernel "
set K = $1
echo "input type of reference model: e.g., prem, iasp, iasp.tc, prem.tc, ...: "
set mdl = $2
echo "input type of G matrix, e.g., G, Gw, Gdif, Gwdif: "
set GW = $3
echo "input frequency of tt measurements, e.g., h, hml,...: "
set fq = $4

# set grid interval in degree
set COLOR = purple

set stsym = ( t0.10 s0.10 d0.10 t0.10 i0.10 a0.10 )
set stcol = ( yellow brown 245/150/0 0/200/0 cyan purple )

set secor = se

set TBTDIR = ..
set PDIR = .
set BDIR = $TBTDIR/BUILDG

\cp $BDIR/mesh.config .
set PROG = preplot_gtg_map

#---------------------------------------
# set hot color scale for kernel plot
#---------------------------------------
set CPT = hot2.cpt
set DDD = `echo $cmin $cmax | awk ' { print ($2-$1)/10 } '`
makecpt -Chot -T$cmin/$cmax/$DDD -I > cpt
set r = ( `tail -4 cpt | head -1 | awk ' { print $6 } '` )
set g = ( `tail -4 cpt | head -1 | awk ' { print $7 } '` )
set b = ( `tail -4 cpt | head -1 | awk ' { print $8 } '` )
cat cpt |  awk ' { if ( $1 == "F" ) { printf("%s\t%d\t%d\t%d\n", $1, '$r', '$g', '$b') } else { print $0 }}' > $CPT

set width = 2.4
set PROJ = M${width}i

set im = 2              #AK135-Continent model
set ips = 1             #P(=1)

set DDEP = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set NLON = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $4 } '`
set NLAT = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $5 } '`
set NDEP = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set DDEP0 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3 } '`
set DDEP1 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3+('$NDEP'-1)*'$DDEP' } '`
#echo DDEP = "$DDEP"  NDEP = "$NDEP"

# for P wave model
set weight = `echo $GW | awk '{print substr($1,1,2)}'`
if ( $weight == "Gw" ) then
set dtype = $mdl.$phase.$fq.wgt
else
set dtype = $mdl.$phase.$fq.nowgt
endif
set PS = gtg.$K.$dtype.8m.ps

set FID = $phase.$fq.$area
set F = $phase.$fq.$area
set SF = $TBTDIR/$FID.stations_used
if ( $K == "ray" ) then
set KDIR = $TBTDIR/BUILDG
set FG = $KDIR/"$GW"_"$K"_d_${FID}.$mdl
else
set KDIR = $TBTDIR/KERNEL
set FG = $KDIR/"$GW"_"$K"_d_${FID}.$mdl
endif
echo $FG 

if ( ! -f $FG.gtg.orig ) then
$PDIR/get_GTG_med <<!> log
$FG
G.gtg
G.counts
!
# outout velocity images
mv G.gtg $FG.gtg.orig
endif

# outout velocity images
set FV3D =  image.3D."$GW"_"$K"_d_${FID}.$mdl.gtg.orig
if ( ! ( -f $FV3D ) ) then
\cp /dev/null image.3D
@ layer = 1
while ( $layer <= $NDEP )
set rad = `echo $NDEP $layer $DDEP $DDEP1 | awk ' { printf("%10.5f\n",$4-($1-$2)*$3) } '`
set depth = `echo $rad | awk ' { print 6371.-$1 } ' `
$PDIR/preplot_gtg_map<<! >/dev/null
$FG.gtg.orig
$layer
!
mv gtg.slice image.lyr$layer
awk ' { print $1, $2, '$rad', sqrt($3) } ' image.lyr$layer >> image.3D
@ layer ++
end # of while
mv image.3D $FV3D
endif

# ----- plot mapview -----
set RMAP =  81/91.5/25.5/35.9
set RBOX =  81/91.5/25.5/35.9
set R = `echo $RBOX | awk -F/ '{print $1,$2,$3,$4}'`
set TILL = illum_topo.grd

set tx = `echo $RMAP | awk -F/ ' { print ($1+$2)/2 } '`
set ty = `echo $RMAP | awk -F/ ' { print $4+1.1 } '`

@ i = 0
#foreach layer ( 32 30 28 26 24 22 19 16 )
foreach layer ( 32 30 29 27 26 25 22 18 )

@ i ++
set rad = `echo $NDEP $layer $DDEP $DDEP1 | awk ' { printf("%10.1f\n",$4-($1-$2)*$3) } '`
set depth = `echo $rad | awk ' { print 6371.-$1 } ' `
# get gtg values
$PDIR/$PROG <<!> junk
$FG.gtg.orig
$layer
!
mv gtg.slice image.lyr$layer
awk ' { print $1, $2, sqrt($3) } ' image.lyr$layer > gtg.slice
triangulate gtg.slice `minmax -I1 gtg.slice` -I0.1/0.1 -Gslice.grd -E > tria.out
# plot in log scale
grdmath slice.grd LOG10 = logslice.grd
mv logslice.grd slice.grd

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
grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -O -K -Ba2f1/a2f1::wens >> $PS
if ( $i == 1 ) then
more +2 $SF | awk '{print $4, $3} ' | psxy -R -J -St0.1 -W1p,0 -O -K >> $PS
endif
#psxy $BDIR/mesh.xy1 -J -R -O -K -M -W1p,black >>$PS
psxy $SUTURE -J -R -M -H1 -A -W.5/0/0/0p -O -K >> $PS
pstext -R -J -O -K -N <<EOF>> $PS
81.5 26 12 0 1 ML $depth km
EOF
end # end of each lyr
psscale -D-3/-0.9/3./0.17h -C$CPT -B1::/:"log(diag(G@+T@+G)@+1/2@+)": -V -E -O -N >> $PS                                                                                                     
#psscale -D-3/-0.9/3./0.17h -C$CPT -B1::/:"log(diag(G@+T@+G)@+1/2@+)": -V -E -O -N -UBL/-8i/-1i/"gtg.$K.$dtype" >> $PS                                                                                                     
/bin/rm image.lyr* junk gtg.slice cpt tria.out slice.grd

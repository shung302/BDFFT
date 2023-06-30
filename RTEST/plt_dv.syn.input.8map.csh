#!/bin/csh
#gmtset BASEMAP_TYPE plain
gmtset ANOT_FONT_SIZE 12p LABEL_FONT_SIZE 12p
gmtset LABEL_FONT Helvetica ANNOT_FONT_PRIMARY Helvetica ANNOT_FONT_SECONDARY Helvetica HEADER_FONT Helvetica
set cmax = 4
set pattern = $1
set ips = $2
# set grid interval in degree
set COLOR = purple

set stsym = ( t0.10 s0.10 d0.10 t0.10 i0.10 a0.10 )
set stcol = ( yellow brown 245/150/0 0/200/0 cyan purple )

set TBTDIR = ..
set PDIR = ../PLOTS
set BDIR = $TBTDIR/BUILDG

set RMAP =  81/91.5/25.5/35.9
\cp $BDIR/mesh.config .
set PROG = preplot_ray_map100
set PROG1 = preplot_ray_map100_nan

# create color table for dvp
set CPT = dv.cpt
set ddd = `echo $cmax | awk '{ if ( $1 >= 2. ) {printf("%6.1f\n",$1/2) } \
else if ( $1 < 2 && $1 > 1 )  {printf("%6.1f\n",$1/3)} else {printf("%6.1f\n",$1/2)}} '`
#goto makectbl
svcpt12_table_cont <<!>! junk
-$cmax $cmax
!
mv svel12.cpt $CPT

set width = 2.4
set PROJ = M${width}i

set im = 2              #AK135-Continent model
#set im = -4              #iasp model
#P(=1) S(=2)
if ( $ips == 1 ) then
set phase = p 
else
set phase = s 
endif
set idemean = 0

set DDEP = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set NLON = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $4 } '`
set NLAT = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $5 } '`
set NDEP = `more +1 $BDIR/mesh.config | head -1 | awk ' { print $6 } '`
set DDEP0 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3 } '`
set DDEP1 = `more +2 $BDIR/mesh.config | head -1 | awk ' { print $3+('$NDEP'-1)*'$DDEP' } '`
#echo DDEP = "$DDEP"  NDEP = "$NDEP"

set TRYDIR = $TBTDIR/RTEST
set PS = dv${phase}.syn.input.$pattern.8m.ps
echo $PS

set FF = try.slowness.syn.$pattern

set F = $TRYDIR/$FF
echo $F

# outout velocity images

# ----- plot mapview -----

set tx = `echo $RMAP | awk -F/ ' { print $1+0.5 } '`
set ty = `echo $RMAP | awk -F/ ' { print $3+0.8 } '`

@ i = 0
foreach layer ( 32 30 28 27 26 23 20 16 )                                                                                                                           
@ i ++
set rad = `echo $NDEP $layer $DDEP $DDEP1 | awk ' { printf("%10.1f\n",$4-($1-$2)*$3) } '`
set depth = `echo $rad | awk ' { print 6371.-$1 } ' `
$PDIR/preplot_ray_map100 <<!>! junk
$idemean
$im
$ips
$F
$layer
!
blockmean slow.slice `minmax -I1 slow.slice` -I0.1/0.1 > slow.tmp
surface slow.tmp `minmax -I1 slow.slice` -I0.1/0.1 -Gslow.grd
grd2xyz slow.grd > slow.xyz
paste slow.xyz | awk '{print $1, $2, $3} ' > slow.slice
#triangulate slow.slice `minmax -I1 slow.slice` -I0.1/0.1 -Gslice.grd -E > tria.out
blockmean slow.slice `minmax -I0.1 slow.slice` -I0.1/0.1 > slow.tmp                                                                 
surface slow.tmp `minmax -I0.1 slow.slice` -I0.1/0.1 -Gslice.grd

if ( $i <= 4 ) then
if ( $i == 1 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::WeNs -X0.6 -Y5 -K > $PS
else if ( $i == 4 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wENs -X2.7 -O -K >> $PS
else
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::weNs -X2.7 -O -K >> $PS
endif
else
if ( $i == 5 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::WenS -X-8.1 -Y-3.0 -O -K >> $PS
else if ( $i == 8 ) then
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wEnS -X2.7 -O -K >> $PS
else
psbasemap -R$RMAP -J$PROJ -Ba2f1/a2f1::wenS -X2.7 -O -K >> $PS
endif
endif

psclip $BDIR/hiclimb_clip.xy -J$PROJ -R$RMAP -O -K >> $PS
grdimage slice.grd -J$PROJ -R$RMAP -C$CPT -O -K -Ba2f1/a2f1::wens >> $PS
#grdcontour slice.grd -J -R -O -K -A1 -C1 -W1p/255 >> $PS 
psclip -C -O -K >> $PS

pstext -R -J -O -K -N <<EOF>> $PS
$tx $ty 14 0 1 ML $depth km
EOF
end # of foreach layer
psscale -D-3/-0.7/2.4/0.17h -C$CPT -B"$cmax"::/:"@~d@~lnV${phase} (%)": -V -E -O -N -UBL/-8/-0.5/"$PS" >> $PS
/bin/rm junk slow.* slice.grd

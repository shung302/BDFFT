#!/bin/csh
#before run mesh_car.f, redefine the dimension of lati, long, depth
# for Taiwan's Chainan region
gmtset BASEMAP_TYPE fancy LABEL_FONT_SIZE 18p ANOT_FONT_SIZE 16p
gmtset PLOT_DEGREE_FORMAT ddd:mm:ss

set ngrid = ( 65 33 65 )
set ngrid = ( 65 33 33 )
set RANGE =  80/93/23/38
set gcp = ( 88.2 26.4 84.5 35.3 )

# make a oblique mesh and plot topo+crust+mesh together
set lon1 = `echo $RANGE | awk -F/ ' { print $1 } '`
set lon2 = `echo $RANGE | awk -F/ ' { print $2 } '`
set lonc = `echo $lon1 $lon2 | awk ' { print ($1+$2)*0.5 } '`
set lat1 = `echo $RANGE | awk -F/ ' { print $3 } '`
set lat2 = `echo $RANGE | awk -F/ ' { print $4 } '`
set latc = `echo $lat1 $lat2 | awk ' { print ($1+$2)*0.5 } '`

set ddep = 900
set ddep0 = 0.
set PS = oblique_mesh.ps
/bin/rm $PS
./oblique_mesh<<!
$gcp[1] $gcp[2] $gcp[3] $gcp[4]
$ngrid[1] $ngrid[2] $ngrid[3]
$ddep0 $ddep
!
set SCALE = B"$lonc"/$latc/20/40/6i

pscoast -J$SCALE -R$RANGE -N1/0.25tap -W0.25p/255/255/255 \
-P -K -X1 -Y3 -B5f1/a5f1 > $PS
psxy mesh.xy -R -J$SCALE -M -O -K -W1p/255/0/0 >> $PS
more +2 stations_tibet | awk ' { print $4, $3 } ' | \
psxy -R -J$SCALE -O -St0.05 -G0 -W1p/0 >>$PS
mv oblique_mesh.config mesh.config


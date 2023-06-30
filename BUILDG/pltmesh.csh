#!/bin/csh
#before run mesh_car.f, redefine the dimension of lati, long, depth
# for Taiwan's Chainan region
set RANGE = 85/93/26.7/34.7
set ngrid = ( 33 33 33 )
set RANGE =  83/91/27/35
set RMAP = 81/93/25/37
set lonc = `echo $RANGE | awk -F/ ' { print ($1+$2)*0.5 } '`
set latc = `echo $RANGE | awk -F/ ' { print ($3+$4)*0.5 } '`
# half-width in longitude and latitud direction
set dlon = `echo $RANGE | awk -F/ ' { print ($2-$1)*0.5 } '`
set dlat = `echo $RANGE | awk -F/ ' { print ($4-$3)*0.5 } '`
set ddep = 900
set ddep0 = 0.
set PS = mesh.ps
./mesh_car<<!
$ngrid[1] $ngrid[2] $ngrid[3]
$lonc $latc
$dlon $dlat $ddep $ddep0
!
set SCALE = C"$lonc"/$latc/6.5
set SCALE = M6i
pscoast -J$SCALE -R$RMAP -N1/0.25tap -W0.25p/255/255/255 \
-P -K -X1 -Y3 -B5f1/a5f1 > $PS
psxy mesh.plot -R -J$SCALE -M -O -K -W1p/255/0/0 >> $PS
more +2 stations_tibet | awk ' { print $4, $3 } ' | \
psxy -R -J$SCALE -O -St0.05 -G0 -W1p/0 >>$PS

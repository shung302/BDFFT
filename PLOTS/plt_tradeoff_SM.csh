#!/bin/csh
# modified by plttbt_tradeoff_kern_SM.cc.cmd in HICLMB
gmtset ANOT_FONT_SIZE 14p LABEL_FONT_SIZE 16p HEADER_FONT_SIZE 18p
echo "enter type of data, e.g., prem.tc, ak1c.HiCLIMB3, ..."
set refmdl = $1
set mdl = `echo $refmdl | awk -F. '{print $1 }'`
echo "enter KS, KM, KM.tc, RM, RM.tc ... "
set qtype = $2
set area = hic
set INVDIR = ../INV
set KDIR = ../KERNEL
set RDIR = ../BUILDG
set PD = .
echo " ===> set up bound for damping factor ??  (y or n) "
set flag_damp = $3
#if ( $flag_damp == "y" ) then
echo " ===> enter minimum damping factor: "
set damp_min = $4
echo " ===> enter maximum damping factor: "
set damp_max = $5
#endif
#echo "====> fix range in y direction:  (y or n) "
set flag_yaxis = $6
echo "====> phase?"
set pha = $7
echo "====> frequency? "
set freq = $8
echo "====> weight of sta and evt corr.? if both, type n-n, where n=weight "
set weight = $9
set wse = ${weight}-${weight}
set secor = se
foreach phase ( $pha )
foreach fq ( $freq )
foreach WDIR ( WGT )
set wg = `echo $WDIR | tr '[A-Z]' '[a-z]' | awk ' { print substr($1,1,3) }'`
foreach PARA ( `echo ${qtype} ` )
set RK = `echo $PARA | awk ' { if  ( substr($1,1,1) == "R" ) { printf "ray" } else {print "kern"  } } '` 
set par = `echo $PARA | awk ' { print substr($1,2,1) } '`
set RG = `echo $par | awk '{ if ( $1 == "S" ) { printf "simple" } \
else if ( $1 == "Q" ) { printf "quell" } else if ( $1 == "H" ) { printf "hybrid" } \
else { printf "multis" } } '`
set FDIR = $INVDIR/$PARA/$WDIR
echo $phase $fq $RG $WDIR
if ( $RK == "ray" ) then
set WGMX = $RDIR/Gw_ray_d_$phase.$fq.$area.$refmdl
set GMX = $RDIR/G_ray_d_$phase.$fq.$area.$refmdl
else
set WGMX = $KDIR/Gw_kernel_d_$phase.$fq.$area.$refmdl
set GMX = $KDIR/G_kernel_d_$phase.$fq.$area.$refmdl
endif
endif
\ls $FDIR/try."$RG"_"$RK"_"$secor""$wse"_1000.$phase.$fq.$area.$refmdl
#\ls $FDIR/try."$RG"_"$RK"_"$secor""$wse"_*.$phase.$fq.$area.$refmdl

cp /dev/null damp.tmp
foreach damp ( `\ls $FDIR/try."$RG"_"$RK"_"$secor""$wse"_*.$phase.$fq.$area.$refmdl | \
awk -F$secor ' { print $2 }' | awk -F_ '{ print $2 }' | awk -F. ' { print $1 } '` )
if ( $flag_damp == "y" ) then
echo $damp | awk ' { if ( $1 >= '$damp_min' && $1 <= '$damp_max' ) { print $1 }} ' >> damp.tmp
else
echo $damp >> damp.tmp
endif
end
set nsize = `wc -l damp.tmp | awk ' { print $1 } '`
if ( $nsize != 0 ) then
sort -n damp.tmp > junk
mv junk damp.tmp
else
echo "no files exist!"
goto skip
endif
set chisq = "NaN"
set TFILE = misfit."$RG"_"$RK"_"$secor""$wse".$phase.$fq.$wg.$refmdl.gmt
cat /dev/null > $TFILE
foreach damp ( `cat damp.tmp` )
set TRY = $FDIR/try."$RG"_"$RK"_"$secor""$wse"_"$damp".$phase.$fq.$area.$refmdl
set SUM = $FDIR/sum."$RG"_"$RK"_"$secor""$wse"_"$damp".$phase.$fq.$area.$refmdl
set vd = `more +6 $SUM | head -1 | awk ' { print $1 } '`
set mvar = `more +6 $SUM | head -1 | awk ' { print $2 } '`
set mnorm = `more +6 $SUM | head -1 | awk ' { print $3 } '`
echo $damp $vd $chisq $mvar $mnorm >> $TFILE
end
echo $TFILE
set PS = vd."$RG"_"$RK"_"$secor""$wse".$phase.$fq.$wg.$refmdl.ps
echo $PS
# plot variance reduction vs. model variance
set PROJ = X4.5i/4.5i
if ( $RK == "ray" ) then
set x = 100000
else
set x = 10000
endif
set logx = `awk ' BEGIN { print int(log('$x')/log(10)+0.01) } '`
set top = `head -1 $TFILE | awk ' { print ($4+$4*0.1)*'$x' }'`
set bottom = `head -1 $TFILE | awk ' { print (-$4*0.1)*'$x' }'`
echo "model norm in y-axis multiplied by $x = $bottom $top"
#if ( $flag_yaxis == "n" ) then
#set RANGE = 0/100/$bottom/$top
#echo "unfixed range = $RANGE "
#else
echo "y axis multiplied by $x : "
set ymax = $10
set RANGE = 0/100/-0.005/$ymax
#endif
echo "tick mark in y axis: "
set ytick = $11
set yanot = `echo $ytick | awk '{print $1*2 }'`
awk '{ print $2, $4*'$x' }' $TFILE > junk
psxy junk -R$RANGE -J$PROJ -Ba10f5:"variance reduction (%)":/a${yanot}f${ytick}:"model variance  @~\264@~ 10@+"$logx"@+"::."$phase.$fq $RK $RG ws-$wse $wg ":WeSn \
-Sc0.08 -G0 -W2 -Y2 -X1 -K -U/0/-1/"$INVDIR/$PARA" > $PS
psxy junk -R -J -W2 -O -K >> $PS
awk ' { printf("%10.2e %10.2e 12 0 4 LB %3d\n", $2, $4*'$x', int($1)) } ' $TFILE | \
pstext -R -J -O -K -N >> $PS

# plot varianc reduction vs. model norm
if ( $RK == "ray" ) then
set x = 1000
else
set x = 100
endif
set logx = `awk ' BEGIN { print int(log('$x')/log(10)+0.01) } '`
set top = `head -1 $TFILE | awk ' { print ($5+$5*0.1)*'$x' }'`
set bottom = `head -1 $TFILE | awk ' { print (-$5*0.1)*'$x' }'`
echo "model variance in y-axis multiplied by $x = $bottom $top"
#if ( $flag_yaxis == "n" ) then
#set RANGE = 0/100/$bottom/$top
#echo "unfixed range = $RANGE "
#else
echo "y axis multiplied by $x :  "
set ymax = $12
set RANGE = 0/100/-0.005/$ymax
#endif
echo "tick mark in y axis: "
set ytick = $13
set yanot = `echo $ytick | awk '{print $1*2 }'`
awk '{ print $2, $5*'$x' }' $TFILE > junk
psxy junk -R$RANGE -J$PROJ -Ba10f5:"variance reduction (%)":/a${yanot}f${ytick}:"model norm  @~\264@~ 10@+"$logx"@+"::."$phase.$fq $RK $RG ws-$wse $wg ":WeSn -Sc0.08 -G0 -W2 -X5.3 -O -K >> $PS
psxy junk -R -J -W2 -O -K >> $PS
awk ' { printf("%10.2e %10.2e 12 0 4 LB %3d\n", $2, $5*'$x', int($1)) } ' $TFILE | \
pstext -R -J -O -N >> $PS

end #foreach PARA
#end #foreach WS
end #foreach WDIR
end #foreach fq
end #foreach phase
# plot tradeoff curve
/bin/rm junk damp.tmp

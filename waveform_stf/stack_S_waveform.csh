#!/bin/csh
# sum of all seismograms 
set DIR = S.2004.0728.0356
set PICKDIR = .
set events = ( S.2004.0728.0356 S.2004.0728.0356 )
set evid = `echo $DIR | awk -FS. ' { print $2 } '`
echo $evid
set wbeg = ( -8 -5 )
set wend = ( 12  30 )
set twidth = ( 0.1 0.15 )
foreach phase ( S )
@ iw = 0
foreach fq ( h m )
@ iw++
set DIR = $events[$iw]
#cp $PICKDIR/$DIR.$fq.t.o.x .
set evid = `echo $DIR | awk -FS. ' { print $2 } '`
set wb = $wbeg[$iw]
set we = $wend[$iw]
set tw = $twidth[$iw]
set TTFILE = $phase.$evid.$fq.t.o.x
echo $DIR $TTFILE
@ i = 0
set ns = `tail +2 $TTFILE | head -1 | awk ' { print $1 } '`
set TT = `tail +3 $TTFILE | head -$ns | awk ' { print -$2 } '`
set FILE1 = `tail +3 $TTFILE | head -1 | awk ' { print $8 ".mccc" }'`
@ ns1 = $ns - 1
set FILE = `tail +4 $TTFILE | head -$ns1 | awk '{ print $8 ".mccc" }'`
echo "r $DIR/$FILE1" > add_seis.m
foreach STN ( `echo $FILE` )
echo "addf $DIR/$STN" >> add_seis.m
echo "w $phase.stack.$fq.t" >> add_seis.m
echo "r $phase.stack.$fq.t" >> add_seis.m
end
echo "cut $wb $we" >> add_seis.m
echo "r $phase.stack.$fq.t" >> add_seis.m
echo "taper t cosine w $tw" >> add_seis.m
echo "fft" >> add_seis.m
echo "wsp $phase.stack.$fq.t" >> add_seis.m
echo "r $phase.stack.$fq.t.am" >> add_seis.m
echo "w alpha $phase.stack.$fq.t.am.ascii" >> add_seis.m
sac<<!
m add_seis.m
quit
!

end
end  # end of foreach phase
/bin/rm S.??.*.t $phase.*.$fq.*1 j?

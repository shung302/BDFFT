#!/bin/csh
set TBTDIR = ..
set PD = .
# raytracing P, PKPdf, PcP etc.
# 1st line: insflag (=1 interactively; =2 from file)
# 2nd line: ivelin (1.$phaREM (prem); -4=iasp91(iasp), 2=AK135-continent(ak1c))
# 3rd line: phase= 1 for P; =18 for PKPdf; =12 for Pdiff; 4 for PcP
#                 =21 for S; =36 for SKSac; =32 for Sdiff; 24 for ScS 
foreach freq ( m )
set insflag = 1
set ivelin = 2
set phase = 21 
set pha = S
# for S
$PD/shootray_evt<<! >run.out
1
$ivelin
$phase
1000. 0.
10. 5. .5
500	
1, 1	
1.0e-10, 0.5	
1.
1
$TBTDIR/S.$freq.hic.evtlst
!
mv shootrays.bp.$pha shootrays.bp.$pha.$freq.hic.ak1c
mv shootrays.$pha  shootrays.$pha.$freq.hic.ak1c
end
exit

# for SKSac
set phase = 36
$PD/shootray_evt<<!>> run.out
1
$ivelin
$phase
1000. 0.
10. 5. .5
500
1, 1
1.0e-10, 0.5
1.
1
$TBTDIR/$pha.$freq.hic.evtlst
!
mv shootrays.bp.$pha shootrays.bp$pha.$freq.hic.ak1c
mv shootrays.$pha  shootrays.$pha.$freq.hic.ak1c

set phase = 24
set pha = ScS
$PD/shootray_evt<<!> run.out
1
$ivelin
$phase
1000. 0.
10. 5. .5
500
1, 1
1.0e-10, 0.5
1.
1
$TBTDIR/$pha.$freq.hic.evtlst
!
mv shootrays.bp.$pha shootrays.bp.$pha.$freq.hic.ak1c
mv shootrays.$pha  shootrays.$pha.$freq.hic.ak1c

end

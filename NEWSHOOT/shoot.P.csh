#!/bin/csh
set PD = .
set TBTDIR = .
set area = tbt
# raytracing $phase, for P and PKIKP temporary
# 1st line: insflag (=1 interactively; =2 from file)
# 2nd line: ivelin (1=PREM (prem); -4=iasp91(iasp), 2=AK135-continent($mdl) 3=AK135-global average ($mdl))
# 3rd line: pha  = 1 for P; =18 for PKPdf; =12 for Pdiff; = 4 for PcP
#                =21 for S; =36 for SKSac; =32 for Sdiff; =24 for ScS
# two output files
foreach phase ( P )
foreach freq ( h )
# ivelin =  1  for PREM
# ivelin = -4 for IASP91
# ivelin = 2 for ak135-continent model
# ivelin = 3 for ak135 global average model
#foreach ivelin ( -4 1 )
foreach ivelin ( -4 )
if ( $ivelin == 1 ) then
set mdl = prem
else if ( $ivelin == 3 ) then
set mdl = ak1o
else if ( $ivelin == 2 ) then
set mdl = ak1c
else
set mdl = iasp
endif
echo $mdl
if ( $phase == "P" ) set pha = 1
if ( $phase == "PcP" ) set pha = 4
if ( $phase == "PKPdf" || $phase == "PKIKP" ) set pha = 18
#if ( -f shootrays.$phase.$freq.$area.$mdl ) continue
echo $pha
echo "Now raytracing $phase phase, $freq band"
$PD/shootray_evt<<!> run.out
1
$ivelin 
$pha
1000. 0.    
10. 5. .5
500          
1, 1          
1.0e-10, 0.5    
1.               
1
$TBTDIR/$phase.$freq.$area.evtlst
!
if ( $pha == 18 ) then
mv shootrays.bp.PKPdf shootrays.bp.$phase.$freq.$area.$mdl
mv shootrays.PKPdf  shootrays.$phase.$freq.$area.$mdl
mv topocorr.PKPdf topocorr.$phase.$freq.$area.$mdl
else
mv shootrays.bp.$phase shootrays.bp.$phase.$freq.$area.$mdl
mv shootrays.$phase  shootrays.$phase.$freq.$area.$mdl
mv topocorr.$phase topocorr.$phase.$freq.$area.$mdl
endif
end
end
end

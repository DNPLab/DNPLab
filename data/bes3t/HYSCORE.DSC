#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	CPLX
XTYP	IDX
YTYP	IDX
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
IIFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	175
XMIN	40.000000
XWID	3480.000000
YPTS	175
YMIN	40.000000
YWID	3480.000000
*
*	Documentational Text:
*
TITL	'HYSCORE310'
IRNAM	'Intensity'
IINAM	'Intensity'
XNAM	'Time'
YNAM	'Time'
IRUNI	''
IIUNI	''
XUNI	'ns'
YUNI	'ns'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    tomc
DATE    02/14/12
TIME    14:46:04
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    PLS
OXS1    TADC
AXS1    B0VL
AXS2    NONE
AXS3    
MWPW    0.0006284
A1CT    0.31
B0VL    0.31
A1SW    0.0
MWFQ    9.683401e+09
AVGS    1
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     cwBridge, 1.0

AcqFineTuning      Never
Power              0.6284 mW
PowerAtten         25.0 dB

.DVC     endor, 1.0

EIEENDORFreq       14.9021800 MHz/3.5 kG
EIEIsotope         H1
EIERFSweepDir      Same
EIEStaticField     3480.000 G
EIEStaticRF        1.000 MHz
RF1Atten           60.000 dB
RF1FreqPos         1.000 MHz
RF1StartFreq       1.000 MHz
RF1SweepWidth      30.000 MHz
RF2Atten           60.000 dB
RF2FreqPos         1.000 MHz
RF2StartFreq       1.000 MHz
RF2SweepWidth      30.000 MHz
RFSrcMixing        Add
SumAtten           10.000 dB
SumAttenStart      0.000 dB
SumAttenWidth      80.000 dB

.DVC     fieldCtrl, 1.0

CenterField        3100.00 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
SweepDirection     Up
SweepWidth         0.0 G

.DVC     freqCounter, 1.0

FrequencyMon       9.683401 GHz
QMonitBridge       On

.DVC     ftBridge, 1.0

Attenuation        7.00 dB
ELDORAtt           30 dB
FrequencyA         9.500000 GHz
VideoBW            200 MHz
VideoGain          48 dB

.DVC     ftEpr, 1.0

AveragesPerScan    50
ELDORFreqStart     9.300000 GHz
ELDORFreqWidth     0.700000 GHz
FTAcqModeSlct      Run from PulseSPEL
PPExtTrg           False
PPExtTrgSlope      Rising
PlsSPELEXPSlct     HYSCORE
PlsSPELGlbTxt      \
;\n\
;       HYSCORE package\n\
;     PulseSPEL general variables, definitions & conventions\n\
begin defs\n\
\n\
; Variables    ; Comments:\n\
\n\
p0   = 16      ; 90 pulse length\n\
p1   = 32      ; 180 pulse length\n\
p2   = 28      ; Inversion Pulse\n\
p3   = 72      ; Matched Pulse\n\
\n\
pg   = 52      ; Integrator Gate Width (RESERVED)\n\
d0   = 452     ; Detection Trigger offset time\n\
\n\
d1   = 152     ; tau\n\
d2   = 40      ; t1\n\
d3   = 40      ; t2\n\
d4   = 800     ; tau2\n\
d5   = 400     ; Remote Detection storage time (< T1)\n\
d6   = 8       ; offset for alignment of Store pulse in Remote Detection\n\
\n\
d9   = 0       ; DAF (Delay After Flash) delay pulse sequence until\n\
               ; after external event (eg. laser pulse)\n\
d30  = 20      ; 1st time scale increment\n\
d31  = 20      ; 2nd time scale increment\n\
\n\
a    = 50      ; (RESERVED) number of transient averages (TRANS)\n\
               ; if greater than 1024, a must be a multiple of 1024\n\
h    = 50      ; (CONVENTION) number of shots per point (INTG and SPT)\n\
n    = 1       ; (CONVENTION) number of sweeps to accumulate\n\
\n\
srt  = 1000 * srtu ; (RESERVED) SRT - Shot Repetition Time (srtu=1.02 us)\n\
\n\
b    = 1       ;\n\
c    = 1       ;\n\
r    = 1       ;\n\
w    = 800     ; (CONVENTION) Split point of sweeps (points)\n\
               ; i.e. pulse separation where TWT gate pulse can be split.\n\
dx   = 0       ; (RESERVED) 1st time scale starting value\n\
dy   = 0       ; (RESERVED) 2nd time scale starting value\n\
end defs\n\
\n\
; Note:\n\
;       Variables Marked 'RESERVED' = variable may be used for other times,\n\
;                                       however only this variable can change\n\
;                                        the indicated parameter.\n\
;                                       eg. only 'a' will change the number of\n\
;                                           transients averaged in SpecJet\n\
\n\
\n
PlsSPELLISTSlct    4-step
PlsSPELPhPrgEx     Normal
PlsSPELPrg         PulseSPEL/4pHYSCORE.exp
PlsSPELPrgTxt      \
;\n\
;\n\
;HYSCORE\n\
;\n\
;P. Hofer, A. Grupp, G. Nebenfuhr and M. Mehring\n\
;		Chem. Phys. Lett. 1986, v. 132, p. 279\n\
;C. Gemperle, G. Aebli, A. Schweiger, and R. R. Ernst\n\
;			J. Man. Reson. 1990, v. 88, p. 241\n\
;\n\
begin defs		; dimension of data array [sx,sy]\n\
 dim s[512]		; 	for HYSCORE set-up\n\
 dim1 s[32]		;	for Inversion Optimization\n\
 dim2 s[175,175]	; 	for HYSCORE\n\
end defs\n\
;\n\
;\n\
begin lists "none"\n\
 ph1 +x\n\
 ph2 +x\n\
 ph3 +x\n\
 ph4 +x\n\
 asg1 +a\n\
 bsg1 +b\n\
end lists\n\
\n\
begin lists1 "2-step"\n\
 ph1 +x +x\n\
 ph2 +x +x\n\
 ph3 +x +x\n\
 ph4 +x -x\n\
 asg1 +a -a\n\
 bsg1 +b -a\n\
end lists1\n\
\n\
begin lists2 "4-step"\n\
 ph1 +x +x +x +x\n\
 ph2 +x +x +x +x\n\
 ph3 +x +x -x -x\n\
 ph4 +x -x +x -x\n\
 asg1 +a -a +a -a\n\
 bsg1 +b -b +b -b\n\
end lists2\n\
\n\
begin lists3 "8-step"\n\
 ph1 +x -x +x -x  +x -x +x -x\n\
 ph2 +x +x -x -x  +x +x -x -x\n\
 ph3 +x +x +x +x  +x +x +x +x\n\
 ph4 +x +x +x +x  -x -x -x -x\n\
 asg1 +a -a -a +a  -a +a +a -a\n\
 bsg1 +b -b -b +b  -b +b +b -b\n\
end lists3\n\
\n\
begin lists4 "16-step"\n\
 ph1 +x +x +x +x   +x +x +x +x   +x +x +x +x   +x +x +x +x\n\
 ph2 +x +y -x -y   +x +y -x -y   +x +y -x -y   +x +y -x -y\n\
 ph3 +x +y -x -y   +x +y -x -y   -x -y +x +y   -x -y +x +y\n\
 ph4 +x +y -x -y   -x -y +x +y   +x +y -x -y   -x -y +x +y\n\
 asg1 +a -a +a -a   -a +a -a +a   +a -a +a -a   -a +a -a +a\n\
 bsg1 +b -b +b -b   -b +b -b +b   +b -b +b -b   -b +b -b +b\n\
end lists4\n\
\n\
;\n\
; HYSCORE setup\n\
;\n\
begin exp "HYSCORE Setup" [TRANS QUAD lists2 lists1 lists3 lists4 lists]\n\
	; QUAD detection of Transient - 4-step phase default\n\
\n\
    shot i=1 to h		; accumulation loop\n\
      d9			; DAF (Delay After Flash)\n\
      p0 [ph1]			; 1st pulse and phase program\n\
      d1			; tau\n\
      p0 [ph2]			; 2nd pulse and phase program\n\
      d2			; t1\n\
      p2 [ph3]  		; inversion pulse\n\
      d3 			; t2\n\
      p0 [ph4]			; 3rd pulse\n\
      d1			; tau\n\
      d0			; constant acquisition delay\n\
      dig [sg1]			; acquisition\n\
    next i 			; end of accumulation loop\n\
\n\
end exp\n\
\n\
;\n\
; Inversion Optimization\n\
;\n\
begin exp1 "Inversion Optimization" [INTG QUAD lists2 lists1 lists3 lists4 lists]\n\
	; QUAD detection Integration - 4-step phase default\n\
p2=4\n\
d30=4				; 4 ns inversion pulse increment\n\
dx=p2				; define x-axis\n\
\n\
sweep x=1 to sx			; inversion pulse sweep loop\n\
  shot i=1 to h			; accumulation loop\n\
   d9				; DAF (Delay After Flash)\n\
   p0 [ph1]			; 1st pulse and phase program\n\
   d1				; tau\n\
   p0 [ph2]			; 2nd pulse and phase program\n\
   d2				; t1\n\
   p2 [ph3]  			; inversion pulse\n\
   d3 				; t2\n\
   p0 [ph4]			; 3rd pulse\n\
   d1				; tau\n\
   d0				; constant acquisition delay\n\
   acq [sg1]			; acquisition\n\
  next i 			; end of accumulation loop\n\
\n\
  p2=p2+d30			; increment inversion pulse\n\
  d3=d3+d30			; increment t2 to keep constant time between pulses\n\
  dx=dx+d30			; increment x-axis\n\
next x				; end of inversion pulse sweep loop\n\
\n\
end exp1\n\
\n\
;\n\
; HYSCORE\n\
;\n\
begin exp2 "HYSCORE" [INTG QUAD lists2 lists1 lists3 lists4 lists]\n\
	; QUAD detection Integration - 4-step phase cycle default\n\
\n\
for k=1 to n			; averaging loop\n\
 totscans(n)			; output of total number of scans\n\
\n\
dy=d2				; starting t1 time and assignment of y-axis\n\
\n\
  for y=1 to sy			; t1 loop\n\
\n\
   dx=d3			; starting t2 time and assignment of x-axis\n\
\n\
   sweep x=1 to sx		; t2 sweep loop\n\
    shot i=1 to h		; accumulation loop\n\
      d9			; DAF (Delay After Flash)\n\
      p0 [ph1]			; 1st pulse and phase program\n\
      d1			; tau\n\
      p0 [ph2]			; 2nd pulse and phase program\n\
      dy			; t1\n\
      p2 [ph3]			; inversion pulse and phase program\n\
      dx			; t2\n\
      p0 [ph4]			; 3rd pulse and phase program\n\
      d1			; tau\n\
      d0			; constant acquisition delay\n\
      acq [sg1]			; acquisition\n\
    next i			; end of accumulation loop\n\
\n\
    dx=dx+d30			; increment t2\n\
   next x			; end of t2 sweep loop\n\
\n\
   dy=dy+d31			; increment t1\n\
  next y			; end of t1 loop\n\
\n\
  scansdone(k)			; output of scans done\n\
 next k				; end of averaging loop\n\
\n\
end exp2\n\
\n\
\n
Psd1               0,0,308,460,3980,4020,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20,20,32,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd10              0,308,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd11              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd12              0,308,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd13              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd14              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd15              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd16              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd17              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd18              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd19              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd2               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd20              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd21              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd22              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd23              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd24              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd25              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd26              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd3               0,0,140,3812,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,440,328,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd4               0,0,8,3680,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,380,268,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd5               0,8,4324,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,52,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd6               0,308,308,460,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,16,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd7               0,300,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd8               0,308,3980,4020,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,28,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
Psd9               0,304,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
QuadDetect         True
RF1Prg             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
RF2Prg             0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000
ReplaceMode        Off
ShotRepTime        1020.00 us
ShotsPLoop         50
SweepsPExp         1
TriggerTimeOut     10 s
XAxisQuant         Time
XSpecRes           175
YAxisQuant         Magnetic Field
YSpecRes           175

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         1
NbScansDone        1
NbScansToDo        1
ReplaceMode        Off

.DVC     routeTrg, 1.0


.DVC     sigChanSmall, 1.0


.DVC     transRec, 1.0


*
************************************************************

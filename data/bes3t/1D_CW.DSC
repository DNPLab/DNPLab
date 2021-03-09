#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	2250
XMIN	3420.550000
XWID	149.933333
*
*	Documentational Text:
*
TITL	'1D_CW'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    11/13/20
TIME    11:12:15
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.349555
A1SW    0.015
MWFQ    9.804448e+09
MWPW    0.001002
AVGS    2
SPTP    0.02
RCAG    30
RCHM    1
B0MA    0.0001
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    2250
RCTC    0.01024
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     ezcw, 1.0

MetalIonsDef       False
OrgRadicalsDef     False

.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        3495.55 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     3420.550 G
SweepDirection     Up
SweepWidth         150.0 G
WidthTM            200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.804448 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        120.0
Power              1.002 mW
PowerAtten         23.0 dB
QValue             6500

.DVC     recorder, 1.0

BaselineCorr       Off
NbScansAcc         2
NbScansDone        2
NbScansToDo        2
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       0

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           9.271e-09
ConvTime           20.00 ms
DModDetectSCT      First
DualDetect         OFF
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               30 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             1.000 G
ModFreq            100.00 kHz
ModPhase           0.0
ModResol           15.000
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0.99652,0,0,0.00737177,0,0,-0.00559614,0,0,-2.88221e-05,0,0,1.00404e-05,0,0,3.43695e-08,0,0,-5.0404e-09,0,0,-1.4783e-11,0,0,-1.29132e-12,0
QuadMode           False
ResCenter          62.5 mm
ResLength          40 mm
Resolution         2250
Resonator          1
SctNorm            True
SctRevision        Allegro
SpuExtension       False
SweepTime          45.00 s
TimeConst          10.24 ms
TimeExp            False
TuneCaps           40

.DVC     tempMonit, 1.0

AcqWaitTime        5.0 s
Temperature        295.00 K
Tolerance          1.00 K

*
************************************************************

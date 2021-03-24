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
XPTS	2500
XMIN	3352.200000
XWID	149.940000
*
*	Documentational Text:
*
TITL	'apple_F091-1_before'
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
DATE    09/10/18
TIME    20:20:14
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.34272
A1SW    0.015
MWFQ    9.627173e+09
MWPW    0.0003162
AVGS    20
SPTP    0.00601
RCAG    40
RCHM    1
B0MA    0.0003
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    2500
RCTC    0.00256
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
CenterField        3427.20 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait stable
GFactor            2.007000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     3352.200 G
SweepDirection     Up
SweepWidth         150.0 G
WidthTM            200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.627173 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     On
AcqSliceFTuning    On
BridgeCalib        113.9
Power              0.3162 mW
PowerAtten         25.00 dB
QValue             3068

.DVC     recorder, 1.0

BaselineCorr       On
NbScansAcc         20
NbScansDone        20
NbScansToDo        20
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       0

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           1.21e-08
ConvTime           6.01 ms
DModDetectSCT      First
DualDetect         OFF
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               40 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             3.000 G
ModFreq            100.00 kHz
ModPhase           0.0
ModResol           50.000
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0.993287,0,0,0.00240032,0,0,-0.00984044,0,0,8.8985e-05,0,0,-4.75288e-06,0,0,-1.05279e-06,0,0,2.9393e-07,0,0,2.70257e-09,0,0,-7.36288e-10,0
QuadMode           False
ResCenter          125 mm
ResLength          25 mm
Resolution         2500
Resonator          1
SctNorm            True
SctRevision        Allegro
SpuExtension       False
SweepTime          15.02 s
TimeConst          2.56 ms
TimeExp            False
TuneCaps           1

*
************************************************************

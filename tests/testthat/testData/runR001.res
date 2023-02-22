Run Start Date:  04/02/2021 
Run Start Time: 17:29:17.86 
Run Stop Date:   04/02/2021 
Run Stop Time:  17:29:35.21 
  
************************************************************************************************************************ 
********************                            CONTROL STREAM                                      ******************** 
************************************************************************************************************************ 
  
$PROB    BOLUS_1CPT_VCL SINGLE DOSE FOCE INTER (120 Ind/2276 Obs) runR001
$INPUT   ID TIME DV=DROP LNDV=DV MDV AMT EVID DOSE VI CLI SSX IIX SD CMT
$DATA    BOLUS_1CPT.csv IGNORE=@ IGNORE (SD.EQ.0)
$SUBR    ADVAN1,TRANS2
$PK
         CL=THETA(1)*EXP(ETA(1))
         V=THETA(2)*EXP(ETA(2))
         S1=V
$ERROR   (ONLY OBSERVATIONS)
         IPRED=LOG(0.0025)
         W=THETA(3)
         IF (F.GT.0) IPRED=LOG(F)
         IRES=IPRED-DV
         IWRES=IRES/W
         Y=IPRED+W*ERR(1)
$THETA   (0,5,15) (0,90,200) (0,0.3,1)
;#       CL V RSV
$OMEGA   0.15 0.15
;#       IIV_CL IIV_V
$SIGMA   1 FIX
;#       RSV_FIX
$EST     NSIG=3 PRINT=5 MAX=9999 NOABORT POSTHOC METHOD=COND INTER NOOBT
$COV
$TABLE  ID TIME LNDV MDV AMT EVID DOSE VI CLI CL V ETA1 ETA2
        IPRED IRES IWRES CWRESI
        ONEHEADER NOPRINT FILE=runR001.csv   
   
************************************************************************************************************************ 
********************                            SYSTEM INFORMATION                                  ******************** 
************************************************************************************************************************ 
   
Operating system: 

Microsoft Windows [Version 10.0.19042.746]
   
Compiler: Intel(R) Parallel Studio XE 2015 Update 3 Composer Edition (package 208)  
   
Compiler settings: /Gs /nologo /nbs /w /fp:strict   
   
************************************************************************************************************************ 
********************                            NMTRAN MESSAGES                                     ******************** 
************************************************************************************************************************ 
   
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.
   
************************************************************************************************************************ 
********************                            RUNNING MPI                                         ******************** 
************************************************************************************************************************ 
   
mpiexec -hosts 1 localhost 1 -noprompt nonmem.exe runR001.ctl -licfile=c:\nm750\license\nonmem.lic -parafile=C:\MPI\node6.pnm : ^
-wdir c:\MPI\worker1 -hosts 1 localhost 1 -noprompt nonmem.exe runR001.ctl -licfile=c:\nm750\license\nonmem.lic -worker : ^
-wdir c:\MPI\worker2 -hosts 1 localhost 1 -noprompt nonmem.exe runR001.ctl -licfile=c:\nm750\license\nonmem.lic -worker : ^
-wdir c:\MPI\worker3 -hosts 1 localhost 1 -noprompt nonmem.exe runR001.ctl -licfile=c:\nm750\license\nonmem.lic -worker : ^
-wdir c:\MPI\worker4 -hosts 1 localhost 1 -noprompt nonmem.exe runR001.ctl -licfile=c:\nm750\license\nonmem.lic -worker : ^
-wdir c:\MPI\worker5 -hosts 1 localhost 1 -noprompt nonmem.exe runR001.ctl -licfile=c:\nm750\license\nonmem.lic -worker
   
   
************************************************************************************************************************ 
********************                            NONMEM EXECUTION                                    ******************** 
************************************************************************************************************************ 
   
License Registered to: Occams
Expiration Date:    14 JUN 2021
Current Date:        4 FEB 2021
Days until program expires : 130
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 BOLUS_1CPT_VCL SINGLE DOSE FOCE INTER (120 Ind/2276 Obs) runR001
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:     2400
 NO. OF DATA ITEMS IN DATA SET:  13
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  4
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   5   0   0   0  13   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME LNDV MDV AMT EVID DOSE VI CLI SSX IIX SD CMT
0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:
 CL V IPRED IRES IWRES
0FORMAT FOR DATA:
 (E4.0,2E7.0,E2.0,E7.0,E2.0,E7.0,E6.0,E5.0,2E3.0,2E2.0)

 TOT. NO. OF OBS RECS:     2276
 TOT. NO. OF INDIVIDUALS:      120
0LENGTH OF THETA:   3
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   2
0DEFAULT OMEGA BOUNDARY TEST OMITTED:   YES
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.5000E+01     0.1500E+02
  0.0000E+00     0.9000E+02     0.2000E+03
  0.0000E+00     0.3000E+00     0.1000E+01
0INITIAL ESTIMATE OF OMEGA:
 0.1500E+00
 0.0000E+00   0.1500E+00
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
0TABLES STEP OMITTED:    NO
 NO. OF TABLES:           1
 SEED NUMBER (SEED):    11456
 RANMETHOD:             3U
 MC SAMPLES (ESAMPLE):    300
 WRES SQUARE ROOT TYPE (WRESCHOL): EIGENVALUE
0-- TABLE   1 --
0RECORDS ONLY:    ALL
04 COLUMNS APPENDED:    YES
 PRINTED:                NO
 HEADERS:               ONE
 FILE TO BE FORWARDED:   NO
 FORMAT:                S1PE11.4
 IDFORMAT:
 LFORMAT:
 RFORMAT:
 FIXED_EFFECT_ETAS:
0USER-CHOSEN ITEMS:
 ID TIME LNDV MDV AMT EVID DOSE VI CLI CL V ETA1 ETA2 IPRED IRES IWRES CWRESI
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 ONE COMPARTMENT MODEL (ADVAN1)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   2
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         CENTRAL      ON         NO         YES        YES        YES
    2         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            3           *           *           *           *
    2            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   5
   COMPT. NO. DATA ITEM IS DATA ITEM NO.:   13

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 OTHERWISE, ERROR SUBROUTINE CALLED ONLY WITH OBSERVATION EVENTS.
 
 #PARA: PARAFILE=C:\MPI\node6.pnm, PROTOCOL=MPI, NODES= 6
 
1
 
 
 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction
 
 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            9999
 NO. OF SIG. FIGURES REQUIRED:            3
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): runR001.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE

 
 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI
 
 MONITORING OF SEARCH:

 
0ITERATION NO.:    0    OBJECTIVE VALUE:  -3465.72885682371        NO. OF FUNC. EVALS.:   5
 CUMULATIVE NO. OF FUNC. EVALS.:        5
 NPARAMETR:  5.0000E+00  9.0000E+01  3.0000E-01  1.5000E-01  1.5000E-01
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   2.3055E+02  2.2706E+02  1.5860E+03  4.8390E+01 -1.4221E+01
 
0ITERATION NO.:    5    OBJECTIVE VALUE:  -4111.52222574939        NO. OF FUNC. EVALS.:  32
 CUMULATIVE NO. OF FUNC. EVALS.:       37
 NPARAMETR:  3.8359E+00  6.7264E+01  1.9172E-01  7.1326E-02  1.2917E-01
 PARAMETER: -2.7516E-01 -3.7906E-01 -4.9154E-01 -2.7169E-01  2.5233E-02
 GRADIENT:  -5.8220E+01  6.2223E+00 -2.6800E+02 -6.0305E+00  6.4533E+01
 
0ITERATION NO.:   10    OBJECTIVE VALUE:  -4125.77900367604        NO. OF FUNC. EVALS.:  37
 CUMULATIVE NO. OF FUNC. EVALS.:       74
 NPARAMETR:  3.9417E+00  6.7383E+01  1.9938E-01  7.1528E-02  9.1955E-02
 PARAMETER: -2.3841E-01 -3.7641E-01 -4.4287E-01 -2.7028E-01 -1.4467E-01
 GRADIENT:  -2.2855E+01 -1.9924E+01 -2.8036E+00 -2.3039E-01 -9.7916E-01
 
0ITERATION NO.:   14    OBJECTIVE VALUE:  -4126.11481640313        NO. OF FUNC. EVALS.:  40
 CUMULATIVE NO. OF FUNC. EVALS.:      114
 NPARAMETR:  3.9801E+00  6.8219E+01  1.9947E-01  7.1556E-02  9.2158E-02
 PARAMETER: -2.2524E-01 -3.5774E-01 -4.4230E-01 -2.7008E-01 -1.4356E-01
 GRADIENT:  -2.2459E-01  2.5575E-01  1.3847E-01  1.1535E-01 -5.3138E-02
 
 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      114
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.
 
 ETABAR:         4.7657E-04 -1.5163E-03
 SE:             2.3984E-02  2.7193E-02
 N:                     120         120
 
 P VAL.:         9.8415E-01  9.5553E-01
 
 ETASHRINKSD(%)  1.3722E+00  1.4618E+00
 ETASHRINKVR(%)  2.7255E+00  2.9021E+00
 EBVSHRINKSD(%)  1.5800E+00  1.7629E+00
 EBVSHRINKVR(%)  3.1350E+00  3.4947E+00
 RELATIVEINF(%)  9.6812E+01  9.6453E+01
 EPSSHRINKSD(%)  5.2203E+00
 EPSSHRINKVR(%)  1.0168E+01
 
  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):         2276
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    4183.00820314767     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -4126.11481640313     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       56.8933867445448     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           240
  
 #TERE:
 Elapsed estimation  time in seconds:     0.83
 Elapsed covariance  time in seconds:     0.34
 Elapsed postprocess time in seconds:     0.08
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -4126.115       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3     
 
         3.98E+00  6.82E+01  1.99E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2     
 
 ETA1
+        7.16E-02
 
 ETA2
+        0.00E+00  9.22E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2     
 
 ETA1
+        2.67E-01
 
 ETA2
+        0.00E+00  3.04E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3     
 
         9.90E-02  1.93E+00  2.80E-03
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2     
 
 ETA1
+        9.13E-03
 
 ETA2
+       .........  1.22E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2     
 
 ETA1
+        1.71E-02
 
 ETA2
+       .........  2.01E-02
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11      OM12      OM22      SG11  
 
 TH 1
+        9.80E-03
 
 TH 2
+        2.08E-02  3.72E+00
 
 TH 3
+        3.86E-07  6.66E-04  7.87E-06
 
 OM11
+        1.81E-05  5.51E-04  1.68E-06  8.33E-05
 
 OM12
+       ......... ......... ......... ......... .........
 
 OM22
+       -9.96E-05  8.12E-04 -2.61E-06 -2.55E-06 .........  1.49E-04
 
 SG11
+       ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11      OM12      OM22      SG11  
 
 TH 1
+        9.90E-02
 
 TH 2
+        1.09E-01  1.93E+00
 
 TH 3
+        1.39E-03  1.23E-01  2.80E-03
 
 OM11
+        2.01E-02  3.13E-02  6.56E-02  9.13E-03
 
 OM12
+       ......... ......... ......... ......... .........
 
 OM22
+       -8.25E-02  3.45E-02 -7.62E-02 -2.29E-02 .........  1.22E-02
 
 SG11
+       ......... ......... ......... ......... ......... ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      OM11      OM12      OM22      SG11  
 
 TH 1
+        1.04E+02
 
 TH 2
+       -6.10E-01  2.77E-01
 
 TH 3
+        7.49E+01 -2.39E+01  1.30E+05
 
 OM11
+       -1.79E+01 -1.29E+00 -2.41E+03  1.21E+04
 
 OM12
+       ......... ......... ......... ......... .........
 
 OM22
+        7.40E+01 -2.36E+00  2.42E+03  1.60E+02 .........  6.83E+03
 
 SG11
+       ......... ......... ......... ......... ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.20
 #CPUT: Total CPU Time in Seconds,        7.266

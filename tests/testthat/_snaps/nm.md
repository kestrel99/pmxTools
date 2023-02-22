# read_nm

    Code
      read_nm("testData/runR001.xml")
    Message <simpleMessage>
      Reading testData/runR001.xml
    Output
      $start_datetime
      $start_datetime[[1]]
      [1] "2021-02-04T17:29:28.815"
      
      
      $control_stream
      $control_stream[[1]]
      [1] "\n$PROB    BOLUS_1CPT_VCL SINGLE DOSE FOCE INTER (120 Ind/2276 Obs) runR001\n$INPUT   ID TIME DV=DROP LNDV=DV MDV AMT EVID DOSE VI CLI SSX IIX SD CMT\n$DATA    BOLUS_1CPT.csv IGNORE=@ IGNORE (SD.EQ.0)\n$SUBR    ADVAN1,TRANS2\n$PK\n         CL=THETA(1)*EXP(ETA(1))\n         V=THETA(2)*EXP(ETA(2))\n         S1=V\n$ERROR   (ONLY OBSERVATIONS)\n         IPRED=LOG(0.0025)\n         W=THETA(3)\n         IF (F.GT.0) IPRED=LOG(F)\n         IRES=IPRED-DV\n         IWRES=IRES/W\n         Y=IPRED+W*ERR(1)\n$THETA   (0,5,15) (0,90,200) (0,0.3,1)\n;#       CL V RSV\n$OMEGA   0.15 0.15\n;#       IIV_CL IIV_V\n$SIGMA   1 FIX\n;#       RSV_FIX\n$EST     NSIG=3 PRINT=5 MAX=9999 NOABORT POSTHOC METHOD=COND INTER NOOBT\n$COV\n$TABLE  ID TIME LNDV MDV AMT EVID DOSE VI CLI CL V ETA1 ETA2\n        IPRED IRES IWRES CWRESI\n        ONEHEADER NOPRINT FILE=runR001.csv\n"
      
      
      $nmtran
      $nmtran[[1]]
      [1] "\n\n WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1\n\n (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.\n"
      
      
      $nonmem
      $nonmem$license_information
      $nonmem$license_information[[1]]
      [1] "\nLicense Registered to: Occams\nExpiration Date:    14 JUN 2021\nCurrent Date:        4 FEB 2021\nDays until program expires : 130\n"
      
      
      $nonmem$program_information
      $nonmem$program_information[[1]]
      [1] "\n1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0\n ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN\n CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,\n AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION\n PERFORMED BY NOUS INFOSYSTEMS.\n"
      
      
      $nonmem$problem
      $nonmem$problem$problem_title
      $nonmem$problem$problem_title[[1]]
      [1] "BOLUS_1CPT_VCL SINGLE DOSE FOCE INTER (120 Ind/2276 Obs) runR001"
      
      
      $nonmem$problem$problem_information
      $nonmem$problem$problem_information[[1]]
      [1] "\n\n PROBLEM NO.:         1\n BOLUS_1CPT_VCL SINGLE DOSE FOCE INTER (120 Ind/2276 Obs) runR001\n0DATA CHECKOUT RUN:              NO\n DATA SET LOCATED ON UNIT NO.:    2\n THIS UNIT TO BE REWOUND:        NO\n NO. OF DATA RECS IN DATA SET:     2400\n NO. OF DATA ITEMS IN DATA SET:  13\n ID DATA ITEM IS DATA ITEM NO.:   1\n DEP VARIABLE IS DATA ITEM NO.:   3\n MDV DATA ITEM IS DATA ITEM NO.:  4\n0INDICES PASSED TO SUBROUTINE PRED:\n   6   2   5   0   0   0  13   0   0   0   0\n0LABELS FOR DATA ITEMS:\n ID TIME LNDV MDV AMT EVID DOSE VI CLI SSX IIX SD CMT\n0(NONBLANK) LABELS FOR PRED-DEFINED ITEMS:\n CL V IPRED IRES IWRES\n0FORMAT FOR DATA:\n (E4.0,2E7.0,E2.0,E7.0,E2.0,E7.0,E6.0,E5.0,2E3.0,2E2.0)\n\n TOT. NO. OF OBS RECS:     2276\n TOT. NO. OF INDIVIDUALS:      120\n0LENGTH OF THETA:   3\n0DEFAULT THETA BOUNDARY TEST OMITTED:    NO\n0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   2\n0DEFAULT OMEGA BOUNDARY TEST OMITTED:   YES\n0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1\n0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO\n0INITIAL ESTIMATE OF THETA:\n LOWER BOUND    INITIAL EST    UPPER BOUND\n  0.0000E+00     0.5000E+01     0.1500E+02\n  0.0000E+00     0.9000E+02     0.2000E+03\n  0.0000E+00     0.3000E+00     0.1000E+01\n0INITIAL ESTIMATE OF OMEGA:\n 0.1500E+00\n 0.0000E+00   0.1500E+00\n0INITIAL ESTIMATE OF SIGMA:\n 0.1000E+01\n0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE\n0COVARIANCE STEP OMITTED:        NO\n EIGENVLS. PRINTED:              NO\n SPECIAL COMPUTATION:            NO\n COMPRESSED FORMAT:              NO\n GRADIENT METHOD USED:     NOSLOW\n SIGDIGITS ETAHAT (SIGLO):                  -1\n SIGDIGITS GRADIENTS (SIGL):                -1\n EXCLUDE COV FOR FOCE (NOFCOV):              NO\n Cholesky Transposition of R Matrix (CHOLROFF):0\n KNUTHSUMOFF:                                -1\n RESUME COV ANALYSIS (RESUME):               NO\n SIR SAMPLE SIZE (SIRSAMPLE):\n NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1\n PRECONDTIONING CYCLES (PRECOND):        0\n PRECONDTIONING TYPES (PRECONDS):        TOS\n FORCED PRECONDTIONING CYCLES (PFCOND):0\n PRECONDTIONING TYPE (PRETYPE):        0\n FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0\n SIMPLE POS. DEFINITE SETTING: (POSDEF):-1\n0TABLES STEP OMITTED:    NO\n NO. OF TABLES:           1\n SEED NUMBER (SEED):    11456\n RANMETHOD:             3U\n MC SAMPLES (ESAMPLE):    300\n WRES SQUARE ROOT TYPE (WRESCHOL): EIGENVALUE\n0-- TABLE   1 --\n0RECORDS ONLY:    ALL\n04 COLUMNS APPENDED:    YES\n PRINTED:                NO\n HEADERS:               ONE\n FILE TO BE FORWARDED:   NO\n FORMAT:                S1PE11.4\n IDFORMAT:\n LFORMAT:\n RFORMAT:\n FIXED_EFFECT_ETAS:\n0USER-CHOSEN ITEMS:\n ID TIME LNDV MDV AMT EVID DOSE VI CLI CL V ETA1 ETA2 IPRED IRES IWRES CWRESI\n1DOUBLE PRECISION PREDPP VERSION 7.5.0\n\n ONE COMPARTMENT MODEL (ADVAN1)\n0MAXIMUM NO. OF BASIC PK PARAMETERS:   2\n0BASIC PK PARAMETERS (AFTER TRANSLATION):\n   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1\n\n TRANSLATOR WILL CONVERT PARAMETERS\n CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)\n0COMPARTMENT ATTRIBUTES\n COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT\n                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.\n    1         CENTRAL      ON         NO         YES        YES        YES\n    2         OUTPUT       OFF        YES        NO         NO         NO\n1\n ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG\n COMPT. NO.                             INDICES\n              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB\n                         FRACTION    RATE        DURATION    LAG\n    1            3           *           *           *           *\n    2            *           -           -           -           -\n             - PARAMETER IS NOT ALLOWED FOR THIS MODEL\n             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;\n               WILL DEFAULT TO ONE IF APPLICABLE\n0DATA ITEM INDICES USED BY PRED ARE:\n   EVENT ID DATA ITEM IS DATA ITEM NO.:      6\n   TIME DATA ITEM IS DATA ITEM NO.:          2\n   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   5\n   COMPT. NO. DATA ITEM IS DATA ITEM NO.:   13\n\n0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.\n PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.\n0DURING SIMULATION, ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.\n OTHERWISE, ERROR SUBROUTINE CALLED ONLY WITH OBSERVATION EVENTS.\n"
      
      
      $nonmem$problem$problem_options
      list()
      attr(,"data_checkout_run")
      [1] "no"
      attr(,"data_unit")
      [1] "2"
      attr(,"data_rewind")
      [1] "no"
      attr(,"data_nrec")
      [1] "2400"
      attr(,"data_nitems")
      [1] "13"
      attr(,"data_id")
      [1] "1"
      attr(,"data_l2")
      [1] "0"
      attr(,"data_dv")
      [1] "3"
      attr(,"data_mdv")
      [1] "4"
      attr(,"data_mrg")
      [1] "0"
      attr(,"data_raw")
      [1] "0"
      attr(,"data_rpt")
      [1] "0"
      attr(,"data_sub_array1")
      [1] "0"
      attr(,"data_sub_array2")
      [1] "0"
      attr(,"data_sub_array3")
      [1] "0"
      attr(,"data_pred_indices")
      [1] "6,2,5,0,0,0,13,0,0,0,0"
      attr(,"data_format")
      [1] "(e4.0,2e7.0,e2.0,e7.0,e2.0,e7.0,e6.0,e5.0,2e3.0,2e2.0)"
      attr(,"data_nobs")
      [1] "2276"
      attr(,"data_nind")
      [1] "120"
      attr(,"data_mdv100")
      [1] "0"
      attr(,"nthetat")
      [1] "3"
      attr(,"theta_bound_test_omitted")
      [1] "no"
      attr(,"omega_diagdim")
      [1] "2"
      attr(,"omega_bound_test_omitted")
      [1] "no"
      attr(,"sigma_diagdim")
      [1] "1"
      attr(,"sigma_bound_test_omitted")
      [1] "no"
      attr(,"cov_omitted")
      [1] "no"
      attr(,"cov_matrix")
      [1] "rsr"
      attr(,"cov_eigen_print")
      [1] "no"
      attr(,"cov_special")
      [1] "no"
      attr(,"cov_compressed")
      [1] "no"
      attr(,"cov_slow_gradient")
      [1] "noslow"
      attr(,"cov_siglocov")
      [1] "-1"
      attr(,"cov_siglcov")
      [1] "-1"
      attr(,"cov_tol")
      [1] "-1"
      attr(,"cov_atol")
      [1] "-1"
      attr(,"cov_nofcov")
      [1] "no"
      attr(,"cov_cholroff")
      [1] "0"
      attr(,"cov_knuthsumoff")
      [1] "-1"
      attr(,"cov_resume")
      [1] "no"
      attr(,"cov_sirsample")
      [1] "BLANK"
      attr(,"cov_thbnd")
      [1] "1"
      attr(,"cov_precond")
      [1] "0"
      attr(,"cov_preconds")
      [1] "tos"
      attr(,"cov_pfcond")
      [1] "0"
      attr(,"cov_pretype")
      [1] "0"
      attr(,"cov_fposdef")
      [1] "0"
      attr(,"cov_posdef")
      [1] "-1"
      attr(,"tab_omitted")
      [1] "no"
      attr(,"tab_number")
      [1] "1"
      attr(,"tab_seed")
      [1] "11456"
      attr(,"tab_clockseed")
      [1] "11456"
      attr(,"tab_ranmethod")
      [1] "3u"
      attr(,"tab_esample")
      [1] "300"
      attr(,"tab_wres")
      [1] "eigenvalue"
      attr(,"tab_fixedetas_01")
      [1] "BLANK"
      attr(,"pred_advan")
      [1] "1"
      attr(,"pred_evid")
      [1] "6"
      attr(,"pred_time")
      [1] "2"
      attr(,"pred_amt")
      [1] "5"
      attr(,"pred_cmt")
      [1] "13"
      
      $nonmem$problem$theta_lb
      $nonmem$problem$theta_lb$val
      $nonmem$problem$theta_lb$val[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"name")
      [1] "1"
      
      $nonmem$problem$theta_lb$val
      $nonmem$problem$theta_lb$val[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"name")
      [1] "2"
      
      $nonmem$problem$theta_lb$val
      $nonmem$problem$theta_lb$val[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"name")
      [1] "3"
      
      
      $nonmem$problem$theta_in
      $nonmem$problem$theta_in$val
      $nonmem$problem$theta_in$val[[1]]
      [1] "5.00000000000000"
      
      attr(,"name")
      [1] "1"
      
      $nonmem$problem$theta_in$val
      $nonmem$problem$theta_in$val[[1]]
      [1] "90.0000000000000"
      
      attr(,"name")
      [1] "2"
      
      $nonmem$problem$theta_in$val
      $nonmem$problem$theta_in$val[[1]]
      [1] "0.300000000000000"
      
      attr(,"name")
      [1] "3"
      
      
      $nonmem$problem$theta_ub
      $nonmem$problem$theta_ub$val
      $nonmem$problem$theta_ub$val[[1]]
      [1] "15.0000000000000"
      
      attr(,"name")
      [1] "1"
      
      $nonmem$problem$theta_ub$val
      $nonmem$problem$theta_ub$val[[1]]
      [1] "200.000000000000"
      
      attr(,"name")
      [1] "2"
      
      $nonmem$problem$theta_ub$val
      $nonmem$problem$theta_ub$val[[1]]
      [1] "1.00000000000000"
      
      attr(,"name")
      [1] "3"
      
      
      $nonmem$problem$estimation
      $nonmem$problem$estimation$parallel_est
      list()
      attr(,"parafile")
      [1] "C:\\MPI\\node6.pnm"
      attr(,"protocol")
      [1] "MPI"
      attr(,"nodes")
      [1] "6"
      
      $nonmem$problem$estimation$table_series
      $nonmem$problem$estimation$table_series[[1]]
      [1] "1"
      
      
      $nonmem$problem$estimation$estimation_method
      $nonmem$problem$estimation$estimation_method[[1]]
      [1] "focei"
      
      
      $nonmem$problem$estimation$estimation_title
      $nonmem$problem$estimation$estimation_title[[1]]
      [1] "First Order Conditional Estimation with Interaction"
      
      
      $nonmem$problem$estimation$estimation_information
      $nonmem$problem$estimation$estimation_information[[1]]
      [1] "\n ESTIMATION STEP OMITTED:                 NO\n ANALYSIS TYPE:                           POPULATION\n NUMBER OF SADDLE POINT RESET ITERATIONS:      0\n GRADIENT METHOD USED:               NOSLOW\n CONDITIONAL ESTIMATES USED:              YES\n CENTERED ETA:                            NO\n EPS-ETA INTERACTION:                     YES\n LAPLACIAN OBJ. FUNC.:                    NO\n NO. OF FUNCT. EVALS. ALLOWED:            9999\n NO. OF SIG. FIGURES REQUIRED:            3\n INTERMEDIATE PRINTOUT:                   YES\n ESTIMATE OUTPUT TO MSF:                  NO\n ABORT WITH PRED EXIT CODE 1:             NO\n IND. OBJ. FUNC. VALUES SORTED:           NO\n NUMERICAL DERIVATIVE\n       FILE REQUEST (NUMDER):               NONE\n MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0\n ETA HESSIAN EVALUATION METHOD (ETADER):    0\n INITIAL ETA FOR MAP ESTIMATION (MCETA):    0\n SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100\n GRADIENT SIGDIGITS OF\n       FIXED EFFECTS PARAMETERS (SIGL):     100\n NOPRIOR SETTING (NOPRIOR):                 0\n NOCOV SETTING (NOCOV):                     OFF\n DERCONT SETTING (DERCONT):                 OFF\n FINAL ETA RE-EVALUATION (FNLETA):          1\n EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS\n       IN SHRINKAGE (ETASTYPE):             NO\n NON-INFL. ETA CORRECTION (NONINFETA):      0\n RAW OUTPUT FILE (FILE): runR001.ext\n EXCLUDE TITLE (NOTITLE):                   NO\n EXCLUDE COLUMN LABELS (NOLABEL):           NO\n FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5\n PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL\n KNUTHSUMOFF:                               0\n INCLUDE LNTWOPI:                           NO\n INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO\n INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO\n ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO\n EM OR BAYESIAN METHOD USED:                 NONE\n\n"
      
      
      $nonmem$problem$estimation$estimation_options
      list()
      attr(,"estim_omitted")
      [1] "no"
      attr(,"analysis_type")
      [1] "pop"
      attr(,"saddle_reset")
      [1] "0"
      attr(,"saddle_hess")
      [1] "0"
      attr(,"slow_gradient")
      [1] "noslow"
      attr(,"cond_estim")
      [1] "yes"
      attr(,"centered_eta")
      [1] "no"
      attr(,"epseta_interaction")
      [1] "yes"
      attr(,"laplace")
      [1] "no"
      attr(,"predflag")
      [1] "0"
      attr(,"maxfn")
      [1] "9999"
      attr(,"nsig")
      [1] "3"
      attr(,"msfo")
      [1] "no"
      attr(,"abort")
      [1] "no"
      attr(,"objsort")
      [1] "no"
      attr(,"numder")
      [1] "0"
      attr(,"optmap")
      [1] "0"
      attr(,"etader")
      [1] "0"
      attr(,"mceta")
      [1] "0"
      attr(,"siglo")
      [1] "100"
      attr(,"sigl")
      [1] "100"
      attr(,"noprior")
      [1] "0"
      attr(,"nocov")
      [1] "0"
      attr(,"dercont")
      [1] "0"
      attr(,"atol")
      [1] "0"
      attr(,"fnleta")
      [1] "1"
      attr(,"etastype")
      [1] "0"
      attr(,"noninfeta")
      [1] "0"
      attr(,"file")
      [1] "runr001.ext"
      attr(,"notitle")
      [1] "0"
      attr(,"nolabel")
      [1] "0"
      attr(,"format")
      [1] "s1pe12.5"
      attr(,"order")
      [1] "tsol"
      attr(,"knuthsumoff")
      [1] "0"
      attr(,"lntwopi")
      [1] "0"
      attr(,"priorc")
      [1] "0"
      attr(,"olntwopi")
      [1] "0"
      attr(,"ctype")
      [1] "0"
      
      $nonmem$problem$estimation$monitor
      $nonmem$problem$estimation$monitor$obj
      $nonmem$problem$estimation$monitor$obj[[1]]
      [1] "-3465.72885682371"
      
      attr(,"iteration")
      [1] "0"
      
      $nonmem$problem$estimation$monitor$obj
      $nonmem$problem$estimation$monitor$obj[[1]]
      [1] "-4111.52222574939"
      
      attr(,"iteration")
      [1] "5"
      
      $nonmem$problem$estimation$monitor$obj
      $nonmem$problem$estimation$monitor$obj[[1]]
      [1] "-4125.77900367604"
      
      attr(,"iteration")
      [1] "10"
      
      $nonmem$problem$estimation$monitor$obj
      $nonmem$problem$estimation$monitor$obj[[1]]
      [1] "-4126.11481640313"
      
      attr(,"iteration")
      [1] "14"
      
      
      $nonmem$problem$estimation$termination_status
      $nonmem$problem$estimation$termination_status[[1]]
      [1] "0"
      
      
      $nonmem$problem$estimation$termination_nfuncevals
      $nonmem$problem$estimation$termination_nfuncevals[[1]]
      [1] "114"
      
      
      $nonmem$problem$estimation$termination_sigdigits
      $nonmem$problem$estimation$termination_sigdigits[[1]]
      [1] "3.02298445301969"
      
      
      $nonmem$problem$estimation$termination_information
      $nonmem$problem$estimation$termination_information[[1]]
      [1] "\n0MINIMIZATION SUCCESSFUL\n NO. OF FUNCTION EVALUATIONS USED:      114\n NO. OF SIG. DIGITS IN FINAL EST.:  3.0\n"
      
      
      $nonmem$problem$estimation$termination_txtmsg
      $nonmem$problem$estimation$termination_txtmsg$val
      $nonmem$problem$estimation$termination_txtmsg$val[[1]]
      [1] "37"
      
      attr(,"name")
      [1] "1"
      
      
      $nonmem$problem$estimation$etabar
      $nonmem$problem$estimation$etabar$row
      $nonmem$problem$estimation$etabar$row$col
      $nonmem$problem$estimation$etabar$row$col[[1]]
      [1] "4.765711648433239E-004"
      
      attr(,"cname")
      [1] "ETA1"
      
      $nonmem$problem$estimation$etabar$row$col
      $nonmem$problem$estimation$etabar$row$col[[1]]
      [1] "-1.516261194686842E-003"
      
      attr(,"cname")
      [1] "ETA2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$etabarse
      $nonmem$problem$estimation$etabarse$row
      $nonmem$problem$estimation$etabarse$row$col
      $nonmem$problem$estimation$etabarse$row$col[[1]]
      [1] "2.398354777032987E-002"
      
      attr(,"cname")
      [1] "ETA1"
      
      $nonmem$problem$estimation$etabarse$row$col
      $nonmem$problem$estimation$etabarse$row$col[[1]]
      [1] "2.719347789006172E-002"
      
      attr(,"cname")
      [1] "ETA2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$etabarn
      $nonmem$problem$estimation$etabarn$row
      $nonmem$problem$estimation$etabarn$row$col
      $nonmem$problem$estimation$etabarn$row$col[[1]]
      [1] "120"
      
      attr(,"cname")
      [1] "ETA1"
      
      $nonmem$problem$estimation$etabarn$row$col
      $nonmem$problem$estimation$etabarn$row$col[[1]]
      [1] "120"
      
      attr(,"cname")
      [1] "ETA2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$etabarpval
      $nonmem$problem$estimation$etabarpval$row
      $nonmem$problem$estimation$etabarpval$row$col
      $nonmem$problem$estimation$etabarpval$row$col[[1]]
      [1] "0.984146393573823"
      
      attr(,"cname")
      [1] "ETA1"
      
      $nonmem$problem$estimation$etabarpval$row$col
      $nonmem$problem$estimation$etabarpval$row$col[[1]]
      [1] "0.955534240611886"
      
      attr(,"cname")
      [1] "ETA2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$etashrinksd
      $nonmem$problem$estimation$etashrinksd$row
      $nonmem$problem$estimation$etashrinksd$row$col
      $nonmem$problem$estimation$etashrinksd$row$col[[1]]
      [1] "1.37215788183357"
      
      attr(,"cname")
      [1] "ETA1"
      
      $nonmem$problem$estimation$etashrinksd$row$col
      $nonmem$problem$estimation$etashrinksd$row$col[[1]]
      [1] "1.46175485792922"
      
      attr(,"cname")
      [1] "ETA2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$etashrinkvr
      $nonmem$problem$estimation$etashrinkvr$row
      $nonmem$problem$estimation$etashrinkvr$row$col
      $nonmem$problem$estimation$etashrinkvr$row$col[[1]]
      [1] "2.72548759114038"
      
      attr(,"cname")
      [1] "ETA1"
      
      $nonmem$problem$estimation$etashrinkvr$row$col
      $nonmem$problem$estimation$etashrinkvr$row$col[[1]]
      [1] "2.90214244321166"
      
      attr(,"cname")
      [1] "ETA2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$ebvshrinksd
      $nonmem$problem$estimation$ebvshrinksd$row
      $nonmem$problem$estimation$ebvshrinksd$row$col
      $nonmem$problem$estimation$ebvshrinksd$row$col[[1]]
      [1] "1.58000136432566"
      
      attr(,"cname")
      [1] "EBV1"
      
      $nonmem$problem$estimation$ebvshrinksd$row$col
      $nonmem$problem$estimation$ebvshrinksd$row$col[[1]]
      [1] "1.76289476258970"
      
      attr(,"cname")
      [1] "EBV2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$ebvshrinkvr
      $nonmem$problem$estimation$ebvshrinkvr$row
      $nonmem$problem$estimation$ebvshrinkvr$row$col
      $nonmem$problem$estimation$ebvshrinkvr$row$col[[1]]
      [1] "3.13503868553859"
      
      attr(,"cname")
      [1] "EBV1"
      
      $nonmem$problem$estimation$ebvshrinkvr$row$col
      $nonmem$problem$estimation$ebvshrinkvr$row$col[[1]]
      [1] "3.49471154573972"
      
      attr(,"cname")
      [1] "EBV2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$relativeinf
      $nonmem$problem$estimation$relativeinf$row
      $nonmem$problem$estimation$relativeinf$row$col
      $nonmem$problem$estimation$relativeinf$row$col[[1]]
      [1] "96.8120802973034"
      
      attr(,"cname")
      [1] "REL1"
      
      $nonmem$problem$estimation$relativeinf$row$col
      $nonmem$problem$estimation$relativeinf$row$col[[1]]
      [1] "96.4526037915573"
      
      attr(,"cname")
      [1] "REL2"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$epsshrinksd
      $nonmem$problem$estimation$epsshrinksd$row
      $nonmem$problem$estimation$epsshrinksd$row$col
      $nonmem$problem$estimation$epsshrinksd$row$col[[1]]
      [1] "5.22027080681049"
      
      attr(,"cname")
      [1] "EPS1"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$epsshrinkvr
      $nonmem$problem$estimation$epsshrinkvr$row
      $nonmem$problem$estimation$epsshrinkvr$row$col
      $nonmem$problem$estimation$epsshrinkvr$row$col[[1]]
      [1] "10.1680293406566"
      
      attr(,"cname")
      [1] "EPS1"
      
      attr(,"rname")
      [1] "SUBPOP1"
      
      
      $nonmem$problem$estimation$estimation_elapsed_time
      $nonmem$problem$estimation$estimation_elapsed_time[[1]]
      [1] "0.83"
      
      
      $nonmem$problem$estimation$parallel_cov
      list()
      attr(,"parafile")
      [1] "C:\\MPI\\node6.pnm"
      attr(,"protocol")
      [1] "MPI"
      attr(,"nodes")
      [1] "6"
      
      $nonmem$problem$estimation$covariance_information
      $nonmem$problem$estimation$covariance_information[[1]]
      [1] "\n"
      
      
      $nonmem$problem$estimation$covariance_status
      list()
      attr(,"error")
      [1] "0"
      attr(,"numnegeigenvalues")
      [1] "-1"
      attr(,"mineigenvalue")
      [1] "0.000000000000000E+000"
      attr(,"maxeigenvalue")
      [1] "0.000000000000000E+000"
      attr(,"rms")
      [1] "0.000000000000000E+000"
      
      $nonmem$problem$estimation$covariance_elapsed_time
      $nonmem$problem$estimation$covariance_elapsed_time[[1]]
      [1] "0.34"
      
      
      $nonmem$problem$estimation$parallel_fnleta
      list()
      attr(,"parafile")
      [1] "C:\\MPI\\node6.pnm"
      attr(,"protocol")
      [1] "MPI"
      attr(,"nodes")
      [1] "6"
      
      $nonmem$problem$estimation$parallel_wres
      list()
      attr(,"parafile")
      [1] "C:\\MPI\\node6.pnm"
      attr(,"protocol")
      [1] "MPI"
      attr(,"nodes")
      [1] "6"
      
      $nonmem$problem$estimation$final_objective_function_text
      $nonmem$problem$estimation$final_objective_function_text[[1]]
      [1] "MINIMUM VALUE OF OBJECTIVE FUNCTION"
      
      
      $nonmem$problem$estimation$final_objective_function
      $nonmem$problem$estimation$final_objective_function[[1]]
      [1] "-4126.11481640313"
      
      
      $nonmem$problem$estimation$theta
      $nonmem$problem$estimation$theta$val
      $nonmem$problem$estimation$theta$val[[1]]
      [1] "3.98013225640529"
      
      attr(,"name")
      [1] "1"
      
      $nonmem$problem$estimation$theta$val
      $nonmem$problem$estimation$theta$val[[1]]
      [1] "68.2193506477679"
      
      attr(,"name")
      [1] "2"
      
      $nonmem$problem$estimation$theta$val
      $nonmem$problem$estimation$theta$val[[1]]
      [1] "0.199472175944961"
      
      attr(,"name")
      [1] "3"
      
      
      $nonmem$problem$estimation$omega
      $nonmem$problem$estimation$omega$row
      $nonmem$problem$estimation$omega$row$col
      $nonmem$problem$estimation$omega$row$col[[1]]
      [1] "7.155554952720695E-002"
      
      attr(,"cname")
      [1] "1"
      
      attr(,"rname")
      [1] "1"
      
      $nonmem$problem$estimation$omega$row
      $nonmem$problem$estimation$omega$row$col
      $nonmem$problem$estimation$omega$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "1"
      
      $nonmem$problem$estimation$omega$row$col
      $nonmem$problem$estimation$omega$row$col[[1]]
      [1] "9.215849908940954E-002"
      
      attr(,"cname")
      [1] "2"
      
      attr(,"rname")
      [1] "2"
      
      
      $nonmem$problem$estimation$sigma
      $nonmem$problem$estimation$sigma$row
      $nonmem$problem$estimation$sigma$row$col
      $nonmem$problem$estimation$sigma$row$col[[1]]
      [1] "1.00000000000000"
      
      attr(,"cname")
      [1] "1"
      
      attr(,"rname")
      [1] "1"
      
      
      $nonmem$problem$estimation$omegac
      $nonmem$problem$estimation$omegac$row
      $nonmem$problem$estimation$omegac$row$col
      $nonmem$problem$estimation$omegac$row$col[[1]]
      [1] "0.267498690701856"
      
      attr(,"cname")
      [1] "1"
      
      attr(,"rname")
      [1] "1"
      
      $nonmem$problem$estimation$omegac$row
      $nonmem$problem$estimation$omegac$row$col
      $nonmem$problem$estimation$omegac$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "1"
      
      $nonmem$problem$estimation$omegac$row$col
      $nonmem$problem$estimation$omegac$row$col[[1]]
      [1] "0.303576183336917"
      
      attr(,"cname")
      [1] "2"
      
      attr(,"rname")
      [1] "2"
      
      
      $nonmem$problem$estimation$sigmac
      $nonmem$problem$estimation$sigmac$row
      $nonmem$problem$estimation$sigmac$row$col
      $nonmem$problem$estimation$sigmac$row$col[[1]]
      [1] "1.00000000000000"
      
      attr(,"cname")
      [1] "1"
      
      attr(,"rname")
      [1] "1"
      
      
      $nonmem$problem$estimation$thetase
      $nonmem$problem$estimation$thetase$val
      $nonmem$problem$estimation$thetase$val[[1]]
      [1] "9.897592474768040E-002"
      
      attr(,"name")
      [1] "1"
      
      $nonmem$problem$estimation$thetase$val
      $nonmem$problem$estimation$thetase$val[[1]]
      [1] "1.92948784532265"
      
      attr(,"name")
      [1] "2"
      
      $nonmem$problem$estimation$thetase$val
      $nonmem$problem$estimation$thetase$val[[1]]
      [1] "2.804748734272427E-003"
      
      attr(,"name")
      [1] "3"
      
      
      $nonmem$problem$estimation$omegase
      $nonmem$problem$estimation$omegase$row
      $nonmem$problem$estimation$omegase$row$col
      $nonmem$problem$estimation$omegase$row$col[[1]]
      [1] "9.128170885017182E-003"
      
      attr(,"cname")
      [1] "1"
      
      attr(,"rname")
      [1] "1"
      
      $nonmem$problem$estimation$omegase$row
      $nonmem$problem$estimation$omegase$row$col
      $nonmem$problem$estimation$omegase$row$col[[1]]
      [1] "10000000000.0000"
      
      attr(,"cname")
      [1] "1"
      
      $nonmem$problem$estimation$omegase$row$col
      $nonmem$problem$estimation$omegase$row$col[[1]]
      [1] "1.220086738410604E-002"
      
      attr(,"cname")
      [1] "2"
      
      attr(,"rname")
      [1] "2"
      
      
      $nonmem$problem$estimation$sigmase
      $nonmem$problem$estimation$sigmase$row
      $nonmem$problem$estimation$sigmase$row$col
      $nonmem$problem$estimation$sigmase$row$col[[1]]
      [1] "10000000000.0000"
      
      attr(,"cname")
      [1] "1"
      
      attr(,"rname")
      [1] "1"
      
      
      $nonmem$problem$estimation$omegacse
      $nonmem$problem$estimation$omegacse$row
      $nonmem$problem$estimation$omegacse$row$col
      $nonmem$problem$estimation$omegacse$row$col[[1]]
      [1] "1.706208516585062E-002"
      
      attr(,"cname")
      [1] "1"
      
      attr(,"rname")
      [1] "1"
      
      $nonmem$problem$estimation$omegacse$row
      $nonmem$problem$estimation$omegacse$row$col
      $nonmem$problem$estimation$omegacse$row$col[[1]]
      [1] "10000000000.0000"
      
      attr(,"cname")
      [1] "1"
      
      $nonmem$problem$estimation$omegacse$row$col
      $nonmem$problem$estimation$omegacse$row$col[[1]]
      [1] "2.009523153297765E-002"
      
      attr(,"cname")
      [1] "2"
      
      attr(,"rname")
      [1] "2"
      
      
      $nonmem$problem$estimation$sigmacse
      $nonmem$problem$estimation$sigmacse$row
      $nonmem$problem$estimation$sigmacse$row$col
      $nonmem$problem$estimation$sigmacse$row$col[[1]]
      [1] "10000000000.0000"
      
      attr(,"cname")
      [1] "1"
      
      attr(,"rname")
      [1] "1"
      
      
      $nonmem$problem$estimation$covariance
      $nonmem$problem$estimation$covariance$row
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "9.796233679658492E-003"
      
      attr(,"cname")
      [1] "THETA1"
      
      attr(,"rname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$covariance$row
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "2.082984482340075E-002"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "3.72292334524783"
      
      attr(,"cname")
      [1] "THETA2"
      
      attr(,"rname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$covariance$row
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "3.857680188527142E-007"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "6.662205065331063E-004"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "7.866615462402783E-006"
      
      attr(,"cname")
      [1] "THETA3"
      
      attr(,"rname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$covariance$row
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "1.813035529535678E-005"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "5.510382599179932E-004"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "1.678782439694377E-006"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "8.332350370607537E-005"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      attr(,"rname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$covariance$row
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      attr(,"rname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$covariance$row
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "-9.964252462000372E-005"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "8.115456604809203E-004"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "-2.606638906468023E-006"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "-2.550978499775127E-006"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "1.488611649245426E-004"
      
      attr(,"cname")
      [1] "OMEGA(2,2)"
      
      attr(,"rname")
      [1] "OMEGA(2,2)"
      
      $nonmem$problem$estimation$covariance$row
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(2,2)"
      
      $nonmem$problem$estimation$covariance$row$col
      $nonmem$problem$estimation$covariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "SIGMA(1,1)"
      
      attr(,"rname")
      [1] "SIGMA(1,1)"
      
      
      $nonmem$problem$estimation$correlation
      $nonmem$problem$estimation$correlation$row
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "9.897592474768040E-002"
      
      attr(,"cname")
      [1] "THETA1"
      
      attr(,"rname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$correlation$row
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.109072286986378"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "1.92948784532265"
      
      attr(,"cname")
      [1] "THETA2"
      
      attr(,"rname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$correlation$row
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "1.389641232920324E-003"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.123106784682911"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "2.804748734272427E-003"
      
      attr(,"cname")
      [1] "THETA3"
      
      attr(,"rname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$correlation$row
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "2.006748675708277E-002"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "3.128642533299852E-002"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "6.557173805929989E-002"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "9.128170885017182E-003"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      attr(,"rname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$correlation$row
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      attr(,"rname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$correlation$row
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "-8.251339336697662E-002"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "3.447308943484225E-002"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "-7.617214858069581E-002"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "-2.290510726407609E-002"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "1.220086738410604E-002"
      
      attr(,"cname")
      [1] "OMEGA(2,2)"
      
      attr(,"rname")
      [1] "OMEGA(2,2)"
      
      $nonmem$problem$estimation$correlation$row
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(2,2)"
      
      $nonmem$problem$estimation$correlation$row$col
      $nonmem$problem$estimation$correlation$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "SIGMA(1,1)"
      
      attr(,"rname")
      [1] "SIGMA(1,1)"
      
      
      $nonmem$problem$estimation$invcovariance
      $nonmem$problem$estimation$invcovariance$row
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "104.159725745341"
      
      attr(,"cname")
      [1] "THETA1"
      
      attr(,"rname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$invcovariance$row
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "-0.609671502471118"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.277005825584178"
      
      attr(,"cname")
      [1] "THETA2"
      
      attr(,"rname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$invcovariance$row
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "74.8759828021359"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "-23.9363102159962"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "130460.861422742"
      
      attr(,"cname")
      [1] "THETA3"
      
      attr(,"rname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$invcovariance$row
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "-17.8737363625336"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "-1.28922038476296"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "-2412.28773609416"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "12067.3196432856"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      attr(,"rname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$invcovariance$row
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      attr(,"rname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$invcovariance$row
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "74.0494905940748"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "-2.35947484707118"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "2423.71409515852"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "159.617096229120"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "6825.27392229574"
      
      attr(,"cname")
      [1] "OMEGA(2,2)"
      
      attr(,"rname")
      [1] "OMEGA(2,2)"
      
      $nonmem$problem$estimation$invcovariance$row
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA1"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA2"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "THETA3"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,1)"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(1,2)"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "OMEGA(2,2)"
      
      $nonmem$problem$estimation$invcovariance$row$col
      $nonmem$problem$estimation$invcovariance$row$col[[1]]
      [1] "0.000000000000000E+000"
      
      attr(,"cname")
      [1] "SIGMA(1,1)"
      
      attr(,"rname")
      [1] "SIGMA(1,1)"
      
      
      attr(,"number")
      [1] "1"
      attr(,"type")
      [1] "0"
      
      $nonmem$problem$post_process_times
      $nonmem$problem$post_process_times$post_elapsed_time
      $nonmem$problem$post_process_times$post_elapsed_time[[1]]
      [1] "0.08"
      
      
      $nonmem$problem$post_process_times$finaloutput_elapsed_time
      $nonmem$problem$post_process_times$finaloutput_elapsed_time[[1]]
      [1] "0.20"
      
      
      
      attr(,"number")
      [1] "1"
      attr(,"subproblem")
      [1] "0"
      attr(,"superproblem1")
      [1] "0"
      attr(,"iteration1")
      [1] "0"
      attr(,"superproblem2")
      [1] "0"
      attr(,"iteration2")
      [1] "0"
      
      attr(,"version")
      [1] "7.5.0"
      
      $stop_datetime
      $stop_datetime[[1]]
      [1] "2021-02-04T17:29:31.574"
      
      
      $total_cputime
      $total_cputime[[1]]
      [1] "7.266"
      
      
      attr(,"schemaLocation")
      [1] "http://namespaces.oreilly.com/xmlnut/address output.xsd"
      attr(,"xmlns:xsi")
      [1] "http://www.w3.org/2001/XMLSchema-instance"
      attr(,"xmlns:nm")
      [1] "http://namespaces.oreilly.com/xmlnut/address"

# rnm

    Code
      rnm(index = "001", prefix = "runR", pathNM = "testData", ext = ".res", extmod = ".ctl")
    Output
      $Theta
        Parameter  Estimate         SE       CV        CIL        CIU P
      1    THETA1  3.980130 0.09897590 2.486750  3.7861408  4.1741192 0
      2    THETA2 68.219400 1.92949000 2.828360 64.4376691 72.0011309 0
      3    THETA3  0.199472 0.00280475 1.406087  0.1939748  0.2049692 0
      
      $Eta
        Eta  Estimate         SE       CV    EtaCV Shrinkage
      1   1 0.0715555 0.00912817 12.75677 26.74986   1.37216
      2   2 0.0921585 0.01220090 13.23904 30.35762   1.46175
      
      $Epsilon
        Epsilon Estimate SE CV EpsCV Shrinkage
      1       1        1  0  0   100   5.22027
      
      $CorTheta
            THETA1   THETA2     THETA3
      1 1.00000000 0.109072 0.00138964
      2 0.10907200 1.000000 0.12310700
      3 0.00138964 0.123107 1.00000000
      
      $CorOmega
            Eta1  Eta2
      Eta1 1.000 0.000
      Eta2 0.000 1.000
      
      $CorSigma
            Eps1
      Eps1 1.000
      
      $OmegaMatrix
                [,1]      [,2]
      [1,] 0.0715555 0.0000000
      [2,] 0.0000000 0.0921585
      
      $SigmaMatrix
           [,1]
      [1,]    1
      
      $CovMatrixTheta
             THETA1      THETA2      THETA3
      1 9.79623e-03 0.020829800 3.85768e-07
      2 2.08298e-02 3.722920000 6.66221e-04
      3 3.85768e-07 0.000666221 7.86662e-06
      
      $CovMatrix
              THETA1      THETA2       THETA3 SIGMA.1.1.   OMEGA.1.1. OMEGA.2.1.
      1  9.79623e-03 0.020829800  3.85768e-07          0  1.81304e-05          0
      2  2.08298e-02 3.722920000  6.66221e-04          0  5.51038e-04          0
      3  3.85768e-07 0.000666221  7.86662e-06          0  1.67878e-06          0
      4  0.00000e+00 0.000000000  0.00000e+00          0  0.00000e+00          0
      5  1.81304e-05 0.000551038  1.67878e-06          0  8.33235e-05          0
      6  0.00000e+00 0.000000000  0.00000e+00          0  0.00000e+00          0
      7 -9.96425e-05 0.000811546 -2.60664e-06          0 -2.55098e-06          0
          OMEGA.2.2.
      1 -9.96425e-05
      2  8.11546e-04
      3 -2.60664e-06
      4  0.00000e+00
      5 -2.55098e-06
      6  0.00000e+00
      7  1.48861e-04
      
      $OFV
      [1] -4126.115
      
      $ThetaString
        Parameter Estimate      SE   CV   Estimate (95% CI)         P
      1    THETA1     3.98  0.0990 2.5%    3.98 (3.79/4.17) < 0.00001
      2    THETA2     68.2    1.93 2.8%    68.2 (64.4/72.0) < 0.00001
      3    THETA3    0.199 0.00280 1.4% 0.199 (0.194/0.205) < 0.00001
      
      $EtaString
        Eta Estimate      SE    CV 100*sqrt(Omega) Shrinkage
      1   1   0.0716 0.00913 12.8%           26.7%      1.4%
      2   2   0.0922  0.0122 13.2%           30.4%      1.5%
      
      $EpsString
        Epsilon Estimate      SE CV 100*sqrt(Omega) Shrinkage
      1       1     1.00 (Fixed)                         5.2%
      
      $RunTime
      [1] 17.35
      
      $conditionN
      NULL
      

# read_nmcov

    Code
      read_nmcov("testData/runR001.cov")
    Message <simpleMessage>
      Reading testData/runR001.cov
    Output
                 THETA1      THETA2       THETA3 SIGMA.1.1.   OMEGA.1.1. OMEGA.2.1.
      [1,]  9.79623e-03 0.020829800  3.85768e-07          0  1.81304e-05          0
      [2,]  2.08298e-02 3.722920000  6.66221e-04          0  5.51038e-04          0
      [3,]  3.85768e-07 0.000666221  7.86662e-06          0  1.67878e-06          0
      [4,]  0.00000e+00 0.000000000  0.00000e+00          0  0.00000e+00          0
      [5,]  1.81304e-05 0.000551038  1.67878e-06          0  8.33235e-05          0
      [6,]  0.00000e+00 0.000000000  0.00000e+00          0  0.00000e+00          0
      [7,] -9.96425e-05 0.000811546 -2.60664e-06          0 -2.55098e-06          0
             OMEGA.2.2.
      [1,] -9.96425e-05
      [2,]  8.11546e-04
      [3,] -2.60664e-06
      [4,]  0.00000e+00
      [5,] -2.55098e-06
      [6,]  0.00000e+00
      [7,]  1.48861e-04

# read_nmext

    Code
      read_nmext("testData/runR001.ext")
    Message <simpleMessage>
      Reading testData/runR001.ext
    Output
      $Termination
      NULL
      
      $OFV
      [1] -4126.115
      
      $Thetas
      [1]  3.980130 68.219400  0.199472
      
      $Omega
      $Omega[[1]]
      [1] 0.0715555
      
      $Omega[[2]]
      [1] 0.0000000 0.0921585
      
      
      $Sigma
      $Sigma[[1]]
      [1] 1
      
      
      $seThetas
      [1] 0.09897590 1.92949000 0.00280475
      
      $seOmega
      $seOmega[[1]]
      [1] 0.00912817
      
      $seOmega[[2]]
      [1]        NA 0.0122009
      
      
      $seSigma
      $seSigma[[1]]
      [1] NA
      
      
      $raw_lst
       [1] "TABLE NO.     1: First Order Conditional Estimation with Interaction: Goal Function=MINIMUM VALUE OF OBJECTIVE FUNCTION: Problem=1 Subproblem=0 Superproblem1=0 Iteration1=0 Superproblem2=0 Iteration2=0"
       [2] " ITERATION    THETA1       THETA2       THETA3       SIGMA(1,1)   OMEGA(1,1)   OMEGA(2,1)   OMEGA(2,2)   OBJ"                                                                                             
       [3] "            0  5.00000E+00  9.00000E+01  3.00000E-01  1.00000E+00  1.50000E-01  0.00000E+00  1.50000E-01   -3465.72885682371"                                                                             
       [4] "            5  3.83588E+00  6.72641E+01  1.91725E-01  1.00000E+00  7.13261E-02  0.00000E+00  1.29166E-01   -4111.52222574939"                                                                             
       [5] "           10  3.94172E+00  6.73826E+01  1.99381E-01  1.00000E+00  7.15277E-02  0.00000E+00  9.19553E-02   -4125.77900367604"                                                                             
       [6] "           14  3.98013E+00  6.82194E+01  1.99472E-01  1.00000E+00  7.15555E-02  0.00000E+00  9.21585E-02   -4126.11481640313"                                                                             
       [7] "  -1000000000  3.98013E+00  6.82194E+01  1.99472E-01  1.00000E+00  7.15555E-02  0.00000E+00  9.21585E-02   -4126.11481640313"                                                                             
       [8] "  -1000000001  9.89759E-02  1.92949E+00  2.80475E-03  1.00000E+10  9.12817E-03  1.00000E+10  1.22009E-02   0.000000000000000E+000"                                                                        
       [9] "  -1000000004  0.00000E+00  0.00000E+00  0.00000E+00  1.00000E+00  2.67499E-01  0.00000E+00  3.03576E-01   0.000000000000000E+000"                                                                        
      [10] "  -1000000005  0.00000E+00  0.00000E+00  0.00000E+00  1.00000E+10  1.70621E-02  1.00000E+10  2.00952E-02   0.000000000000000E+000"                                                                        
      [11] "  -1000000006  0.00000E+00  0.00000E+00  0.00000E+00  1.00000E+00  0.00000E+00  1.00000E+00  0.00000E+00   0.000000000000000E+000"                                                                        
      [12] "  -1000000007  0.00000E+00  3.70000E+01  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00   0.000000000000000E+000"                                                                        
      [13] "  -1000000008  3.88104E-02 -2.92122E-03 -4.33821E-01  0.00000E+00 -4.03024E-01  0.00000E+00  1.44149E-01   0.000000000000000E+000"                                                                        
      

# get_est_table

    Code
      get_est_table(nm1, thetaLabels = c("CL", "V", "KA"), omegaLabels = c("CL", "V"),
      sigmaLabels = "Residual")
    Output
        Parameter Estimate   RSE          CI95 Shrinkage
      1        CL   3.9800  2.49     3.79-4.17         -
      2         V  68.2000  2.83     64.4-72.0         -
      3        KA   0.1990  1.41   0.194-0.205         -
      4        CL   0.0716 12.80 0.0537-0.0894      1.37
      5     OM2,1   0.0000    NA          <NA>         -
      6         V   0.0922 13.20  0.0682-0.116      1.46
      7  Residual   1.0000    NA           ---      5.22

# get_probinfo

    Code
      get_probinfo(nm1)
    Output
      $prob_title
      [1] "BOLUS_1CPT_VCL SINGLE DOSE FOCE INTER (120 Ind/2276 Obs) runR001"
      
      $prog_info
      [1] "\n1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0\n ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN\n CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,\n AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION\n PERFORMED BY NOUS INFOSYSTEMS.\n"
      
      $lic_info
      [1] "\nLicense Registered to: Occams\nExpiration Date:    14 JUN 2021\nCurrent Date:        4 FEB 2021\nDays until program expires : 130\n"
      
      $prob_times
                      startTime                stopTime elapsedTime
      1 2021-02-04T17:29:31.574 2021-02-04T17:29:31.574       7.266
      
      $prob_info
        nRecords nOBS nIND nTHETA nOMEGA nSIGMA COV_STEP_OMITTED COV_MATRIX
      1     2400 2276  120      3     NA      1               NO        rsr
      
      $est_info
        final_objective_function estimation_method
      1                 -4126.11             focei
                                           estimation_title estimation_elapsed_time
      1 First Order Conditional Estimation with Interaction                    0.83
        termination_status condition_number
      1         SUCCESSFUL               NA
      

# sample_uncertainty

    Code
      sample_uncert("testData/runR001.xml", n = 50, seed = 740727)
    Message <simpleMessage>
      Reading testData/runR001.xml
      Reading testData/runR001.cov
    Output
               V1       V2        V3 V4         V5 V6         V7
      1  4.165433 69.39547 0.2027448  1 0.06848604  0 0.09347681
      2  3.980225 69.50231 0.1994999  1 0.07953186  0 0.08515790
      3  4.008908 67.31330 0.2020217  1 0.04367090  0 0.08225640
      4  4.001296 68.48560 0.2047266  1 0.07295125  0 0.09326988
      5  3.978739 67.59651 0.1971732  1 0.08556432  0 0.08930405
      6  4.042743 71.34557 0.1996294  1 0.05363680  0 0.10454367
      7  3.969076 69.34515 0.1982914  1 0.07396419  0 0.09076399
      8  3.869667 69.34122 0.1993334  1 0.06957219  0 0.08570421
      9  4.083545 68.47258 0.2010389  1 0.07686870  0 0.08509724
      10 4.153517 68.51932 0.2009742  1 0.06672599  0 0.11383226
      11 4.028990 67.55515 0.2025752  1 0.06291851  0 0.08666884
      12 4.095689 67.97826 0.2024989  1 0.07236835  0 0.09128803
      13 4.028240 71.82350 0.2024290  1 0.07231574  0 0.09965155
      14 3.983857 65.11927 0.1980753  1 0.07718695  0 0.09164850
      15 3.942173 67.80190 0.2045056  1 0.06681079  0 0.08288682
      16 3.945592 68.78126 0.1997871  1 0.06646189  0 0.10900627
      17 3.842798 68.74961 0.2024411  1 0.06104701  0 0.09405646
      18 3.843747 67.80288 0.2056648  1 0.07270782  0 0.07868634
      19 4.005418 67.09663 0.2037665  1 0.07063172  0 0.10041499
      20 3.988539 67.15033 0.2015780  1 0.06780923  0 0.09632218
      21 4.315006 68.16353 0.1958835  1 0.07113135  0 0.10648345
      22 3.975249 68.53082 0.2040278  1 0.08818818  0 0.08984173
      23 3.921943 69.60221 0.1979640  1 0.07267085  0 0.12511286
      24 3.983742 67.76555 0.2031760  1 0.06717891  0 0.10367067
      25 3.895660 69.47939 0.1993874  1 0.06817941  0 0.08973712
      26 3.958638 69.74054 0.1949952  1 0.07217236  0 0.09218806
      27 3.985585 68.95600 0.1997755  1 0.06356015  0 0.10439573
      28 3.920165 70.13811 0.1952001  1 0.07006469  0 0.06607858
      29 4.135315 69.01056 0.1978472  1 0.07398465  0 0.07513122
      30 3.827007 68.56629 0.2015357  1 0.05506670  0 0.09303089
      31 3.893408 67.49121 0.1975671  1 0.04499286  0 0.08351562
      32 4.025988 70.45182 0.1990793  1 0.07242308  0 0.09782781
      33 4.019942 70.96364 0.2014000  1 0.08447572  0 0.08777964
      34 3.884504 66.55993 0.1971079  1 0.07450460  0 0.08856721
      35 3.996054 66.93821 0.1985125  1 0.06129866  0 0.08274537
      36 4.081618 68.97580 0.2001548  1 0.07202699  0 0.09638920
      37 3.993148 65.10763 0.1985139  1 0.06485045  0 0.07416839
      38 3.908453 65.67802 0.1990315  1 0.06266032  0 0.09186744
      39 3.935636 66.32350 0.2013148  1 0.07732190  0 0.08450629
      40 3.803549 66.89265 0.1998202  1 0.06222067  0 0.09201371
      41 3.794152 65.91168 0.2045937  1 0.08406271  0 0.09105954
      42 3.966106 69.23118 0.1967659  1 0.06171561  0 0.09058125
      43 3.889720 71.25744 0.1962289  1 0.07186124  0 0.10326919
      44 4.154339 67.57126 0.1984505  1 0.06970990  0 0.09025829
      45 3.930934 65.24055 0.1965228  1 0.07108583  0 0.08314434
      46 3.846555 65.76844 0.1990567  1 0.07902672  0 0.07322697
      47 3.758647 65.37698 0.1983753  1 0.08090015  0 0.07478860
      48 4.209726 67.88323 0.1984188  1 0.07493575  0 0.08783786
      49 4.077548 68.53444 0.1997585  1 0.06783113  0 0.07974986
      50 3.901412 65.80751 0.2001533  1 0.06834075  0 0.08115616

# sample_omega

    Code
      sample_omega("testData/runR001.xml", n = 50, seed = 740727)
    Message <simpleMessage>
      Reading testData/runR001.xml
    Output
               OMEGA1       OMEGA2
      1  -0.485774017 -0.185195155
      2   0.019057330 -0.201841117
      3  -0.092220115  0.142531846
      4  -0.053452246 -0.041899840
      5  -0.005860621  0.097999946
      6  -0.122117141 -0.491895640
      7   0.047158306 -0.177104223
      8   0.317216216 -0.176399096
      9  -0.277511640 -0.039923984
      10 -0.466124693 -0.047340588
      11 -0.143078090  0.104464520
      12 -0.317857393  0.037836851
      13 -0.075626875 -0.567076755
      14 -0.057431324  0.487740577
      15  0.096573312  0.065719009
      16  0.102980688 -0.088367559
      17  0.381530989 -0.083297510
      18  0.364015055  0.065651205
      19 -0.085608934  0.176622932
      20 -0.039014555  0.168189713
      21 -0.910806987  0.008494193
      22  0.017880132 -0.048991123
      23  0.180236394 -0.217507094
      24 -0.016373515  0.071400401
      25  0.248784681 -0.198159213
      26  0.081610619 -0.239302183
      27 -0.003202234 -0.115894102
      28  0.191494306 -0.301816487
      29 -0.410335617 -0.124609535
      30  0.421659420 -0.054440752
      31  0.224514080  0.114643305
      32 -0.090473167 -0.351267547
      33 -0.066580283 -0.431786741
      34  0.234546310  0.261169214
      35 -0.063058274  0.201557182
      36 -0.264244743 -0.119094904
      37 -0.083333567  0.489565129
      38  0.156144909  0.399898917
      39  0.091808239  0.298321052
      40  0.459855852  0.208893265
      41  0.470330694  0.363237058
      42  0.053542262 -0.159170186
      43  0.292422796 -0.477895354
      44 -0.483526786  0.101818212
      45  0.088062709  0.468706416
      46  0.325173266  0.385727743
      47  0.558226034  0.447393843
      48 -0.629438185  0.052686743
      49 -0.260387425 -0.049650628
      50  0.176918857  0.379532347

# sample_sigma

    Code
      sample_sigma("testData/runR001.xml", n = 50, seed = 740727)
    Message <simpleMessage>
      Reading testData/runR001.xml
    Output
              SIGMA1
      1  -0.61004507
      2  -0.66487797
      3   0.46950931
      4  -0.13802084
      5   0.32281829
      6  -1.62033672
      7  -0.58339300
      8  -0.58107027
      9  -0.13151224
      10 -0.15594302
      11  0.34411303
      12  0.12463708
      13 -1.86798828
      14  1.60664967
      15  0.21648276
      16 -0.29108857
      17 -0.27438750
      18  0.21625940
      19  0.58180760
      20  0.55402802
      21  0.02798043
      22 -0.16137999
      23 -0.71648273
      24  0.23519764
      25 -0.65274953
      26 -0.78827719
      27 -0.38176282
      28 -0.99420344
      29 -0.41047204
      30 -0.17933143
      31  0.37764262
      32 -1.15709850
      33 -1.42233404
      34  0.86030864
      35  0.66394267
      36 -0.39230648
      37  1.61265986
      38  1.31729344
      39  0.98268924
      40  0.68810821
      41  1.19652685
      42 -0.52431711
      43 -1.57421885
      44  0.33539591
      45  1.54394989
      46  1.27061266
      47  1.47374486
      48  0.17355361
      49 -0.16355245
      50  1.25020461


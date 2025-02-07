1

   Abaqus 2021                                  Date 14-Feb-2023   Time 14:21:07
   For use by HARVARD COLLEGE PRESIDENT AND FELLOWS OF under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     D Y N A M I C   A N A L Y S I S - I M P L I C I T  I N T E G R A T I O N


                                                                                          

          DYNAMIC APPLICATION PACKED SCHEME IS                 QUASI-STATIC        
          BACKWARD EULER (BWE) TIME INTEGRATOR IS USED

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                0.100    
          AND A TOTAL TIME PERIOD OF                            2.00    
          THE MINIMUM TIME INCREMENT ALLOWED IS                5.000E-12
          THE MAXIMUM TIME INCREMENT ALLOWED IS                 2.00    

          THE HALF-INCREMENT RESIDUALS WILL NOT BE CALCULATED

          INITIAL ACCELERATIONS WILL NOT BE CALCULATED
          AGGRESSIVE   TIME INCREMENTATION IS USED
          STEP AMPLITUDE IS                                    RAMP

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE FORCE                               1.000E-02
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

 CONVERGENCE TOLERANCE PARAMETERS FOR MOMENT   
     CRITERION FOR RESIDUAL MOMENT    FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR ROTATION CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE MOMENT                              1.000E-02
     AVERAGE MOMENT    IS TIME AVERAGE MOMENT   
     ALTERNATE CRIT. FOR RESIDUAL MOMENT    FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO MOMENT    RELATIVE TO TIME AVRG. MOMENT        1.000E-05
     CRITERION FOR RESIDUAL MOMENT    WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR ROTATION CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL MOMENT    FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO MOMENT    REL. TO TIME AVRG. MAX. MOMENT       1.000E-05

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 TIME INCREMENTATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AT WHICH LOG. CONVERGENCE RATE CHECK BEGINS         8
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     EQUILIBRIUM ITERATION COUNT FOR CUT-BACK IN NEXT INCREMENT               10
     MAXIMUM EQUILIB. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      4
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            50
 *** MAXIMUM ATTEMPTS ALLOWED IN AN INCREMENT                                 25
     MAXIMUM DISCON. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      50
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE   50
     CUT-BACK FACTOR AFTER DIVERGENCE                                 0.2500    
     CUT-BACK FACTOR FOR TOO SLOW CONVERGENCE                         0.5000    
     CUT-BACK FACTOR AFTER TOO MANY EQUILIBRIUM ITERATIONS            0.7500    
     CUT-BACK FACTOR AFTER TOO MANY SEVERE DISCONTINUITY ITERATIONS   0.2500    
     CUT-BACK FACTOR AFTER PROBLEMS IN ELEMENT ASSEMBLY               0.2500    
     INCREASE FACTOR AFTER TWO INCREMENTS THAT CONVERGE QUICKLY        1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED                       1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DYNAMICS)            1.250    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DIFFUSION)           2.000    
     MINIMUM TIME INCREMENT RATIO FOR EXTRAPOLATION TO OCCUR          0.1000    
     MAX. RATIO OF TIME INCREMENT TO STABILITY LIMIT                   1.000    
     FRACTION OF STABILITY LIMIT FOR NEW TIME INCREMENT               0.9500    
     TIME INCREMENT INCREASE FACTOR BEFORE A TIME POINT                1.000    

 *** INDICATES USE OF NON-DEFAULT TIME INCREMENTATION CONTROL PARAMETER
     SELECTED WITH *CONTROLS
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                     500

          LARGE DISPLACEMENT THEORY WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     0.854    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS
          THE INITIAL TIME INCREMENT WAS REDUCED TO   7.400E-02 TO MEET THE NEXT OUTPUT TIME POINT

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.400E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.927E-03   TIME AVG. FORCE       1.927E-03
 LARGEST RESIDUAL FORCE            -0.181       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -3.46       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -3.46       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.167E-04   TIME AVG. MOMENT      2.167E-04
 LARGEST RESIDUAL MOMENT            7.827E-03   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.819E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.819E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 2 NEGATIVE EIGENVALUES.

     NEGATIVE EIGENVALUES MEAN THAT THE SYSTEM MATRIX IS NOT 
     POSITIVE DEFINITE:
     FOR EXAMPLE, A BIFURCATION (BUCKLING) LOAD MAY HAVE BEEN EXCEEDED.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.572E-06   TIME AVG. FORCE       4.572E-06
 LARGEST RESIDUAL FORCE             4.667E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -3.47       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.001E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.700E-06   TIME AVG. MOMENT      1.700E-06
 LARGEST RESIDUAL MOMENT            4.163E-05   AT NODE        399   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.847E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.853E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 2 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.332E-08   TIME AVG. FORCE       1.172E-02
                                                (FROM MOMENT    FIELD)
 LARGEST RESIDUAL FORCE            -6.686E-07   AT NODE        292   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -3.48       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.447E-03   AT NODE        321   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     8.797E-07   TIME AVG. MOMENT      8.797E-07
 LARGEST RESIDUAL MOMENT            1.278E-07   AT NODE        218   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.835E-02   AT NODE        356   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.680E-04   AT NODE        292   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.753E-08   TIME AVG. FORCE       1.172E-02
                                                (FROM MOMENT    FIELD)
 LARGEST RESIDUAL FORCE            -2.058E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -3.46       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.850E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     8.765E-07   TIME AVG. MOMENT      8.765E-07
 LARGEST RESIDUAL MOMENT           -2.547E-09   AT NODE        179   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.813E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.057E-04   AT NODE        237   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.177E-08   TIME AVG. FORCE       1.172E-02
                                                (FROM MOMENT    FIELD)
 LARGEST RESIDUAL FORCE            -7.589E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -3.46       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.965E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE

 AVERAGE MOMENT                     8.765E-07   TIME AVG. MOMENT      8.765E-07
 LARGEST RESIDUAL MOMENT            3.632E-11   AT NODE        234   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.813E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.875E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   7.039E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.400E-02,  FRACTION OF STEP COMPLETED  3.700E-02
 STEP TIME COMPLETED       7.400E-02,  TOTAL TIME COMPLETED        7.400E-02


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.039E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 7 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.660       TIME AVG. FORCE       0.660    
 LARGEST RESIDUAL FORCE             -52.3       AT NODE        264   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          18.7       AT NODE        360   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         21.6       AT NODE        364   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      1.43       TIME AVG. MOMENT      0.717    
 LARGEST RESIDUAL MOMENT            -93.7       AT NODE        263   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      0.847       AT NODE        256   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     0.899       AT NODE        256   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 231 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       1.05       TIME AVG. FORCE        1.05    
 LARGEST RESIDUAL FORCE             -294.       AT NODE        345   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          21.5       AT NODE        341   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -7.29       AT NODE        329   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      2.27       TIME AVG. MOMENT       1.14    
 LARGEST RESIDUAL MOMENT            -119.       AT NODE        278   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      -3.84       AT NODE        301   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION      8.66       AT NODE        341   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 394 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       16.1       TIME AVG. FORCE        16.1    
 LARGEST RESIDUAL FORCE             4.154E+03   AT NODE         60   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -144.       AT NODE         76   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -143.       AT NODE         76   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      165.       TIME AVG. MOMENT       82.3    
 LARGEST RESIDUAL MOMENT           -2.945E+04   AT NODE         59   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      -6.06       AT NODE         17   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION      41.7       AT NODE        376   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 960 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       191.       TIME AVG. FORCE        191.    
 LARGEST RESIDUAL FORCE             2.437E+05   AT NODE         33   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.083E+03   AT NODE         48   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.083E+03   AT NODE         48   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.359E+05   TIME AVG. MOMENT      6.795E+04
 LARGEST RESIDUAL MOMENT           -8.337E+07   AT NODE         38   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      -5.02       AT NODE         69   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.749E+03   AT NODE         24   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     2 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.760E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.253E-04   TIME AVG. FORCE       2.253E-04
 LARGEST RESIDUAL FORCE            -5.137E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.61       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.784       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.523E-05   TIME AVG. MOMENT      1.306E-05
 LARGEST RESIDUAL MOMENT           -6.515E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.590E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.208E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.051E-05   TIME AVG. FORCE       1.051E-05
 LARGEST RESIDUAL FORCE            -4.847E-05   AT NODE         52   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.59       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.728E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.117E-06   TIME AVG. MOMENT      2.497E-06
 LARGEST RESIDUAL MOMENT           -1.868E-05   AT NODE         67   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.576E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.262E-03   AT NODE         52   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.154E-04   TIME AVG. FORCE       2.154E-04
 LARGEST RESIDUAL FORCE            -6.293E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.948       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.642       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.282E-05   TIME AVG. MOMENT      1.685E-05
 LARGEST RESIDUAL MOMENT           -1.682E-03   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.736E-02   AT NODE        208   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.497E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.048E-06   TIME AVG. FORCE       4.048E-06
 LARGEST RESIDUAL FORCE            -3.012E-04   AT NODE        264   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.946       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.874E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.633E-06   TIME AVG. MOMENT      1.755E-06
 LARGEST RESIDUAL MOMENT            7.830E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.700E-02   AT NODE        208   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.622E-03   AT NODE        264   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.373E-05   TIME AVG. FORCE       1.373E-05
 LARGEST RESIDUAL FORCE            -4.636E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.08       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.139       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.100E-06   TIME AVG. MOMENT      1.988E-06
 LARGEST RESIDUAL MOMENT            3.108E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.036E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.529E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.765E-06   TIME AVG. FORCE       3.765E-06
 LARGEST RESIDUAL FORCE             3.492E-04   AT NODE        269   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.08       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.942E-03   AT NODE        284   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.577E-06   TIME AVG. MOMENT      1.727E-06
 LARGEST RESIDUAL MOMENT            9.635E-06   AT NODE        278   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.988E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.041E-03   AT NODE        269   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      1.171E-06   TIME AVG. FORCE       1.171E-06
 LARGEST RESIDUAL FORCE            -1.614E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.06       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.479E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.371E-06   TIME AVG. MOMENT      1.624E-06
 LARGEST RESIDUAL MOMENT           -4.042E-07   AT NODE        158   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.812E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.765E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      1.188E-06   TIME AVG. FORCE       1.188E-06
 LARGEST RESIDUAL FORCE            -1.172E-04   AT NODE        228   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.06       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.675E-03   AT NODE        224   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.379E-06   TIME AVG. MOMENT      1.628E-06
 LARGEST RESIDUAL MOMENT           -4.238E-06   AT NODE        242   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.800E-02   AT NODE        389   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.501E-03   AT NODE        228   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      1.845E-07   TIME AVG. FORCE       1.845E-07
 LARGEST RESIDUAL FORCE             5.931E-07   AT NODE        153   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.06       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.443E-04   AT NODE        169   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.324E-06   TIME AVG. MOMENT      1.600E-06
 LARGEST RESIDUAL MOMENT           -7.582E-08   AT NODE        191   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.799E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.494E-04   AT NODE        153   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                      6.095E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.809E-06   AT NODE        204   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.06       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.396E-04   AT NODE        204   DOF  3
   INSTANCE: PART-1-1                                                                        
          ALL DISP.    CORRECTIONS ARE ZERO

 AVERAGE MOMENT                     2.367E-06   TIME AVG. MOMENT      1.622E-06
 LARGEST RESIDUAL MOMENT           -1.409E-07   AT NODE        227   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.797E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.313E-04   AT NODE        204   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     2 STARTS. ATTEMPT NUMBER  3, TIME INCREMENT  8.799E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.432E-06   TIME AVG. FORCE       4.432E-06
 LARGEST RESIDUAL FORCE            -5.668E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.505       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.377E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.260E-06   TIME AVG. MOMENT      3.068E-06
 LARGEST RESIDUAL MOMENT           -9.002E-05   AT NODE        391   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.045E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.536E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.999E-07   TIME AVG. FORCE       7.999E-07
 LARGEST RESIDUAL FORCE            -3.425E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.486       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.926E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.143E-06   TIME AVG. MOMENT      1.510E-06
 LARGEST RESIDUAL MOMENT            2.481E-06   AT NODE        266   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.624E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.230E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.483E-07   TIME AVG. FORCE       8.483E-07
 LARGEST RESIDUAL FORCE            -1.116E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.467       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.855E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.028E-06   TIME AVG. MOMENT      1.452E-06
 LARGEST RESIDUAL MOMENT           -1.226E-06   AT NODE        263   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.036E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.378E-03   AT NODE        237   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.564E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.096E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.466       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.822E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.018E-06   TIME AVG. MOMENT      1.447E-06
 LARGEST RESIDUAL MOMENT            2.630E-08   AT NODE        175   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.889E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.473E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.467E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -5.000E-07   AT NODE        232   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.465       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.985E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.065E-06   TIME AVG. MOMENT      1.471E-06
 LARGEST RESIDUAL MOMENT           -1.638E-08   AT NODE        238   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.833E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.288E-04   AT NODE        225   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.007E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE             5.645E-12   AT NODE        197   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.465       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.458E-06   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     2.071E-06   TIME AVG. MOMENT      1.474E-06
 LARGEST RESIDUAL MOMENT            2.757E-12   AT NODE        210   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.833E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.698E-07   AT NODE        204   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.799E-03,  FRACTION OF STEP COMPLETED  4.140E-02
 STEP TIME COMPLETED       8.280E-02,  TOTAL TIME COMPLETED        8.280E-02


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.799E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.381E-07   TIME AVG. FORCE       7.381E-07
 LARGEST RESIDUAL FORCE            -6.207E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.526       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.048E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.578E-06   TIME AVG. MOMENT      1.842E-06
 LARGEST RESIDUAL MOMENT           -8.029E-06   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.940E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.107E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.107E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.122E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.525       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.717E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.340E-06   TIME AVG. MOMENT      1.762E-06
 LARGEST RESIDUAL MOMENT            3.355E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.893E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.662E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.080E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -6.511E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.525       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.460E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     2.340E-06   TIME AVG. MOMENT      1.762E-06
 LARGEST RESIDUAL MOMENT            1.859E-12   AT NODE        314   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.882E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.116E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.799E-03,  FRACTION OF STEP COMPLETED  4.580E-02
 STEP TIME COMPLETED       9.160E-02,  TOTAL TIME COMPLETED        9.160E-02


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.799E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.295E-07   TIME AVG. FORCE       9.295E-07
 LARGEST RESIDUAL FORCE            -8.263E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.591       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.609E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.031E-06   TIME AVG. MOMENT      2.080E-06
 LARGEST RESIDUAL MOMENT           -1.198E-05   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.017E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.286E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.675E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.522E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.589       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.391E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.637E-06   TIME AVG. MOMENT      1.981E-06
 LARGEST RESIDUAL MOMENT            7.997E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.008E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.256E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.582E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.653E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.589       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.446E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     2.642E-06   TIME AVG. MOMENT      1.982E-06
 LARGEST RESIDUAL MOMENT            1.032E-11   AT NODE        318   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.006E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.644E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.320E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.799E-03,  FRACTION OF STEP COMPLETED  5.020E-02
 STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED        0.100    


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.320E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.692E-06   TIME AVG. FORCE       7.692E-06
 LARGEST RESIDUAL FORCE            -9.818E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.07       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.184       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.046E-06   TIME AVG. MOMENT      2.995E-06
 LARGEST RESIDUAL MOMENT           -9.673E-05   AT NODE        386   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.952E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.436E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.023E-07   TIME AVG. FORCE       1.023E-07
 LARGEST RESIDUAL FORCE            -2.373E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.06       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.635E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.182E-06   TIME AVG. MOMENT      2.222E-06
 LARGEST RESIDUAL MOMENT           -5.296E-07   AT NODE        199   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.930E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.617E-04   AT NODE        131   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.015E-07   TIME AVG. FORCE       3.015E-07
 LARGEST RESIDUAL FORCE            -8.007E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.04       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.034E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.131E-06   TIME AVG. MOMENT      2.212E-06
 LARGEST RESIDUAL MOMENT           -8.329E-08   AT NODE        374   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.806E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.243E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.425E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -7.068E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.04       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.792E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     3.174E-06   TIME AVG. MOMENT      2.221E-06
 LARGEST RESIDUAL MOMENT            1.052E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.804E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.514E-05   AT NODE        253   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.980E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.320E-02,  FRACTION OF STEP COMPLETED  5.680E-02
 STEP TIME COMPLETED       0.114    ,  TOTAL TIME COMPLETED        0.114    


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.980E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.312E-04   TIME AVG. FORCE       4.312E-04
 LARGEST RESIDUAL FORCE            -8.749E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -2.62       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -1.06       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.817E-05   TIME AVG. MOMENT      9.879E-06
 LARGEST RESIDUAL MOMENT           -1.724E-03   AT NODE        302   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.009E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.303E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.269E-05   TIME AVG. FORCE       2.269E-05
 LARGEST RESIDUAL FORCE            -6.634E-04   AT NODE         56   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -2.62       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.654E-02   AT NODE        321   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.552E-06   TIME AVG. MOMENT      3.276E-06
 LARGEST RESIDUAL MOMENT           -2.244E-04   AT NODE         71   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.998E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.351E-03   AT NODE         56   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.706E-03   TIME AVG. FORCE       3.706E-03
 LARGEST RESIDUAL FORCE            -0.725       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -3.89       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -1.28       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.531E-03   TIME AVG. MOMENT      4.237E-04
 LARGEST RESIDUAL MOMENT            9.458E-02   AT NODE        294   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.192       AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.133       AT NODE        293   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.939E-04   TIME AVG. FORCE       4.939E-04
 LARGEST RESIDUAL FORCE            -2.215E-02   AT NODE        258   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -3.57       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.439       AT NODE        272   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.413E-04   TIME AVG. MOMENT      2.539E-05
 LARGEST RESIDUAL MOMENT           -8.199E-03   AT NODE        291   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.196       AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.968E-02   AT NODE        248   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 9 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.478E-02   TIME AVG. FORCE       2.478E-02
 LARGEST RESIDUAL FORCE              2.68       AT NODE        273   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -6.15       AT NODE        397   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -2.75       AT NODE        393   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.442E-02   TIME AVG. MOMENT      5.738E-03
 LARGEST RESIDUAL MOMENT             1.72       AT NODE        278   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.616       AT NODE        381   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.447       AT NODE        381   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     6 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  4.949E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.107E-07   TIME AVG. FORCE       2.107E-07
 LARGEST RESIDUAL FORCE             2.136E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.413       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.318E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.490E-06   TIME AVG. MOMENT      2.432E-06
 LARGEST RESIDUAL MOMENT           -5.005E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.269E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.030E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.145E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.050E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.412       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.549E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     3.391E-06   TIME AVG. MOMENT      2.416E-06
 LARGEST RESIDUAL MOMENT            1.095E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.193E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.667E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.677E-08
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.088E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -2.076E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.412       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.907E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     3.395E-06   TIME AVG. MOMENT      2.416E-06
 LARGEST RESIDUAL MOMENT            7.081E-13   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.191E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.010E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.424E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.949E-03,  FRACTION OF STEP COMPLETED  5.927E-02
 STEP TIME COMPLETED       0.119    ,  TOTAL TIME COMPLETED        0.119    


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.424E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.311E-07   TIME AVG. FORCE       7.311E-07
 LARGEST RESIDUAL FORCE            -6.029E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.672       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.386E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.139E-06   TIME AVG. MOMENT      2.662E-06
 LARGEST RESIDUAL MOMENT           -1.299E-05   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.192E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.132E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.629E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -7.887E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.670       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.864E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     3.744E-06   TIME AVG. MOMENT      2.606E-06
 LARGEST RESIDUAL MOMENT            9.197E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.179E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.231E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -8.710E-08
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.433E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.341E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.670       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.986E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     3.746E-06   TIME AVG. MOMENT      2.606E-06
 LARGEST RESIDUAL MOMENT            2.156E-11   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.177E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.538E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.114E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.424E-03,  FRACTION OF STEP COMPLETED  6.298E-02
 STEP TIME COMPLETED       0.126    ,  TOTAL TIME COMPLETED        0.126    


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.114E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.908E-06   TIME AVG. FORCE       7.908E-06
 LARGEST RESIDUAL FORCE            -1.193E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.17       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.164       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.474E-06   TIME AVG. MOMENT      3.465E-06
 LARGEST RESIDUAL MOMENT           -1.489E-04   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.266E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.002E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.914E-07   TIME AVG. FORCE       1.914E-07
 LARGEST RESIDUAL FORCE            -7.627E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.16       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.173E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.387E-06   TIME AVG. MOMENT      2.829E-06
 LARGEST RESIDUAL MOMENT            1.249E-06   AT NODE        231   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.226E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.342E-04   AT NODE        141   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.761E-07   TIME AVG. FORCE       8.761E-07
 LARGEST RESIDUAL FORCE            -2.663E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.12       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.338E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.304E-06   TIME AVG. MOMENT      2.819E-06
 LARGEST RESIDUAL MOMENT           -4.491E-07   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.998E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.280E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.497E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -8.977E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.12       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.365E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     4.341E-06   TIME AVG. MOMENT      2.823E-06
 LARGEST RESIDUAL MOMENT            1.244E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.994E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.379E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.670E-02
          TIME INCREMENT REDUCED TO   7.286E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.114E-02,  FRACTION OF STEP COMPLETED  6.855E-02
 STEP TIME COMPLETED       0.137    ,  TOTAL TIME COMPLETED        0.137    


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.286E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.567E-06   TIME AVG. FORCE       1.567E-06
 LARGEST RESIDUAL FORCE             2.422E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.801       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.584E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.253E-06   TIME AVG. MOMENT      3.204E-06
 LARGEST RESIDUAL MOMENT           -5.245E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.512E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.075E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.602E-07   TIME AVG. FORCE       1.602E-07
 LARGEST RESIDUAL FORCE            -2.626E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.789       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.150E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.754E-06   TIME AVG. MOMENT      3.038E-06
 LARGEST RESIDUAL MOMENT            1.820E-07   AT NODE        298   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.441E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.144E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      9.409E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -8.920E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.785       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.650E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     4.756E-06   TIME AVG. MOMENT      3.038E-06
 LARGEST RESIDUAL MOMENT            9.460E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.399E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.173E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    2.169E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.434E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.704E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.785       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.768E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     4.752E-06   TIME AVG. MOMENT      3.037E-06
 LARGEST RESIDUAL MOMENT            2.655E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.396E-02   AT NODE        389   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.672E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.506E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.286E-03,  FRACTION OF STEP COMPLETED  7.219E-02
 STEP TIME COMPLETED       0.144    ,  TOTAL TIME COMPLETED        0.144    


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.506E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 2 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.199       TIME AVG. FORCE       0.199    
 LARGEST RESIDUAL FORCE             -33.9       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -17.0       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -14.3       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     0.376       TIME AVG. MOMENT      3.765E-02
 LARGEST RESIDUAL MOMENT             12.9       AT NODE        334   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.673       AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.625       AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 175 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.849       TIME AVG. FORCE       0.849    
 LARGEST RESIDUAL FORCE             -611.       AT NODE        241   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -15.7       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -2.93       AT NODE        395   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      2.05       TIME AVG. MOMENT      0.205    
 LARGEST RESIDUAL MOMENT            -289.       AT NODE        235   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      -4.48       AT NODE        360   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     -10.4       AT NODE        401   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 280 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       24.7       TIME AVG. FORCE        24.7    
 LARGEST RESIDUAL FORCE             6.652E+03   AT NODE        311   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -119.       AT NODE        336   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -110.       AT NODE        196   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      279.       TIME AVG. MOMENT       27.9    
 LARGEST RESIDUAL MOMENT            4.233E+04   AT NODE        311   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      -5.98       AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     -62.8       AT NODE        260   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 944 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       190.       TIME AVG. FORCE        190.    
 LARGEST RESIDUAL FORCE             1.026E+05   AT NODE        386   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -566.       AT NODE        313   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -551.       AT NODE        313   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.568E+04   TIME AVG. MOMENT      5.568E+03
 LARGEST RESIDUAL MOMENT           -1.353E+07   AT NODE        313   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION       5.10       AT NODE        322   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION      651.       AT NODE        400   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    10 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  6.264E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.559E-07   TIME AVG. FORCE       6.559E-07
 LARGEST RESIDUAL FORCE             9.240E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.717       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.213E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.705E-06   TIME AVG. MOMENT      3.304E-06
 LARGEST RESIDUAL MOMENT           -2.064E-05   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.312E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.121E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.889E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -6.964E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.711       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.083E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     5.128E-06   TIME AVG. MOMENT      3.246E-06
 LARGEST RESIDUAL MOMENT            2.460E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.276E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.666E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.125E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.661E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.711       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.892E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          ALL DISP.    CORRECTIONS ARE ZERO

 AVERAGE MOMENT                     5.132E-06   TIME AVG. MOMENT      3.247E-06
 LARGEST RESIDUAL MOMENT            4.354E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.270E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.704E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    9.396E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.264E-03,  FRACTION OF STEP COMPLETED  7.533E-02
 STEP TIME COMPLETED       0.151    ,  TOTAL TIME COMPLETED        0.151    


  INCREMENT    11 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.396E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.638E-06   TIME AVG. FORCE       6.638E-06
 LARGEST RESIDUAL FORCE            -1.086E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.20       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.133       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.091E-05   TIME AVG. MOMENT      3.943E-06
 LARGEST RESIDUAL MOMENT           -1.705E-04   AT NODE        395   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.396E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.910E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.663E-07   TIME AVG. FORCE       2.663E-07
 LARGEST RESIDUAL FORCE            -1.724E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.18       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.602E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.830E-06   TIME AVG. MOMENT      3.482E-06
 LARGEST RESIDUAL MOMENT            1.567E-06   AT NODE        250   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.337E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.205E-04   AT NODE        304   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.067E-06   TIME AVG. FORCE       1.067E-06
 LARGEST RESIDUAL FORCE            -3.456E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.15       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.397E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.739E-06   TIME AVG. MOMENT      3.473E-06
 LARGEST RESIDUAL MOMENT           -6.401E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.076E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.609E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      9.278E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.854E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.15       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.210E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     5.752E-06   TIME AVG. MOMENT      3.475E-06
 LARGEST RESIDUAL MOMENT            3.342E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.067E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.781E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.256E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -4.330E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.15       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.050E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          ALL DISP.    CORRECTIONS ARE ZERO

 AVERAGE MOMENT                     5.745E-06   TIME AVG. MOMENT      3.474E-06
 LARGEST RESIDUAL MOMENT            2.707E-11   AT NODE        147   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.058E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.176E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.396E-03,  FRACTION OF STEP COMPLETED  8.002E-02
 STEP TIME COMPLETED       0.160    ,  TOTAL TIME COMPLETED        0.160    


  INCREMENT    12 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.396E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.536E-05   TIME AVG. FORCE       2.536E-05
 LARGEST RESIDUAL FORCE            -4.565E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.36       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.209       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.239E-05   TIME AVG. MOMENT      5.050E-06
 LARGEST RESIDUAL MOMENT           -4.775E-04   AT NODE        375   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.083E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.025E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.541E-06   TIME AVG. FORCE       1.541E-06
 LARGEST RESIDUAL FORCE            -2.780E-05   AT NODE         88   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.34       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.083E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.037E-06   TIME AVG. MOMENT      3.771E-06
 LARGEST RESIDUAL MOMENT            1.417E-05   AT NODE        215   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.044E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.708E-03   AT NODE         88   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.137E-05   TIME AVG. FORCE       1.137E-05
 LARGEST RESIDUAL FORCE            -2.963E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.23       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.106       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.861E-06   TIME AVG. MOMENT      3.756E-06
 LARGEST RESIDUAL MOMENT            1.770E-05   AT NODE        234   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.253E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.912E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.961E-07   TIME AVG. FORCE       1.961E-07
 LARGEST RESIDUAL FORCE            -6.690E-07   AT NODE        156   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.23       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.871E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.411E-06   TIME AVG. MOMENT      3.719E-06
 LARGEST RESIDUAL MOMENT           -3.040E-07   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.249E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.635E-04   AT NODE        156   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.728E-07   TIME AVG. FORCE       1.728E-07
 LARGEST RESIDUAL FORCE            -7.059E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.22       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.826E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.399E-06   TIME AVG. MOMENT      3.718E-06
 LARGEST RESIDUAL MOMENT            3.634E-07   AT NODE        194   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.207E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.866E-04   AT NODE        265   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.062E-07   TIME AVG. FORCE       1.062E-07
 LARGEST RESIDUAL FORCE            -4.345E-10   AT NODE        208   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.22       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.836E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     6.403E-06   TIME AVG. MOMENT      3.718E-06
 LARGEST RESIDUAL MOMENT            2.859E-10   AT NODE        243   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.207E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.691E-06   AT NODE        208   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.396E-03,  FRACTION OF STEP COMPLETED  8.472E-02
 STEP TIME COMPLETED       0.169    ,  TOTAL TIME COMPLETED        0.169    


  INCREMENT    13 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.396E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.578E-04   TIME AVG. FORCE       7.898E-05
 LARGEST RESIDUAL FORCE            -9.986E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.51       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.283       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.217E-05   TIME AVG. MOMENT      7.445E-06
 LARGEST RESIDUAL MOMENT           -1.194E-03   AT NODE        307   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.788E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.634E-02   AT NODE        277   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.333E-05   TIME AVG. FORCE       6.718E-06
 LARGEST RESIDUAL FORCE            -2.462E-04   AT NODE         76   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.49       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.960E-02   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.330E-06   TIME AVG. MOMENT      4.150E-06
 LARGEST RESIDUAL MOMENT           -1.121E-04   AT NODE         99   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.757E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.087E-03   AT NODE         76   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.439E-05   TIME AVG. FORCE       3.225E-05
 LARGEST RESIDUAL FORCE            -6.492E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.31       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.175       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.077E-05   TIME AVG. MOMENT      4.261E-06
 LARGEST RESIDUAL MOMENT            1.042E-04   AT NODE        255   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.547E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.210E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.404E-06   TIME AVG. FORCE       1.255E-06
 LARGEST RESIDUAL FORCE            -1.207E-05   AT NODE         68   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.32       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.210E-03   AT NODE        232   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.268E-06   TIME AVG. MOMENT      3.991E-06
 LARGEST RESIDUAL MOMENT           -2.257E-06   AT NODE         83   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.538E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.129E-03   AT NODE         68   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.066E-06   TIME AVG. FORCE       4.586E-06
 LARGEST RESIDUAL FORCE            -4.416E-04   AT NODE        264   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.30       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.651E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.771E-06   TIME AVG. MOMENT      4.030E-06
 LARGEST RESIDUAL MOMENT            4.092E-05   AT NODE        198   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.362E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.629E-03   AT NODE        261   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.623E-07   TIME AVG. FORCE       2.343E-07
 LARGEST RESIDUAL FORCE            -4.146E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.30       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.573E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.114E-06   TIME AVG. MOMENT      3.979E-06
 LARGEST RESIDUAL MOMENT            1.318E-07   AT NODE        187   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.353E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.544E-04   AT NODE        156   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      2.848E-07   TIME AVG. FORCE       1.955E-07
 LARGEST RESIDUAL FORCE            -1.699E-06   AT NODE        220   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.30       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.217E-04   AT NODE        224   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.112E-06   TIME AVG. MOMENT      3.979E-06
 LARGEST RESIDUAL MOMENT            1.891E-07   AT NODE        174   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.349E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.119E-04   AT NODE        220   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      2.613E-07   TIME AVG. FORCE       1.837E-07
 LARGEST RESIDUAL FORCE            -1.685E-11   AT NODE        173   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.30       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.594E-06   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.109E-06   TIME AVG. MOMENT      3.979E-06
 LARGEST RESIDUAL MOMENT            1.288E-11   AT NODE        178   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.349E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.438E-07   AT NODE        176   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.396E-03,  FRACTION OF STEP COMPLETED  8.942E-02
 STEP TIME COMPLETED       0.179    ,  TOTAL TIME COMPLETED        0.179    


  INCREMENT    14 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.396E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      0.106       TIME AVG. FORCE       3.523E-02
 LARGEST RESIDUAL FORCE             -4.07       AT NODE        284   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -2.84       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -1.54       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     0.121       TIME AVG. MOMENT      8.612E-03
 LARGEST RESIDUAL MOMENT             3.90       AT NODE        282   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.463       AT NODE        281   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.462       AT NODE        284   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 124 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                       9.51       TIME AVG. FORCE        3.17    
 LARGEST RESIDUAL FORCE             1.321E+03   AT NODE        151   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          39.0       AT NODE        165   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         39.3       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      22.7       TIME AVG. MOMENT       1.62    
 LARGEST RESIDUAL MOMENT            5.046E+03   AT NODE        154   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      -5.06       AT NODE        320   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     -11.9       AT NODE        312   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 592 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       24.8       TIME AVG. FORCE        8.27    
 LARGEST RESIDUAL FORCE             5.016E+03   AT NODE        374   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -46.6       AT NODE        373   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -82.6       AT NODE        373   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      227.       TIME AVG. MOMENT       16.2    
 LARGEST RESIDUAL MOMENT            1.043E+05   AT NODE        374   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION       5.04       AT NODE        151   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     -130.       AT NODE        384   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 861 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       163.       TIME AVG. FORCE        54.3    
 LARGEST RESIDUAL FORCE             3.490E+04   AT NODE        135   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          679.       AT NODE        368   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         691.       AT NODE        372   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.604E+03   TIME AVG. MOMENT       686.    
 LARGEST RESIDUAL MOMENT            2.030E+06   AT NODE        140   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION       4.45       AT NODE        156   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION      145.       AT NODE        135   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    14 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  2.349E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.246E-07   TIME AVG. FORCE       2.307E-07
 LARGEST RESIDUAL FORCE             5.240E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.329       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.125E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.290E-06   TIME AVG. MOMENT      4.215E-06
 LARGEST RESIDUAL MOMENT           -1.310E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.001E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.292E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.889E-07   TIME AVG. FORCE       2.188E-07
 LARGEST RESIDUAL FORCE            -9.107E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.328       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.608E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.291E-06   TIME AVG. MOMENT      4.216E-06
 LARGEST RESIDUAL MOMENT            9.166E-11   AT NODE        362   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.959E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.246E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.887E-07   TIME AVG. FORCE       2.187E-07
 LARGEST RESIDUAL FORCE             2.966E-12   AT NODE        387   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.328       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.458E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.291E-06   TIME AVG. MOMENT      4.216E-06
 LARGEST RESIDUAL MOMENT           -4.527E-13   AT NODE        390   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.959E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.799E-08   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.349E-03,  FRACTION OF STEP COMPLETED  9.060E-02
 STEP TIME COMPLETED       0.181    ,  TOTAL TIME COMPLETED        0.181    


  INCREMENT    15 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.349E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.018E-07   TIME AVG. FORCE       2.395E-07
 LARGEST RESIDUAL FORCE             1.566E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.333       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.391E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.454E-06   TIME AVG. MOMENT      4.431E-06
 LARGEST RESIDUAL MOMENT           -3.802E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.049E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.022E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.827E-07   TIME AVG. FORCE       2.347E-07
 LARGEST RESIDUAL FORCE            -6.691E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.333       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.101E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.475E-06   TIME AVG. MOMENT      4.433E-06
 LARGEST RESIDUAL MOMENT            1.497E-12   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.046E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.635E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.523E-03

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.349E-03,  FRACTION OF STEP COMPLETED  9.177E-02
 STEP TIME COMPLETED       0.184    ,  TOTAL TIME COMPLETED        0.184    


  INCREMENT    16 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.523E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.852E-07   TIME AVG. FORCE       2.648E-07
 LARGEST RESIDUAL FORCE             7.366E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.509       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.970E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.775E-06   TIME AVG. MOMENT      4.642E-06
 LARGEST RESIDUAL MOMENT           -1.756E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.284E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.159E-04   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.937E-07   TIME AVG. FORCE       2.465E-07
 LARGEST RESIDUAL FORCE            -3.993E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.509       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.875E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.753E-06   TIME AVG. MOMENT      4.640E-06
 LARGEST RESIDUAL MOMENT            6.834E-11   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.256E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.797E-05   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.936E-07   TIME AVG. FORCE       2.465E-07
 LARGEST RESIDUAL FORCE            -1.634E-12   AT NODE        299   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.509       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.931E-07   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.753E-06   TIME AVG. MOMENT      4.640E-06
 LARGEST RESIDUAL MOMENT            4.024E-13   AT NODE        397   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.256E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.996E-08   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.285E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.523E-03,  FRACTION OF STEP COMPLETED  9.353E-02
 STEP TIME COMPLETED       0.187    ,  TOTAL TIME COMPLETED        0.187    


  INCREMENT    17 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.285E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.469E-07   TIME AVG. FORCE       3.632E-07
 LARGEST RESIDUAL FORCE             6.439E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.788       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.530E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.550E-06   TIME AVG. MOMENT      4.870E-06
 LARGEST RESIDUAL MOMENT           -1.462E-05   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.465E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.688E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.448E-07   TIME AVG. FORCE       2.629E-07
 LARGEST RESIDUAL FORCE            -6.040E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.784       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.239E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.189E-06   TIME AVG. MOMENT      4.849E-06
 LARGEST RESIDUAL MOMENT            1.157E-08   AT NODE        334   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.431E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.443E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.726E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.129E-07   TIME AVG. FORCE       2.575E-07
 LARGEST RESIDUAL FORCE            -2.051E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.784       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.849E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.184E-06   TIME AVG. MOMENT      4.849E-06
 LARGEST RESIDUAL MOMENT            1.390E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.429E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.032E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.128E-07   TIME AVG. FORCE       2.575E-07
 LARGEST RESIDUAL FORCE             2.408E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.784       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.830E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     8.184E-06   TIME AVG. MOMENT      4.849E-06
 LARGEST RESIDUAL MOMENT           -4.889E-13   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.429E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.160E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.928E-03

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.285E-03,  FRACTION OF STEP COMPLETED  9.617E-02
 STEP TIME COMPLETED       0.192    ,  TOTAL TIME COMPLETED        0.192    


  INCREMENT    18 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.928E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.856E-05   TIME AVG. FORCE       2.872E-06
 LARGEST RESIDUAL FORCE             1.547E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.29       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.116       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.573E-05   TIME AVG. MOMENT      5.453E-06
 LARGEST RESIDUAL MOMENT           -2.626E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.745E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.021E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.121E-06   TIME AVG. FORCE       3.809E-07
 LARGEST RESIDUAL FORCE            -2.694E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.27       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.161E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.027E-06   TIME AVG. MOMENT      5.081E-06
 LARGEST RESIDUAL MOMENT            3.206E-06   AT NODE        258   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.672E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.584E-04   AT NODE        284   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.061E-06   TIME AVG. FORCE       9.437E-07
 LARGEST RESIDUAL FORCE            -7.793E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.22       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.556E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.867E-06   TIME AVG. MOMENT      5.072E-06
 LARGEST RESIDUAL MOMENT           -2.261E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.274E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.981E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.521E-07   TIME AVG. FORCE       2.710E-07
 LARGEST RESIDUAL FORCE            -7.970E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.22       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.797E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.870E-06   TIME AVG. MOMENT      5.072E-06
 LARGEST RESIDUAL MOMENT            1.094E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.261E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.248E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.638E-07   TIME AVG. FORCE       2.727E-07
 LARGEST RESIDUAL FORCE            -3.512E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.22       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.748E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.858E-06   TIME AVG. MOMENT      5.072E-06
 LARGEST RESIDUAL MOMENT           -5.519E-10   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.235E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.653E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.339E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.462E-07   TIME AVG. FORCE       2.702E-07
 LARGEST RESIDUAL FORCE            -1.669E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.22       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.089E-06   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     8.858E-06   TIME AVG. MOMENT      5.072E-06
 LARGEST RESIDUAL MOMENT            2.565E-12   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.235E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.753E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.928E-03,  FRACTION OF STEP COMPLETED  0.100    
 STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED        0.200    


  INCREMENT    19 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.928E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.488E-05   TIME AVG. FORCE       8.346E-06
 LARGEST RESIDUAL FORCE            -5.750E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.41       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.190       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.975E-05   TIME AVG. MOMENT      6.370E-06
 LARGEST RESIDUAL MOMENT           -6.504E-04   AT NODE        375   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.435E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.200E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.517E-06   TIME AVG. FORCE       6.761E-07
 LARGEST RESIDUAL FORCE             3.385E-05   AT NODE        228   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.38       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.534E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.047E-05   TIME AVG. MOMENT      5.356E-06
 LARGEST RESIDUAL MOMENT            2.534E-05   AT NODE        227   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.392E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.632E-03   AT NODE        104   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.752E-05   TIME AVG. FORCE       4.927E-06
 LARGEST RESIDUAL FORCE            -4.637E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.27       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.119       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.063E-05   TIME AVG. MOMENT      5.364E-06
 LARGEST RESIDUAL MOMENT            3.673E-05   AT NODE        254   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.385E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.007E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.632E-07   TIME AVG. FORCE       3.068E-07
 LARGEST RESIDUAL FORCE            -7.300E-07   AT NODE        164   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.26       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.199E-03   AT NODE        228   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.569E-06   TIME AVG. MOMENT      5.308E-06
 LARGEST RESIDUAL MOMENT           -4.515E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.381E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.731E-04   AT NODE        164   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.129E-07   TIME AVG. FORCE       3.130E-07
 LARGEST RESIDUAL FORCE            -1.360E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.26       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.322E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.556E-06   TIME AVG. MOMENT      5.308E-06
 LARGEST RESIDUAL MOMENT            7.512E-07   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.322E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.131E-04   AT NODE        277   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.848E-07   TIME AVG. FORCE       2.845E-07
 LARGEST RESIDUAL FORCE            -6.076E-10   AT NODE        220   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.26       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.851E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     9.560E-06   TIME AVG. MOMENT      5.308E-06
 LARGEST RESIDUAL MOMENT            4.765E-10   AT NODE        259   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.322E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.785E-06   AT NODE        220   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   3.301E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.928E-03,  FRACTION OF STEP COMPLETED  0.104    
 STEP TIME COMPLETED       0.208    ,  TOTAL TIME COMPLETED        0.208    


  INCREMENT    20 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.301E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.686E-07   TIME AVG. FORCE       3.161E-07
 LARGEST RESIDUAL FORCE             2.999E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.531       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.026E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.965E-06   TIME AVG. MOMENT      5.541E-06
 LARGEST RESIDUAL MOMENT           -7.314E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.006E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.914E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.234E-07   TIME AVG. FORCE       3.000E-07
 LARGEST RESIDUAL FORCE            -3.024E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.529       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.394E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.855E-06   TIME AVG. MOMENT      5.535E-06
 LARGEST RESIDUAL MOMENT            2.980E-09   AT NODE        354   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.810E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.469E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.006E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.096E-07   TIME AVG. FORCE       2.984E-07
 LARGEST RESIDUAL FORCE            -1.696E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.529       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.238E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     9.856E-06   TIME AVG. MOMENT      5.535E-06
 LARGEST RESIDUAL MOMENT            6.280E-12   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.808E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.889E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.301E-03,  FRACTION OF STEP COMPLETED  0.106    
 STEP TIME COMPLETED       0.212    ,  TOTAL TIME COMPLETED        0.212    


  INCREMENT    21 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.928E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.790E-05   TIME AVG. FORCE       2.059E-06
 LARGEST RESIDUAL FORCE             1.678E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.38       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.107       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.753E-05   TIME AVG. MOMENT      6.106E-06
 LARGEST RESIDUAL MOMENT           -2.724E-04   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.962E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.066E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.282E-06   TIME AVG. FORCE       3.968E-07
 LARGEST RESIDUAL FORCE            -2.862E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.36       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.196E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.075E-05   TIME AVG. MOMENT      5.784E-06
 LARGEST RESIDUAL MOMENT            2.836E-06   AT NODE        262   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.887E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.010E-03   AT NODE        288   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.682E-06   TIME AVG. FORCE       8.367E-07
 LARGEST RESIDUAL FORCE            -8.288E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.31       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.605E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.059E-05   TIME AVG. MOMENT      5.776E-06
 LARGEST RESIDUAL MOMENT           -2.559E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.473E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.136E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.489E-07   TIME AVG. FORCE       3.135E-07
 LARGEST RESIDUAL FORCE            -9.045E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.31       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.905E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.060E-05   TIME AVG. MOMENT      5.776E-06
 LARGEST RESIDUAL MOMENT            1.254E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.460E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.308E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.674E-07   TIME AVG. FORCE       3.153E-07
 LARGEST RESIDUAL FORCE            -4.584E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.30       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.103E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.059E-05   TIME AVG. MOMENT      5.776E-06
 LARGEST RESIDUAL MOMENT           -8.515E-10   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.429E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.052E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.073E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.429E-07   TIME AVG. FORCE       3.129E-07
 LARGEST RESIDUAL FORCE            -2.756E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.30       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.907E-06   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.059E-05   TIME AVG. MOMENT      5.776E-06
 LARGEST RESIDUAL MOMENT            4.070E-12   AT NODE        383   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.429E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.503E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.928E-03,  FRACTION OF STEP COMPLETED  0.110    
 STEP TIME COMPLETED       0.219    ,  TOTAL TIME COMPLETED        0.219    


  INCREMENT    22 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.928E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.333E-04   TIME AVG. FORCE       2.150E-05
 LARGEST RESIDUAL FORCE            -1.576E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.59       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.286       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.232E-05   TIME AVG. MOMENT      9.255E-06
 LARGEST RESIDUAL MOMENT           -2.072E-03   AT NODE        314   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.520E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.091E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.704E-05   TIME AVG. FORCE       1.833E-06
 LARGEST RESIDUAL FORCE            -2.700E-04   AT NODE         92   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.57       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.563E-02   AT NODE          1   DOF  1
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.529E-05   TIME AVG. MOMENT      6.208E-06
 LARGEST RESIDUAL MOMENT           -2.055E-04   AT NODE        119   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.487E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.313E-03   AT NODE         89   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.148E-04   TIME AVG. FORCE       1.072E-05
 LARGEST RESIDUAL FORCE            -1.231E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.35       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.221       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.991E-05   TIME AVG. MOMENT      6.418E-06
 LARGEST RESIDUAL MOMENT            2.381E-04   AT NODE        287   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.785E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.702E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.964E-06   TIME AVG. FORCE       5.539E-07
 LARGEST RESIDUAL FORCE            -2.134E-05   AT NODE        172   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.35       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.211E-03   AT NODE        224   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.158E-05   TIME AVG. MOMENT      6.040E-06
 LARGEST RESIDUAL MOMENT           -3.752E-06   AT NODE        207   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.776E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.427E-03   AT NODE        172   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.006E-05   TIME AVG. FORCE       2.108E-06
 LARGEST RESIDUAL FORCE            -1.046E-03   AT NODE        284   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.33       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.996E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.443E-05   TIME AVG. MOMENT      6.169E-06
 LARGEST RESIDUAL MOMENT            1.311E-04   AT NODE        230   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.520E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.760E-03   AT NODE        281   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.582E-07   TIME AVG. FORCE       3.534E-07
 LARGEST RESIDUAL FORCE            -9.924E-07   AT NODE        176   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.33       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.203E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.135E-05   TIME AVG. MOMENT      6.029E-06
 LARGEST RESIDUAL MOMENT            4.969E-07   AT NODE        207   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.506E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.121E-04   AT NODE        176   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      7.314E-07   TIME AVG. FORCE       3.509E-07
 LARGEST RESIDUAL FORCE            -1.786E-05   AT NODE        232   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.33       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.972E-03   AT NODE        292   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.134E-05   TIME AVG. MOMENT      6.029E-06
 LARGEST RESIDUAL MOMENT            3.158E-06   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.494E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.280E-03   AT NODE        232   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      4.844E-07   TIME AVG. FORCE       3.285E-07
 LARGEST RESIDUAL FORCE             7.705E-10   AT NODE        208   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.33       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.734E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.134E-05   TIME AVG. MOMENT      6.029E-06
 LARGEST RESIDUAL MOMENT            6.286E-10   AT NODE        207   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.494E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.127E-06   AT NODE        188   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.928E-03,  FRACTION OF STEP COMPLETED  0.114    
 STEP TIME COMPLETED       0.227    ,  TOTAL TIME COMPLETED        0.227    


  INCREMENT    23 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.928E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.784E-04   TIME AVG. FORCE       4.017E-05
 LARGEST RESIDUAL FORCE            -2.513E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.68       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.352       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.772E-04   TIME AVG. MOMENT      1.347E-05
 LARGEST RESIDUAL MOMENT           -5.260E-03   AT NODE        302   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.222E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.466E-02   AT NODE        293   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.975E-05   TIME AVG. FORCE       3.613E-06
 LARGEST RESIDUAL FORCE             6.397E-04   AT NODE        120   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.69       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.120E-02   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.480E-05   TIME AVG. MOMENT      6.845E-06
 LARGEST RESIDUAL MOMENT           -5.635E-04   AT NODE        115   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.183E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.055E-03   AT NODE         93   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.587E-04   TIME AVG. FORCE       1.352E-05
 LARGEST RESIDUAL FORCE            -1.228E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.43       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.264       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.211E-05   TIME AVG. MOMENT      6.728E-06
 LARGEST RESIDUAL MOMENT           -3.070E-04   AT NODE        315   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.418E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.765E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.978E-06   TIME AVG. FORCE       8.826E-07
 LARGEST RESIDUAL FORCE            -6.792E-05   AT NODE        176   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.42       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.166E-03   AT NODE        200   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.288E-05   TIME AVG. MOMENT      6.327E-06
 LARGEST RESIDUAL MOMENT           -2.575E-05   AT NODE         95   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.393E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.489E-03   AT NODE        175   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.765E-05   TIME AVG. FORCE       5.938E-06
 LARGEST RESIDUAL FORCE            -3.559E-03   AT NODE        292   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.35       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.983E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.987E-05   TIME AVG. MOMENT      7.065E-06
 LARGEST RESIDUAL MOMENT            5.880E-04   AT NODE        263   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.648E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.779E-02   AT NODE        289   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      2.391E-06   TIME AVG. FORCE       5.004E-07
 LARGEST RESIDUAL FORCE            -1.135E-05   AT NODE        180   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.35       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.526E-03   AT NODE        196   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.225E-05   TIME AVG. MOMENT      6.299E-06
 LARGEST RESIDUAL MOMENT           -4.196E-06   AT NODE        203   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.635E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.046E-03   AT NODE        180   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      8.256E-06   TIME AVG. FORCE       9.891E-07
 LARGEST RESIDUAL FORCE            -4.445E-04   AT NODE        236   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.35       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.021E-02   AT NODE        296   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.325E-05   TIME AVG. MOMENT      6.343E-06
 LARGEST RESIDUAL MOMENT            8.235E-05   AT NODE        202   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.559E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.322E-03   AT NODE        236   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      5.767E-07   TIME AVG. FORCE       3.492E-07
 LARGEST RESIDUAL FORCE            -3.277E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.35       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.562E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.211E-05   TIME AVG. MOMENT      6.293E-06
 LARGEST RESIDUAL MOMENT            1.841E-07   AT NODE        199   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.551E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.382E-04   AT NODE        176   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      5.375E-07   TIME AVG. FORCE       3.459E-07
 LARGEST RESIDUAL FORCE            -1.088E-06   AT NODE        204   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.35       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.710E-04   AT NODE        252   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.211E-05   TIME AVG. MOMENT      6.293E-06
 LARGEST RESIDUAL MOMENT            2.161E-07   AT NODE        174   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.550E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.149E-04   AT NODE        204   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    23 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  3.964E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.027E-06   TIME AVG. FORCE       3.867E-07
 LARGEST RESIDUAL FORCE             1.081E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.677       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.657E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.214E-05   TIME AVG. MOMENT      6.294E-06
 LARGEST RESIDUAL MOMENT           -2.501E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.338E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.124E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.417E-07   TIME AVG. FORCE       3.546E-07
 LARGEST RESIDUAL FORCE            -2.553E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.669       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.525E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.173E-05   TIME AVG. MOMENT      6.276E-06
 LARGEST RESIDUAL MOMENT            3.126E-08   AT NODE        346   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.266E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.253E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -9.068E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.115E-07   TIME AVG. FORCE       3.437E-07
 LARGEST RESIDUAL FORCE            -6.935E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.669       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.703E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.172E-05   TIME AVG. MOMENT      6.276E-06
 LARGEST RESIDUAL MOMENT            1.138E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.262E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.822E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.113E-07   TIME AVG. FORCE       3.437E-07
 LARGEST RESIDUAL FORCE            -8.367E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.669       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.132E-06   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.172E-05   TIME AVG. MOMENT      6.276E-06
 LARGEST RESIDUAL MOMENT           -6.631E-13   AT NODE        393   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.262E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.307E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.964E-03,  FRACTION OF STEP COMPLETED  0.116    
 STEP TIME COMPLETED       0.231    ,  TOTAL TIME COMPLETED        0.231    


  INCREMENT    24 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.964E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.346E-07   TIME AVG. FORCE       3.738E-07
 LARGEST RESIDUAL FORCE             3.094E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.676       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.112E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.224E-05   TIME AVG. MOMENT      6.525E-06
 LARGEST RESIDUAL MOMENT           -7.601E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.302E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.052E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.550E-07   TIME AVG. FORCE       3.600E-07
 LARGEST RESIDUAL FORCE            -3.454E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.674       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.722E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.211E-05   TIME AVG. MOMENT      6.519E-06
 LARGEST RESIDUAL MOMENT            3.284E-09   AT NODE        350   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.276E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.670E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.154E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.377E-07   TIME AVG. FORCE       3.586E-07
 LARGEST RESIDUAL FORCE            -3.335E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.674       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.877E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.211E-05   TIME AVG. MOMENT      6.519E-06
 LARGEST RESIDUAL MOMENT            1.101E-11   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.275E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.657E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.946E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.964E-03,  FRACTION OF STEP COMPLETED  0.118    
 STEP TIME COMPLETED       0.235    ,  TOTAL TIME COMPLETED        0.235    


  INCREMENT    25 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.542E-06   TIME AVG. FORCE       5.860E-07
 LARGEST RESIDUAL FORCE             4.191E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.04       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.395E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.442E-05   TIME AVG. MOMENT      6.835E-06
 LARGEST RESIDUAL MOMENT           -8.349E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.146E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.332E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.230E-06   TIME AVG. FORCE       4.209E-07
 LARGEST RESIDUAL FORCE            -8.222E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.791E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.271E-05   TIME AVG. MOMENT      6.767E-06
 LARGEST RESIDUAL MOMENT            3.126E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.016E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.305E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.725E-07   TIME AVG. FORCE       3.810E-07
 LARGEST RESIDUAL FORCE            -1.988E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.700E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.271E-05   TIME AVG. MOMENT      6.767E-06
 LARGEST RESIDUAL MOMENT            5.357E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.951E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.478E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.745E-07   TIME AVG. FORCE       3.741E-07
 LARGEST RESIDUAL FORCE            -5.958E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.555E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.270E-05   TIME AVG. MOMENT      6.766E-06
 LARGEST RESIDUAL MOMENT            2.605E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.940E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.111E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.717E-07   TIME AVG. FORCE       3.738E-07
 LARGEST RESIDUAL FORCE            -8.634E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.465E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.269E-05   TIME AVG. MOMENT      6.766E-06
 LARGEST RESIDUAL MOMENT            3.463E-12   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.939E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.233E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.121    
 STEP TIME COMPLETED       0.241    ,  TOTAL TIME COMPLETED        0.241    


  INCREMENT    26 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.139E-06   TIME AVG. FORCE       8.249E-07
 LARGEST RESIDUAL FORCE             9.360E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.07       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.830E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.679E-05   TIME AVG. MOMENT      7.152E-06
 LARGEST RESIDUAL MOMENT           -1.687E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.325E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.857E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.619E-06   TIME AVG. FORCE       4.568E-07
 LARGEST RESIDUAL FORCE            -8.925E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.05       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.274E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.336E-05   TIME AVG. MOMENT      7.020E-06
 LARGEST RESIDUAL MOMENT            1.045E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.189E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.366E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.474E-06   TIME AVG. FORCE       4.472E-07
 LARGEST RESIDUAL FORCE            -1.615E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.636E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.327E-05   TIME AVG. MOMENT      7.016E-06
 LARGEST RESIDUAL MOMENT           -5.985E-08   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.003E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.860E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.260E-07   TIME AVG. FORCE       3.907E-07
 LARGEST RESIDUAL FORCE            -4.229E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.654E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.330E-05   TIME AVG. MOMENT      7.017E-06
 LARGEST RESIDUAL MOMENT            5.733E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.973E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.963E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.057E-07   TIME AVG. FORCE       3.893E-07
 LARGEST RESIDUAL FORCE            -4.216E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.681E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.329E-05   TIME AVG. MOMENT      7.017E-06
 LARGEST RESIDUAL MOMENT            4.442E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.964E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.421E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.040E-07   TIME AVG. FORCE       3.892E-07
 LARGEST RESIDUAL FORCE            -2.696E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.198E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.329E-05   TIME AVG. MOMENT      7.017E-06
 LARGEST RESIDUAL MOMENT            8.902E-13   AT NODE        383   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.963E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.616E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.124    
 STEP TIME COMPLETED       0.247    ,  TOTAL TIME COMPLETED        0.247    


  INCREMENT    27 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.678E-06   TIME AVG. FORCE       8.447E-07
 LARGEST RESIDUAL FORCE             1.028E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.07       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.792E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.764E-05   TIME AVG. MOMENT      7.411E-06
 LARGEST RESIDUAL MOMENT           -1.819E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.370E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.063E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.749E-06   TIME AVG. FORCE       4.742E-07
 LARGEST RESIDUAL FORCE            -8.751E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.05       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.307E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.397E-05   TIME AVG. MOMENT      7.274E-06
 LARGEST RESIDUAL MOMENT            1.164E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.234E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.358E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.725E-06   TIME AVG. FORCE       4.727E-07
 LARGEST RESIDUAL FORCE            -1.958E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.791E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.387E-05   TIME AVG. MOMENT      7.271E-06
 LARGEST RESIDUAL MOMENT           -1.280E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.028E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.058E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.686E-07   TIME AVG. FORCE       4.067E-07
 LARGEST RESIDUAL FORCE            -4.631E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.790E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.390E-05   TIME AVG. MOMENT      7.272E-06
 LARGEST RESIDUAL MOMENT            7.382E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.997E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.111E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.477E-07   TIME AVG. FORCE       4.053E-07
 LARGEST RESIDUAL FORCE            -6.936E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.810E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.388E-05   TIME AVG. MOMENT      7.271E-06
 LARGEST RESIDUAL MOMENT            6.084E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.985E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.214E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.442E-07   TIME AVG. FORCE       4.051E-07
 LARGEST RESIDUAL FORCE            -5.747E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.793E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.388E-05   TIME AVG. MOMENT      7.271E-06
 LARGEST RESIDUAL MOMENT            1.112E-12   AT NODE        386   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.985E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.089E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.127    
 STEP TIME COMPLETED       0.253    ,  TOTAL TIME COMPLETED        0.253    


  INCREMENT    28 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.151E-06   TIME AVG. FORCE       8.607E-07
 LARGEST RESIDUAL FORCE             1.110E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.08       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.722E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.847E-05   TIME AVG. MOMENT      7.671E-06
 LARGEST RESIDUAL MOMENT           -1.934E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.408E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.231E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.950E-06   TIME AVG. FORCE       4.960E-07
 LARGEST RESIDUAL FORCE            -8.683E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.05       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.335E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.458E-05   TIME AVG. MOMENT      7.532E-06
 LARGEST RESIDUAL MOMENT            1.271E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.272E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.359E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.974E-06   TIME AVG. FORCE       4.974E-07
 LARGEST RESIDUAL FORCE            -2.279E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.920E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.447E-05   TIME AVG. MOMENT      7.528E-06
 LARGEST RESIDUAL MOMENT           -2.019E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.049E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.231E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.115E-07   TIME AVG. FORCE       4.232E-07
 LARGEST RESIDUAL FORCE            -5.016E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.913E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.450E-05   TIME AVG. MOMENT      7.529E-06
 LARGEST RESIDUAL MOMENT            9.078E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.016E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.249E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.909E-07   TIME AVG. FORCE       4.219E-07
 LARGEST RESIDUAL FORCE            -1.040E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.196E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.449E-05   TIME AVG. MOMENT      7.529E-06
 LARGEST RESIDUAL MOMENT            7.792E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.002E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.495E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.848E-07   TIME AVG. FORCE       4.216E-07
 LARGEST RESIDUAL FORCE            -1.042E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.177E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.449E-05   TIME AVG. MOMENT      7.529E-06
 LARGEST RESIDUAL MOMENT            2.505E-12   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.001E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.470E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.130    
 STEP TIME COMPLETED       0.259    ,  TOTAL TIME COMPLETED        0.259    


  INCREMENT    29 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.495E-06   TIME AVG. FORCE       8.701E-07
 LARGEST RESIDUAL FORCE             1.173E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.07       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.850E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.922E-05   TIME AVG. MOMENT      7.932E-06
 LARGEST RESIDUAL MOMENT           -2.018E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.436E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.346E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.091E-06   TIME AVG. FORCE       5.143E-07
 LARGEST RESIDUAL FORCE            -8.702E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.05       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.356E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.519E-05   TIME AVG. MOMENT      7.793E-06
 LARGEST RESIDUAL MOMENT            1.352E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.299E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.366E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.193E-06   TIME AVG. FORCE       5.200E-07
 LARGEST RESIDUAL FORCE            -2.520E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.004E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.508E-05   TIME AVG. MOMENT      7.789E-06
 LARGEST RESIDUAL MOMENT           -2.632E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.063E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.359E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.581E-07   TIME AVG. FORCE       4.403E-07
 LARGEST RESIDUAL FORCE            -5.363E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.014E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.511E-05   TIME AVG. MOMENT      7.790E-06
 LARGEST RESIDUAL MOMENT            1.054E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.030E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.373E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.374E-07   TIME AVG. FORCE       4.391E-07
 LARGEST RESIDUAL FORCE            -1.389E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.374E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.509E-05   TIME AVG. MOMENT      7.790E-06
 LARGEST RESIDUAL MOMENT            9.344E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.012E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.737E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.290E-07   TIME AVG. FORCE       4.387E-07
 LARGEST RESIDUAL FORCE            -1.629E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.466E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.509E-05   TIME AVG. MOMENT      7.790E-06
 LARGEST RESIDUAL MOMENT            3.460E-12   AT NODE        390   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.012E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.848E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.133    
 STEP TIME COMPLETED       0.265    ,  TOTAL TIME COMPLETED        0.265    


  INCREMENT    30 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.681E-06   TIME AVG. FORCE       8.725E-07
 LARGEST RESIDUAL FORCE             1.212E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.07       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.928E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.992E-05   TIME AVG. MOMENT      8.194E-06
 LARGEST RESIDUAL MOMENT           -2.067E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.453E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.405E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.193E-06   TIME AVG. FORCE       5.310E-07
 LARGEST RESIDUAL FORCE            -8.774E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.05       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.366E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.580E-05   TIME AVG. MOMENT      8.057E-06
 LARGEST RESIDUAL MOMENT            1.401E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.315E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.378E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.432E-06   TIME AVG. FORCE       5.436E-07
 LARGEST RESIDUAL FORCE            -2.649E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.039E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.569E-05   TIME AVG. MOMENT      8.053E-06
 LARGEST RESIDUAL MOMENT           -2.991E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.072E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.433E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      8.118E-07   TIME AVG. FORCE       4.583E-07
 LARGEST RESIDUAL FORCE            -5.658E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.086E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.570E-05   TIME AVG. MOMENT      8.054E-06
 LARGEST RESIDUAL MOMENT            1.157E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.037E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.480E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.885E-07   TIME AVG. FORCE       4.571E-07
 LARGEST RESIDUAL FORCE            -1.653E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.490E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.570E-05   TIME AVG. MOMENT      8.053E-06
 LARGEST RESIDUAL MOMENT            1.060E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.018E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.905E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.776E-07   TIME AVG. FORCE       4.565E-07
 LARGEST RESIDUAL FORCE            -2.208E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.698E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.570E-05   TIME AVG. MOMENT      8.053E-06
 LARGEST RESIDUAL MOMENT            4.662E-12   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.017E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.164E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.672E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.135    
 STEP TIME COMPLETED       0.271    ,  TOTAL TIME COMPLETED        0.271    


  INCREMENT    31 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.672E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.875E-06   TIME AVG. FORCE       5.274E-07
 LARGEST RESIDUAL FORCE             2.160E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.811       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.158E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.693E-05   TIME AVG. MOMENT      8.340E-06
 LARGEST RESIDUAL MOMENT           -4.778E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.728E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.425E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.391E-06   TIME AVG. FORCE       5.032E-07
 LARGEST RESIDUAL FORCE            -6.809E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.798       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.292E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.620E-05   TIME AVG. MOMENT      8.316E-06
 LARGEST RESIDUAL MOMENT           -9.188E-08   AT NODE        342   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.605E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.228E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.846E-07   TIME AVG. FORCE       4.779E-07
 LARGEST RESIDUAL FORCE            -1.326E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.797       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.219E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.619E-05   TIME AVG. MOMENT      8.316E-06
 LARGEST RESIDUAL MOMENT            1.467E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.588E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.723E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.750E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      8.756E-07   TIME AVG. FORCE       4.774E-07
 LARGEST RESIDUAL FORCE            -1.839E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.797       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.419E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.618E-05   TIME AVG. MOMENT      8.315E-06
 LARGEST RESIDUAL MOMENT            3.357E-11   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.586E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.015E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      8.756E-07   TIME AVG. FORCE       4.774E-07
 LARGEST RESIDUAL FORCE            -3.695E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.797       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.811E-07   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.618E-05   TIME AVG. MOMENT      8.315E-06
 LARGEST RESIDUAL MOMENT            9.462E-13   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.586E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.139E-08   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.672E-03,  FRACTION OF STEP COMPLETED  0.138    
 STEP TIME COMPLETED       0.276    ,  TOTAL TIME COMPLETED        0.276    


  INCREMENT    32 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.257E-06   TIME AVG. FORCE       7.050E-07
 LARGEST RESIDUAL FORCE             7.422E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.04       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.907E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.943E-05   TIME AVG. MOMENT      8.663E-06
 LARGEST RESIDUAL MOMENT           -1.364E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.339E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.211E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.222E-06   TIME AVG. FORCE       5.605E-07
 LARGEST RESIDUAL FORCE            -1.076E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.175E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.684E-05   TIME AVG. MOMENT      8.582E-06
 LARGEST RESIDUAL MOMENT            7.495E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.185E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.545E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.513E-06   TIME AVG. FORCE       5.267E-07
 LARGEST RESIDUAL FORCE            -8.001E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.01       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.087E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.682E-05   TIME AVG. MOMENT      8.581E-06
 LARGEST RESIDUAL MOMENT            9.356E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.050E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.348E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.010E-06   TIME AVG. FORCE       5.028E-07
 LARGEST RESIDUAL FORCE            -4.076E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.01       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.448E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.682E-05   TIME AVG. MOMENT      8.581E-06
 LARGEST RESIDUAL MOMENT            3.011E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.020E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.004E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    9.668E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.801E-07   TIME AVG. FORCE       5.014E-07
 LARGEST RESIDUAL FORCE            -1.090E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.01       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.744E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.679E-05   TIME AVG. MOMENT      8.580E-06
 LARGEST RESIDUAL MOMENT            2.535E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.015E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.935E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      9.794E-07   TIME AVG. FORCE       5.014E-07
 LARGEST RESIDUAL FORCE             6.165E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.01       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.053E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.679E-05   TIME AVG. MOMENT      8.580E-06
 LARGEST RESIDUAL MOMENT            6.850E-13   AT NODE        381   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.015E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.034E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.141    
 STEP TIME COMPLETED       0.282    ,  TOTAL TIME COMPLETED        0.282    


  INCREMENT    33 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.380E-06   TIME AVG. FORCE       8.595E-07
 LARGEST RESIDUAL FORCE             1.193E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.04       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.854E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.147E-05   TIME AVG. MOMENT      8.971E-06
 LARGEST RESIDUAL MOMENT           -2.022E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.444E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.293E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.381E-06   TIME AVG. FORCE       5.868E-07
 LARGEST RESIDUAL FORCE            -9.320E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.346E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.750E-05   TIME AVG. MOMENT      8.851E-06
 LARGEST RESIDUAL MOMENT            1.377E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.300E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.443E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.640E-06   TIME AVG. FORCE       5.986E-07
 LARGEST RESIDUAL FORCE            -2.401E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.997       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.890E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.740E-05   TIME AVG. MOMENT      8.847E-06
 LARGEST RESIDUAL MOMENT           -2.316E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.064E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.355E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.042E-06   TIME AVG. FORCE       5.259E-07
 LARGEST RESIDUAL FORCE            -6.394E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.994       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.195E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.741E-05   TIME AVG. MOMENT      8.848E-06
 LARGEST RESIDUAL MOMENT            1.169E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.027E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.760E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.002E-06   TIME AVG. FORCE       5.241E-07
 LARGEST RESIDUAL FORCE            -1.617E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.992       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.443E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.740E-05   TIME AVG. MOMENT      8.848E-06
 LARGEST RESIDUAL MOMENT            1.290E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.007E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.916E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      9.883E-07   TIME AVG. FORCE       5.235E-07
 LARGEST RESIDUAL FORCE            -2.922E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.992       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.913E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.740E-05   TIME AVG. MOMENT      8.848E-06
 LARGEST RESIDUAL MOMENT            5.419E-12   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.007E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.531E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.144    
 STEP TIME COMPLETED       0.288    ,  TOTAL TIME COMPLETED        0.288    


  INCREMENT    34 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.973E-06   TIME AVG. FORCE       8.474E-07
 LARGEST RESIDUAL FORCE             1.144E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.720E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.187E-05   TIME AVG. MOMENT      9.231E-06
 LARGEST RESIDUAL MOMENT           -1.948E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.424E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.164E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.450E-06   TIME AVG. FORCE       6.072E-07
 LARGEST RESIDUAL FORCE            -9.672E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.999       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.323E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.809E-05   TIME AVG. MOMENT      9.119E-06
 LARGEST RESIDUAL MOMENT            1.314E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.276E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.480E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.403E-06   TIME AVG. FORCE       6.052E-07
 LARGEST RESIDUAL FORCE            -2.139E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.982       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.763E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.799E-05   TIME AVG. MOMENT      9.117E-06
 LARGEST RESIDUAL MOMENT           -1.603E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.052E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.237E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.049E-06   TIME AVG. FORCE       5.463E-07
 LARGEST RESIDUAL FORCE            -6.633E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.978       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.205E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.803E-05   TIME AVG. MOMENT      9.118E-06
 LARGEST RESIDUAL MOMENT            1.079E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.013E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.858E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.019E-06   TIME AVG. FORCE       5.450E-07
 LARGEST RESIDUAL FORCE            -1.356E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.977       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.310E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.800E-05   TIME AVG. MOMENT      9.117E-06
 LARGEST RESIDUAL MOMENT            1.305E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.996E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.766E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.008E-06   TIME AVG. FORCE       5.445E-07
 LARGEST RESIDUAL FORCE            -2.648E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.977       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.818E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.800E-05   TIME AVG. MOMENT      9.117E-06
 LARGEST RESIDUAL MOMENT            4.503E-12   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.995E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.444E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.147    
 STEP TIME COMPLETED       0.293    ,  TOTAL TIME COMPLETED        0.293    


  INCREMENT    35 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.447E-06   TIME AVG. FORCE       8.321E-07
 LARGEST RESIDUAL FORCE             1.071E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.00       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.521E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.221E-05   TIME AVG. MOMENT      9.491E-06
 LARGEST RESIDUAL MOMENT           -1.843E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.394E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.983E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.519E-06   TIME AVG. FORCE       6.268E-07
 LARGEST RESIDUAL FORCE            -1.018E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.980       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.296E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.868E-05   TIME AVG. MOMENT      9.390E-06
 LARGEST RESIDUAL MOMENT            1.221E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.241E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.529E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.219E-06   TIME AVG. FORCE       6.143E-07
 LARGEST RESIDUAL FORCE            -1.810E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.964       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.601E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.859E-05   TIME AVG. MOMENT      9.387E-06
 LARGEST RESIDUAL MOMENT           -7.453E-08   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.034E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.071E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.084E-06   TIME AVG. FORCE       5.670E-07
 LARGEST RESIDUAL FORCE            -6.851E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.960       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.201E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.862E-05   TIME AVG. MOMENT      9.388E-06
 LARGEST RESIDUAL MOMENT            9.438E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.994E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.952E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.047E-06   TIME AVG. FORCE       5.655E-07
 LARGEST RESIDUAL FORCE            -1.019E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.959       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.125E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.860E-05   TIME AVG. MOMENT      9.388E-06
 LARGEST RESIDUAL MOMENT            1.249E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.979E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.541E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.039E-06   TIME AVG. FORCE       5.651E-07
 LARGEST RESIDUAL FORCE            -2.205E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.959       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.627E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.860E-05   TIME AVG. MOMENT      9.388E-06
 LARGEST RESIDUAL MOMENT            3.515E-12   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.978E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.224E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.150    
 STEP TIME COMPLETED       0.299    ,  TOTAL TIME COMPLETED        0.299    


  INCREMENT    36 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.807E-06   TIME AVG. FORCE       8.148E-07
 LARGEST RESIDUAL FORCE             9.762E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.980       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.256E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.250E-05   TIME AVG. MOMENT      9.752E-06
 LARGEST RESIDUAL MOMENT           -1.708E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.353E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.748E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.595E-06   TIME AVG. FORCE       6.463E-07
 LARGEST RESIDUAL FORCE            -1.082E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.957       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.265E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.925E-05   TIME AVG. MOMENT      9.662E-06
 LARGEST RESIDUAL MOMENT            1.093E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.194E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.589E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.018E-06   TIME AVG. FORCE       6.232E-07
 LARGEST RESIDUAL FORCE            -1.437E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.943       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.407E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.919E-05   TIME AVG. MOMENT      9.660E-06
 LARGEST RESIDUAL MOMENT            5.241E-08   AT NODE        351   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.009E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.858E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.123E-06   TIME AVG. FORCE       5.874E-07
 LARGEST RESIDUAL FORCE            -6.876E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.940       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.145E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.922E-05   TIME AVG. MOMENT      9.661E-06
 LARGEST RESIDUAL MOMENT            7.622E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.969E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.994E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.084E-06   TIME AVG. FORCE       5.859E-07
 LARGEST RESIDUAL FORCE            -6.555E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.939       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.927E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.919E-05   TIME AVG. MOMENT      9.660E-06
 LARGEST RESIDUAL MOMENT            1.074E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.956E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.244E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.078E-06   TIME AVG. FORCE       5.857E-07
 LARGEST RESIDUAL FORCE            -1.416E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.939       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.308E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.919E-05   TIME AVG. MOMENT      9.660E-06
 LARGEST RESIDUAL MOMENT            2.037E-12   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.956E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.820E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.153    
 STEP TIME COMPLETED       0.305    ,  TOTAL TIME COMPLETED        0.305    


  INCREMENT    37 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.109E-06   TIME AVG. FORCE       7.981E-07
 LARGEST RESIDUAL FORCE             8.653E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.954       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.929E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.272E-05   TIME AVG. MOMENT      1.001E-05
 LARGEST RESIDUAL MOMENT           -1.547E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.303E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.470E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.668E-06   TIME AVG. FORCE       6.657E-07
 LARGEST RESIDUAL FORCE            -1.151E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.932       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.223E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.984E-05   TIME AVG. MOMENT      9.935E-06
 LARGEST RESIDUAL MOMENT            9.348E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.138E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.653E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.820E-06   TIME AVG. FORCE       6.331E-07
 LARGEST RESIDUAL FORCE            -1.064E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.920       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.194E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.980E-05   TIME AVG. MOMENT      9.934E-06
 LARGEST RESIDUAL MOMENT            9.640E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.977E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.610E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.163E-06   TIME AVG. FORCE       6.079E-07
 LARGEST RESIDUAL FORCE            -6.489E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.917       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.995E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.980E-05   TIME AVG. MOMENT      9.934E-06
 LARGEST RESIDUAL MOMENT            5.554E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.938E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.915E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.123E-06   TIME AVG. FORCE       6.063E-07
 LARGEST RESIDUAL FORCE            -3.451E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.917       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.408E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.978E-05   TIME AVG. MOMENT      9.933E-06
 LARGEST RESIDUAL MOMENT            7.772E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.929E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.085E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.121E-06   TIME AVG. FORCE       6.062E-07
 LARGEST RESIDUAL FORCE            -6.955E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.917       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.933E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.978E-05   TIME AVG. MOMENT      9.933E-06
 LARGEST RESIDUAL MOMENT            8.879E-13   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.928E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.266E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.156    
 STEP TIME COMPLETED       0.311    ,  TOTAL TIME COMPLETED        0.311    


  INCREMENT    38 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.401E-06   TIME AVG. FORCE       7.838E-07
 LARGEST RESIDUAL FORCE             7.457E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.927       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.550E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.291E-05   TIME AVG. MOMENT      1.027E-05
 LARGEST RESIDUAL MOMENT           -1.371E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.244E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.159E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.722E-06   TIME AVG. FORCE       6.846E-07
 LARGEST RESIDUAL FORCE            -1.217E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.905       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.170E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.040E-05   TIME AVG. MOMENT      1.021E-05
 LARGEST RESIDUAL MOMENT            7.551E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.073E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.714E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.653E-06   TIME AVG. FORCE       6.450E-07
 LARGEST RESIDUAL FORCE            -7.310E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.895       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.740E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.037E-05   TIME AVG. MOMENT      1.021E-05
 LARGEST RESIDUAL MOMENT            1.363E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.939E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.343E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.203E-06   TIME AVG. FORCE       6.283E-07
 LARGEST RESIDUAL FORCE            -5.570E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.892       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.721E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.037E-05   TIME AVG. MOMENT      1.021E-05
 LARGEST RESIDUAL MOMENT            3.574E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.902E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.660E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    9.595E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.168E-06   TIME AVG. FORCE       6.270E-07
 LARGEST RESIDUAL FORCE            -1.425E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.892       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.073E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.035E-05   TIME AVG. MOMENT      1.021E-05
 LARGEST RESIDUAL MOMENT            4.468E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.896E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.876E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.167E-06   TIME AVG. FORCE       6.270E-07
 LARGEST RESIDUAL FORCE            -2.751E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.892       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.916E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.035E-05   TIME AVG. MOMENT      1.021E-05
 LARGEST RESIDUAL MOMENT            1.230E-12   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.896E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.099E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.159    
 STEP TIME COMPLETED       0.317    ,  TOTAL TIME COMPLETED        0.317    


  INCREMENT    39 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.724E-06   TIME AVG. FORCE       7.733E-07
 LARGEST RESIDUAL FORCE             6.241E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.897       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.126E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.308E-05   TIME AVG. MOMENT      1.054E-05
 LARGEST RESIDUAL MOMENT           -1.187E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.179E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.826E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.778E-06   TIME AVG. FORCE       7.038E-07
 LARGEST RESIDUAL FORCE            -1.269E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.876       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.102E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.094E-05   TIME AVG. MOMENT      1.048E-05
 LARGEST RESIDUAL MOMENT            5.653E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.002E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.766E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.531E-06   TIME AVG. FORCE       6.593E-07
 LARGEST RESIDUAL FORCE            -4.618E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.868       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.614E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.094E-05   TIME AVG. MOMENT      1.048E-05
 LARGEST RESIDUAL MOMENT            1.465E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.894E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.075E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.248E-06   TIME AVG. FORCE       6.492E-07
 LARGEST RESIDUAL FORCE            -4.204E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.866       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.320E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.092E-05   TIME AVG. MOMENT      1.048E-05
 LARGEST RESIDUAL MOMENT            1.981E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.862E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.208E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    4.340E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.218E-06   TIME AVG. FORCE       6.481E-07
 LARGEST RESIDUAL FORCE            -4.379E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.866       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.232E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.091E-05   TIME AVG. MOMENT      1.048E-05
 LARGEST RESIDUAL MOMENT            1.939E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.859E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.278E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.218E-06   TIME AVG. FORCE       6.481E-07
 LARGEST RESIDUAL FORCE             3.852E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.866       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.052E-06   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.091E-05   TIME AVG. MOMENT      1.048E-05
 LARGEST RESIDUAL MOMENT            7.085E-13   AT NODE        375   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.859E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.019E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.162    
 STEP TIME COMPLETED       0.323    ,  TOTAL TIME COMPLETED        0.323    


  INCREMENT    40 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.049E-06   TIME AVG. FORCE       7.654E-07
 LARGEST RESIDUAL FORCE             5.056E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.865       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.666E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.325E-05   TIME AVG. MOMENT      1.080E-05
 LARGEST RESIDUAL MOMENT           -1.001E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.107E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.477E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.841E-06   TIME AVG. FORCE       7.237E-07
 LARGEST RESIDUAL FORCE            -1.293E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.845       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.014E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.149E-05   TIME AVG. MOMENT      1.076E-05
 LARGEST RESIDUAL MOMENT            3.791E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.927E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.799E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.468E-06   TIME AVG. FORCE       6.764E-07
 LARGEST RESIDUAL FORCE            -2.643E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.839       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.657E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.148E-05   TIME AVG. MOMENT      1.076E-05
 LARGEST RESIDUAL MOMENT            1.319E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.845E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.184E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.305E-06   TIME AVG. FORCE       6.707E-07
 LARGEST RESIDUAL FORCE            -2.681E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.837       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.821E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.147E-05   TIME AVG. MOMENT      1.076E-05
 LARGEST RESIDUAL MOMENT            9.164E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.819E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.585E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    1.796E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.283E-06   TIME AVG. FORCE       6.700E-07
 LARGEST RESIDUAL FORCE            -9.329E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.837       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.018E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.146E-05   TIME AVG. MOMENT      1.076E-05
 LARGEST RESIDUAL MOMENT            5.902E-11   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.818E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.522E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.284E-06   TIME AVG. FORCE       6.700E-07
 LARGEST RESIDUAL FORCE             3.067E-12   AT NODE        395   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.837       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.029E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.146E-05   TIME AVG. MOMENT      1.076E-05
 LARGEST RESIDUAL MOMENT           -7.541E-13   AT NODE        235   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.818E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.034E-08   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.165    
 STEP TIME COMPLETED       0.329    ,  TOTAL TIME COMPLETED        0.329    


  INCREMENT    41 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.473E-06   TIME AVG. FORCE       7.635E-07
 LARGEST RESIDUAL FORCE             3.948E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.832       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.179E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.341E-05   TIME AVG. MOMENT      1.107E-05
 LARGEST RESIDUAL MOMENT           -8.201E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.030E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.121E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.859E-06   TIME AVG. FORCE       7.430E-07
 LARGEST RESIDUAL FORCE            -1.272E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.813       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.903E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.202E-05   TIME AVG. MOMENT      1.103E-05
 LARGEST RESIDUAL MOMENT            2.901E-07   AT NODE        119   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.850E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.802E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.477E-06   TIME AVG. FORCE       6.969E-07
 LARGEST RESIDUAL FORCE            -1.341E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.809       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.948E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.201E-05   TIME AVG. MOMENT      1.103E-05
 LARGEST RESIDUAL MOMENT            1.017E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.791E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.861E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.406E-06   TIME AVG. FORCE       6.945E-07
 LARGEST RESIDUAL FORCE            -1.377E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.808       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.285E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.200E-05   TIME AVG. MOMENT      1.103E-05
 LARGEST RESIDUAL MOMENT            3.386E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.772E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.869E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    6.223E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.392E-06   TIME AVG. FORCE       6.941E-07
 LARGEST RESIDUAL FORCE            -1.250E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.807       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.681E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.200E-05   TIME AVG. MOMENT      1.103E-05
 LARGEST RESIDUAL MOMENT            1.109E-11   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.772E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.604E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   1.943E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.168    
 STEP TIME COMPLETED       0.335    ,  TOTAL TIME COMPLETED        0.335    


  INCREMENT    42 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.943E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.298E-06   TIME AVG. FORCE       7.136E-07
 LARGEST RESIDUAL FORCE            -9.473E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.260       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.074E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.215E-05   TIME AVG. MOMENT      1.130E-05
 LARGEST RESIDUAL MOMENT            4.420E-07   AT NODE        283   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.755E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.006E-05   AT NODE        204   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.294E-06   TIME AVG. FORCE       7.134E-07
 LARGEST RESIDUAL FORCE            -1.266E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.260       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.054E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.217E-05   TIME AVG. MOMENT      1.130E-05
 LARGEST RESIDUAL MOMENT           -1.602E-11   AT NODE        154   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.737E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.823E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.294E-06   TIME AVG. FORCE       7.134E-07
 LARGEST RESIDUAL FORCE             3.594E-12   AT NODE        366   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.260       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.035E-08   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.217E-05   TIME AVG. MOMENT      1.130E-05
 LARGEST RESIDUAL MOMENT           -8.359E-13   AT NODE        361   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.737E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.631E-09   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.943E-03,  FRACTION OF STEP COMPLETED  0.169    
 STEP TIME COMPLETED       0.337    ,  TOTAL TIME COMPLETED        0.337    


  INCREMENT    43 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.946E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.593E-06   TIME AVG. FORCE       7.409E-07
 LARGEST RESIDUAL FORCE             1.537E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.773       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.209E-02   AT NODE        384   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.275E-05   TIME AVG. MOMENT      1.156E-05
 LARGEST RESIDUAL MOMENT            8.389E-06   AT NODE        235   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.778E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.272E-04   AT NODE        204   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.625E-06   TIME AVG. FORCE       7.419E-07
 LARGEST RESIDUAL FORCE            -2.088E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.766       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.005E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.269E-05   TIME AVG. MOMENT      1.156E-05
 LARGEST RESIDUAL MOMENT           -4.239E-08   AT NODE        131   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.704E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.418E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -3.749E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.428E-06   TIME AVG. FORCE       7.358E-07
 LARGEST RESIDUAL FORCE             2.080E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.766       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.353E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.268E-05   TIME AVG. MOMENT      1.156E-05
 LARGEST RESIDUAL MOMENT           -5.218E-11   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.704E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.670E-06   AT NODE        380   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    8.919E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.946E-03,  FRACTION OF STEP COMPLETED  0.171    
 STEP TIME COMPLETED       0.343    ,  TOTAL TIME COMPLETED        0.343    


  INCREMENT    44 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.795E-05   TIME AVG. FORCE       2.773E-06
 LARGEST RESIDUAL FORCE             8.020E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.22       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.138       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.322E-05   TIME AVG. MOMENT      1.228E-05
 LARGEST RESIDUAL MOMENT           -7.501E-04   AT NODE        390   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.835E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.279E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.871E-06   TIME AVG. FORCE       8.308E-07
 LARGEST RESIDUAL FORCE             2.236E-05   AT NODE        236   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.20       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.836E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.413E-05   TIME AVG. MOMENT      1.185E-05
 LARGEST RESIDUAL MOMENT           -1.954E-05   AT NODE        239   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.796E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.162E-03   AT NODE        252   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.808E-05   TIME AVG. FORCE       2.777E-06
 LARGEST RESIDUAL FORCE             6.171E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.09       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.123       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.448E-05   TIME AVG. MOMENT      1.185E-05
 LARGEST RESIDUAL MOMENT           -7.904E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.512E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.284E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.496E-06   TIME AVG. FORCE       7.588E-07
 LARGEST RESIDUAL FORCE            -1.006E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.09       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.857E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.347E-05   TIME AVG. MOMENT      1.183E-05
 LARGEST RESIDUAL MOMENT           -2.494E-07   AT NODE        267   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.509E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.525E-04   AT NODE        284   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.802E-06   TIME AVG. FORCE       7.681E-07
 LARGEST RESIDUAL FORCE            -2.382E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.08       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.911E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.342E-05   TIME AVG. MOMENT      1.183E-05
 LARGEST RESIDUAL MOMENT           -8.339E-08   AT NODE        159   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.429E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.096E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.471E-06   TIME AVG. FORCE       7.580E-07
 LARGEST RESIDUAL FORCE             1.315E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.08       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.101E-06   AT NODE        337   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.342E-05   TIME AVG. MOMENT      1.183E-05
 LARGEST RESIDUAL MOMENT           -4.193E-11   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.429E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.228E-06   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.176    
 STEP TIME COMPLETED       0.352    ,  TOTAL TIME COMPLETED        0.352    


  INCREMENT    45 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.140E-04   TIME AVG. FORCE       4.088E-06
 LARGEST RESIDUAL FORCE             1.307E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.18       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.185       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.303E-05   TIME AVG. MOMENT      1.275E-05
 LARGEST RESIDUAL MOMENT           -1.025E-03   AT NODE        382   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.107E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.678E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.494E-06   TIME AVG. FORCE       8.679E-07
 LARGEST RESIDUAL FORCE             4.518E-05   AT NODE        228   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.16       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.930E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.527E-05   TIME AVG. MOMENT      1.213E-05
 LARGEST RESIDUAL MOMENT           -4.098E-05   AT NODE        227   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.078E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.464E-03   AT NODE        244   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.207E-04   TIME AVG. FORCE       4.285E-06
 LARGEST RESIDUAL FORCE             1.132E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.01       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.175       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.696E-05   TIME AVG. MOMENT      1.217E-05
 LARGEST RESIDUAL MOMENT           -1.867E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.361E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.717E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.541E-06   TIME AVG. FORCE       7.811E-07
 LARGEST RESIDUAL FORCE            -2.559E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.01       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.950E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.413E-05   TIME AVG. MOMENT      1.210E-05
 LARGEST RESIDUAL MOMENT           -7.055E-07   AT NODE        263   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.358E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.678E-04   AT NODE        272   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.882E-06   TIME AVG. FORCE       7.911E-07
 LARGEST RESIDUAL FORCE             1.716E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.999       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.173E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.412E-05   TIME AVG. MOMENT      1.210E-05
 LARGEST RESIDUAL MOMENT            5.503E-07   AT NODE        258   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.287E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.179E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.520E-06   TIME AVG. FORCE       7.804E-07
 LARGEST RESIDUAL FORCE             9.873E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.999       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.268E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.411E-05   TIME AVG. MOMENT      1.210E-05
 LARGEST RESIDUAL MOMENT           -1.139E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.287E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.700E-06   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.180    
 STEP TIME COMPLETED       0.361    ,  TOTAL TIME COMPLETED        0.361    


  INCREMENT    46 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.325E-05   TIME AVG. FORCE       2.280E-06
 LARGEST RESIDUAL FORCE             6.615E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.05       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.118       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.117E-05   TIME AVG. MOMENT      1.273E-05
 LARGEST RESIDUAL MOMENT           -6.493E-04   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.422E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.135E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.016E-06   TIME AVG. FORCE       8.729E-07
 LARGEST RESIDUAL FORCE             1.781E-05   AT NODE        248   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.058E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.537E-05   TIME AVG. MOMENT      1.239E-05
 LARGEST RESIDUAL MOMENT           -1.452E-05   AT NODE        251   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.376E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.151E-03   AT NODE        260   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.036E-05   TIME AVG. FORCE       2.197E-06
 LARGEST RESIDUAL FORCE             5.129E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.929       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.112       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.547E-05   TIME AVG. MOMENT      1.239E-05
 LARGEST RESIDUAL MOMENT           -5.668E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.228E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.147E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.554E-06   TIME AVG. FORCE       8.026E-07
 LARGEST RESIDUAL FORCE            -8.225E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.927       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.993E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.479E-05   TIME AVG. MOMENT      1.238E-05
 LARGEST RESIDUAL MOMENT            2.811E-07   AT NODE        339   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.225E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.739E-04   AT NODE        288   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.952E-06   TIME AVG. FORCE       8.139E-07
 LARGEST RESIDUAL FORCE             3.403E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.918       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.031E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.475E-05   TIME AVG. MOMENT      1.238E-05
 LARGEST RESIDUAL MOMENT           -5.967E-08   AT NODE        163   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.130E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.481E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.012E-04
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.553E-06   TIME AVG. FORCE       8.025E-07
 LARGEST RESIDUAL FORCE             3.816E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.918       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.171E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.475E-05   TIME AVG. MOMENT      1.238E-05
 LARGEST RESIDUAL MOMENT            2.344E-11   AT NODE        335   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.130E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.291E-06   AT NODE        304   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.185    
 STEP TIME COMPLETED       0.370    ,  TOTAL TIME COMPLETED        0.370    


  INCREMENT    47 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.316E-05   TIME AVG. FORCE       1.423E-06
 LARGEST RESIDUAL FORCE             2.947E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.922       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.416E-02   AT NODE        260   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.389E-05   TIME AVG. MOMENT      1.284E-05
 LARGEST RESIDUAL MOMENT           -3.690E-04   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.846E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.158E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.910E-06   TIME AVG. FORCE       8.888E-07
 LARGEST RESIDUAL FORCE            -3.089E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.898       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.185E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.564E-05   TIME AVG. MOMENT      1.266E-05
 LARGEST RESIDUAL MOMENT           -4.334E-06   AT NODE        279   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.760E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.868E-03   AT NODE        280   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.647E-05   TIME AVG. FORCE       1.238E-06
 LARGEST RESIDUAL FORCE             1.755E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.848       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.154E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.546E-05   TIME AVG. MOMENT      1.266E-05
 LARGEST RESIDUAL MOMENT           -1.102E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.093E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.673E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.612E-06   TIME AVG. FORCE       8.250E-07
 LARGEST RESIDUAL FORCE            -1.421E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.846       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.170E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.539E-05   TIME AVG. MOMENT      1.265E-05
 LARGEST RESIDUAL MOMENT            3.596E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.078E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.855E-04   AT NODE        316   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.036E-06   TIME AVG. FORCE       8.368E-07
 LARGEST RESIDUAL FORCE             5.201E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.837       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.071E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.535E-05   TIME AVG. MOMENT      1.265E-05
 LARGEST RESIDUAL MOMENT           -5.812E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.963E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.152E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.862E-04
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.591E-06   TIME AVG. FORCE       8.244E-07
 LARGEST RESIDUAL FORCE            -2.167E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.837       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.104E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.534E-05   TIME AVG. MOMENT      1.265E-05
 LARGEST RESIDUAL MOMENT            4.502E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.963E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.867E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.189    
 STEP TIME COMPLETED       0.379    ,  TOTAL TIME COMPLETED        0.379    


  INCREMENT    48 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.753E-06   TIME AVG. FORCE       1.066E-06
 LARGEST RESIDUAL FORCE             1.102E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.810       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.811E-02   AT NODE        280   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.981E-05   TIME AVG. MOMENT      1.301E-05
 LARGEST RESIDUAL MOMENT           -1.808E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.368E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.050E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.236E-06   TIME AVG. FORCE       9.166E-07
 LARGEST RESIDUAL FORCE             1.077E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.784       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.447E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.601E-05   TIME AVG. MOMENT      1.293E-05
 LARGEST RESIDUAL MOMENT           -3.044E-06   AT NODE        103   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.201E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.813E-03   AT NODE        328   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.874E-06   TIME AVG. FORCE       9.068E-07
 LARGEST RESIDUAL FORCE             2.948E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.766       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.373E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.590E-05   TIME AVG. MOMENT      1.293E-05
 LARGEST RESIDUAL MOMENT           -3.580E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.927E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.739E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.802E-06   TIME AVG. FORCE       8.508E-07
 LARGEST RESIDUAL FORCE             1.802E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.760       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.150E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.592E-05   TIME AVG. MOMENT      1.293E-05
 LARGEST RESIDUAL MOMENT            2.819E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.858E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.853E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.749E-06   TIME AVG. FORCE       8.494E-07
 LARGEST RESIDUAL FORCE             1.678E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.756       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.585E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.589E-05   TIME AVG. MOMENT      1.293E-05
 LARGEST RESIDUAL MOMENT            7.069E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.793E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.502E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    1.630E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.635E-06   TIME AVG. FORCE       8.463E-07
 LARGEST RESIDUAL FORCE             6.785E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.755       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.660E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.589E-05   TIME AVG. MOMENT      1.293E-05
 LARGEST RESIDUAL MOMENT            1.029E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.789E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.192E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      1.634E-06   TIME AVG. FORCE       8.463E-07
 LARGEST RESIDUAL FORCE             2.042E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.755       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.928E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.589E-05   TIME AVG. MOMENT      1.293E-05
 LARGEST RESIDUAL MOMENT            9.140E-13   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.789E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.278E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   8.355E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.194    
 STEP TIME COMPLETED       0.388    ,  TOTAL TIME COMPLETED        0.388    


  INCREMENT    49 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.355E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.801E-06   TIME AVG. FORCE       9.241E-07
 LARGEST RESIDUAL FORCE             1.388E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.661       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.080E-02   AT NODE        340   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.731E-05   TIME AVG. MOMENT      1.322E-05
 LARGEST RESIDUAL MOMENT           -3.742E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.788E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.494E-03   AT NODE        168   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.091E-06   TIME AVG. FORCE       9.317E-07
 LARGEST RESIDUAL FORCE             2.275E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.638       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.006E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.632E-05   TIME AVG. MOMENT      1.320E-05
 LARGEST RESIDUAL MOMENT           -1.401E-06   AT NODE        119   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.549E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.389E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.671E-06   TIME AVG. FORCE       8.680E-07
 LARGEST RESIDUAL FORCE            -9.641E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.637       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.283E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.637E-05   TIME AVG. MOMENT      1.320E-05
 LARGEST RESIDUAL MOMENT            5.196E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.532E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.630E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -6.046E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.679E-06   TIME AVG. FORCE       8.682E-07
 LARGEST RESIDUAL FORCE             5.517E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.637       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.006E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.635E-05   TIME AVG. MOMENT      1.320E-05
 LARGEST RESIDUAL MOMENT           -1.649E-10   AT NODE        175   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.521E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.173E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.674E-06   TIME AVG. FORCE       8.681E-07
 LARGEST RESIDUAL FORCE            -4.985E-12   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.637       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.656E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.635E-05   TIME AVG. MOMENT      1.320E-05
 LARGEST RESIDUAL MOMENT           -8.794E-13   AT NODE        384   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.521E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.322E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.355E-03,  FRACTION OF STEP COMPLETED  0.198    
 STEP TIME COMPLETED       0.396    ,  TOTAL TIME COMPLETED        0.396    


  INCREMENT    50 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.038E-06   TIME AVG. FORCE       9.494E-07
 LARGEST RESIDUAL FORCE            -3.041E-05   AT NODE        240   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.619       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.046E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.749E-05   TIME AVG. MOMENT      1.349E-05
 LARGEST RESIDUAL MOMENT            3.600E-05   AT NODE        227   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.649E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.621E-03   AT NODE        180   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.384E-06   TIME AVG. FORCE       9.326E-07
 LARGEST RESIDUAL FORCE             1.656E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.601       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.518E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.674E-05   TIME AVG. MOMENT      1.347E-05
 LARGEST RESIDUAL MOMENT           -1.265E-06   AT NODE        115   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.446E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.036E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.708E-06   TIME AVG. FORCE       8.896E-07
 LARGEST RESIDUAL FORCE            -1.267E-08   AT NODE        244   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.601       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.454E-04   AT NODE        204   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.679E-05   TIME AVG. MOMENT      1.348E-05
 LARGEST RESIDUAL MOMENT           -8.451E-09   AT NODE        123   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.445E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.619E-05   AT NODE        280   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.712E-06   TIME AVG. FORCE       8.897E-07
 LARGEST RESIDUAL FORCE             1.548E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.601       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.639E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.679E-05   TIME AVG. MOMENT      1.348E-05
 LARGEST RESIDUAL MOMENT            2.361E-11   AT NODE        250   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.445E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.346E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.202    
 STEP TIME COMPLETED       0.405    ,  TOTAL TIME COMPLETED        0.405    


  INCREMENT    51 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.465E-06   TIME AVG. FORCE       9.791E-07
 LARGEST RESIDUAL FORCE            -5.494E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.535       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.602E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.795E-05   TIME AVG. MOMENT      1.376E-05
 LARGEST RESIDUAL MOMENT            3.168E-05   AT NODE        251   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.390E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.584E-03   AT NODE        200   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.266E-06   TIME AVG. FORCE       9.241E-07
 LARGEST RESIDUAL FORCE             5.112E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.525       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.420E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.717E-05   TIME AVG. MOMENT      1.374E-05
 LARGEST RESIDUAL MOMENT            8.456E-07   AT NODE        234   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.276E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.139E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.751E-06   TIME AVG. FORCE       9.113E-07
 LARGEST RESIDUAL FORCE             3.537E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.524       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.342E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.717E-05   TIME AVG. MOMENT      1.374E-05
 LARGEST RESIDUAL MOMENT            6.757E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.267E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.405E-05   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.750E-06   TIME AVG. FORCE       9.112E-07
 LARGEST RESIDUAL FORCE             1.481E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.524       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.810E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.717E-05   TIME AVG. MOMENT      1.374E-05
 LARGEST RESIDUAL MOMENT           -1.784E-11   AT NODE        155   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.265E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.914E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.751E-06   TIME AVG. FORCE       9.113E-07
 LARGEST RESIDUAL FORCE            -3.555E-12   AT NODE        374   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.524       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.654E-07   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.717E-05   TIME AVG. MOMENT      1.374E-05
 LARGEST RESIDUAL MOMENT           -7.234E-13   AT NODE        357   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.265E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.803E-08   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.207    
 STEP TIME COMPLETED       0.414    ,  TOTAL TIME COMPLETED        0.414    


  INCREMENT    52 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.966E-06   TIME AVG. FORCE       1.010E-06
 LARGEST RESIDUAL FORCE            -1.143E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.455       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.880E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.821E-05   TIME AVG. MOMENT      1.402E-05
 LARGEST RESIDUAL MOMENT            2.673E-05   AT NODE        286   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.150E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.579E-03   AT NODE        228   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.879E-06   TIME AVG. FORCE       9.349E-07
 LARGEST RESIDUAL FORCE             9.711E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.451       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.479E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.753E-05   TIME AVG. MOMENT      1.401E-05
 LARGEST RESIDUAL MOMENT            4.834E-07   AT NODE        230   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.100E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.058E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.792E-06   TIME AVG. FORCE       9.327E-07
 LARGEST RESIDUAL FORCE             8.188E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.450       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.636E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.750E-05   TIME AVG. MOMENT      1.401E-05
 LARGEST RESIDUAL MOMENT           -3.087E-09   AT NODE        131   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.086E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.406E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -8.977E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.785E-06   TIME AVG. FORCE       9.326E-07
 LARGEST RESIDUAL FORCE             1.073E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.449       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.808E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.749E-05   TIME AVG. MOMENT      1.401E-05
 LARGEST RESIDUAL MOMENT            5.861E-12   AT NODE        391   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.085E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.187E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.211    
 STEP TIME COMPLETED       0.423    ,  TOTAL TIME COMPLETED        0.423    


  INCREMENT    53 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.919E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.324E-06   TIME AVG. FORCE       1.037E-06
 LARGEST RESIDUAL FORCE            -1.364E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.381       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.864E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.831E-05   TIME AVG. MOMENT      1.428E-05
 LARGEST RESIDUAL MOMENT            2.143E-05   AT NODE        318   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.388E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.605E-03   AT NODE        281   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.827E-06   TIME AVG. FORCE       9.539E-07
 LARGEST RESIDUAL FORCE             1.671E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.379       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.808E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.779E-05   TIME AVG. MOMENT      1.427E-05
 LARGEST RESIDUAL MOMENT            2.384E-07   AT NODE        226   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.170E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.177E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.820E-06   TIME AVG. FORCE       9.537E-07
 LARGEST RESIDUAL FORCE             3.430E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.378       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.020E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.778E-05   TIME AVG. MOMENT      1.427E-05
 LARGEST RESIDUAL MOMENT           -1.010E-09   AT NODE        135   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.079E-03   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.076E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.820E-06   TIME AVG. FORCE       9.537E-07
 LARGEST RESIDUAL FORCE             1.483E-11   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.378       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.451E-06   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.778E-05   TIME AVG. MOMENT      1.427E-05
 LARGEST RESIDUAL MOMENT            7.211E-13   AT NODE        331   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.079E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.966E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.338E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.919E-03,  FRACTION OF STEP COMPLETED  0.216    
 STEP TIME COMPLETED       0.432    ,  TOTAL TIME COMPLETED        0.432    


  INCREMENT    54 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.338E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.792E-05   TIME AVG. FORCE       1.348E-06
 LARGEST RESIDUAL FORCE            -4.512E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.430       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.138       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.983E-05   TIME AVG. MOMENT      1.456E-05
 LARGEST RESIDUAL MOMENT            5.492E-05   AT NODE        322   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.071E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.299E-03   AT NODE        248   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.917E-06   TIME AVG. FORCE       9.761E-07
 LARGEST RESIDUAL FORCE            -3.547E-06   AT NODE        272   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.426       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.963E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.811E-05   TIME AVG. MOMENT      1.452E-05
 LARGEST RESIDUAL MOMENT            1.930E-06   AT NODE        214   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.044E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.965E-04   AT NODE        272   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.007E-06   TIME AVG. FORCE       9.782E-07
 LARGEST RESIDUAL FORCE             1.508E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.421       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.403E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.808E-05   TIME AVG. MOMENT      1.452E-05
 LARGEST RESIDUAL MOMENT            1.924E-07   AT NODE        230   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.836E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.069E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.857E-06   TIME AVG. FORCE       9.747E-07
 LARGEST RESIDUAL FORCE            -7.566E-09   AT NODE        284   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.421       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.436E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.807E-05   TIME AVG. MOMENT      1.452E-05
 LARGEST RESIDUAL MOMENT            1.987E-09   AT NODE        226   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.827E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.245E-05   AT NODE        281   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.861E-06   TIME AVG. FORCE       9.748E-07
 LARGEST RESIDUAL FORCE             1.233E-11   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.421       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.578E-06   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.807E-05   TIME AVG. MOMENT      1.452E-05
 LARGEST RESIDUAL MOMENT            1.044E-12   AT NODE         31   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.827E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.830E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   7.530E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.338E-02,  FRACTION OF STEP COMPLETED  0.222    
 STEP TIME COMPLETED       0.445    ,  TOTAL TIME COMPLETED        0.445    


  INCREMENT    55 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.530E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.494E-06   TIME AVG. FORCE       1.032E-06
 LARGEST RESIDUAL FORCE            -5.928E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.192       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.422E-02   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.824E-05   TIME AVG. MOMENT      1.477E-05
 LARGEST RESIDUAL MOMENT            5.696E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.396E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.142E-03   AT NODE        365   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.884E-06   TIME AVG. FORCE       9.955E-07
 LARGEST RESIDUAL FORCE            -1.714E-08   AT NODE        276   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.192       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.190E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.820E-05   TIME AVG. MOMENT      1.477E-05
 LARGEST RESIDUAL MOMENT            7.709E-09   AT NODE        222   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.376E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.460E-05   AT NODE        276   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.887E-06   TIME AVG. FORCE       9.955E-07
 LARGEST RESIDUAL FORCE             5.759E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.192       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.843E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.820E-05   TIME AVG. MOMENT      1.477E-05
 LARGEST RESIDUAL MOMENT           -2.410E-12   AT NODE        150   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.374E-03   AT NODE        321   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.708E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.530E-03,  FRACTION OF STEP COMPLETED  0.226    
 STEP TIME COMPLETED       0.453    ,  TOTAL TIME COMPLETED        0.453    


  INCREMENT    56 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.338E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.618E-05   TIME AVG. FORCE       1.333E-06
 LARGEST RESIDUAL FORCE             4.083E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.220       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.126       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.869E-05   TIME AVG. MOMENT      1.502E-05
 LARGEST RESIDUAL MOMENT            1.680E-05   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.858E-03   AT NODE        181   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.247E-03   AT NODE        385   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.907E-06   TIME AVG. FORCE       1.016E-06
 LARGEST RESIDUAL FORCE            -3.261E-07   AT NODE        268   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.220       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.794E-04   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.834E-05   TIME AVG. MOMENT      1.501E-05
 LARGEST RESIDUAL MOMENT            1.593E-07   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.843E-03   AT NODE        181   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.507E-04   AT NODE        268   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.916E-06   TIME AVG. FORCE       1.016E-06
 LARGEST RESIDUAL FORCE             4.903E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.219       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.940E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.834E-05   TIME AVG. MOMENT      1.501E-05
 LARGEST RESIDUAL MOMENT            1.568E-09   AT NODE        214   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.844E-03   AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.481E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.917E-06   TIME AVG. FORCE       1.016E-06
 LARGEST RESIDUAL FORCE             3.720E-12   AT NODE        350   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.219       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.323E-07   AT NODE        345   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.834E-05   TIME AVG. MOMENT      1.501E-05
 LARGEST RESIDUAL MOMENT            9.948E-13   AT NODE        375   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.844E-03   AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.643E-07   AT NODE        273   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.007E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.338E-02,  FRACTION OF STEP COMPLETED  0.233    
 STEP TIME COMPLETED       0.466    ,  TOTAL TIME COMPLETED        0.466    


  INCREMENT    57 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.007E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.762E-05   TIME AVG. FORCE       2.029E-06
 LARGEST RESIDUAL FORCE             1.336E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.118       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.224       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.870E-05   TIME AVG. MOMENT      1.525E-05
 LARGEST RESIDUAL MOMENT            1.827E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.423E-03   AT NODE        105   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.698E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.964E-06   TIME AVG. FORCE       1.037E-06
 LARGEST RESIDUAL FORCE            -1.052E-06   AT NODE        260   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.117       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.699E-03   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.839E-05   TIME AVG. MOMENT      1.525E-05
 LARGEST RESIDUAL MOMENT            5.474E-07   AT NODE        198   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.412E-03   AT NODE        105   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.703E-04   AT NODE        260   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.958E-06   TIME AVG. FORCE       1.036E-06
 LARGEST RESIDUAL FORCE            -9.609E-08   AT NODE        252   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.117       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.231E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.838E-05   TIME AVG. MOMENT      1.525E-05
 LARGEST RESIDUAL MOMENT            5.069E-08   AT NODE        198   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.415E-03   AT NODE        105   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.093E-05   AT NODE        252   DOF  6
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -7.495E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.967E-06   TIME AVG. FORCE       1.037E-06
 LARGEST RESIDUAL FORCE             1.711E-10   AT NODE        269   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.117       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.503E-05   AT NODE        353   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.839E-05   TIME AVG. MOMENT      1.525E-05
 LARGEST RESIDUAL MOMENT            5.657E-11   AT NODE        214   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.415E-03   AT NODE        108   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.402E-06   AT NODE        272   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.100E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.007E-02,  FRACTION OF STEP COMPLETED  0.243    
 STEP TIME COMPLETED       0.486    ,  TOTAL TIME COMPLETED        0.486    


  INCREMENT    58 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.100E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.613E-05   TIME AVG. FORCE       1.570E-06
 LARGEST RESIDUAL FORCE             7.176E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         6.214E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.162       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.827E-05   TIME AVG. MOMENT      1.547E-05
 LARGEST RESIDUAL MOMENT            2.204E-06   AT NODE        114   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      3.007E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.084E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.048E-06   TIME AVG. FORCE       1.058E-06
 LARGEST RESIDUAL FORCE             8.906E-09   AT NODE        233   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         6.182E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.143E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.829E-05   TIME AVG. MOMENT      1.547E-05
 LARGEST RESIDUAL MOMENT           -5.425E-09   AT NODE         83   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      3.005E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.452E-05   AT NODE        233   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.041E-06   TIME AVG. FORCE       1.058E-06
 LARGEST RESIDUAL FORCE             3.667E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         6.174E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.617E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.829E-05   TIME AVG. MOMENT      1.547E-05
 LARGEST RESIDUAL MOMENT           -5.693E-11   AT NODE         91   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      3.003E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.074E-06   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.150E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.100E-02,  FRACTION OF STEP COMPLETED  0.253    
 STEP TIME COMPLETED       0.507    ,  TOTAL TIME COMPLETED        0.507    


  INCREMENT    59 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.150E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.994E-06   TIME AVG. FORCE       1.223E-06
 LARGEST RESIDUAL FORCE             1.500E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         0.170       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.708E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.813E-05   TIME AVG. MOMENT      1.569E-05
 LARGEST RESIDUAL MOMENT            3.601E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      6.435E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.932E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.530E-06   TIME AVG. FORCE       1.089E-06
 LARGEST RESIDUAL FORCE             1.783E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         0.169       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.190E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.810E-05   TIME AVG. MOMENT      1.569E-05
 LARGEST RESIDUAL MOMENT           -7.240E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      6.432E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.386E-06   AT NODE        205   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.726E-02
          TIME INCREMENT REDUCED TO   2.096E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.150E-02,  FRACTION OF STEP COMPLETED  0.269    
 STEP TIME COMPLETED       0.538    ,  TOTAL TIME COMPLETED        0.538    


  INCREMENT    60 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.096E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.448E-06   TIME AVG. FORCE       1.157E-06
 LARGEST RESIDUAL FORCE             6.041E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         6.879E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.397E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.804E-05   TIME AVG. MOMENT      1.589E-05
 LARGEST RESIDUAL MOMENT           -2.999E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      3.161E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.120E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.512E-06   TIME AVG. FORCE       1.118E-06
 LARGEST RESIDUAL FORCE             1.715E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         6.872E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.792E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.802E-05   TIME AVG. MOMENT      1.589E-05
 LARGEST RESIDUAL MOMENT           -3.854E-10   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      3.156E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.383E-06   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.063E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.096E-02,  FRACTION OF STEP COMPLETED  0.280    
 STEP TIME COMPLETED       0.559    ,  TOTAL TIME COMPLETED        0.559    


  INCREMENT    61 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.063E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.080E-04   TIME AVG. FORCE       1.526E-05
 LARGEST RESIDUAL FORCE             1.960E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.730       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.857       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.384E-05   TIME AVG. MOMENT      1.635E-05
 LARGEST RESIDUAL MOMENT            3.054E-04   AT NODE        318   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.610E-02   AT NODE        161   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.178E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.072E-06   TIME AVG. FORCE       1.137E-06
 LARGEST RESIDUAL FORCE            -6.795E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.721       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.124E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.843E-05   TIME AVG. MOMENT      1.610E-05
 LARGEST RESIDUAL MOMENT           -2.056E-07   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.610E-02   AT NODE        161   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.213E-05   AT NODE        172   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.255E-06   TIME AVG. FORCE       1.140E-06
 LARGEST RESIDUAL FORCE             3.544E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.712       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.127E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.844E-05   TIME AVG. MOMENT      1.610E-05
 LARGEST RESIDUAL MOMENT            1.874E-07   AT NODE        198   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.589E-02   AT NODE        160   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.995E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.181E-06   TIME AVG. FORCE       1.139E-06
 LARGEST RESIDUAL FORCE             6.421E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.712       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.680E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.845E-05   TIME AVG. MOMENT      1.610E-05
 LARGEST RESIDUAL MOMENT           -1.866E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.590E-02   AT NODE        160   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.863E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.180E-06   TIME AVG. FORCE       1.139E-06
 LARGEST RESIDUAL FORCE             2.510E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.712       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.807E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.845E-05   TIME AVG. MOMENT      1.610E-05
 LARGEST RESIDUAL MOMENT            1.077E-12   AT NODE        239   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.590E-02   AT NODE        159   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.438E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.891E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.063E-02,  FRACTION OF STEP COMPLETED  0.305    
 STEP TIME COMPLETED       0.610    ,  TOTAL TIME COMPLETED        0.610    


  INCREMENT    62 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.891E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.534E-03   TIME AVG. FORCE       3.119E-05
 LARGEST RESIDUAL FORCE             4.111E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.99       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -1.32       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.764E-04   TIME AVG. MOMENT      1.868E-05
 LARGEST RESIDUAL MOMENT           -3.894E-03   AT NODE        115   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.542E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.198E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.754E-06   TIME AVG. FORCE       1.229E-06
 LARGEST RESIDUAL FORCE            -2.902E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -2.00       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.065E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.018E-05   TIME AVG. MOMENT      1.633E-05
 LARGEST RESIDUAL MOMENT            1.016E-05   AT NODE        171   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.545E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.651E-04   AT NODE         40   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.105E-03   TIME AVG. FORCE       2.278E-05
 LARGEST RESIDUAL FORCE             4.655E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.05       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         1.02       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.032E-04   TIME AVG. MOMENT      1.912E-05
 LARGEST RESIDUAL MOMENT           -2.991E-03   AT NODE        111   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.784E-02   AT NODE        149   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.347E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      8.252E-06   TIME AVG. FORCE       1.278E-06
 LARGEST RESIDUAL FORCE            -3.436E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.068E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.933E-05   TIME AVG. MOMENT      1.631E-05
 LARGEST RESIDUAL MOMENT            2.629E-05   AT NODE        287   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.775E-02   AT NODE        149   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.917E-04   AT NODE         33   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.365E-03   TIME AVG. FORCE       2.788E-05
 LARGEST RESIDUAL FORCE             7.716E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -2.05       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -1.07       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.448E-04   TIME AVG. MOMENT      1.817E-05
 LARGEST RESIDUAL MOMENT            3.762E-03   AT NODE        366   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.667E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.474E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    62 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.223E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.803E-06   TIME AVG. FORCE       1.211E-06
 LARGEST RESIDUAL FORCE             9.596E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.224       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.796E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.879E-05   TIME AVG. MOMENT      1.630E-05
 LARGEST RESIDUAL MOMENT           -5.590E-06   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.959E-03   AT NODE        196   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.410E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.172E-06   TIME AVG. FORCE       1.159E-06
 LARGEST RESIDUAL FORCE             7.469E-08   AT NODE        253   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.224       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.233E-04   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.860E-05   TIME AVG. MOMENT      1.630E-05
 LARGEST RESIDUAL MOMENT            3.441E-08   AT NODE        198   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.947E-03   AT NODE        196   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.218E-05   AT NODE        253   DOF  6
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -4.443E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.182E-06   TIME AVG. FORCE       1.159E-06
 LARGEST RESIDUAL FORCE             9.138E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.224       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.047E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.860E-05   TIME AVG. MOMENT      1.630E-05
 LARGEST RESIDUAL MOMENT           -1.255E-12   AT NODE         75   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.947E-03   AT NODE        195   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.619E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.223E-02,  FRACTION OF STEP COMPLETED  0.311    
 STEP TIME COMPLETED       0.622    ,  TOTAL TIME COMPLETED        0.622    


  INCREMENT    63 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.223E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.271E-06   TIME AVG. FORCE       1.239E-06
 LARGEST RESIDUAL FORCE             1.072E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.286       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.221E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.891E-05   TIME AVG. MOMENT      1.650E-05
 LARGEST RESIDUAL MOMENT           -4.034E-06   AT NODE        239   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.296E-03   AT NODE        252   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.510E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.204E-06   TIME AVG. FORCE       1.180E-06
 LARGEST RESIDUAL FORCE             1.545E-08   AT NODE        257   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.286       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.884E-04   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.879E-05   TIME AVG. MOMENT      1.650E-05
 LARGEST RESIDUAL MOMENT            6.995E-09   AT NODE        202   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.292E-03   AT NODE        249   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.309E-05   AT NODE        257   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.204E-06   TIME AVG. FORCE       1.180E-06
 LARGEST RESIDUAL FORCE             3.262E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.286       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.174E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.879E-05   TIME AVG. MOMENT      1.650E-05
 LARGEST RESIDUAL MOMENT           -8.965E-13   AT NODE         31   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.291E-03   AT NODE        249   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.748E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.834E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.223E-02,  FRACTION OF STEP COMPLETED  0.317    
 STEP TIME COMPLETED       0.635    ,  TOTAL TIME COMPLETED        0.635    


  INCREMENT    64 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.834E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.046E-05   TIME AVG. FORCE       1.543E-06
 LARGEST RESIDUAL FORCE             6.531E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.578       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.149       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.066E-05   TIME AVG. MOMENT      1.672E-05
 LARGEST RESIDUAL MOMENT           -3.143E-05   AT NODE        211   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.307E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.692E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.254E-06   TIME AVG. FORCE       1.200E-06
 LARGEST RESIDUAL FORCE             1.879E-06   AT NODE        257   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.577       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.185E-03   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.920E-05   TIME AVG. MOMENT      1.670E-05
 LARGEST RESIDUAL MOMENT            8.899E-07   AT NODE        202   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.301E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.677E-04   AT NODE        257   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.245E-06   TIME AVG. FORCE       1.200E-06
 LARGEST RESIDUAL FORCE             9.472E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.574       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.838E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.918E-05   TIME AVG. MOMENT      1.670E-05
 LARGEST RESIDUAL MOMENT           -1.945E-09   AT NODE         98   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.286E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.489E-04   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -3.255E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.259E-06   TIME AVG. FORCE       1.200E-06
 LARGEST RESIDUAL FORCE            -5.621E-11   AT NODE        276   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.574       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.588E-06   AT NODE        352   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.918E-05   TIME AVG. MOMENT      1.670E-05
 LARGEST RESIDUAL MOMENT            1.596E-11   AT NODE        214   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.286E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.985E-06   AT NODE        276   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.751E-02
          TIME INCREMENT REDUCED TO   6.114E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.834E-02,  FRACTION OF STEP COMPLETED  0.326    
 STEP TIME COMPLETED       0.653    ,  TOTAL TIME COMPLETED        0.653    


  INCREMENT    65 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.114E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.422E-06   TIME AVG. FORCE       1.223E-06
 LARGEST RESIDUAL FORCE             1.053E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.208       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.706E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.933E-05   TIME AVG. MOMENT      1.689E-05
 LARGEST RESIDUAL MOMENT           -1.714E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.681E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.928E-04   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.264E-06   TIME AVG. FORCE       1.220E-06
 LARGEST RESIDUAL FORCE             6.057E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.208       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.561E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.932E-05   TIME AVG. MOMENT      1.689E-05
 LARGEST RESIDUAL MOMENT            3.722E-10   AT NODE        222   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.668E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.268E-05   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.731E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.114E-03,  FRACTION OF STEP COMPLETED  0.329    
 STEP TIME COMPLETED       0.659    ,  TOTAL TIME COMPLETED        0.659    


  INCREMENT    66 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.731E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.550E-02   TIME AVG. FORCE       2.830E-04
 LARGEST RESIDUAL FORCE             0.402       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.          2.89       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         4.50       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.235E-03   TIME AVG. MOMENT      9.595E-05
 LARGEST RESIDUAL MOMENT           -0.102       AT NODE        223   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      8.020E-02   AT NODE        220   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     0.115       AT NODE        224   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 3 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.974E-02   TIME AVG. FORCE       5.419E-04
 LARGEST RESIDUAL FORCE            -0.854       AT NODE        220   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          2.88       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.883       AT NODE        164   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.332E-02   TIME AVG. MOMENT      3.699E-04
 LARGEST RESIDUAL MOMENT            0.898       AT NODE        215   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      0.206       AT NODE         32   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.208       AT NODE        136   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 32 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       5.45       TIME AVG. FORCE       9.916E-02
 LARGEST RESIDUAL FORCE             -702.       AT NODE        196   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -8.92       AT NODE        293   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -10.7       AT NODE        293   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      8.81       TIME AVG. MOMENT      0.134    
 LARGEST RESIDUAL MOMENT            1.037E+03   AT NODE        262   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION       5.82       AT NODE        261   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION      10.6       AT NODE        196   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 374 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       634.       TIME AVG. FORCE        11.5    
 LARGEST RESIDUAL FORCE             9.214E+04   AT NODE         67   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         3.285E+03   AT NODE        149   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.286E+03   AT NODE        149   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.125E+05   TIME AVG. MOMENT      1.705E+03
 LARGEST RESIDUAL MOMENT            1.561E+07   AT NODE         66   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION       5.01       AT NODE        340   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.372E+03   AT NODE        149   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    66 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.183E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.029E-06   TIME AVG. FORCE       1.289E-06
 LARGEST RESIDUAL FORCE             1.150E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.464       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.097E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.002E-05   TIME AVG. MOMENT      1.709E-05
 LARGEST RESIDUAL MOMENT           -1.053E-05   AT NODE        390   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.048E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.450E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.306E-06   TIME AVG. FORCE       1.239E-06
 LARGEST RESIDUAL FORCE             1.119E-07   AT NODE        265   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.463       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.765E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.965E-05   TIME AVG. MOMENT      1.709E-05
 LARGEST RESIDUAL MOMENT            4.523E-08   AT NODE        206   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.043E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.087E-05   AT NODE        265   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.329E-06   TIME AVG. FORCE       1.240E-06
 LARGEST RESIDUAL FORCE             2.834E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.463       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.427E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.964E-05   TIME AVG. MOMENT      1.708E-05
 LARGEST RESIDUAL MOMENT           -2.004E-11   AT NODE        119   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.040E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.544E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.329E-06   TIME AVG. FORCE       1.240E-06
 LARGEST RESIDUAL FORCE             3.859E-12   AT NODE        347   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.463       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.359E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     2.964E-05   TIME AVG. MOMENT      1.708E-05
 LARGEST RESIDUAL MOMENT           -9.407E-13   AT NODE        357   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.040E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.229E-08   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.774E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.183E-02,  FRACTION OF STEP COMPLETED  0.335    
 STEP TIME COMPLETED       0.671    ,  TOTAL TIME COMPLETED        0.671    


  INCREMENT    67 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.774E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.403E-05   TIME AVG. FORCE       1.647E-06
 LARGEST RESIDUAL FORCE             1.278E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.857       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.163       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.474E-05   TIME AVG. MOMENT      1.735E-05
 LARGEST RESIDUAL MOMENT           -1.130E-04   AT NODE        394   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.049E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.891E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.554E-06   TIME AVG. FORCE       1.263E-06
 LARGEST RESIDUAL FORCE             5.733E-06   AT NODE        265   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.848       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.471E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.044E-05   TIME AVG. MOMENT      1.728E-05
 LARGEST RESIDUAL MOMENT            2.273E-06   AT NODE        222   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.022E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.646E-04   AT NODE        265   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.141E-06   TIME AVG. FORCE       1.292E-06
 LARGEST RESIDUAL FORCE             1.504E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.816       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.173E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.016E-05   TIME AVG. MOMENT      1.728E-05
 LARGEST RESIDUAL MOMENT           -2.407E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.838E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.843E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.412E-06   TIME AVG. FORCE       1.261E-06
 LARGEST RESIDUAL FORCE             2.759E-08   AT NODE        285   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.816       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.621E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.020E-05   TIME AVG. MOMENT      1.728E-05
 LARGEST RESIDUAL MOMENT            5.845E-09   AT NODE        222   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.836E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.500E-05   AT NODE        288   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.419E-06   TIME AVG. FORCE       1.261E-06
 LARGEST RESIDUAL FORCE             1.013E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.816       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.119E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.020E-05   TIME AVG. MOMENT      1.728E-05
 LARGEST RESIDUAL MOMENT            5.214E-12   AT NODE        134   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.834E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.512E-05   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      2.420E-06   TIME AVG. FORCE       1.261E-06
 LARGEST RESIDUAL FORCE             4.033E-12   AT NODE        391   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.816       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.370E-08   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.020E-05   TIME AVG. MOMENT      1.728E-05
 LARGEST RESIDUAL MOMENT           -9.833E-13   AT NODE        351   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.834E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.167E-09   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.774E-02,  FRACTION OF STEP COMPLETED  0.344    
 STEP TIME COMPLETED       0.689    ,  TOTAL TIME COMPLETED        0.689    


  INCREMENT    68 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.774E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.038E-05   TIME AVG. FORCE       2.474E-06
 LARGEST RESIDUAL FORCE             5.164E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.07       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.253       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.391E-05   TIME AVG. MOMENT      1.767E-05
 LARGEST RESIDUAL MOMENT           -3.640E-04   AT NODE        390   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.833E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.984E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.180E-06   TIME AVG. FORCE       1.295E-06
 LARGEST RESIDUAL FORCE            -1.267E-05   AT NODE        288   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.05       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.741E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.113E-05   TIME AVG. MOMENT      1.748E-05
 LARGEST RESIDUAL MOMENT           -7.215E-06   AT NODE        147   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.803E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.026E-03   AT NODE        286   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.232E-05   TIME AVG. FORCE       1.806E-06
 LARGEST RESIDUAL FORCE             2.508E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.921       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.131       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.100E-05   TIME AVG. MOMENT      1.748E-05
 LARGEST RESIDUAL MOMENT           -1.605E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.064E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.408E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.529E-06   TIME AVG. FORCE       1.283E-06
 LARGEST RESIDUAL FORCE            -2.085E-07   AT NODE        296   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.921       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.148E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.084E-05   TIME AVG. MOMENT      1.748E-05
 LARGEST RESIDUAL MOMENT            4.488E-08   AT NODE        271   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.062E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.255E-04   AT NODE        293   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.531E-06   TIME AVG. FORCE       1.283E-06
 LARGEST RESIDUAL FORCE             8.932E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.918       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.464E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.082E-05   TIME AVG. MOMENT      1.748E-05
 LARGEST RESIDUAL MOMENT           -1.992E-10   AT NODE        279   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.048E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.408E-04   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      2.535E-06   TIME AVG. FORCE       1.283E-06
 LARGEST RESIDUAL FORCE             3.177E-12   AT NODE        278   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.918       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.104E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.082E-05   TIME AVG. MOMENT      1.748E-05
 LARGEST RESIDUAL MOMENT            1.106E-12   AT NODE        324   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.048E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.920E-08   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.774E-02,  FRACTION OF STEP COMPLETED  0.353    
 STEP TIME COMPLETED       0.706    ,  TOTAL TIME COMPLETED        0.706    


  INCREMENT    69 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.774E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.534E-04   TIME AVG. FORCE       3.906E-06
 LARGEST RESIDUAL FORCE             1.191E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.26       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.345       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.494E-05   TIME AVG. MOMENT      1.802E-05
 LARGEST RESIDUAL MOMENT           -6.477E-04   AT NODE        374   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.580E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.532E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.292E-06   TIME AVG. FORCE       1.387E-06
 LARGEST RESIDUAL FORCE            -1.241E-04   AT NODE        284   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.24       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.258E-02   AT NODE        328   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.311E-05   TIME AVG. MOMENT      1.771E-05
 LARGEST RESIDUAL MOMENT            5.608E-05   AT NODE        250   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.575E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.220E-03   AT NODE        282   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.089E-04   TIME AVG. FORCE       3.139E-06
 LARGEST RESIDUAL FORCE             8.538E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.996       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.245       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.292E-05   TIME AVG. MOMENT      1.770E-05
 LARGEST RESIDUAL MOMENT           -9.651E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.218E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.357E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.128E-06   TIME AVG. FORCE       1.315E-06
 LARGEST RESIDUAL FORCE             2.892E-05   AT NODE        273   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.997       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.792E-03   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.151E-05   TIME AVG. MOMENT      1.768E-05
 LARGEST RESIDUAL MOMENT            9.115E-06   AT NODE        218   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.206E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.499E-03   AT NODE        273   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.716E-06   TIME AVG. FORCE       1.308E-06
 LARGEST RESIDUAL FORCE             1.172E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.993       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.016E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.149E-05   TIME AVG. MOMENT      1.768E-05
 LARGEST RESIDUAL MOMENT            1.156E-07   AT NODE        202   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.190E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.713E-04   AT NODE        292   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      2.720E-06   TIME AVG. FORCE       1.308E-06
 LARGEST RESIDUAL FORCE             1.747E-07   AT NODE        273   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.993       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.401E-04   AT NODE        356   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.149E-05   TIME AVG. MOMENT      1.768E-05
 LARGEST RESIDUAL MOMENT            5.087E-08   AT NODE        218   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.189E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.164E-04   AT NODE        276   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      2.731E-06   TIME AVG. FORCE       1.308E-06
 LARGEST RESIDUAL FORCE             5.700E-12   AT NODE        269   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.993       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.706E-06   AT NODE        337   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.149E-05   TIME AVG. MOMENT      1.768E-05
 LARGEST RESIDUAL MOMENT            3.309E-12   AT NODE         27   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.189E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.645E-07   AT NODE        272   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.774E-02,  FRACTION OF STEP COMPLETED  0.362    
 STEP TIME COMPLETED       0.724    ,  TOTAL TIME COMPLETED        0.724    


  INCREMENT    70 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.774E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.575E-04   TIME AVG. FORCE       5.650E-06
 LARGEST RESIDUAL FORCE             2.129E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.43       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.433       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.680E-05   TIME AVG. MOMENT      1.838E-05
 LARGEST RESIDUAL MOMENT           -9.282E-04   AT NODE        358   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.245E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.056E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.684E-04   TIME AVG. FORCE       5.836E-06
 LARGEST RESIDUAL FORCE            -9.728E-03   AT NODE        272   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.39       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.124       AT NODE          1   DOF  1
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.691E-04   TIME AVG. MOMENT      1.985E-05
 LARGEST RESIDUAL MOMENT            4.258E-03   AT NODE        230   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.458E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.794E-02   AT NODE        272   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.812E-03   TIME AVG. FORCE       8.284E-05
 LARGEST RESIDUAL FORCE             0.273       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -3.51       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -2.12       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.076E-04   TIME AVG. MOMENT      2.611E-05
 LARGEST RESIDUAL MOMENT            1.939E-02   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.124       AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.233E-02   AT NODE        316   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.544E-04   TIME AVG. FORCE       3.902E-06
 LARGEST RESIDUAL FORCE             6.156E-03   AT NODE        325   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -3.54       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.108       AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.135E-04   TIME AVG. MOMENT      1.905E-05
 LARGEST RESIDUAL MOMENT            2.481E-03   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.128       AT NODE        397   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.331E-02   AT NODE        325   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 2 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.807E-04   TIME AVG. FORCE       1.452E-05
 LARGEST RESIDUAL FORCE             7.130E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -3.70       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.326       AT NODE        361   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.714E-04   TIME AVG. MOMENT      1.988E-05
 LARGEST RESIDUAL MOMENT            7.574E-03   AT NODE        314   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.170       AT NODE        317   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.036E-02   AT NODE        313   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.074E-05   TIME AVG. FORCE       2.316E-06
 LARGEST RESIDUAL FORCE             2.000E-03   AT NODE        185   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -3.67       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.593E-02   AT NODE        289   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.685E-05   TIME AVG. MOMENT      1.839E-05
 LARGEST RESIDUAL MOMENT            1.678E-03   AT NODE        170   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.170       AT NODE        317   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.235E-02   AT NODE        185   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      1.143E-02   TIME AVG. FORCE       1.951E-04
 LARGEST RESIDUAL FORCE             0.226       AT NODE        320   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.818       AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         2.96       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.807E-03   TIME AVG. MOMENT      7.182E-05
 LARGEST RESIDUAL MOMENT           -0.219       AT NODE        319   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.788E-02   AT NODE        316   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     0.163       AT NODE        317   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    70 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  4.435E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.823E-06   TIME AVG. FORCE       1.334E-06
 LARGEST RESIDUAL FORCE             3.608E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.254       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.789E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.161E-05   TIME AVG. MOMENT      1.788E-05
 LARGEST RESIDUAL MOMENT           -1.106E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.598E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.252E-04   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.831E-06   TIME AVG. FORCE       1.334E-06
 LARGEST RESIDUAL FORCE             1.460E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.254       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.093E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.166E-05   TIME AVG. MOMENT      1.788E-05
 LARGEST RESIDUAL MOMENT           -1.077E-10   AT NODE        350   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.541E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.652E-05   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    5.504E-09
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.834E-06   TIME AVG. FORCE       1.334E-06
 LARGEST RESIDUAL FORCE             4.769E-12   AT NODE        383   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.254       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.988E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.166E-05   TIME AVG. MOMENT      1.788E-05
 LARGEST RESIDUAL MOMENT           -7.915E-13   AT NODE        328   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.541E-03   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.486E-08   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.435E-03,  FRACTION OF STEP COMPLETED  0.364    
 STEP TIME COMPLETED       0.728    ,  TOTAL TIME COMPLETED        0.728    


  INCREMENT    71 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.435E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.906E-06   TIME AVG. FORCE       1.360E-06
 LARGEST RESIDUAL FORCE             6.602E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.257       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.795E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.179E-05   TIME AVG. MOMENT      1.808E-05
 LARGEST RESIDUAL MOMENT           -2.064E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.598E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.402E-05   AT NODE        152   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.905E-06   TIME AVG. FORCE       1.360E-06
 LARGEST RESIDUAL FORCE             1.070E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.257       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.573E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.183E-05   TIME AVG. MOMENT      1.808E-05
 LARGEST RESIDUAL MOMENT           -1.920E-11   AT NODE         66   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.593E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.876E-06   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    6.652E-03

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.435E-03,  FRACTION OF STEP COMPLETED  0.366    
 STEP TIME COMPLETED       0.733    ,  TOTAL TIME COMPLETED        0.733    


  INCREMENT    72 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.652E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.036E-06   TIME AVG. FORCE       1.388E-06
 LARGEST RESIDUAL FORCE             2.581E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.393       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.914E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.204E-05   TIME AVG. MOMENT      1.827E-05
 LARGEST RESIDUAL MOMENT           -8.220E-07   AT NODE        239   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.506E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.334E-04   AT NODE        216   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.037E-06   TIME AVG. FORCE       1.388E-06
 LARGEST RESIDUAL FORCE             6.036E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.392       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.801E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.209E-05   TIME AVG. MOMENT      1.827E-05
 LARGEST RESIDUAL MOMENT            1.075E-10   AT NODE        198   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.469E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.648E-05   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.039E-06   TIME AVG. FORCE       1.388E-06
 LARGEST RESIDUAL FORCE             4.435E-12   AT NODE        391   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.392       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.944E-07   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.209E-05   TIME AVG. MOMENT      1.827E-05
 LARGEST RESIDUAL MOMENT           -9.187E-13   AT NODE        290   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.469E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.941E-08   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    9.979E-03

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.652E-03,  FRACTION OF STEP COMPLETED  0.370    
 STEP TIME COMPLETED       0.740    ,  TOTAL TIME COMPLETED        0.740    


  INCREMENT    73 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.979E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.547E-06   TIME AVG. FORCE       1.423E-06
 LARGEST RESIDUAL FORCE             3.128E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.604       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.544E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.260E-05   TIME AVG. MOMENT      1.847E-05
 LARGEST RESIDUAL MOMENT           -9.055E-06   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.327E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.687E-04   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.347E-06   TIME AVG. FORCE       1.419E-06
 LARGEST RESIDUAL FORCE             7.837E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.597       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.462E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.248E-05   TIME AVG. MOMENT      1.847E-05
 LARGEST RESIDUAL MOMENT           -9.678E-09   AT NODE        298   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.286E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.116E-04   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    4.400E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.308E-06   TIME AVG. FORCE       1.419E-06
 LARGEST RESIDUAL FORCE             2.061E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.597       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.384E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.247E-05   TIME AVG. MOMENT      1.847E-05
 LARGEST RESIDUAL MOMENT            2.106E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.284E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.133E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.307E-06   TIME AVG. FORCE       1.419E-06
 LARGEST RESIDUAL FORCE             1.520E-11   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.597       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.266E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.247E-05   TIME AVG. MOMENT      1.847E-05
 LARGEST RESIDUAL MOMENT            1.010E-12   AT NODE        357   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.284E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.815E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.497E-02
          TIME INCREMENT REDUCED TO   2.562E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.979E-03,  FRACTION OF STEP COMPLETED  0.375    
 STEP TIME COMPLETED       0.750    ,  TOTAL TIME COMPLETED        0.750    


  INCREMENT    74 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.562E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.257E-06   TIME AVG. FORCE       1.448E-06
 LARGEST RESIDUAL FORCE             1.192E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.154       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.925E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.257E-05   TIME AVG. MOMENT      1.866E-05
 LARGEST RESIDUAL MOMENT            2.984E-08   AT NODE        175   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.306E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.116E-05   AT NODE         35   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.253E-06   TIME AVG. FORCE       1.448E-06
 LARGEST RESIDUAL FORCE             2.643E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.154       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.119E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.257E-05   TIME AVG. MOMENT      1.866E-05
 LARGEST RESIDUAL MOMENT           -1.157E-11   AT NODE         86   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.304E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.408E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.245E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.562E-03,  FRACTION OF STEP COMPLETED  0.376    
 STEP TIME COMPLETED       0.752    ,  TOTAL TIME COMPLETED        0.752    


  INCREMENT    75 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.245E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.516E-04   TIME AVG. FORCE       5.356E-06
 LARGEST RESIDUAL FORCE             2.014E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.76       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.411       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.511E-05   TIME AVG. MOMENT      1.928E-05
 LARGEST RESIDUAL MOMENT           -8.680E-04   AT NODE        358   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.883E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.988E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.388E-03   TIME AVG. FORCE       5.437E-05
 LARGEST RESIDUAL FORCE             0.151       AT NODE        273   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.77       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.516       AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.813E-03   TIME AVG. MOMENT      4.258E-05
 LARGEST RESIDUAL MOMENT           -4.488E-02   AT NODE        290   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      0.110       AT NODE        273   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     0.112       AT NODE        276   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 2 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.197E-03   TIME AVG. FORCE       2.013E-05
 LARGEST RESIDUAL FORCE            -1.904E-02   AT NODE        332   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.96       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.384       AT NODE        300   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.878E-04   TIME AVG. MOMENT      2.091E-05
 LARGEST RESIDUAL MOMENT            1.409E-02   AT NODE        327   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      0.114       AT NODE        273   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.687E-02   AT NODE        328   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 3 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.441E-04   TIME AVG. FORCE       3.677E-06
 LARGEST RESIDUAL FORCE            -4.498E-03   AT NODE        208   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.96       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.371E-02   AT NODE        324   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.001E-04   TIME AVG. MOMENT      1.975E-05
 LARGEST RESIDUAL MOMENT           -3.831E-03   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      0.124       AT NODE        269   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.652E-02   AT NODE        208   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.946E-03   TIME AVG. FORCE       1.568E-04
 LARGEST RESIDUAL FORCE             0.730       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.736       AT NODE        297   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         2.31       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.920E-03   TIME AVG. MOMENT      4.401E-05
 LARGEST RESIDUAL MOMENT           -9.131E-02   AT NODE        326   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      0.151       AT NODE        297   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     0.167       AT NODE        333   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    75 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  5.613E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.483E-06   TIME AVG. FORCE       1.480E-06
 LARGEST RESIDUAL FORCE             3.383E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.339       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.078E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.278E-05   TIME AVG. MOMENT      1.885E-05
 LARGEST RESIDUAL MOMENT           -8.686E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.257E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.099E-05   AT NODE         32   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.485E-06   TIME AVG. FORCE       1.480E-06
 LARGEST RESIDUAL FORCE             4.924E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.339       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.318E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.279E-05   TIME AVG. MOMENT      1.885E-05
 LARGEST RESIDUAL MOMENT           -1.760E-11   AT NODE         78   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.246E-03   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.030E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    8.419E-03

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.613E-03,  FRACTION OF STEP COMPLETED  0.379    
 STEP TIME COMPLETED       0.758    ,  TOTAL TIME COMPLETED        0.758    


  INCREMENT    76 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.419E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.570E-06   TIME AVG. FORCE       1.512E-06
 LARGEST RESIDUAL FORCE             2.998E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.508       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.037E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.303E-05   TIME AVG. MOMENT      1.903E-05
 LARGEST RESIDUAL MOMENT            8.606E-07   AT NODE        175   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.094E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.233E-04   AT NODE        104   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.662E-06   TIME AVG. FORCE       1.513E-06
 LARGEST RESIDUAL FORCE             1.332E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.506       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.446E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.312E-05   TIME AVG. MOMENT      1.903E-05
 LARGEST RESIDUAL MOMENT           -2.278E-09   AT NODE        110   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.077E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.692E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -4.478E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.670E-06   TIME AVG. FORCE       1.513E-06
 LARGEST RESIDUAL FORCE             5.286E-12   AT NODE        315   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.506       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.650E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.312E-05   TIME AVG. MOMENT      1.903E-05
 LARGEST RESIDUAL MOMENT            9.536E-13   AT NODE        320   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.077E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.163E-07   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.263E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.419E-03,  FRACTION OF STEP COMPLETED  0.383    
 STEP TIME COMPLETED       0.766    ,  TOTAL TIME COMPLETED        0.766    


  INCREMENT    77 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.263E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.614E-06   TIME AVG. FORCE       1.560E-06
 LARGEST RESIDUAL FORCE             1.272E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.769       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.868E-02   AT NODE        296   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.406E-05   TIME AVG. MOMENT      1.923E-05
 LARGEST RESIDUAL MOMENT           -3.205E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.741E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.253E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.828E-06   TIME AVG. FORCE       1.564E-06
 LARGEST RESIDUAL FORCE             9.322E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.744       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.574E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.363E-05   TIME AVG. MOMENT      1.922E-05
 LARGEST RESIDUAL MOMENT           -4.269E-07   AT NODE         99   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.600E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.407E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.599E-06   TIME AVG. FORCE       1.545E-06
 LARGEST RESIDUAL FORCE             3.765E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.740       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.380E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.363E-05   TIME AVG. MOMENT      1.922E-05
 LARGEST RESIDUAL MOMENT            3.201E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.572E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.857E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.142E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.607E-06   TIME AVG. FORCE       1.545E-06
 LARGEST RESIDUAL FORCE             5.726E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.739       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.358E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.360E-05   TIME AVG. MOMENT      1.922E-05
 LARGEST RESIDUAL MOMENT            4.676E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.561E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.101E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.614E-06   TIME AVG. FORCE       1.545E-06
 LARGEST RESIDUAL FORCE             1.121E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.739       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.816E-05   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.359E-05   TIME AVG. MOMENT      1.922E-05
 LARGEST RESIDUAL MOMENT            1.632E-12   AT NODE        382   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.560E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.543E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.263E-02,  FRACTION OF STEP COMPLETED  0.389    
 STEP TIME COMPLETED       0.779    ,  TOTAL TIME COMPLETED        0.779    


  INCREMENT    78 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.263E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.229E-06   TIME AVG. FORCE       1.600E-06
 LARGEST RESIDUAL FORCE             1.913E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.748       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.498E-02   AT NODE        300   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.479E-05   TIME AVG. MOMENT      1.942E-05
 LARGEST RESIDUAL MOMENT           -4.526E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.719E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.584E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.208E-06   TIME AVG. FORCE       1.600E-06
 LARGEST RESIDUAL FORCE             1.369E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.715       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.257E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.408E-05   TIME AVG. MOMENT      1.941E-05
 LARGEST RESIDUAL MOMENT           -9.036E-07   AT NODE         99   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.549E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.701E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.456E-06   TIME AVG. FORCE       1.574E-06
 LARGEST RESIDUAL FORCE             1.076E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.710       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.647E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.415E-05   TIME AVG. MOMENT      1.942E-05
 LARGEST RESIDUAL MOMENT            7.648E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.501E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.815E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -4.075E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.431E-06   TIME AVG. FORCE       1.573E-06
 LARGEST RESIDUAL FORCE             3.388E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.706       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.332E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.405E-05   TIME AVG. MOMENT      1.941E-05
 LARGEST RESIDUAL MOMENT            2.817E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.474E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.670E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    9.836E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.424E-06   TIME AVG. FORCE       1.573E-06
 LARGEST RESIDUAL FORCE             5.071E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.706       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.228E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.405E-05   TIME AVG. MOMENT      1.941E-05
 LARGEST RESIDUAL MOMENT            4.372E-11   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.473E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.039E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   5.104E-03 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.263E-02,  FRACTION OF STEP COMPLETED  0.396    
 STEP TIME COMPLETED       0.791    ,  TOTAL TIME COMPLETED        0.791    


  INCREMENT    79 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.104E-03

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.348E-06   TIME AVG. FORCE       1.599E-06
 LARGEST RESIDUAL FORCE            -3.846E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.280       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.111E-03   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.424E-05   TIME AVG. MOMENT      1.960E-05
 LARGEST RESIDUAL MOMENT            1.131E-06   AT NODE        290   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.821E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.618E-04   AT NODE        260   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.360E-06   TIME AVG. FORCE       1.599E-06
 LARGEST RESIDUAL FORCE             6.169E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.279       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.316E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.422E-05   TIME AVG. MOMENT      1.960E-05
 LARGEST RESIDUAL MOMENT            2.879E-10   AT NODE        250   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.784E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.658E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.358E-06   TIME AVG. FORCE       1.599E-06
 LARGEST RESIDUAL FORCE             4.651E-12   AT NODE        338   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.279       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.080E-07   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.422E-05   TIME AVG. MOMENT      1.960E-05
 LARGEST RESIDUAL MOMENT           -8.978E-13   AT NODE        333   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.784E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.300E-08   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.104E-03,  FRACTION OF STEP COMPLETED  0.398    
 STEP TIME COMPLETED       0.796    ,  TOTAL TIME COMPLETED        0.796    


  INCREMENT    80 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.263E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.425E-06   TIME AVG. FORCE       1.640E-06
 LARGEST RESIDUAL FORCE            -4.460E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.651       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.002E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.486E-05   TIME AVG. MOMENT      1.979E-05
 LARGEST RESIDUAL MOMENT            1.404E-05   AT NODE        263   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.373E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.440E-04   AT NODE        208   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.345E-06   TIME AVG. FORCE       1.625E-06
 LARGEST RESIDUAL FORCE             1.515E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.642       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.693E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.462E-05   TIME AVG. MOMENT      1.979E-05
 LARGEST RESIDUAL MOMENT            3.078E-07   AT NODE        222   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.316E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.752E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.224E-06   TIME AVG. FORCE       1.623E-06
 LARGEST RESIDUAL FORCE             1.704E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.641       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.661E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.463E-05   TIME AVG. MOMENT      1.979E-05
 LARGEST RESIDUAL MOMENT            1.094E-09   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.310E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.974E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.223E-06   TIME AVG. FORCE       1.623E-06
 LARGEST RESIDUAL FORCE             7.467E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.641       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.041E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.463E-05   TIME AVG. MOMENT      1.979E-05
 LARGEST RESIDUAL MOMENT            1.424E-12   AT NODE        346   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.309E-02   AT NODE        403   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.977E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    1.894E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.263E-02,  FRACTION OF STEP COMPLETED  0.405    
 STEP TIME COMPLETED       0.809    ,  TOTAL TIME COMPLETED        0.809    


  INCREMENT    81 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.894E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.563E-05   TIME AVG. FORCE       1.823E-06
 LARGEST RESIDUAL FORCE             8.477E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.958       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.916E-02   AT NODE        284   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.826E-05   TIME AVG. MOMENT      2.002E-05
 LARGEST RESIDUAL MOMENT           -1.372E-04   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.332E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.686E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.275E-06   TIME AVG. FORCE       1.689E-06
 LARGEST RESIDUAL FORCE             1.042E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.904       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.373E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.530E-05   TIME AVG. MOMENT      1.998E-05
 LARGEST RESIDUAL MOMENT            5.455E-06   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.185E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.951E-03   AT NODE        300   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.410E-06   TIME AVG. FORCE       1.691E-06
 LARGEST RESIDUAL FORCE             4.151E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.863       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.160E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.511E-05   TIME AVG. MOMENT      1.998E-05
 LARGEST RESIDUAL MOMENT           -1.977E-06   AT NODE         99   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.889E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.962E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.509E-06   TIME AVG. FORCE       1.650E-06
 LARGEST RESIDUAL FORCE             2.676E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.847       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.581E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.520E-05   TIME AVG. MOMENT      1.998E-05
 LARGEST RESIDUAL MOMENT            3.598E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.814E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.533E-04   AT NODE        368   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.944E-06   TIME AVG. FORCE       1.656E-06
 LARGEST RESIDUAL FORCE             9.528E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.826       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.066E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.511E-05   TIME AVG. MOMENT      1.998E-05
 LARGEST RESIDUAL MOMENT           -6.148E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.673E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.406E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          ROTATION CORRECTION TOO LARGE COMPARED TO ROTATION INCREMENT

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.190E-06   TIME AVG. FORCE       1.645E-06
 LARGEST RESIDUAL FORCE             1.263E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.824       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.601E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.514E-05   TIME AVG. MOMENT      1.998E-05
 LARGEST RESIDUAL MOMENT            1.867E-08   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.657E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.629E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      3.192E-06   TIME AVG. FORCE       1.645E-06
 LARGEST RESIDUAL FORCE             2.903E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.823       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.137E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.513E-05   TIME AVG. MOMENT      1.998E-05
 LARGEST RESIDUAL MOMENT            1.020E-10   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.649E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.742E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      3.192E-06   TIME AVG. FORCE       1.645E-06
 LARGEST RESIDUAL FORCE             1.353E-11   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.823       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.604E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.513E-05   TIME AVG. MOMENT      1.998E-05
 LARGEST RESIDUAL MOMENT            1.049E-12   AT NODE        367   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.649E-02   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.148E-07   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   1.150E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.894E-02,  FRACTION OF STEP COMPLETED  0.414    
 STEP TIME COMPLETED       0.828    ,  TOTAL TIME COMPLETED        0.828    


  INCREMENT    82 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.150E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.672E-06   TIME AVG. FORCE       1.688E-06
 LARGEST RESIDUAL FORCE            -5.373E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.453       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.603E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.566E-05   TIME AVG. MOMENT      2.017E-05
 LARGEST RESIDUAL MOMENT            1.506E-05   AT NODE        263   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.291E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.072E-03   AT NODE        216   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.267E-06   TIME AVG. FORCE       1.668E-06
 LARGEST RESIDUAL FORCE             9.704E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.445       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.678E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.540E-05   TIME AVG. MOMENT      2.017E-05
 LARGEST RESIDUAL MOMENT            3.557E-07   AT NODE        218   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.828E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.626E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.212E-06   TIME AVG. FORCE       1.667E-06
 LARGEST RESIDUAL FORCE             2.031E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.444       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.908E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.539E-05   TIME AVG. MOMENT      2.017E-05
 LARGEST RESIDUAL MOMENT           -1.103E-09   AT NODE        119   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.763E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.493E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.207E-06   TIME AVG. FORCE       1.667E-06
 LARGEST RESIDUAL FORCE             2.990E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.444       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.047E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.539E-05   TIME AVG. MOMENT      2.017E-05
 LARGEST RESIDUAL MOMENT            1.058E-12   AT NODE        383   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.760E-03   AT NODE        402   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.512E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.150E-02,  FRACTION OF STEP COMPLETED  0.420    
 STEP TIME COMPLETED       0.840    ,  TOTAL TIME COMPLETED        0.840    


  INCREMENT    83 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.894E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.462E-05   TIME AVG. FORCE       1.847E-06
 LARGEST RESIDUAL FORCE            -2.421E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.603       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.129       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.706E-05   TIME AVG. MOMENT      2.037E-05
 LARGEST RESIDUAL MOMENT            4.704E-05   AT NODE        282   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.241E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.791E-03   AT NODE        204   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.386E-06   TIME AVG. FORCE       1.691E-06
 LARGEST RESIDUAL FORCE             1.052E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.593       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.012E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.585E-05   TIME AVG. MOMENT      2.036E-05
 LARGEST RESIDUAL MOMENT            2.889E-06   AT NODE        202   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.190E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.167E-04   AT NODE        260   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.539E-06   TIME AVG. FORCE       1.693E-06
 LARGEST RESIDUAL FORCE             2.995E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.578       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.486E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.573E-05   TIME AVG. MOMENT      2.035E-05
 LARGEST RESIDUAL MOMENT            5.335E-07   AT NODE        214   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.110E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.956E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.222E-06   TIME AVG. FORCE       1.689E-06
 LARGEST RESIDUAL FORCE             8.852E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.577       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.895E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.573E-05   TIME AVG. MOMENT      2.035E-05
 LARGEST RESIDUAL MOMENT            1.073E-08   AT NODE        218   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.105E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.275E-05   AT NODE        273   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.237E-06   TIME AVG. FORCE       1.689E-06
 LARGEST RESIDUAL FORCE             2.617E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.577       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.136E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.573E-05   TIME AVG. MOMENT      2.035E-05
 LARGEST RESIDUAL MOMENT            1.652E-11   AT NODE        234   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.105E-02   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.456E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.288E-02

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.894E-02,  FRACTION OF STEP COMPLETED  0.429    
 STEP TIME COMPLETED       0.858    ,  TOTAL TIME COMPLETED        0.858    


  INCREMENT    84 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.288E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.459E-05   TIME AVG. FORCE       2.003E-06
 LARGEST RESIDUAL FORCE            -4.002E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.520       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.177       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.798E-05   TIME AVG. MOMENT      2.056E-05
 LARGEST RESIDUAL MOMENT            6.131E-05   AT NODE        298   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.038E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.705E-03   AT NODE        212   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.417E-06   TIME AVG. FORCE       1.713E-06
 LARGEST RESIDUAL FORCE            -4.523E-06   AT NODE        252   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.513       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.586E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.609E-05   TIME AVG. MOMENT      2.054E-05
 LARGEST RESIDUAL MOMENT            3.568E-06   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.006E-02   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.469E-04   AT NODE        252   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.083E-06   TIME AVG. FORCE       1.722E-06
 LARGEST RESIDUAL FORCE             6.089E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.491       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.259E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.608E-05   TIME AVG. MOMENT      2.054E-05
 LARGEST RESIDUAL MOMENT            2.562E-06   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.144E-03   AT NODE        193   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.147E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.324E-06   TIME AVG. FORCE       1.712E-06
 LARGEST RESIDUAL FORCE            -1.407E-07   AT NODE        272   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.490       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.843E-04   AT NODE        353   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.600E-05   TIME AVG. MOMENT      2.054E-05
 LARGEST RESIDUAL MOMENT            5.196E-08   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.135E-03   AT NODE        193   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.124E-04   AT NODE        269   DOF  6
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.280E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.337E-06   TIME AVG. FORCE       1.712E-06
 LARGEST RESIDUAL FORCE            -4.436E-09   AT NODE        276   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.490       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.719E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.601E-05   TIME AVG. MOMENT      2.054E-05
 LARGEST RESIDUAL MOMENT            1.599E-09   AT NODE        218   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.134E-03   AT NODE        193   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.998E-05   AT NODE        273   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.288E-02,  FRACTION OF STEP COMPLETED  0.441    
 STEP TIME COMPLETED       0.881    ,  TOTAL TIME COMPLETED        0.881    


  INCREMENT    85 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.288E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.085E-05   TIME AVG. FORCE       1.970E-06
 LARGEST RESIDUAL FORCE             4.986E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.330       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.160       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.678E-05   TIME AVG. MOMENT      2.073E-05
 LARGEST RESIDUAL MOMENT            2.965E-05   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.451E-03   AT NODE        141   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.364E-03   AT NODE        349   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.509E-06   TIME AVG. FORCE       1.736E-06
 LARGEST RESIDUAL FORCE            -2.569E-06   AT NODE        256   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.328       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.105E-03   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.620E-05   TIME AVG. MOMENT      2.072E-05
 LARGEST RESIDUAL MOMENT            1.638E-06   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.410E-03   AT NODE        137   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.854E-04   AT NODE        256   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.488E-06   TIME AVG. FORCE       1.736E-06
 LARGEST RESIDUAL FORCE             1.863E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.325       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.922E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.616E-05   TIME AVG. MOMENT      2.072E-05
 LARGEST RESIDUAL MOMENT            1.700E-07   AT NODE        202   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.415E-03   AT NODE        137   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.039E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.500E-06   TIME AVG. FORCE       1.736E-06
 LARGEST RESIDUAL FORCE             3.597E-09   AT NODE        265   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.325       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.816E-05   AT NODE        353   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.618E-05   TIME AVG. MOMENT      2.072E-05
 LARGEST RESIDUAL MOMENT            1.340E-09   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.414E-03   AT NODE        139   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.802E-05   AT NODE        268   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   1.776E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.288E-02,  FRACTION OF STEP COMPLETED  0.452    
 STEP TIME COMPLETED       0.904    ,  TOTAL TIME COMPLETED        0.904    


  INCREMENT    86 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.776E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.395E-06   TIME AVG. FORCE       1.811E-06
 LARGEST RESIDUAL FORCE             1.184E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.178       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.435E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.621E-05   TIME AVG. MOMENT      2.090E-05
 LARGEST RESIDUAL MOMENT            3.608E-06   AT NODE        250   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.938E-03   AT NODE        117   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.585E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.719E-06   TIME AVG. FORCE       1.762E-06
 LARGEST RESIDUAL FORCE            -4.335E-08   AT NODE        252   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.178       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.920E-04   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.626E-05   TIME AVG. MOMENT      2.091E-05
 LARGEST RESIDUAL MOMENT            2.798E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.933E-03   AT NODE        113   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.302E-05   AT NODE        252   DOF  6
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -4.887E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.728E-06   TIME AVG. FORCE       1.762E-06
 LARGEST RESIDUAL FORCE             2.909E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.178       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.027E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.626E-05   TIME AVG. MOMENT      2.091E-05
 LARGEST RESIDUAL MOMENT            3.394E-11   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.933E-03   AT NODE        116   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.531E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.953E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.776E-02,  FRACTION OF STEP COMPLETED  0.461    
 STEP TIME COMPLETED       0.922    ,  TOTAL TIME COMPLETED        0.922    


  INCREMENT    87 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.953E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.386E-06   TIME AVG. FORCE       1.850E-06
 LARGEST RESIDUAL FORCE             1.129E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.322       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.280E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.635E-05   TIME AVG. MOMENT      2.108E-05
 LARGEST RESIDUAL MOMENT            5.327E-06   AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.815E-03   AT NODE        105   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.560E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.321E-06   TIME AVG. FORCE       1.796E-06
 LARGEST RESIDUAL FORCE            -2.998E-07   AT NODE        244   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.322       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.084E-03   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.639E-05   TIME AVG. MOMENT      2.108E-05
 LARGEST RESIDUAL MOMENT            2.208E-07   AT NODE        178   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.791E-03   AT NODE        105   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.656E-04   AT NODE        244   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.343E-06   TIME AVG. FORCE       1.796E-06
 LARGEST RESIDUAL FORCE             3.813E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.321       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.117E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.639E-05   TIME AVG. MOMENT      2.108E-05
 LARGEST RESIDUAL MOMENT            1.003E-09   AT NODE        178   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.788E-03   AT NODE        108   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.843E-06   AT NODE        244   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.929E-02
          TIME INCREMENT REDUCED TO   3.847E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.953E-02,  FRACTION OF STEP COMPLETED  0.481    
 STEP TIME COMPLETED       0.962    ,  TOTAL TIME COMPLETED        0.962    


  INCREMENT    88 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.847E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.706E-05   TIME AVG. FORCE       1.995E-06
 LARGEST RESIDUAL FORCE             3.793E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.445       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.132       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.703E-05   TIME AVG. MOMENT      2.126E-05
 LARGEST RESIDUAL MOMENT           -1.843E-05   AT NODE        131   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.314E-03   AT NODE        124   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.758E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.329E-06   TIME AVG. FORCE       1.829E-06
 LARGEST RESIDUAL FORCE             9.153E-07   AT NODE        233   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.444       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.937E-03   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.659E-05   TIME AVG. MOMENT      2.126E-05
 LARGEST RESIDUAL MOMENT            5.941E-07   AT NODE        178   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.252E-03   AT NODE        124   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.898E-04   AT NODE        233   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.296E-06   TIME AVG. FORCE       1.829E-06
 LARGEST RESIDUAL FORCE             5.377E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.443       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.663E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.660E-05   TIME AVG. MOMENT      2.126E-05
 LARGEST RESIDUAL MOMENT           -6.311E-09   AT NODE        119   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.227E-03   AT NODE        124   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.556E-05   AT NODE        208   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.299E-06   TIME AVG. FORCE       1.829E-06
 LARGEST RESIDUAL FORCE             1.431E-11   AT NODE         24   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.443       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.194E-06   AT NODE        360   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.659E-05   TIME AVG. MOMENT      2.126E-05
 LARGEST RESIDUAL MOMENT           -9.108E-12   AT NODE         23   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.227E-03   AT NODE        123   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.080E-06   AT NODE        260   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    8.894E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.847E-02,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         1.00    


  INCREMENT    89 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.894E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.449E-02   TIME AVG. FORCE       5.722E-04
 LARGEST RESIDUAL FORCE             0.541       AT NODE        220   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          6.27       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         7.29       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.943E-02   TIME AVG. MOMENT      3.517E-04
 LARGEST RESIDUAL MOMENT           -0.514       AT NODE        215   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      0.163       AT NODE        242   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     0.181       AT NODE        148   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 5 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.196E-02   TIME AVG. FORCE       1.053E-03
 LARGEST RESIDUAL FORCE             -10.3       AT NODE        132   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          8.27       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         2.01       AT NODE        381   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.787E-02   TIME AVG. MOMENT      8.960E-04
 LARGEST RESIDUAL MOMENT            -6.38       AT NODE        110   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      0.733       AT NODE        132   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.754       AT NODE        132   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 103 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                       1.15       TIME AVG. FORCE       1.477E-02
 LARGEST RESIDUAL FORCE              437.       AT NODE        155   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          4.62       AT NODE        349   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         4.32       AT NODE        160   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      1.86       TIME AVG. MOMENT      2.096E-02
 LARGEST RESIDUAL MOMENT            -537.       AT NODE        155   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION       4.19       AT NODE        142   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION      8.84       AT NODE        160   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 191 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                       19.9       TIME AVG. FORCE       0.255    
 LARGEST RESIDUAL FORCE             7.197E+03   AT NODE        172   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -43.2       AT NODE        172   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -44.8       AT NODE        172   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      159.       TIME AVG. MOMENT       1.78    
 LARGEST RESIDUAL MOMENT           -1.146E+05   AT NODE        172   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION       5.81       AT NODE        220   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     -198.       AT NODE        172   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    89 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  2.224E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.265E-06   TIME AVG. FORCE       1.899E-06
 LARGEST RESIDUAL FORCE             1.235E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.288       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.247E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.671E-05   TIME AVG. MOMENT      2.143E-05
 LARGEST RESIDUAL MOMENT            4.312E-06   AT NODE         87   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.437E-03   AT NODE        288   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.544E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.417E-06   TIME AVG. FORCE       1.862E-06
 LARGEST RESIDUAL FORCE             2.562E-08   AT NODE        137   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.288       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.221E-04   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.677E-05   TIME AVG. MOMENT      2.143E-05
 LARGEST RESIDUAL MOMENT            2.637E-08   AT NODE        138   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.430E-03   AT NODE        288   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.316E-05   AT NODE        176   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.423E-06   TIME AVG. FORCE       1.862E-06
 LARGEST RESIDUAL FORCE             1.335E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.288       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.100E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.677E-05   TIME AVG. MOMENT      2.143E-05
 LARGEST RESIDUAL MOMENT            2.673E-11   AT NODE         19   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.427E-03   AT NODE        285   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.198E-06   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.335E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.224E-02,  FRACTION OF STEP COMPLETED  0.511    
 STEP TIME COMPLETED        1.02    ,  TOTAL TIME COMPLETED         1.02    


  INCREMENT    90 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.335E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.028E-05   TIME AVG. FORCE       1.969E-06
 LARGEST RESIDUAL FORCE            -9.950E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.358       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.415E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.697E-05   TIME AVG. MOMENT      2.161E-05
 LARGEST RESIDUAL MOMENT            8.430E-06   AT NODE        195   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.803E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.462E-03   AT NODE        112   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.454E-06   TIME AVG. FORCE       1.895E-06
 LARGEST RESIDUAL FORCE             6.513E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.355       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.823E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.701E-05   TIME AVG. MOMENT      2.161E-05
 LARGEST RESIDUAL MOMENT            3.492E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.677E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.145E-04   AT NODE        244   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.477E-06   TIME AVG. FORCE       1.895E-06
 LARGEST RESIDUAL FORCE            -1.408E-08   AT NODE        232   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.354       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.881E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.700E-05   TIME AVG. MOMENT      2.161E-05
 LARGEST RESIDUAL MOMENT            1.288E-08   AT NODE        174   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.660E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.603E-05   AT NODE        232   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.474E-06   TIME AVG. FORCE       1.895E-06
 LARGEST RESIDUAL FORCE            -1.111E-11   AT NODE         24   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.354       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.082E-06   AT NODE        361   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.700E-05   TIME AVG. MOMENT      2.161E-05
 LARGEST RESIDUAL MOMENT            7.018E-12   AT NODE         23   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.660E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.462E-07   AT NODE        257   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.003E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.335E-02,  FRACTION OF STEP COMPLETED  0.528    
 STEP TIME COMPLETED        1.06    ,  TOTAL TIME COMPLETED         1.06    


  INCREMENT    91 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.003E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.230E-05   TIME AVG. FORCE       3.025E-06
 LARGEST RESIDUAL FORCE             2.282E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.190       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.341       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.938E-05   TIME AVG. MOMENT      2.180E-05
 LARGEST RESIDUAL MOMENT            7.715E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.651E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.002E-03   AT NODE        345   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.330E-06   TIME AVG. FORCE       1.926E-06
 LARGEST RESIDUAL FORCE            -2.413E-06   AT NODE        276   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.190       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.929E-03   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.716E-05   TIME AVG. MOMENT      2.178E-05
 LARGEST RESIDUAL MOMENT            1.264E-06   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.592E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.803E-04   AT NODE        276   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.712E-06   TIME AVG. FORCE       1.930E-06
 LARGEST RESIDUAL FORCE             2.611E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.172       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.856E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.713E-05   TIME AVG. MOMENT      2.178E-05
 LARGEST RESIDUAL MOMENT            1.730E-06   AT NODE        198   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.840E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.524E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.290E-06   TIME AVG. FORCE       1.925E-06
 LARGEST RESIDUAL FORCE             2.672E-08   AT NODE        273   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.172       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.090E-04   AT NODE        357   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.710E-05   TIME AVG. MOMENT      2.178E-05
 LARGEST RESIDUAL MOMENT            8.879E-09   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.837E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.969E-05   AT NODE        276   DOF  6
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -3.491E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.298E-06   TIME AVG. FORCE       1.925E-06
 LARGEST RESIDUAL FORCE             3.537E-10   AT NODE        269   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.172       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.513E-05   AT NODE        365   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.711E-05   TIME AVG. MOMENT      2.178E-05
 LARGEST RESIDUAL MOMENT            1.344E-10   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.837E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.719E-06   AT NODE        269   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.003E-02,  FRACTION OF STEP COMPLETED  0.553    
 STEP TIME COMPLETED        1.11    ,  TOTAL TIME COMPLETED         1.11    


  INCREMENT    92 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.003E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.964E-05   TIME AVG. FORCE       2.638E-06
 LARGEST RESIDUAL FORCE             1.703E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         9.805E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.270       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.711E-05   TIME AVG. MOMENT      2.194E-05
 LARGEST RESIDUAL MOMENT            2.308E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.014E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.850E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.335E-06   TIME AVG. FORCE       1.955E-06
 LARGEST RESIDUAL FORCE            -4.399E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         9.774E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.650E-04   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.706E-05   TIME AVG. MOMENT      2.194E-05
 LARGEST RESIDUAL MOMENT            6.657E-09   AT NODE        210   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.011E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.481E-05   AT NODE        273   DOF  6
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.004E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.353E-06   TIME AVG. FORCE       1.955E-06
 LARGEST RESIDUAL FORCE             3.002E-09   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         9.794E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.962E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.705E-05   TIME AVG. MOMENT      2.194E-05
 LARGEST RESIDUAL MOMENT            6.435E-10   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.019E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.992E-06   AT NODE        244   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.003E-02,  FRACTION OF STEP COMPLETED  0.578    
 STEP TIME COMPLETED        1.16    ,  TOTAL TIME COMPLETED         1.16    


  INCREMENT    93 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.003E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.241E-05   TIME AVG. FORCE       2.083E-06
 LARGEST RESIDUAL FORCE             2.354E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         0.201       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.103       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.696E-05   TIME AVG. MOMENT      2.210E-05
 LARGEST RESIDUAL MOMENT            6.468E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      4.197E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.178E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.465E-06   TIME AVG. FORCE       1.986E-06
 LARGEST RESIDUAL FORCE            -9.040E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         0.201       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.023E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.692E-05   TIME AVG. MOMENT      2.210E-05
 LARGEST RESIDUAL MOMENT            1.728E-10   AT NODE        219   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      4.197E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.702E-06   AT NODE        261   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.504E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.003E-02,  FRACTION OF STEP COMPLETED  0.603    
 STEP TIME COMPLETED        1.21    ,  TOTAL TIME COMPLETED         1.21    


  INCREMENT    94 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.504E-02

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.164E-05   TIME AVG. FORCE       2.102E-06
 LARGEST RESIDUAL FORCE             2.218E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         0.203       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.899E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.692E-05   TIME AVG. MOMENT      2.226E-05
 LARGEST RESIDUAL MOMENT           -7.754E-06   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      4.242E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.064E-03   AT NODE        392   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.434E-06   TIME AVG. FORCE       2.015E-06
 LARGEST RESIDUAL FORCE            -7.121E-08   AT NODE        264   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         0.203       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.640E-04   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.680E-05   TIME AVG. MOMENT      2.226E-05
 LARGEST RESIDUAL MOMENT            3.375E-08   AT NODE        198   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      4.227E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.082E-05   AT NODE        264   DOF  6
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    3.518E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.453E-06   TIME AVG. FORCE       2.016E-06
 LARGEST RESIDUAL FORCE             5.841E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         0.201       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.046E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.681E-05   TIME AVG. MOMENT      2.226E-05
 LARGEST RESIDUAL MOMENT           -1.531E-10   AT NODE        171   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      4.192E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.432E-05   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.451E-06   TIME AVG. FORCE       2.016E-06
 LARGEST RESIDUAL FORCE            -5.318E-12   AT NODE        303   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         0.201       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.320E-08   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.681E-05   TIME AVG. MOMENT      2.226E-05
 LARGEST RESIDUAL MOMENT           -1.113E-12   AT NODE        299   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      4.192E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.075E-09   AT NODE        265   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.113    

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.504E-02,  FRACTION OF STEP COMPLETED  0.640    
 STEP TIME COMPLETED        1.28    ,  TOTAL TIME COMPLETED         1.28    


  INCREMENT    95 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.113    

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.347E-05   TIME AVG. FORCE       2.747E-06
 LARGEST RESIDUAL FORCE             1.636E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         2.802E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.274       AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.701E-05   TIME AVG. MOMENT      2.242E-05
 LARGEST RESIDUAL MOMENT           -1.503E-05   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      6.060E-04   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.713E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.414E-06   TIME AVG. FORCE       2.044E-06
 LARGEST RESIDUAL FORCE             3.582E-07   AT NODE        265   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.761E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.002E-03   AT NODE        401   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.680E-05   TIME AVG. MOMENT      2.241E-05
 LARGEST RESIDUAL MOMENT            1.742E-07   AT NODE        202   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      5.920E-04   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.813E-04   AT NODE        268   DOF  6
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.459E-06   TIME AVG. FORCE       2.045E-06
 LARGEST RESIDUAL FORCE             1.127E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         2.610E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.516E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.679E-05   TIME AVG. MOMENT      2.241E-05
 LARGEST RESIDUAL MOMENT           -1.745E-08   AT NODE         83   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      5.417E-04   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.725E-05   AT NODE        240   DOF  6
   INSTANCE: PART-1-1                                                                        
          ROTATION CORRECTION TOO LARGE COMPARED TO ROTATION INCREMENT

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.450E-06   TIME AVG. FORCE       2.045E-06
 LARGEST RESIDUAL FORCE            -1.222E-11   AT NODE         24   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.610E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.524E-06   AT NODE        365   DOF  1
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.679E-05   TIME AVG. MOMENT      2.241E-05
 LARGEST RESIDUAL MOMENT            7.706E-12   AT NODE         23   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      5.416E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.700E-07   AT NODE        265   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.169    

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.113    ,  FRACTION OF STEP COMPLETED  0.697    
 STEP TIME COMPLETED        1.39    ,  TOTAL TIME COMPLETED         1.39    


  INCREMENT    96 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.169    

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.150E-05   TIME AVG. FORCE       2.156E-06
 LARGEST RESIDUAL FORCE             2.035E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -5.816E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.731E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.687E-05   TIME AVG. MOMENT      2.256E-05
 LARGEST RESIDUAL MOMENT            1.960E-07   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.211E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.024E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.455E-06   TIME AVG. FORCE       2.073E-06
 LARGEST RESIDUAL FORCE            -6.219E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -5.821E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.672E-05   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.683E-05   TIME AVG. MOMENT      2.256E-05
 LARGEST RESIDUAL MOMENT           -5.126E-11   AT NODE        148   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.211E-03   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.957E-07   AT NODE        257   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.253    

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.169    ,  FRACTION OF STEP COMPLETED  0.781    
 STEP TIME COMPLETED        1.56    ,  TOTAL TIME COMPLETED         1.56    


  INCREMENT    97 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.253    

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.098E-06   TIME AVG. FORCE       2.143E-06
 LARGEST RESIDUAL FORCE             1.077E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -1.643E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.088E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.687E-05   TIME AVG. MOMENT      2.271E-05
 LARGEST RESIDUAL MOMENT            3.553E-07   AT NODE        399   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.422E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.475E-03   AT NODE        401   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.457E-06   TIME AVG. FORCE       2.101E-06
 LARGEST RESIDUAL FORCE            -1.831E-10   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -1.645E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.264E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.684E-05   TIME AVG. MOMENT      2.271E-05
 LARGEST RESIDUAL MOMENT           -3.281E-11   AT NODE        132   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.421E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.242E-06   AT NODE        264   DOF  6
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.253    ,  FRACTION OF STEP COMPLETED  0.908    
 STEP TIME COMPLETED        1.82    ,  TOTAL TIME COMPLETED         1.82    


  INCREMENT    98 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.185    

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.566E-06   TIME AVG. FORCE       2.129E-06
 LARGEST RESIDUAL FORCE             4.184E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         1.951E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.394E-02   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.684E-05   TIME AVG. MOMENT      2.285E-05
 LARGEST RESIDUAL MOMENT            1.572E-09   AT NODE        398   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      4.066E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.900E-04   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
          ROTATION CORRECTION TOO LARGE COMPARED TO ROTATION INCREMENT

  NUMBER OF EQUATIONS =      2415     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  5.10E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.458E-06   TIME AVG. FORCE       2.128E-06
 LARGEST RESIDUAL FORCE             6.935E-12   AT NODE        367   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.950E-03   AT NODE        401   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.681E-06   AT NODE        401   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     3.684E-05   TIME AVG. MOMENT      2.285E-05
 LARGEST RESIDUAL MOMENT            1.340E-12   AT NODE        141   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      4.066E-05   AT NODE        404   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.813E-09   AT NODE         88   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.185    ,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        2.00    ,  TOTAL TIME COMPLETED         2.00    


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF         98  INCREMENTS
                      11  CUTBACKS IN AUTOMATIC INCREMENTATION
                     490  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     490  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     490  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       2  WARNING MESSAGES DURING USER INPUT PROCESSING
                      40  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                      32  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   23.290    
       SYSTEM TIME (SEC)    =   1.1000    
       TOTAL CPU TIME (SEC) =   24.390    
       WALLCLOCK TIME (SEC) =         28

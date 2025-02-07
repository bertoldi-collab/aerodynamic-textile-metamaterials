1

   Abaqus 2021                                  Date 14-Feb-2023   Time 14:19:51
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

     CHARACTERISTIC ELEMENT LENGTH     0.851    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS
          THE INITIAL TIME INCREMENT WAS REDUCED TO   7.400E-02 TO MEET THE NEXT OUTPUT TIME POINT

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.400E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.372E-04   TIME AVG. FORCE       2.372E-04
 LARGEST RESIDUAL FORCE            -2.323E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.751       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.751       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.241E-05   TIME AVG. MOMENT      2.241E-05
 LARGEST RESIDUAL MOMENT            9.092E-04   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.847E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.847E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.919E-07   TIME AVG. FORCE       1.919E-07
 LARGEST RESIDUAL FORCE             1.442E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.752       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.217E-03   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.333E-07   TIME AVG. MOMENT      5.333E-07
 LARGEST RESIDUAL MOMENT           -2.409E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.853E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.486E-05   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.063E-08   TIME AVG. FORCE       1.176E-02
                                                (FROM MOMENT    FIELD)
 LARGEST RESIDUAL FORCE            -2.020E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.750       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.770E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     4.978E-07   TIME AVG. MOMENT      4.978E-07
 LARGEST RESIDUAL MOMENT           -3.937E-09   AT NODE        143   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.844E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.268E-05   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.998E-08   TIME AVG. FORCE       1.176E-02
                                                (FROM MOMENT    FIELD)
 LARGEST RESIDUAL FORCE            -1.175E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.751       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.225E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE

 AVERAGE MOMENT                     4.981E-07   TIME AVG. MOMENT      4.981E-07
 LARGEST RESIDUAL MOMENT           -8.011E-12   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.846E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.032E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   7.039E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.400E-02,  FRACTION OF STEP COMPLETED  3.700E-02
 STEP TIME COMPLETED       7.400E-02,  TOTAL TIME COMPLETED        7.400E-02


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.039E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.

     NEGATIVE EIGENVALUES MEAN THAT THE SYSTEM MATRIX IS NOT 
     POSITIVE DEFINITE:
     FOR EXAMPLE, A BIFURCATION (BUCKLING) LOAD MAY HAVE BEEN EXCEEDED.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.029E-02   TIME AVG. FORCE       1.029E-02
 LARGEST RESIDUAL FORCE             0.595       AT NODE        105   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          2.23       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.         2.95       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.041E-03   TIME AVG. MOMENT      2.021E-03
 LARGEST RESIDUAL MOMENT            0.104       AT NODE        102   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.185       AT NODE        105   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.185       AT NODE        105   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.332E-04   TIME AVG. FORCE       4.332E-04
 LARGEST RESIDUAL FORCE            -1.977E-02   AT NODE        109   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          2.29       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.314       AT NODE        189   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.363E-04   TIME AVG. MOMENT      6.839E-05
 LARGEST RESIDUAL MOMENT            2.893E-03   AT NODE        185   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.186       AT NODE        109   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.183E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 10 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.676E-04   TIME AVG. FORCE       6.676E-04
 LARGEST RESIDUAL FORCE             3.535E-02   AT NODE         45   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          2.74       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.549       AT NODE        108   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.413E-04   TIME AVG. MOMENT      2.709E-04
 LARGEST RESIDUAL MOMENT           -2.485E-02   AT NODE         46   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.190       AT NODE        112   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.786E-02   AT NODE         37   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.231E-04   TIME AVG. FORCE       1.231E-04
 LARGEST RESIDUAL FORCE             3.584E-03   AT NODE         45   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          2.69       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.110E-02   AT NODE        189   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.104E-04   TIME AVG. MOMENT      5.544E-05
 LARGEST RESIDUAL MOMENT            2.045E-03   AT NODE         54   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.186       AT NODE        105   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.317E-02   AT NODE         41   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.422E-04   TIME AVG. FORCE       4.422E-04
 LARGEST RESIDUAL FORCE            -1.449E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.          1.93       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.759       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.064E-04   TIME AVG. MOMENT      1.534E-04
 LARGEST RESIDUAL MOMENT           -6.506E-03   AT NODE         79   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.172       AT NODE        109   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.775E-02   AT NODE         41   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      8.244E-05   TIME AVG. FORCE       8.244E-05
 LARGEST RESIDUAL FORCE            -4.781E-03   AT NODE        148   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.          1.91       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.360E-02   AT NODE        156   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.635E-05   TIME AVG. MOMENT      1.842E-05
 LARGEST RESIDUAL MOMENT            3.432E-04   AT NODE        138   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.169       AT NODE        105   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.682E-02   AT NODE        148   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      8.314E-04   TIME AVG. FORCE       8.314E-04
 LARGEST RESIDUAL FORCE             2.370E-02   AT NODE        101   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         0.793       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -1.11       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.934E-04   TIME AVG. MOMENT      2.970E-04
 LARGEST RESIDUAL MOMENT           -1.683E-02   AT NODE        102   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.153       AT NODE        117   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.550E-02   AT NODE         76   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     2 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  1.760E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.964E-06   TIME AVG. FORCE       1.964E-06
 LARGEST RESIDUAL FORCE            -1.455E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.230       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.164E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.524E-06   TIME AVG. MOMENT      1.511E-06
 LARGEST RESIDUAL MOMENT           -3.208E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.080E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.313E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.723E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -4.343E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.223       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.942E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     1.314E-06   TIME AVG. MOMENT      9.062E-07
 LARGEST RESIDUAL MOMENT            1.189E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.687E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.930E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.720E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -4.890E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.217       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.495E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     1.298E-06   TIME AVG. MOMENT      8.983E-07
 LARGEST RESIDUAL MOMENT            9.041E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.273E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.145E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    3.150E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.639E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.487E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.217       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.645E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     1.298E-06   TIME AVG. MOMENT      8.983E-07
 LARGEST RESIDUAL MOMENT            1.111E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.262E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.105E-05   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    2.640E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.760E-02,  FRACTION OF STEP COMPLETED  4.580E-02
 STEP TIME COMPLETED       9.160E-02,  TOTAL TIME COMPLETED        9.160E-02


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.640E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.293E-06   TIME AVG. FORCE       2.293E-06
 LARGEST RESIDUAL FORCE            -1.664E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.381       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.618E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.191E-06   TIME AVG. MOMENT      1.662E-06
 LARGEST RESIDUAL MOMENT           -3.698E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.488E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.486E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.593E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -3.766E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.374       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.228E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     1.809E-06   TIME AVG. MOMENT      1.202E-06
 LARGEST RESIDUAL MOMENT            1.291E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.451E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.660E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      9.291E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.091E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.363       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.123E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     1.782E-06   TIME AVG. MOMENT      1.193E-06
 LARGEST RESIDUAL MOMENT           -5.303E-09   AT NODE        175   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.389E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.189E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.542E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.726E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.073E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.362       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.410E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     1.783E-06   TIME AVG. MOMENT      1.193E-06
 LARGEST RESIDUAL MOMENT            4.349E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.387E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.936E-05   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    3.959E-02
          TIME INCREMENT REDUCED TO   2.640E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.640E-02,  FRACTION OF STEP COMPLETED  5.900E-02
 STEP TIME COMPLETED       0.118    ,  TOTAL TIME COMPLETED        0.118    


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.640E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.483E-06   TIME AVG. FORCE       2.483E-06
 LARGEST RESIDUAL FORCE             2.517E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.410       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.784E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.326E-06   TIME AVG. MOMENT      1.976E-06
 LARGEST RESIDUAL MOMENT           -5.930E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.652E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.644E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.452E-07   TIME AVG. FORCE       2.452E-07
 LARGEST RESIDUAL FORCE            -3.023E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.389       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.141E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.329E-06   TIME AVG. MOMENT      1.477E-06
 LARGEST RESIDUAL MOMENT            3.707E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.548E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.034E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.448E-07   TIME AVG. FORCE       3.448E-07
 LARGEST RESIDUAL FORCE            -5.699E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.364       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.476E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.277E-06   TIME AVG. MOMENT      1.464E-06
 LARGEST RESIDUAL MOMENT           -2.946E-08   AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.407E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.416E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.558E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -2.096E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.359       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.851E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.267E-06   TIME AVG. MOMENT      1.462E-06
 LARGEST RESIDUAL MOMENT            3.109E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.380E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.710E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.889E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -4.483E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.357       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.156E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.267E-06   TIME AVG. MOMENT      1.462E-06
 LARGEST RESIDUAL MOMENT            4.636E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.367E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.255E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.834E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.330E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.357       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.718E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     2.267E-06   TIME AVG. MOMENT      1.462E-06
 LARGEST RESIDUAL MOMENT            2.029E-12   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.367E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.163E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.640E-02,  FRACTION OF STEP COMPLETED  7.219E-02
 STEP TIME COMPLETED       0.144    ,  TOTAL TIME COMPLETED        0.144    


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.959E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.872E-05   TIME AVG. FORCE       4.872E-05
 LARGEST RESIDUAL FORCE            -5.692E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.815       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.279       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.469E-05   TIME AVG. MOMENT      6.107E-06
 LARGEST RESIDUAL MOMENT           -4.973E-04   AT NODE        171   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.502E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.452E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.670E-07   TIME AVG. FORCE       3.670E-07
 LARGEST RESIDUAL FORCE            -9.107E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.808       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.690E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.526E-06   TIME AVG. MOMENT      1.875E-06
 LARGEST RESIDUAL MOMENT            2.680E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.477E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.987E-04   AT NODE         68   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.926E-05   TIME AVG. FORCE       4.926E-05
 LARGEST RESIDUAL FORCE            -7.299E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.488       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.320       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.542E-06   TIME AVG. MOMENT      2.678E-06
 LARGEST RESIDUAL MOMENT           -1.631E-04   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.867E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.610E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.792E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE             1.380E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.487       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.932E-03   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.916E-06   TIME AVG. MOMENT      1.753E-06
 LARGEST RESIDUAL MOMENT            2.190E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.864E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.365E-04   AT NODE        112   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.609E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -6.933E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.484       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.622E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.911E-06   TIME AVG. MOMENT      1.752E-06
 LARGEST RESIDUAL MOMENT           -8.804E-10   AT NODE        159   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.848E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.566E-04   AT NODE        181   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.380E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.489E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.484       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.287E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     2.911E-06   TIME AVG. MOMENT      1.752E-06
 LARGEST RESIDUAL MOMENT           -3.422E-12   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.849E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.301E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   2.752E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-02,  FRACTION OF STEP COMPLETED  9.199E-02
 STEP TIME COMPLETED       0.184    ,  TOTAL TIME COMPLETED        0.184    


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.752E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.410E-06   TIME AVG. FORCE       2.410E-06
 LARGEST RESIDUAL FORCE             2.680E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.381       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.445E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.410E-06   TIME AVG. MOMENT      2.361E-06
 LARGEST RESIDUAL MOMENT           -6.201E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.549E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.642E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.759E-07   TIME AVG. FORCE       3.759E-07
 LARGEST RESIDUAL FORCE            -4.382E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.354       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.632E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.409E-06   TIME AVG. MOMENT      2.028E-06
 LARGEST RESIDUAL MOMENT            4.211E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.425E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.245E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.201E-07   TIME AVG. FORCE       4.201E-07
 LARGEST RESIDUAL FORCE            -6.595E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.328       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.651E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.359E-06   TIME AVG. MOMENT      2.020E-06
 LARGEST RESIDUAL MOMENT           -1.272E-08   AT NODE        171   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.272E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.526E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      9.495E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -4.917E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.320       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.551E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     3.343E-06   TIME AVG. MOMENT      2.017E-06
 LARGEST RESIDUAL MOMENT            5.387E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.231E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.154E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.764E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.416E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.316       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.871E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     3.338E-06   TIME AVG. MOMENT      2.016E-06
 LARGEST RESIDUAL MOMENT            1.828E-09   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.208E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.233E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -3.210E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.422E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -2.476E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.316       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.622E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     3.340E-06   TIME AVG. MOMENT      2.016E-06
 LARGEST RESIDUAL MOMENT            2.730E-11   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.207E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.323E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.752E-02,  FRACTION OF STEP COMPLETED  0.106    
 STEP TIME COMPLETED       0.212    ,  TOTAL TIME COMPLETED        0.212    


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.959E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.338E-05   TIME AVG. FORCE       1.338E-05
 LARGEST RESIDUAL FORCE            -1.262E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.592       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.137       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.188E-05   TIME AVG. MOMENT      3.425E-06
 LARGEST RESIDUAL MOMENT           -1.981E-04   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.433E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.960E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.329E-07   TIME AVG. FORCE       1.329E-07
 LARGEST RESIDUAL FORCE            -4.380E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.581       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.066E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.175E-06   TIME AVG. MOMENT      2.325E-06
 LARGEST RESIDUAL MOMENT            1.088E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.394E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.380E-04   AT NODE        121   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.046E-06   TIME AVG. FORCE       8.046E-06
 LARGEST RESIDUAL FORCE            -1.207E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.453       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.128       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.180E-06   TIME AVG. MOMENT      2.325E-06
 LARGEST RESIDUAL MOMENT           -1.094E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.740E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.539E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      8.802E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -2.150E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.452       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.389E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     3.942E-06   TIME AVG. MOMENT      2.291E-06
 LARGEST RESIDUAL MOMENT            6.183E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.737E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.875E-05   AT NODE        116   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.099E-07   TIME AVG. FORCE       1.099E-07
 LARGEST RESIDUAL FORCE            -4.417E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.445       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.777E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.928E-06   TIME AVG. MOMENT      2.289E-06
 LARGEST RESIDUAL MOMENT           -2.439E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.698E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.948E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.557E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      8.924E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE             9.728E-13   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.445       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.201E-06   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     3.932E-06   TIME AVG. MOMENT      2.290E-06
 LARGEST RESIDUAL MOMENT           -3.754E-13   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.698E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.767E-08   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   2.454E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-02,  FRACTION OF STEP COMPLETED  0.126    
 STEP TIME COMPLETED       0.251    ,  TOTAL TIME COMPLETED        0.251    


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.454E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.868E-07   TIME AVG. FORCE       9.868E-07
 LARGEST RESIDUAL FORCE             1.066E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.302       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.687E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.089E-06   TIME AVG. MOMENT      2.640E-06
 LARGEST RESIDUAL MOMENT           -2.481E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.198E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.461E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.407E-07   TIME AVG. FORCE       2.407E-07
 LARGEST RESIDUAL FORCE            -2.475E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.285       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.775E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.317E-06   TIME AVG. MOMENT      2.543E-06
 LARGEST RESIDUAL MOMENT            8.175E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.105E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.361E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.147E-07   TIME AVG. FORCE       1.147E-07
 LARGEST RESIDUAL FORCE            -4.248E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.278       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.599E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.304E-06   TIME AVG. MOMENT      2.542E-06
 LARGEST RESIDUAL MOMENT            1.759E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.066E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.874E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      9.898E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.692E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.277       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.309E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     4.306E-06   TIME AVG. MOMENT      2.542E-06
 LARGEST RESIDUAL MOMENT            6.614E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.058E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.727E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.884E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.809E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.277       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.191E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     4.305E-06   TIME AVG. MOMENT      2.542E-06
 LARGEST RESIDUAL MOMENT            9.865E-13   AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.058E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.523E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.454E-02,  FRACTION OF STEP COMPLETED  0.138    
 STEP TIME COMPLETED       0.276    ,  TOTAL TIME COMPLETED        0.276    


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.959E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.129E-05   TIME AVG. FORCE       1.129E-05
 LARGEST RESIDUAL FORCE            -9.993E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.573       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.126       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.156E-05   TIME AVG. MOMENT      3.543E-06
 LARGEST RESIDUAL MOMENT           -1.658E-04   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.333E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.261E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.471E-07   TIME AVG. FORCE       1.471E-07
 LARGEST RESIDUAL FORCE            -4.353E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.563       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.003E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.105E-06   TIME AVG. MOMENT      2.827E-06
 LARGEST RESIDUAL MOMENT            8.722E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.295E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.126E-04   AT NODE        125   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.344E-06   TIME AVG. FORCE       5.344E-06
 LARGEST RESIDUAL FORCE            -7.793E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.460       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.103       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.034E-06   TIME AVG. MOMENT      2.819E-06
 LARGEST RESIDUAL MOMENT           -5.696E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.769E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.263E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.141E-07   TIME AVG. FORCE       1.141E-07
 LARGEST RESIDUAL FORCE            -1.622E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.459       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.616E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.918E-06   TIME AVG. MOMENT      2.806E-06
 LARGEST RESIDUAL MOMENT            4.686E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.766E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.031E-05   AT NODE        137   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.236E-07   TIME AVG. FORCE       1.236E-07
 LARGEST RESIDUAL FORCE            -2.346E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.453       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.596E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.913E-06   TIME AVG. MOMENT      2.805E-06
 LARGEST RESIDUAL MOMENT           -9.365E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.737E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.882E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -5.761E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.149E-07   TIME AVG. FORCE       1.149E-07
 LARGEST RESIDUAL FORCE            -1.120E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.453       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.195E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     4.913E-06   TIME AVG. MOMENT      2.805E-06
 LARGEST RESIDUAL MOMENT            3.602E-13   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.737E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.201E-08   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   2.181E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-02,  FRACTION OF STEP COMPLETED  0.158    
 STEP TIME COMPLETED       0.315    ,  TOTAL TIME COMPLETED        0.315    


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  2.181E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.813E-06   TIME AVG. FORCE       9.639E-07
 LARGEST RESIDUAL FORCE             6.681E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.268       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.834E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.716E-06   TIME AVG. MOMENT      3.096E-06
 LARGEST RESIDUAL MOMENT           -1.564E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.060E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.033E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.591E-07   TIME AVG. FORCE       3.870E-07
 LARGEST RESIDUAL FORCE            -1.801E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.254       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.438E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.258E-06   TIME AVG. MOMENT      3.051E-06
 LARGEST RESIDUAL MOMENT            2.807E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.799E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.000E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.753E-07   TIME AVG. FORCE       2.451E-07
 LARGEST RESIDUAL FORCE            -6.913E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.251       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.536E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.253E-06   TIME AVG. MOMENT      3.050E-06
 LARGEST RESIDUAL MOMENT            5.706E-09   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.643E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.564E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    3.179E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.726E-07   TIME AVG. FORCE       2.438E-07
 LARGEST RESIDUAL FORCE            -1.219E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.251       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.339E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.253E-06   TIME AVG. MOMENT      3.050E-06
 LARGEST RESIDUAL MOMENT            2.939E-11   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.622E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.079E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.726E-07   TIME AVG. FORCE       2.438E-07
 LARGEST RESIDUAL FORCE            -1.211E-12   AT NODE        178   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.251       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.560E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     5.253E-06   TIME AVG. MOMENT      3.050E-06
 LARGEST RESIDUAL MOMENT           -1.377E-13   AT NODE        177   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.622E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.900E-08   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  2.181E-02,  FRACTION OF STEP COMPLETED  0.169    
 STEP TIME COMPLETED       0.337    ,  TOTAL TIME COMPLETED        0.337    


  INCREMENT    11 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.959E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.855E-05   TIME AVG. FORCE       9.680E-06
 LARGEST RESIDUAL FORCE            -8.606E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.572       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.117       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.206E-05   TIME AVG. MOMENT      3.869E-06
 LARGEST RESIDUAL MOMENT           -1.577E-04   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.337E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.903E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.167E-07   TIME AVG. FORCE       3.348E-07
 LARGEST RESIDUAL FORCE            -5.530E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.561       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.123E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.046E-06   TIME AVG. MOMENT      3.322E-06
 LARGEST RESIDUAL MOMENT            8.710E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.294E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.615E-04   AT NODE        125   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.545E-05   TIME AVG. FORCE       5.312E-06
 LARGEST RESIDUAL FORCE            -7.349E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.461       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.935E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.970E-06   TIME AVG. MOMENT      3.316E-06
 LARGEST RESIDUAL MOMENT           -5.242E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.782E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.123E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.208E-07   TIME AVG. FORCE       3.028E-07
 LARGEST RESIDUAL FORCE            -1.978E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.460       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.797E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.875E-06   TIME AVG. MOMENT      3.307E-06
 LARGEST RESIDUAL MOMENT            5.797E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.778E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.928E-05   AT NODE        137   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.695E-07   TIME AVG. FORCE       3.190E-07
 LARGEST RESIDUAL FORCE            -3.542E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.453       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.850E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.862E-06   TIME AVG. MOMENT      3.306E-06
 LARGEST RESIDUAL MOMENT           -1.734E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.742E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.549E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.062E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.237E-07   TIME AVG. FORCE       3.038E-07
 LARGEST RESIDUAL FORCE            -4.149E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.453       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.221E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     5.864E-06   TIME AVG. MOMENT      3.306E-06
 LARGEST RESIDUAL MOMENT            1.278E-12   AT NODE        175   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.742E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.506E-07   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   1.930E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-02,  FRACTION OF STEP COMPLETED  0.188    
 STEP TIME COMPLETED       0.377    ,  TOTAL TIME COMPLETED        0.377    


  INCREMENT    12 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.930E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.029E-06   TIME AVG. FORCE       4.851E-07
 LARGEST RESIDUAL FORCE             3.105E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.230       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.049E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.353E-06   TIME AVG. MOMENT      3.560E-06
 LARGEST RESIDUAL MOMENT           -7.418E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.066E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.733E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.814E-07   TIME AVG. FORCE       3.732E-07
 LARGEST RESIDUAL FORCE            -9.387E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.220       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.842E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.156E-06   TIME AVG. MOMENT      3.543E-06
 LARGEST RESIDUAL MOMENT           -8.976E-09   AT NODE        131   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.488E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.788E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    7.003E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.488E-07   TIME AVG. FORCE       3.400E-07
 LARGEST RESIDUAL FORCE            -2.762E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.220       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.718E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.154E-06   TIME AVG. MOMENT      3.543E-06
 LARGEST RESIDUAL MOMENT            6.102E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.456E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.121E-05   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.486E-07   TIME AVG. FORCE       3.400E-07
 LARGEST RESIDUAL FORCE            -1.044E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.220       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.905E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     6.154E-06   TIME AVG. MOMENT      3.543E-06
 LARGEST RESIDUAL MOMENT            3.387E-13   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.455E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.926E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.930E-02,  FRACTION OF STEP COMPLETED  0.198    
 STEP TIME COMPLETED       0.396    ,  TOTAL TIME COMPLETED        0.396    


  INCREMENT    13 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.959E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.143E-05   TIME AVG. FORCE       4.557E-06
 LARGEST RESIDUAL FORCE             6.340E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.547       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.635E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.160E-05   TIME AVG. MOMENT      4.163E-06
 LARGEST RESIDUAL MOMENT           -1.264E-04   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.230E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.951E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.538E-07   TIME AVG. FORCE       4.027E-07
 LARGEST RESIDUAL FORCE            -8.018E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.534       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.305E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.914E-06   TIME AVG. MOMENT      3.803E-06
 LARGEST RESIDUAL MOMENT            7.405E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.177E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.396E-04   AT NODE        137   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.033E-05   TIME AVG. FORCE       2.338E-06
 LARGEST RESIDUAL FORCE            -4.789E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.455       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.958E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.822E-06   TIME AVG. MOMENT      3.795E-06
 LARGEST RESIDUAL MOMENT           -2.742E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.762E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.145E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.991E-07   TIME AVG. FORCE       3.718E-07
 LARGEST RESIDUAL FORCE            -2.084E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.453       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.531E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.773E-06   TIME AVG. MOMENT      3.792E-06
 LARGEST RESIDUAL MOMENT            6.621E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.756E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.568E-05   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.585E-07   TIME AVG. FORCE       3.837E-07
 LARGEST RESIDUAL FORCE            -4.198E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.446       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.375E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.754E-06   TIME AVG. MOMENT      3.790E-06
 LARGEST RESIDUAL MOMENT           -2.180E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.717E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.873E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.276E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      5.018E-07   TIME AVG. FORCE       3.723E-07
 LARGEST RESIDUAL FORCE            -1.565E-11   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.446       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.136E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     6.756E-06   TIME AVG. MOMENT      3.790E-06
 LARGEST RESIDUAL MOMENT            5.196E-12   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.717E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.734E-07   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   1.699E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-02,  FRACTION OF STEP COMPLETED  0.218    
 STEP TIME COMPLETED       0.436    ,  TOTAL TIME COMPLETED        0.436    


  INCREMENT    14 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.699E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.609E-07   TIME AVG. FORCE       4.371E-07
 LARGEST RESIDUAL FORCE             1.357E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.196       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.366E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.095E-06   TIME AVG. MOMENT      4.026E-06
 LARGEST RESIDUAL MOMENT           -3.301E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.677E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.109E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.601E-07   TIME AVG. FORCE       4.036E-07
 LARGEST RESIDUAL FORCE            -3.266E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.190       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.506E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.021E-06   TIME AVG. MOMENT      4.021E-06
 LARGEST RESIDUAL MOMENT           -2.260E-09   AT NODE        139   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.335E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.422E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    2.343E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.249E-07   TIME AVG. FORCE       3.978E-07
 LARGEST RESIDUAL FORCE            -6.074E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.190       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.428E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.021E-06   TIME AVG. MOMENT      4.021E-06
 LARGEST RESIDUAL MOMENT            3.122E-11   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -7.331E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.597E-06   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.699E-02,  FRACTION OF STEP COMPLETED  0.226    
 STEP TIME COMPLETED       0.453    ,  TOTAL TIME COMPLETED        0.453    


  INCREMENT    15 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.959E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.597E-05   TIME AVG. FORCE       2.623E-06
 LARGEST RESIDUAL FORCE             4.966E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.523       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.022E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.133E-05   TIME AVG. MOMENT      4.508E-06
 LARGEST RESIDUAL MOMENT           -9.952E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.123E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.147E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.905E-07   TIME AVG. FORCE       4.539E-07
 LARGEST RESIDUAL FORCE            -1.062E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.509       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.440E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.734E-06   TIME AVG. MOMENT      4.269E-06
 LARGEST RESIDUAL MOMENT            5.975E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.062E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.173E-04   AT NODE        172   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.425E-06   TIME AVG. FORCE       1.259E-06
 LARGEST RESIDUAL FORCE            -2.851E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.448       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.101E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.645E-06   TIME AVG. MOMENT      4.263E-06
 LARGEST RESIDUAL MOMENT           -1.232E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.741E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.205E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.783E-07   TIME AVG. FORCE       4.235E-07
 LARGEST RESIDUAL FORCE            -3.385E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.445       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.192E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.625E-06   TIME AVG. MOMENT      4.261E-06
 LARGEST RESIDUAL MOMENT            6.204E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.731E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.027E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.218E-07   TIME AVG. FORCE       4.298E-07
 LARGEST RESIDUAL FORCE            -3.363E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.439       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.537E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.614E-06   TIME AVG. MOMENT      4.261E-06
 LARGEST RESIDUAL MOMENT           -1.424E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.696E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.476E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -7.980E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      5.789E-07   TIME AVG. FORCE       4.236E-07
 LARGEST RESIDUAL FORCE            -4.717E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.439       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.325E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.613E-06   TIME AVG. MOMENT      4.261E-06
 LARGEST RESIDUAL MOMENT            8.737E-12   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.696E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.214E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   1.485E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-02,  FRACTION OF STEP COMPLETED  0.246    
 STEP TIME COMPLETED       0.492    ,  TOTAL TIME COMPLETED        0.492    


  INCREMENT    16 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.485E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.014E-07   TIME AVG. FORCE       4.584E-07
 LARGEST RESIDUAL FORCE             6.216E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.167       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.193E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.867E-06   TIME AVG. MOMENT      4.486E-06
 LARGEST RESIDUAL MOMENT           -1.535E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.539E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.763E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.075E-07   TIME AVG. FORCE       4.466E-07
 LARGEST RESIDUAL FORCE            -1.065E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.164       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.964E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.834E-06   TIME AVG. MOMENT      4.484E-06
 LARGEST RESIDUAL MOMENT           -5.439E-10   AT NODE        143   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.343E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.959E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    6.940E-08
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.005E-07   TIME AVG. FORCE       4.457E-07
 LARGEST RESIDUAL FORCE            -1.266E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.164       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.209E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.834E-06   TIME AVG. MOMENT      4.484E-06
 LARGEST RESIDUAL MOMENT            1.585E-12   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -6.342E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.448E-07   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.485E-02,  FRACTION OF STEP COMPLETED  0.253    
 STEP TIME COMPLETED       0.507    ,  TOTAL TIME COMPLETED        0.507    


  INCREMENT    17 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.959E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.244E-05   TIME AVG. FORCE       1.779E-06
 LARGEST RESIDUAL FORCE             4.008E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.505       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.819E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.138E-05   TIME AVG. MOMENT      4.890E-06
 LARGEST RESIDUAL MOMENT           -8.076E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.046E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.553E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.344E-07   TIME AVG. FORCE       5.000E-07
 LARGEST RESIDUAL FORCE            -1.334E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.489       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.560E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.530E-06   TIME AVG. MOMENT      4.722E-06
 LARGEST RESIDUAL MOMENT            4.884E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.976E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.944E-04   AT NODE        188   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.259E-06   TIME AVG. FORCE       8.694E-07
 LARGEST RESIDUAL FORCE            -1.765E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.442       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.773E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.444E-06   TIME AVG. MOMENT      4.717E-06
 LARGEST RESIDUAL MOMENT           -5.707E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.724E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.528E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.594E-07   TIME AVG. FORCE       4.695E-07
 LARGEST RESIDUAL FORCE            -6.108E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.439       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.890E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.440E-06   TIME AVG. MOMENT      4.717E-06
 LARGEST RESIDUAL MOMENT            5.364E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.709E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.430E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.842E-07   TIME AVG. FORCE       4.722E-07
 LARGEST RESIDUAL FORCE            -2.360E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.433       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.438E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.429E-06   TIME AVG. MOMENT      4.716E-06
 LARGEST RESIDUAL MOMENT           -6.176E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.680E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.920E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -3.361E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.566E-07   TIME AVG. FORCE       4.692E-07
 LARGEST RESIDUAL FORCE            -1.095E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.433       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.593E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     8.428E-06   TIME AVG. MOMENT      4.716E-06
 LARGEST RESIDUAL MOMENT            9.681E-12   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.680E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.912E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   1.287E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-02,  FRACTION OF STEP COMPLETED  0.273    
 STEP TIME COMPLETED       0.547    ,  TOTAL TIME COMPLETED        0.547    


  INCREMENT    18 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.287E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.159E-07   TIME AVG. FORCE       4.938E-07
 LARGEST RESIDUAL FORCE             2.700E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.142       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.643E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.625E-06   TIME AVG. MOMENT      4.933E-06
 LARGEST RESIDUAL MOMENT           -6.775E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.557E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.516E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.777E-07   TIME AVG. FORCE       4.900E-07
 LARGEST RESIDUAL FORCE            -3.395E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.140       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.564E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.620E-06   TIME AVG. MOMENT      4.933E-06
 LARGEST RESIDUAL MOMENT           -1.221E-10   AT NODE        147   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.447E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.109E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    1.998E-08
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.757E-07   TIME AVG. FORCE       4.898E-07
 LARGEST RESIDUAL FORCE             2.096E-12   AT NODE        182   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.140       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.513E-07   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     8.620E-06   TIME AVG. MOMENT      4.933E-06
 LARGEST RESIDUAL MOMENT            3.133E-13   AT NODE        180   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.446E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.215E-08   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.287E-02,  FRACTION OF STEP COMPLETED  0.280    
 STEP TIME COMPLETED       0.559    ,  TOTAL TIME COMPLETED        0.559    


  INCREMENT    19 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.959E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.763E-06   TIME AVG. FORCE       1.333E-06
 LARGEST RESIDUAL FORCE             3.231E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.489       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.769E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.156E-05   TIME AVG. MOMENT      5.282E-06
 LARGEST RESIDUAL MOMENT           -6.568E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.979E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.042E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.088E-06   TIME AVG. FORCE       5.442E-07
 LARGEST RESIDUAL FORCE            -1.663E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.472       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.691E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.293E-06   TIME AVG. MOMENT      5.162E-06
 LARGEST RESIDUAL MOMENT            3.937E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.902E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.784E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.889E-06   TIME AVG. FORCE       7.079E-07
 LARGEST RESIDUAL FORCE            -1.070E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.435       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.693E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.219E-06   TIME AVG. MOMENT      5.159E-06
 LARGEST RESIDUAL MOMENT           -2.369E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.704E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.973E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.418E-07   TIME AVG. FORCE       5.127E-07
 LARGEST RESIDUAL FORCE            -1.003E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.431       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.645E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.211E-06   TIME AVG. MOMENT      5.158E-06
 LARGEST RESIDUAL MOMENT            4.255E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.686E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.871E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.479E-07   TIME AVG. FORCE       5.133E-07
 LARGEST RESIDUAL FORCE            -1.417E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.427       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.188E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.212E-06   TIME AVG. MOMENT      5.158E-06
 LARGEST RESIDUAL MOMENT            8.821E-11   AT NODE        127   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.663E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.268E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    4.702E-07
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.334E-07   TIME AVG. FORCE       5.120E-07
 LARGEST RESIDUAL FORCE            -1.778E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.427       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.612E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     9.212E-06   TIME AVG. MOMENT      5.158E-06
 LARGEST RESIDUAL MOMENT            7.680E-12   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.663E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.487E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   1.103E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.959E-02,  FRACTION OF STEP COMPLETED  0.300    
 STEP TIME COMPLETED       0.599    ,  TOTAL TIME COMPLETED        0.599    


  INCREMENT    20 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  1.103E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.633E-07   TIME AVG. FORCE       5.329E-07
 LARGEST RESIDUAL FORCE             1.055E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.119       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.222E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.370E-06   TIME AVG. MOMENT      5.369E-06
 LARGEST RESIDUAL MOMENT           -2.692E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.680E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.688E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.507E-07   TIME AVG. FORCE       5.318E-07
 LARGEST RESIDUAL FORCE            -1.022E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.118       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.938E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.372E-06   TIME AVG. MOMENT      5.369E-06
 LARGEST RESIDUAL MOMENT           -2.507E-11   AT NODE        159   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.619E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.099E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    5.679E-09
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      7.501E-07   TIME AVG. FORCE       5.318E-07
 LARGEST RESIDUAL FORCE            -1.867E-12   AT NODE        135   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.118       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.117E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     9.372E-06   TIME AVG. MOMENT      5.369E-06
 LARGEST RESIDUAL MOMENT           -2.711E-13   AT NODE        176   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.619E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.799E-09   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    4.891E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  1.103E-02,  FRACTION OF STEP COMPLETED  0.305    
 STEP TIME COMPLETED       0.610    ,  TOTAL TIME COMPLETED        0.610    


  INCREMENT    21 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.891E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.669E-05   TIME AVG. FORCE       4.083E-06
 LARGEST RESIDUAL FORCE            -1.316E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.674       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.149       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.900E-05   TIME AVG. MOMENT      6.018E-06
 LARGEST RESIDUAL MOMENT           -2.155E-04   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.810E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.628E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.782E-07   TIME AVG. FORCE       5.585E-07
 LARGEST RESIDUAL FORCE            -3.612E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.665       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.809E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.033E-05   TIME AVG. MOMENT      5.605E-06
 LARGEST RESIDUAL MOMENT            1.051E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.776E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.751E-04   AT NODE        109   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.357E-05   TIME AVG. FORCE       3.073E-06
 LARGEST RESIDUAL FORCE            -1.382E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.526       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.139       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.042E-05   TIME AVG. MOMENT      5.609E-06
 LARGEST RESIDUAL MOMENT           -1.431E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.061E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.148E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      8.208E-07   TIME AVG. FORCE       5.540E-07
 LARGEST RESIDUAL FORCE            -1.674E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.525       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.633E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.009E-05   TIME AVG. MOMENT      5.594E-06
 LARGEST RESIDUAL MOMENT            4.504E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.060E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.929E-05   AT NODE        109   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      8.632E-07   TIME AVG. FORCE       5.573E-07
 LARGEST RESIDUAL FORCE            -2.852E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.519       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.353E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.009E-05   TIME AVG. MOMENT      5.593E-06
 LARGEST RESIDUAL MOMENT           -1.346E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.027E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.230E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -9.651E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      8.235E-07   TIME AVG. FORCE       5.542E-07
 LARGEST RESIDUAL FORCE             3.863E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.519       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.294E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.009E-05   TIME AVG. MOMENT      5.593E-06
 LARGEST RESIDUAL MOMENT           -1.091E-12   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.027E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.438E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.731E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.891E-02,  FRACTION OF STEP COMPLETED  0.329    
 STEP TIME COMPLETED       0.659    ,  TOTAL TIME COMPLETED        0.659    


  INCREMENT    22 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.731E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.951E-04   TIME AVG. FORCE       2.874E-05
 LARGEST RESIDUAL FORCE            -1.607E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.996       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.494       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.720E-05   TIME AVG. MOMENT      7.485E-06
 LARGEST RESIDUAL MOMENT           -6.860E-04   AT NODE        139   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.472E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.511E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.260E-06   TIME AVG. FORCE       6.761E-07
 LARGEST RESIDUAL FORCE            -2.059E-05   AT NODE         48   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.991       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.801E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.181E-05   TIME AVG. MOMENT      5.876E-06
 LARGEST RESIDUAL MOMENT            9.719E-06   AT NODE        126   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.451E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.083E-03   AT NODE         48   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.036E-04   TIME AVG. FORCE       3.649E-05
 LARGEST RESIDUAL FORCE            -1.918E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.459       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.532       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.644E-05   TIME AVG. MOMENT      6.995E-06
 LARGEST RESIDUAL MOMENT           -7.919E-04   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.737E-02   AT NODE        153   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.720E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      9.639E-07   TIME AVG. FORCE       5.835E-07
 LARGEST RESIDUAL FORCE             9.469E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.457       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.339E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.069E-05   TIME AVG. MOMENT      5.825E-06
 LARGEST RESIDUAL MOMENT            1.637E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.733E-02   AT NODE        153   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.199E-05   AT NODE        132   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.058E-06   TIME AVG. FORCE       8.045E-07
 LARGEST RESIDUAL FORCE            -1.619E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.500       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.313E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.077E-05   TIME AVG. MOMENT      5.829E-06
 LARGEST RESIDUAL MOMENT            6.322E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.968E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.435E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      8.973E-07   TIME AVG. FORCE       5.787E-07
 LARGEST RESIDUAL FORCE            -3.236E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.498       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.811E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.077E-05   TIME AVG. MOMENT      5.829E-06
 LARGEST RESIDUAL MOMENT            3.567E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.958E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.026E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      9.058E-07   TIME AVG. FORCE       5.794E-07
 LARGEST RESIDUAL FORCE            -9.466E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.495       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.428E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.077E-05   TIME AVG. MOMENT      5.829E-06
 LARGEST RESIDUAL MOMENT           -1.351E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.940E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.866E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      8.966E-07   TIME AVG. FORCE       5.787E-07
 LARGEST RESIDUAL FORCE            -1.082E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.495       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.115E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.077E-05   TIME AVG. MOMENT      5.829E-06
 LARGEST RESIDUAL MOMENT            1.395E-12   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.939E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.014E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.580E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.731E-02,  FRACTION OF STEP COMPLETED  0.353    
 STEP TIME COMPLETED       0.706    ,  TOTAL TIME COMPLETED        0.706    


  INCREMENT    23 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.580E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.429E-04   TIME AVG. FORCE       1.674E-05
 LARGEST RESIDUAL FORCE            -9.124E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.853       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.374       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.055E-05   TIME AVG. MOMENT      7.339E-06
 LARGEST RESIDUAL MOMENT           -5.810E-04   AT NODE        159   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.808E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.930E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.536E-06   TIME AVG. FORCE       6.425E-07
 LARGEST RESIDUAL FORCE            -1.136E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.848       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.129E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.216E-05   TIME AVG. MOMENT      6.104E-06
 LARGEST RESIDUAL MOMENT            6.105E-06   AT NODE        130   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.785E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.652E-04   AT NODE         48   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.626E-04   TIME AVG. FORCE       1.805E-05
 LARGEST RESIDUAL FORCE            -1.006E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.466       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.382       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.054E-05   TIME AVG. MOMENT      6.469E-06
 LARGEST RESIDUAL MOMENT           -2.978E-04   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.819E-02   AT NODE        181   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.966E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      9.556E-07   TIME AVG. FORCE       6.038E-07
 LARGEST RESIDUAL FORCE             2.675E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.465       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.778E-03   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.140E-05   TIME AVG. MOMENT      6.071E-06
 LARGEST RESIDUAL MOMENT            2.981E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.817E-02   AT NODE        181   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.555E-05   AT NODE        132   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.021E-06   TIME AVG. FORCE       6.082E-07
 LARGEST RESIDUAL FORCE            -4.294E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.472       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.513E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.141E-05   TIME AVG. MOMENT      6.072E-06
 LARGEST RESIDUAL MOMENT            5.277E-09   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.857E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.986E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -7.057E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      9.679E-07   TIME AVG. FORCE       6.046E-07
 LARGEST RESIDUAL FORCE            -3.909E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.471       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.112E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.142E-05   TIME AVG. MOMENT      6.072E-06
 LARGEST RESIDUAL MOMENT            1.101E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.856E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.184E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      9.678E-07   TIME AVG. FORCE       6.046E-07
 LARGEST RESIDUAL FORCE            -7.222E-12   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.471       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.856E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.142E-05   TIME AVG. MOMENT      6.072E-06
 LARGEST RESIDUAL MOMENT           -4.426E-13   AT NODE        185   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.856E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.473E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.440E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.580E-02,  FRACTION OF STEP COMPLETED  0.376    
 STEP TIME COMPLETED       0.752    ,  TOTAL TIME COMPLETED        0.752    


  INCREMENT    24 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.440E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.523E-04   TIME AVG. FORCE       1.009E-05
 LARGEST RESIDUAL FORCE            -5.250E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.743       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.287       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.401E-05   TIME AVG. MOMENT      7.236E-06
 LARGEST RESIDUAL MOMENT           -4.660E-04   AT NODE        171   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.295E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.496E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.239E-06   TIME AVG. FORCE       6.443E-07
 LARGEST RESIDUAL FORCE            -7.827E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.738       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.866E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.259E-05   TIME AVG. MOMENT      6.343E-06
 LARGEST RESIDUAL MOMENT            3.621E-06   AT NODE        130   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -3.270E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.420E-04   AT NODE        136   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.495E-04   TIME AVG. FORCE       9.910E-06
 LARGEST RESIDUAL FORCE            -5.752E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.452       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.286       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.571E-05   TIME AVG. MOMENT      6.474E-06
 LARGEST RESIDUAL MOMENT           -1.285E-04   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.784E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.486E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.029E-06   TIME AVG. FORCE       6.312E-07
 LARGEST RESIDUAL FORCE             8.721E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.451       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.553E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.205E-05   TIME AVG. MOMENT      6.321E-06
 LARGEST RESIDUAL MOMENT            7.954E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.783E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.462E-06   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.040E-06   TIME AVG. FORCE       6.318E-07
 LARGEST RESIDUAL FORCE            -1.013E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.450       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.628E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.205E-05   TIME AVG. MOMENT      6.321E-06
 LARGEST RESIDUAL MOMENT           -1.439E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.777E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.900E-05   AT NODE        125   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.038E-06   TIME AVG. FORCE       6.317E-07
 LARGEST RESIDUAL FORCE            -1.405E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.450       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.314E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.205E-05   TIME AVG. MOMENT      6.321E-06
 LARGEST RESIDUAL MOMENT           -5.994E-13   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.777E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.260E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.307E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.440E-02,  FRACTION OF STEP COMPLETED  0.398    
 STEP TIME COMPLETED       0.796    ,  TOTAL TIME COMPLETED        0.796    


  INCREMENT    25 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.307E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.702E-05   TIME AVG. FORCE       6.302E-06
 LARGEST RESIDUAL FORCE            -3.048E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.657       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.221       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.876E-05   TIME AVG. MOMENT      7.219E-06
 LARGEST RESIDUAL MOMENT           -3.603E-04   AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.892E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.168E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.181E-06   TIME AVG. FORCE       6.641E-07
 LARGEST RESIDUAL FORCE            -5.528E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.650       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.956E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.304E-05   TIME AVG. MOMENT      6.590E-06
 LARGEST RESIDUAL MOMENT            2.177E-06   AT NODE        138   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.863E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.099E-04   AT NODE        136   DOF  5
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.735E-05   TIME AVG. FORCE       5.733E-06
 LARGEST RESIDUAL FORCE            -3.376E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.434       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.216       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.417E-05   TIME AVG. MOMENT      6.635E-06
 LARGEST RESIDUAL MOMENT           -5.777E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.725E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.138E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.102E-06   TIME AVG. FORCE       6.594E-07
 LARGEST RESIDUAL FORCE            -2.855E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.434       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.586E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.265E-05   TIME AVG. MOMENT      6.574E-06
 LARGEST RESIDUAL MOMENT            4.833E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.725E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.087E-05   AT NODE         93   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.122E-06   TIME AVG. FORCE       6.606E-07
 LARGEST RESIDUAL FORCE            -1.277E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.429       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.376E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.265E-05   TIME AVG. MOMENT      6.574E-06
 LARGEST RESIDUAL MOMENT           -5.689E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.703E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.189E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.576E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.107E-06   TIME AVG. FORCE       6.597E-07
 LARGEST RESIDUAL FORCE             1.326E-11   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.429       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.336E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.265E-05   TIME AVG. MOMENT      6.574E-06
 LARGEST RESIDUAL MOMENT           -2.374E-12   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.703E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.678E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.182E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.307E-02,  FRACTION OF STEP COMPLETED  0.420    
 STEP TIME COMPLETED       0.840    ,  TOTAL TIME COMPLETED        0.840    


  INCREMENT    26 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.182E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.325E-05   TIME AVG. FORCE       4.137E-06
 LARGEST RESIDUAL FORCE             1.899E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.588       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.171       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.491E-05   TIME AVG. MOMENT      7.279E-06
 LARGEST RESIDUAL MOMENT           -2.733E-04   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.573E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.192E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.236E-06   TIME AVG. FORCE       6.917E-07
 LARGEST RESIDUAL FORCE            -3.808E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.580       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.349E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.355E-05   TIME AVG. MOMENT      6.842E-06
 LARGEST RESIDUAL MOMENT            1.418E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.539E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.863E-04   AT NODE        101   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.139E-05   TIME AVG. FORCE       3.478E-06
 LARGEST RESIDUAL FORCE            -1.989E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.416       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.163       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.386E-05   TIME AVG. MOMENT      6.854E-06
 LARGEST RESIDUAL MOMENT           -2.617E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.665E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.741E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.171E-06   TIME AVG. FORCE       6.880E-07
 LARGEST RESIDUAL FORCE            -1.979E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.416       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.619E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.322E-05   TIME AVG. MOMENT      6.830E-06
 LARGEST RESIDUAL MOMENT            4.792E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.664E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.707E-05   AT NODE         97   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.209E-06   TIME AVG. FORCE       6.902E-07
 LARGEST RESIDUAL FORCE            -2.462E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.410       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.823E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.323E-05   TIME AVG. MOMENT      6.830E-06
 LARGEST RESIDUAL MOMENT           -1.175E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.633E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.048E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -7.475E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.174E-06   TIME AVG. FORCE       6.882E-07
 LARGEST RESIDUAL FORCE             8.160E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.410       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.508E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.323E-05   TIME AVG. MOMENT      6.830E-06
 LARGEST RESIDUAL MOMENT           -2.003E-12   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.633E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.369E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.064E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.182E-02,  FRACTION OF STEP COMPLETED  0.441    
 STEP TIME COMPLETED       0.881    ,  TOTAL TIME COMPLETED        0.881    


  INCREMENT    27 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.064E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.165E-05   TIME AVG. FORCE       2.844E-06
 LARGEST RESIDUAL FORCE             1.333E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.532       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.134       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.222E-05   TIME AVG. MOMENT      7.400E-06
 LARGEST RESIDUAL MOMENT           -2.054E-04   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.315E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.275E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.332E-06   TIME AVG. FORCE       7.221E-07
 LARGEST RESIDUAL FORCE            -4.846E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.522       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.997E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.401E-05   TIME AVG. MOMENT      7.096E-06
 LARGEST RESIDUAL MOMENT            1.100E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.274E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.416E-04   AT NODE        105   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.003E-05   TIME AVG. FORCE       2.233E-06
 LARGEST RESIDUAL FORCE            -1.157E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.399       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.123       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.404E-05   TIME AVG. MOMENT      7.097E-06
 LARGEST RESIDUAL MOMENT           -1.162E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.606E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.674E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.236E-06   TIME AVG. FORCE       7.171E-07
 LARGEST RESIDUAL FORCE            -1.972E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.399       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.552E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.378E-05   TIME AVG. MOMENT      7.087E-06
 LARGEST RESIDUAL MOMENT            5.767E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.604E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.817E-05   AT NODE        109   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.291E-06   TIME AVG. FORCE       7.199E-07
 LARGEST RESIDUAL FORCE            -3.523E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.392       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.769E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.378E-05   TIME AVG. MOMENT      7.087E-06
 LARGEST RESIDUAL MOMENT           -1.921E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.568E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.658E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.218E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.240E-06   TIME AVG. FORCE       7.172E-07
 LARGEST RESIDUAL FORCE            -3.265E-12   AT NODE        178   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.392       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.212E-06   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.378E-05   TIME AVG. MOMENT      7.087E-06
 LARGEST RESIDUAL MOMENT           -7.109E-13   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.568E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.683E-08   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   3.953E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.064E-02,  FRACTION OF STEP COMPLETED  0.461    
 STEP TIME COMPLETED       0.922    ,  TOTAL TIME COMPLETED        0.922    


  INCREMENT    28 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.953E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.807E-05   TIME AVG. FORCE       2.085E-06
 LARGEST RESIDUAL FORCE             9.410E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.486       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.105       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.041E-05   TIME AVG. MOMENT      7.563E-06
 LARGEST RESIDUAL MOMENT           -1.538E-04   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.104E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.790E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.460E-06   TIME AVG. FORCE       7.544E-07
 LARGEST RESIDUAL FORCE            -6.841E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.474       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.184E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.449E-05   TIME AVG. MOMENT      7.352E-06
 LARGEST RESIDUAL MOMENT            8.901E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.054E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.124E-04   AT NODE        117   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.732E-05   TIME AVG. FORCE       1.547E-06
 LARGEST RESIDUAL FORCE            -6.551E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.383       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.099E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.441E-05   TIME AVG. MOMENT      7.349E-06
 LARGEST RESIDUAL MOMENT           -4.939E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.551E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.030E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.301E-06   TIME AVG. FORCE       7.464E-07
 LARGEST RESIDUAL FORCE            -2.042E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.382       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.182E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.431E-05   TIME AVG. MOMENT      7.345E-06
 LARGEST RESIDUAL MOMENT            6.860E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.546E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.854E-05   AT NODE        144   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.367E-06   TIME AVG. FORCE       7.497E-07
 LARGEST RESIDUAL FORCE            -4.189E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.375       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.199E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.431E-05   TIME AVG. MOMENT      7.345E-06
 LARGEST RESIDUAL MOMENT           -2.445E-09   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.506E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.001E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.426E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.303E-06   TIME AVG. FORCE       7.466E-07
 LARGEST RESIDUAL FORCE            -7.476E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.375       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.642E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.431E-05   TIME AVG. MOMENT      7.345E-06
 LARGEST RESIDUAL MOMENT            2.823E-12   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.506E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.741E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   3.847E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.953E-02,  FRACTION OF STEP COMPLETED  0.481    
 STEP TIME COMPLETED       0.962    ,  TOTAL TIME COMPLETED        0.962    


  INCREMENT    29 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.847E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.927E-05   TIME AVG. FORCE       1.628E-06
 LARGEST RESIDUAL FORCE             6.684E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.448       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.271E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.922E-05   TIME AVG. MOMENT      7.755E-06
 LARGEST RESIDUAL MOMENT           -1.152E-04   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.929E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.629E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.614E-06   TIME AVG. FORCE       7.879E-07
 LARGEST RESIDUAL FORCE            -9.924E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.434       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.380E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.496E-05   TIME AVG. MOMENT      7.608E-06
 LARGEST RESIDUAL MOMENT            7.140E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.867E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.168E-04   AT NODE        164   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      9.878E-06   TIME AVG. FORCE       1.181E-06
 LARGEST RESIDUAL FORCE            -3.558E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.368       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.605E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.486E-05   TIME AVG. MOMENT      7.604E-06
 LARGEST RESIDUAL MOMENT           -1.942E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.496E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.714E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.366E-06   TIME AVG. FORCE       7.761E-07
 LARGEST RESIDUAL FORCE            -3.098E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.366       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.915E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.482E-05   TIME AVG. MOMENT      7.603E-06
 LARGEST RESIDUAL MOMENT            7.282E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.487E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.019E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.424E-06   TIME AVG. FORCE       7.788E-07
 LARGEST RESIDUAL FORCE            -3.965E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.359       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.859E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.482E-05   TIME AVG. MOMENT      7.603E-06
 LARGEST RESIDUAL MOMENT           -2.126E-09   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.448E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.905E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.140E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.367E-06   TIME AVG. FORCE       7.761E-07
 LARGEST RESIDUAL FORCE            -3.932E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.359       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.843E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.482E-05   TIME AVG. MOMENT      7.603E-06
 LARGEST RESIDUAL MOMENT            9.150E-12   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.448E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.007E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.847E-02,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         1.00    


  INCREMENT    30 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.891E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.854E-05   TIME AVG. FORCE       2.493E-06
 LARGEST RESIDUAL FORCE            -3.709E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.316       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.141       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.116E-05   TIME AVG. MOMENT      8.055E-06
 LARGEST RESIDUAL MOMENT            1.439E-04   AT NODE        114   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.570E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.995E-03   AT NODE         73   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.809E-06   TIME AVG. FORCE       8.685E-07
 LARGEST RESIDUAL FORCE            -9.106E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.292       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.379E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.528E-05   TIME AVG. MOMENT      7.859E-06
 LARGEST RESIDUAL MOMENT            3.763E-06   AT NODE        115   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.378E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.939E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.036E-05   TIME AVG. FORCE       1.666E-06
 LARGEST RESIDUAL FORCE            -5.603E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.187       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.105       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.524E-05   TIME AVG. MOMENT      7.857E-06
 LARGEST RESIDUAL MOMENT           -3.639E-06   AT NODE        134   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.112E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.701E-03   AT NODE        164   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.408E-06   TIME AVG. FORCE       8.048E-07
 LARGEST RESIDUAL FORCE            -8.932E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.180       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.177E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.511E-05   TIME AVG. MOMENT      7.853E-06
 LARGEST RESIDUAL MOMENT            7.084E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.551E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.657E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.594E-06   TIME AVG. FORCE       1.177E-06
 LARGEST RESIDUAL FORCE            -3.588E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.113       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.703E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.501E-05   TIME AVG. MOMENT      7.850E-06
 LARGEST RESIDUAL MOMENT           -1.976E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.817E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.735E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.318E-06   TIME AVG. FORCE       8.007E-07
 LARGEST RESIDUAL FORCE            -1.437E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.112       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.687E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.497E-05   TIME AVG. MOMENT      7.848E-06
 LARGEST RESIDUAL MOMENT            4.574E-08   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.761E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.622E-05   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      1.332E-06   TIME AVG. FORCE       8.014E-07
 LARGEST RESIDUAL FORCE            -1.366E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.108       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.031E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.498E-05   TIME AVG. MOMENT      7.849E-06
 LARGEST RESIDUAL MOMENT           -4.246E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.531E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.295E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.130E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      1.319E-06   TIME AVG. FORCE       8.008E-07
 LARGEST RESIDUAL FORCE            -2.638E-12   AT NODE        183   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.108       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.289E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.498E-05   TIME AVG. MOMENT      7.849E-06
 LARGEST RESIDUAL MOMENT            8.405E-13   AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.531E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.807E-07   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.891E-02,  FRACTION OF STEP COMPLETED  0.524    
 STEP TIME COMPLETED        1.05    ,  TOTAL TIME COMPLETED         1.05    


  INCREMENT    31 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.891E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.395E-05   TIME AVG. FORCE       2.677E-06
 LARGEST RESIDUAL FORCE            -1.183E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         5.198E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.160       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.490E-05   TIME AVG. MOMENT      8.076E-06
 LARGEST RESIDUAL MOMENT            4.362E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.188E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.719E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.373E-06   TIME AVG. FORCE       8.257E-07
 LARGEST RESIDUAL FORCE             3.813E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         5.217E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.499E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.489E-05   TIME AVG. MOMENT      8.076E-06
 LARGEST RESIDUAL MOMENT            2.516E-09   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.187E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.607E-06   AT NODE        132   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.376E-06   TIME AVG. FORCE       8.258E-07
 LARGEST RESIDUAL FORCE            -6.224E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         5.245E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.785E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.490E-05   TIME AVG. MOMENT      8.076E-06
 LARGEST RESIDUAL MOMENT            1.083E-11   AT NODE         95   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.202E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.552E-05   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.376E-06   TIME AVG. FORCE       8.258E-07
 LARGEST RESIDUAL FORCE            -2.112E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         5.245E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.635E-07   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.490E-05   TIME AVG. MOMENT      8.076E-06
 LARGEST RESIDUAL MOMENT            5.013E-13   AT NODE        187   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.202E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.007E-09   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.891E-02,  FRACTION OF STEP COMPLETED  0.549    
 STEP TIME COMPLETED        1.10    ,  TOTAL TIME COMPLETED         1.10    


  INCREMENT    32 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.891E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.392E-06   TIME AVG. FORCE       8.494E-07
 LARGEST RESIDUAL FORCE             8.253E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         4.779E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.660E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.481E-05   TIME AVG. MOMENT      8.287E-06
 LARGEST RESIDUAL MOMENT           -1.090E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.012E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.906E-04   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.368E-06   TIME AVG. FORCE       8.484E-07
 LARGEST RESIDUAL FORCE            -4.754E-12   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         4.779E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.513E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.483E-05   TIME AVG. MOMENT      8.287E-06
 LARGEST RESIDUAL MOMENT            1.105E-12   AT NODE        139   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.012E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.803E-07   AT NODE        161   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    7.336E-02

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.891E-02,  FRACTION OF STEP COMPLETED  0.573    
 STEP TIME COMPLETED        1.15    ,  TOTAL TIME COMPLETED         1.15    


  INCREMENT    33 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  7.336E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.883E-06   TIME AVG. FORCE       1.130E-06
 LARGEST RESIDUAL FORCE            -1.862E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         7.735E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.396E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.485E-05   TIME AVG. MOMENT      8.486E-06
 LARGEST RESIDUAL MOMENT           -1.645E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      3.268E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.692E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.367E-06   TIME AVG. FORCE       8.691E-07
 LARGEST RESIDUAL FORCE             1.152E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         7.757E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.477E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.482E-05   TIME AVG. MOMENT      8.485E-06
 LARGEST RESIDUAL MOMENT           -1.090E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      3.265E-04   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.994E-07   AT NODE        129   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.110    

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.336E-02,  FRACTION OF STEP COMPLETED  0.610    
 STEP TIME COMPLETED        1.22    ,  TOTAL TIME COMPLETED         1.22    


  INCREMENT    34 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.110    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.791E-06   TIME AVG. FORCE       9.046E-07
 LARGEST RESIDUAL FORCE            -1.602E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -6.991E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.863E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.483E-05   TIME AVG. MOMENT      8.672E-06
 LARGEST RESIDUAL MOMENT            2.509E-09   AT NODE        107   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.942E-04   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.839E-04   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
          ROTATION CORRECTION TOO LARGE COMPARED TO ROTATION INCREMENT

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.368E-06   TIME AVG. FORCE       8.883E-07
 LARGEST RESIDUAL FORCE             7.844E-12   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -6.989E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.256E-06   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.483E-05   TIME AVG. MOMENT      8.672E-06
 LARGEST RESIDUAL MOMENT           -3.802E-12   AT NODE         81   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.942E-04   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.460E-08   AT NODE        144   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.165    

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.110    ,  FRACTION OF STEP COMPLETED  0.665    
 STEP TIME COMPLETED        1.33    ,  TOTAL TIME COMPLETED         1.33    


  INCREMENT    35 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.165    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.447E-06   TIME AVG. FORCE       9.090E-07
 LARGEST RESIDUAL FORCE            -3.861E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -1.351E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.133E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.482E-05   TIME AVG. MOMENT      8.847E-06
 LARGEST RESIDUAL MOMENT            5.276E-09   AT NODE        163   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.684E-05   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.844E-04   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
          ROTATION CORRECTION TOO LARGE COMPARED TO ROTATION INCREMENT

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.369E-06   TIME AVG. FORCE       9.061E-07
 LARGEST RESIDUAL FORCE            -2.863E-12   AT NODE        186   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -1.351E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.547E-06   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.483E-05   TIME AVG. MOMENT      8.848E-06
 LARGEST RESIDUAL MOMENT           -9.216E-13   AT NODE         93   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -5.684E-05   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.759E-09   AT NODE         65   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.248    

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.165    ,  FRACTION OF STEP COMPLETED  0.748    
 STEP TIME COMPLETED        1.50    ,  TOTAL TIME COMPLETED         1.50    


  INCREMENT    36 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.248    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.372E-06   TIME AVG. FORCE       9.227E-07
 LARGEST RESIDUAL FORCE            -2.283E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         1.974E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.224E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.483E-05   TIME AVG. MOMENT      9.014E-06
 LARGEST RESIDUAL MOMENT            5.561E-11   AT NODE        123   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      8.315E-06   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.358E-05   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
          ROTATION CORRECTION TOO LARGE COMPARED TO ROTATION INCREMENT

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.369E-06   TIME AVG. FORCE       9.226E-07
 LARGEST RESIDUAL FORCE            -2.783E-12   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         1.975E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.212E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.483E-05   TIME AVG. MOMENT      9.014E-06
 LARGEST RESIDUAL MOMENT           -3.570E-13   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      8.314E-06   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.987E-10   AT NODE        165   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.257    

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.248    ,  FRACTION OF STEP COMPLETED  0.871    
 STEP TIME COMPLETED        1.74    ,  TOTAL TIME COMPLETED         1.74    


  INCREMENT    37 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.257    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.369E-06   TIME AVG. FORCE       9.380E-07
 LARGEST RESIDUAL FORCE            -1.454E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         2.765E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.774E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE MOMENT                     1.483E-05   TIME AVG. MOMENT      9.171E-06
 LARGEST RESIDUAL MOMENT           -3.846E-13   AT NODE        186   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      1.164E-06   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.472E-06   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
          ROTATION CORRECTION TOO LARGE COMPARED TO ROTATION INCREMENT

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.369E-06   TIME AVG. FORCE       9.380E-07
 LARGEST RESIDUAL FORCE             4.250E-12   AT NODE        187   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         2.765E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.038E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.483E-05   TIME AVG. MOMENT      9.171E-06
 LARGEST RESIDUAL MOMENT            5.600E-13   AT NODE        186   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      1.164E-06   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.973E-12   AT NODE         65   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.257    ,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        2.00    ,  TOTAL TIME COMPLETED         2.00    


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF         37  INCREMENTS
                       1  CUTBACKS IN AUTOMATIC INCREMENTATION
                     186  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     186  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     186  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       2  WARNING MESSAGES DURING USER INPUT PROCESSING
                       9  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       2  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   5.0200    
       SYSTEM TIME (SEC)    =  0.29000    
       TOTAL CPU TIME (SEC) =   5.3100    
       WALLCLOCK TIME (SEC) =          6

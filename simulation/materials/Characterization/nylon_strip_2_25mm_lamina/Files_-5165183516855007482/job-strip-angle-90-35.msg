1

   Abaqus 2021                                  Date 14-Feb-2023   Time 14:20:57
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


 AVERAGE FORCE                      1.316E-04   TIME AVG. FORCE       1.316E-04
 LARGEST RESIDUAL FORCE            -1.427E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.431       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.431       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.562E-06   TIME AVG. MOMENT      4.562E-06
 LARGEST RESIDUAL MOMENT            1.734E-04   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.639E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.639E-02   AT NODE        191   DOF  4
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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.995E-08   TIME AVG. FORCE       1.176E-02
                                                (FROM MOMENT    FIELD)
 LARGEST RESIDUAL FORCE             2.878E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.431       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.034E-03   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     5.336E-07   TIME AVG. MOMENT      5.336E-07
 LARGEST RESIDUAL MOMENT            6.890E-08   AT NODE        187   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.639E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.429E-06   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.895E-08   TIME AVG. FORCE       1.176E-02
                                                (FROM MOMENT    FIELD)
 LARGEST RESIDUAL FORCE            -2.423E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.431       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.442E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE

 AVERAGE MOMENT                     5.327E-07   TIME AVG. MOMENT      5.327E-07
 LARGEST RESIDUAL MOMENT           -2.198E-11   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.639E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.108E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   7.039E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.387E-05   TIME AVG. FORCE       2.387E-05
 LARGEST RESIDUAL FORCE            -2.817E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.577       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.167       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.106E-05   TIME AVG. MOMENT      5.798E-06
 LARGEST RESIDUAL MOMENT           -2.015E-04   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.295E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.362E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      6.359E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.628E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.572       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.318E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.552E-06   TIME AVG. MOMENT      1.542E-06
 LARGEST RESIDUAL MOMENT            5.522E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.276E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.850E-04   AT NODE        148   DOF  4
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


 AVERAGE FORCE                      7.696E-06   TIME AVG. FORCE       7.696E-06
 LARGEST RESIDUAL FORCE            -1.146E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.474       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.836E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.337E-06   TIME AVG. MOMENT      1.435E-06
 LARGEST RESIDUAL MOMENT           -3.941E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.811E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.654E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      4.219E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -2.380E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.474       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.714E-04   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     2.264E-06   TIME AVG. MOMENT      1.398E-06
 LARGEST RESIDUAL MOMENT            7.828E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.810E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.595E-06   AT NODE        109   DOF  5
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


 AVERAGE FORCE                      4.259E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.947E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.473       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.268E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     2.260E-06   TIME AVG. MOMENT      1.396E-06
 LARGEST RESIDUAL MOMENT           -1.008E-11   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.804E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.063E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.223E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.227E-12   AT NODE        175   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.473       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.363E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     2.260E-06   TIME AVG. MOMENT      1.396E-06
 LARGEST RESIDUAL MOMENT           -3.335E-14   AT NODE        188   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.804E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.912E-09   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   6.711E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  7.039E-02,  FRACTION OF STEP COMPLETED  7.219E-02
 STEP TIME COMPLETED       0.144    ,  TOTAL TIME COMPLETED        0.144    


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.711E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.086E-05   TIME AVG. FORCE       2.086E-05
 LARGEST RESIDUAL FORCE            -2.606E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.598       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.148       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.335E-05   TIME AVG. MOMENT      5.382E-06
 LARGEST RESIDUAL MOMENT           -2.414E-04   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.434E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.137E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.256E-07   TIME AVG. FORCE       1.256E-07
 LARGEST RESIDUAL FORCE            -4.795E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.590       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.290E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.811E-06   TIME AVG. MOMENT      2.201E-06
 LARGEST RESIDUAL MOMENT            9.073E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.402E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.403E-04   AT NODE        132   DOF  4
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


 AVERAGE FORCE                      1.846E-05   TIME AVG. FORCE       1.846E-05
 LARGEST RESIDUAL FORCE            -2.744E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.439       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.151       AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.698E-06   TIME AVG. MOMENT      2.164E-06
 LARGEST RESIDUAL MOMENT           -1.461E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.681E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.207E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      6.518E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.330E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.439       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.200E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          ALL DISP.    CORRECTIONS ARE ZERO

 AVERAGE MOMENT                     3.357E-06   TIME AVG. MOMENT      2.050E-06
 LARGEST RESIDUAL MOMENT            4.376E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.680E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.175E-05   AT NODE        113   DOF  5
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


 AVERAGE FORCE                      1.054E-07   TIME AVG. FORCE       1.054E-07
 LARGEST RESIDUAL FORCE            -6.863E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.431       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.609E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.329E-06   TIME AVG. MOMENT      2.041E-06
 LARGEST RESIDUAL MOMENT           -1.875E-09   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.644E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.601E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.543E-05
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


 AVERAGE FORCE                      6.539E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE             5.135E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.431       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.175E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     3.331E-06   TIME AVG. MOMENT      2.041E-06
 LARGEST RESIDUAL MOMENT           -1.626E-12   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.644E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.160E-07   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   6.413E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.711E-02,  FRACTION OF STEP COMPLETED  0.106    
 STEP TIME COMPLETED       0.212    ,  TOTAL TIME COMPLETED        0.212    


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.413E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.293E-06   TIME AVG. FORCE       9.293E-06
 LARGEST RESIDUAL FORCE            -1.074E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.508       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.589E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.631E-06   TIME AVG. MOMENT      3.939E-06
 LARGEST RESIDUAL MOMENT           -1.336E-04   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.032E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.614E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.742E-07   TIME AVG. FORCE       1.742E-07
 LARGEST RESIDUAL FORCE            -1.069E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.496       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.118E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.617E-06   TIME AVG. MOMENT      2.685E-06
 LARGEST RESIDUAL MOMENT            5.521E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.987E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.588E-04   AT NODE        152   DOF  4
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


 AVERAGE FORCE                      5.921E-06   TIME AVG. FORCE       5.921E-06
 LARGEST RESIDUAL FORCE            -8.730E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.412       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.444E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.396E-06   TIME AVG. MOMENT      2.630E-06
 LARGEST RESIDUAL MOMENT           -2.589E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.580E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.066E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      8.847E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -1.414E-07   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.411       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.176E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     IS ZERO EVERYWHERE BUT THE FORCE     RESIDUAL OR
          THE DISP.    CORRECTION IS NON-ZERO

 AVERAGE MOMENT                     4.361E-06   TIME AVG. MOMENT      2.621E-06
 LARGEST RESIDUAL MOMENT            4.803E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.575E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.007E-05   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.295E-07   TIME AVG. FORCE       1.295E-07
 LARGEST RESIDUAL FORCE            -7.396E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.403       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.748E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     4.330E-06   TIME AVG. MOMENT      2.614E-06
 LARGEST RESIDUAL MOMENT           -2.002E-09   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.538E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.739E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.559E-05
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


 AVERAGE FORCE                      8.848E-08   TIME AVG. FORCE       1.000E-02
 LARGEST RESIDUAL FORCE            -9.699E-12   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.403       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.397E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THERE IS ZERO FORCE     EVERYWHERE
 

 ***WARNING: THERE IS ZERO FORCE EVERYWHERE IN THE MODEL BASED ON THE DEFAULT 
             CRITERION. PLEASE CHECK THE VALUE OF THE AVERAGE FORCE DURING THE 
             CURRENT ITERATION TO VERIFY THAT THE FORCE IS SMALL ENOUGH TO BE 
             TREATED AS ZERO. IF NOT, PLEASE USE THE SOLUTION CONTROLS TO RESET 
             THE CRITERION FOR ZERO FORCE.

 AVERAGE MOMENT                     4.331E-06   TIME AVG. MOMENT      2.614E-06
 LARGEST RESIDUAL MOMENT            3.255E-12   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.538E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.963E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   6.140E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.413E-02,  FRACTION OF STEP COMPLETED  0.138    
 STEP TIME COMPLETED       0.276    ,  TOTAL TIME COMPLETED        0.276    


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.140E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.684E-06   TIME AVG. FORCE       5.684E-06
 LARGEST RESIDUAL FORCE            -6.123E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.460       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.367E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.728E-06   TIME AVG. MOMENT      3.837E-06
 LARGEST RESIDUAL MOMENT           -9.005E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.823E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.509E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      2.202E-07   TIME AVG. FORCE       2.202E-07
 LARGEST RESIDUAL FORCE            -1.502E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.447       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.253E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.480E-06   TIME AVG. MOMENT      3.187E-06
 LARGEST RESIDUAL MOMENT            3.740E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.769E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.351E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      2.556E-06   TIME AVG. FORCE       2.556E-06
 LARGEST RESIDUAL FORCE            -3.674E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.393       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.456E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.310E-06   TIME AVG. MOMENT      3.153E-06
 LARGEST RESIDUAL MOMENT           -6.795E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.506E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.639E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.127E-07   TIME AVG. FORCE       1.127E-07
 LARGEST RESIDUAL FORCE            -4.320E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.391       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.957E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.307E-06   TIME AVG. MOMENT      3.152E-06
 LARGEST RESIDUAL MOMENT            3.669E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.497E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.998E-05   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.309E-07   TIME AVG. FORCE       1.309E-07
 LARGEST RESIDUAL FORCE            -3.896E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.385       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.575E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     5.287E-06   TIME AVG. MOMENT      3.148E-06
 LARGEST RESIDUAL MOMENT           -6.576E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.469E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.715E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -4.866E-06
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


 AVERAGE FORCE                      1.118E-07   TIME AVG. FORCE       1.118E-07
 LARGEST RESIDUAL FORCE            -4.748E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.385       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.948E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     5.287E-06   TIME AVG. MOMENT      3.148E-06
 LARGEST RESIDUAL MOMENT            4.086E-12   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.469E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.429E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   5.889E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.140E-02,  FRACTION OF STEP COMPLETED  0.169    
 STEP TIME COMPLETED       0.337    ,  TOTAL TIME COMPLETED        0.337    


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.889E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.111E-05   TIME AVG. FORCE       5.613E-06
 LARGEST RESIDUAL FORCE            -3.840E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.428       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.875E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.656E-06   TIME AVG. MOMENT      4.066E-06
 LARGEST RESIDUAL MOMENT           -6.631E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.690E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.806E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      8.184E-07   TIME AVG. FORCE       4.651E-07
 LARGEST RESIDUAL FORCE            -2.053E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.414       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.417E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.347E-06   TIME AVG. MOMENT      3.681E-06
 LARGEST RESIDUAL MOMENT            2.740E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.628E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.241E-04   AT NODE        185   DOF  4
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


 AVERAGE FORCE                      3.977E-06   TIME AVG. FORCE       2.045E-06
 LARGEST RESIDUAL FORCE            -1.812E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.376       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.815E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.224E-06   TIME AVG. MOMENT      3.661E-06
 LARGEST RESIDUAL MOMENT           -2.100E-07   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.442E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.855E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      4.163E-07   TIME AVG. FORCE       2.641E-07
 LARGEST RESIDUAL FORCE            -9.799E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.373       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.862E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.220E-06   TIME AVG. MOMENT      3.660E-06
 LARGEST RESIDUAL MOMENT            2.732E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.429E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.359E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      4.322E-07   TIME AVG. FORCE       2.720E-07
 LARGEST RESIDUAL FORCE            -2.040E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.369       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.014E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     6.204E-06   TIME AVG. MOMENT      3.658E-06
 LARGEST RESIDUAL MOMENT           -8.673E-11   AT NODE        179   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.409E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.966E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -6.240E-07
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


 AVERAGE FORCE                      4.074E-07   TIME AVG. FORCE       2.596E-07
 LARGEST RESIDUAL FORCE            -1.238E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.369       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.122E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     6.204E-06   TIME AVG. MOMENT      3.658E-06
 LARGEST RESIDUAL MOMENT            3.474E-12   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.409E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.526E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   5.658E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.889E-02,  FRACTION OF STEP COMPLETED  0.198    
 STEP TIME COMPLETED       0.396    ,  TOTAL TIME COMPLETED        0.396    


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.658E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.587E-06   TIME AVG. FORCE       2.702E-06
 LARGEST RESIDUAL FORCE            -2.390E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.401       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.662E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.833E-06   TIME AVG. MOMENT      4.397E-06
 LARGEST RESIDUAL MOMENT           -4.899E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.578E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.241E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.021E-06   TIME AVG. FORCE       5.135E-07
 LARGEST RESIDUAL FORCE            -2.686E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.385       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.581E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.186E-06   TIME AVG. MOMENT      4.162E-06
 LARGEST RESIDUAL MOMENT            1.922E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.506E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.144E-04   AT NODE        185   DOF  4
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


 AVERAGE FORCE                      2.100E-06   TIME AVG. FORCE       8.730E-07
 LARGEST RESIDUAL FORCE            -8.475E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.359       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.598E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.104E-06   TIME AVG. MOMENT      4.150E-06
 LARGEST RESIDUAL MOMENT           -4.382E-08   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.379E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.269E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      4.961E-07   TIME AVG. FORCE       3.384E-07
 LARGEST RESIDUAL FORCE            -1.587E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.356       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.585E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.096E-06   TIME AVG. MOMENT      4.149E-06
 LARGEST RESIDUAL MOMENT            1.617E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.362E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.733E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -6.394E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.836E-07   TIME AVG. FORCE       3.343E-07
 LARGEST RESIDUAL FORCE            -6.925E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.353       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.330E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.093E-06   TIME AVG. MOMENT      4.148E-06
 LARGEST RESIDUAL MOMENT            1.428E-10   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.351E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.146E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      4.780E-07   TIME AVG. FORCE       3.324E-07
 LARGEST RESIDUAL FORCE            -1.075E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.353       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.901E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.093E-06   TIME AVG. MOMENT      4.148E-06
 LARGEST RESIDUAL MOMENT            1.172E-12   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.350E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.426E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   5.445E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.658E-02,  FRACTION OF STEP COMPLETED  0.226    
 STEP TIME COMPLETED       0.453    ,  TOTAL TIME COMPLETED        0.453    


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.445E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.316E-06   TIME AVG. FORCE       1.578E-06
 LARGEST RESIDUAL FORCE             1.657E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.377       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.707E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.190E-06   TIME AVG. MOMENT      4.779E-06
 LARGEST RESIDUAL MOMENT           -3.622E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.479E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.795E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.177E-06   TIME AVG. FORCE       5.435E-07
 LARGEST RESIDUAL FORCE            -3.159E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.360       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.682E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.999E-06   TIME AVG. MOMENT      4.630E-06
 LARGEST RESIDUAL MOMENT            1.249E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.401E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.752E-04   AT NODE        185   DOF  4
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


 AVERAGE FORCE                      1.224E-06   TIME AVG. FORCE       5.553E-07
 LARGEST RESIDUAL FORCE            -3.712E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.343       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.712E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.950E-06   TIME AVG. MOMENT      4.624E-06
 LARGEST RESIDUAL MOMENT            6.339E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.317E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.405E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    4.266E-05
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


 AVERAGE FORCE                      5.668E-07   TIME AVG. FORCE       3.910E-07
 LARGEST RESIDUAL FORCE            -1.527E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.340       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.485E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.934E-06   TIME AVG. MOMENT      4.622E-06
 LARGEST RESIDUAL MOMENT            6.861E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.300E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.702E-04   AT NODE        192   DOF  4
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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.509E-07   TIME AVG. FORCE       3.870E-07
 LARGEST RESIDUAL FORCE            -1.246E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.339       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.845E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     7.934E-06   TIME AVG. MOMENT      4.622E-06
 LARGEST RESIDUAL MOMENT            7.831E-11   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.296E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.865E-05   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      5.500E-07   TIME AVG. FORCE       3.868E-07
 LARGEST RESIDUAL FORCE            -1.633E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.339       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.167E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     7.934E-06   TIME AVG. MOMENT      4.622E-06
 LARGEST RESIDUAL MOMENT            1.323E-13   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.295E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.767E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   5.247E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.445E-02,  FRACTION OF STEP COMPLETED  0.253    
 STEP TIME COMPLETED       0.507    ,  TOTAL TIME COMPLETED        0.507    


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.247E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.828E-06   TIME AVG. FORCE       1.075E-06
 LARGEST RESIDUAL FORCE             1.224E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.356       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.966E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.648E-06   TIME AVG. MOMENT      5.180E-06
 LARGEST RESIDUAL MOMENT           -2.696E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.393E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.446E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.268E-06   TIME AVG. FORCE       5.630E-07
 LARGEST RESIDUAL FORCE            -3.322E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.339       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.701E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.799E-06   TIME AVG. MOMENT      5.086E-06
 LARGEST RESIDUAL MOMENT            7.520E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.313E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.954E-04   AT NODE        189   DOF  4
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


 AVERAGE FORCE                      8.836E-07   TIME AVG. FORCE       4.862E-07
 LARGEST RESIDUAL FORCE            -1.558E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.328       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.105E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.759E-06   TIME AVG. MOMENT      5.081E-06
 LARGEST RESIDUAL MOMENT            1.387E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.259E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.449E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    1.005E-04
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


 AVERAGE FORCE                      6.305E-07   TIME AVG. FORCE       4.355E-07
 LARGEST RESIDUAL FORCE            -9.002E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.325       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.660E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.748E-06   TIME AVG. MOMENT      5.080E-06
 LARGEST RESIDUAL MOMENT            2.190E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.246E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.308E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    2.065E-05
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.219E-07   TIME AVG. FORCE       4.338E-07
 LARGEST RESIDUAL FORCE            -1.275E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.325       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.139E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     8.747E-06   TIME AVG. MOMENT      5.080E-06
 LARGEST RESIDUAL MOMENT            1.554E-11   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.244E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.557E-05   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      6.219E-07   TIME AVG. FORCE       4.338E-07
 LARGEST RESIDUAL FORCE            -1.204E-12   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.325       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.170E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     8.747E-06   TIME AVG. MOMENT      5.080E-06
 LARGEST RESIDUAL MOMENT            7.925E-14   AT NODE        188   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.244E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.077E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   5.063E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.247E-02,  FRACTION OF STEP COMPLETED  0.280    
 STEP TIME COMPLETED       0.559    ,  TOTAL TIME COMPLETED        0.559    


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.063E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  2.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.872E-06   TIME AVG. FORCE       8.401E-07
 LARGEST RESIDUAL FORCE             9.108E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.337       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.384E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.019E-05   TIME AVG. MOMENT      5.590E-06
 LARGEST RESIDUAL MOMENT           -2.021E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.318E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.171E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.300E-06   TIME AVG. FORCE       5.781E-07
 LARGEST RESIDUAL FORCE            -3.169E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.321       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.645E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.558E-06   TIME AVG. MOMENT      5.528E-06
 LARGEST RESIDUAL MOMENT            4.213E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.240E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.774E-04   AT NODE        189   DOF  4
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


 AVERAGE FORCE                      7.887E-07   TIME AVG. FORCE       4.930E-07
 LARGEST RESIDUAL FORCE            -6.308E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.314       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.005E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.535E-06   TIME AVG. MOMENT      5.525E-06
 LARGEST RESIDUAL MOMENT            1.066E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.205E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.470E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    8.783E-05
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


 AVERAGE FORCE                      6.958E-07   TIME AVG. FORCE       4.775E-07
 LARGEST RESIDUAL FORCE            -3.657E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.312       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.687E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.533E-06   TIME AVG. MOMENT      5.525E-06
 LARGEST RESIDUAL MOMENT            5.663E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.197E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.346E-05   AT NODE        192   DOF  4
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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.932E-07   TIME AVG. FORCE       4.771E-07
 LARGEST RESIDUAL FORCE            -8.371E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.312       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.016E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     9.533E-06   TIME AVG. MOMENT      5.525E-06
 LARGEST RESIDUAL MOMENT            1.670E-12   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.197E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.994E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.891E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.063E-02,  FRACTION OF STEP COMPLETED  0.305    
 STEP TIME COMPLETED       0.610    ,  TOTAL TIME COMPLETED        0.610    


  INCREMENT    11 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.891E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.256E-06   TIME AVG. FORCE       7.313E-07
 LARGEST RESIDUAL FORCE             6.813E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.321       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.922E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.076E-05   TIME AVG. MOMENT      6.001E-06
 LARGEST RESIDUAL MOMENT           -1.524E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.251E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.528E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.285E-06   TIME AVG. FORCE       5.925E-07
 LARGEST RESIDUAL FORCE            -2.786E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.305       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.530E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.030E-05   TIME AVG. MOMENT      5.959E-06
 LARGEST RESIDUAL MOMENT            2.213E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.178E-02   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.293E-04   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.059E-06
 MOMENT    EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      7.965E-07   TIME AVG. FORCE       5.227E-07
 LARGEST RESIDUAL FORCE            -2.461E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.301       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.359E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.029E-05   TIME AVG. MOMENT      5.958E-06
 LARGEST RESIDUAL MOMENT            6.370E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.157E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.169E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    6.244E-05
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


 AVERAGE FORCE                      7.656E-07   TIME AVG. FORCE       5.183E-07
 LARGEST RESIDUAL FORCE            -1.132E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.300       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.349E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.029E-05   TIME AVG. MOMENT      5.958E-06
 LARGEST RESIDUAL MOMENT            1.241E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.152E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.648E-05   AT NODE        192   DOF  4
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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.645E-07   TIME AVG. FORCE       5.181E-07
 LARGEST RESIDUAL FORCE            -3.862E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.300       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.718E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.029E-05   TIME AVG. MOMENT      5.958E-06
 LARGEST RESIDUAL MOMENT            1.449E-13   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.152E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.598E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.731E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.891E-02,  FRACTION OF STEP COMPLETED  0.329    
 STEP TIME COMPLETED       0.659    ,  TOTAL TIME COMPLETED        0.659    


  INCREMENT    12 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.731E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.869E-06   TIME AVG. FORCE       6.870E-07
 LARGEST RESIDUAL FORCE             5.123E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.306       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.554E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.136E-05   TIME AVG. MOMENT      6.408E-06
 LARGEST RESIDUAL MOMENT           -1.156E-05   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.192E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -7.782E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.253E-06   TIME AVG. FORCE       6.100E-07
 LARGEST RESIDUAL FORCE            -2.294E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.292       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.380E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.101E-05   TIME AVG. MOMENT      6.379E-06
 LARGEST RESIDUAL MOMENT            1.102E-08   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.126E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.623E-04   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -6.315E-07
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


 AVERAGE FORCE                      8.463E-07   TIME AVG. FORCE       5.591E-07
 LARGEST RESIDUAL FORCE            -9.238E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.289       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.660E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.101E-05   TIME AVG. MOMENT      6.379E-06
 LARGEST RESIDUAL MOMENT            3.332E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.113E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.330E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    4.018E-05
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


 AVERAGE FORCE                      8.367E-07   TIME AVG. FORCE       5.579E-07
 LARGEST RESIDUAL FORCE            -2.869E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.289       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.689E-04   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.101E-05   TIME AVG. MOMENT      6.379E-06
 LARGEST RESIDUAL MOMENT            2.374E-11   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.110E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.342E-05   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      8.365E-07   TIME AVG. FORCE       5.579E-07
 LARGEST RESIDUAL FORCE            -1.347E-12   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.289       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.203E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.101E-05   TIME AVG. MOMENT      6.379E-06
 LARGEST RESIDUAL MOMENT            1.329E-13   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.110E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.610E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.580E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.731E-02,  FRACTION OF STEP COMPLETED  0.353    
 STEP TIME COMPLETED       0.706    ,  TOTAL TIME COMPLETED        0.706    


  INCREMENT    13 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.580E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.630E-06   TIME AVG. FORCE       6.771E-07
 LARGEST RESIDUAL FORCE             3.871E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.292       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.259E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.196E-05   TIME AVG. MOMENT      6.808E-06
 LARGEST RESIDUAL MOMENT           -8.817E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.139E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -6.378E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.226E-06   TIME AVG. FORCE       6.322E-07
 LARGEST RESIDUAL FORCE            -1.800E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.280       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.216E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.171E-05   TIME AVG. MOMENT      6.789E-06
 LARGEST RESIDUAL MOMENT            7.442E-09   AT NODE        103   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.080E-02   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.873E-04   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -4.957E-07
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


 AVERAGE FORCE                      9.101E-07   TIME AVG. FORCE       5.970E-07
 LARGEST RESIDUAL FORCE            -3.340E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.279       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.593E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.170E-05   TIME AVG. MOMENT      6.788E-06
 LARGEST RESIDUAL MOMENT            1.600E-09   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.072E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.999E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      9.072E-07   TIME AVG. FORCE       5.967E-07
 LARGEST RESIDUAL FORCE            -6.283E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.278       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.186E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.170E-05   TIME AVG. MOMENT      6.788E-06
 LARGEST RESIDUAL MOMENT            4.095E-12   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.071E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.097E-05   AT NODE        192   DOF  4
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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.071E-07   TIME AVG. FORCE       5.967E-07
 LARGEST RESIDUAL FORCE            -1.143E-12   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.278       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.349E-07   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.170E-05   TIME AVG. MOMENT      6.788E-06
 LARGEST RESIDUAL MOMENT           -1.150E-13   AT NODE        175   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.071E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.702E-08   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.440E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.580E-02,  FRACTION OF STEP COMPLETED  0.376    
 STEP TIME COMPLETED       0.752    ,  TOTAL TIME COMPLETED        0.752    


  INCREMENT    14 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.440E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.489E-06   TIME AVG. FORCE       6.859E-07
 LARGEST RESIDUAL FORCE             2.940E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.280       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.022E-02   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.257E-05   TIME AVG. MOMENT      7.201E-06
 LARGEST RESIDUAL MOMENT           -6.757E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.090E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.245E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.211E-06   TIME AVG. FORCE       6.582E-07
 LARGEST RESIDUAL FORCE            -1.367E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.269       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.055E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.237E-05   TIME AVG. MOMENT      7.187E-06
 LARGEST RESIDUAL MOMENT            5.234E-09   AT NODE        107   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.039E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.123E-04   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -3.968E-07
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


 AVERAGE FORCE                      9.781E-07   TIME AVG. FORCE       6.349E-07
 LARGEST RESIDUAL FORCE            -1.169E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.268       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.382E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.237E-05   TIME AVG. MOMENT      7.187E-06
 LARGEST RESIDUAL MOMENT            7.258E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.034E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.736E-05   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      9.772E-07   TIME AVG. FORCE       6.348E-07
 LARGEST RESIDUAL FORCE            -1.244E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.268       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.693E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.237E-05   TIME AVG. MOMENT      7.187E-06
 LARGEST RESIDUAL MOMENT            6.777E-13   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.034E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.888E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT REDUCED TO   4.307E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.440E-02,  FRACTION OF STEP COMPLETED  0.398    
 STEP TIME COMPLETED       0.796    ,  TOTAL TIME COMPLETED        0.796    


  INCREMENT    15 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.307E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.413E-06   TIME AVG. FORCE       7.055E-07
 LARGEST RESIDUAL FORCE             2.243E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.269       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.300E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.317E-05   TIME AVG. MOMENT      7.585E-06
 LARGEST RESIDUAL MOMENT           -5.200E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.046E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.325E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.215E-06   TIME AVG. FORCE       6.875E-07
 LARGEST RESIDUAL FORCE            -1.017E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.260       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.066E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.302E-05   TIME AVG. MOMENT      7.576E-06
 LARGEST RESIDUAL MOMENT            3.600E-09   AT NODE        107   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.002E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.423E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -3.061E-07
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


 AVERAGE FORCE                      1.048E-06   TIME AVG. FORCE       6.723E-07
 LARGEST RESIDUAL FORCE            -3.977E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.259       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.445E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.302E-05   TIME AVG. MOMENT      7.576E-06
 LARGEST RESIDUAL MOMENT            3.165E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.993E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.763E-05   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.048E-06   TIME AVG. FORCE       6.723E-07
 LARGEST RESIDUAL FORCE            -2.300E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.259       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.144E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.302E-05   TIME AVG. MOMENT      7.576E-06
 LARGEST RESIDUAL MOMENT            2.353E-13   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.991E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.100E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.150    
          TIME INCREMENT REDUCED TO   4.182E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.307E-02,  FRACTION OF STEP COMPLETED  0.420    
 STEP TIME COMPLETED       0.840    ,  TOTAL TIME COMPLETED        0.840    


  INCREMENT    16 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.182E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.380E-06   TIME AVG. FORCE       7.313E-07
 LARGEST RESIDUAL FORCE             1.718E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.258       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.745E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.376E-05   TIME AVG. MOMENT      7.962E-06
 LARGEST RESIDUAL MOMENT           -4.019E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -1.006E-02   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -3.576E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.236E-06   TIME AVG. FORCE       7.193E-07
 LARGEST RESIDUAL FORCE            -7.477E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.250       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.745E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.365E-05   TIME AVG. MOMENT      7.955E-06
 LARGEST RESIDUAL MOMENT            2.451E-09   AT NODE        111   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.680E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.796E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -2.315E-07
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


 AVERAGE FORCE                      1.117E-06   TIME AVG. FORCE       7.093E-07
 LARGEST RESIDUAL FORCE            -1.323E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.250       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.125E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.365E-05   TIME AVG. MOMENT      7.955E-06
 LARGEST RESIDUAL MOMENT            1.347E-10   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.664E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.594E-05   AT NODE        192   DOF  4
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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.117E-06   TIME AVG. FORCE       7.093E-07
 LARGEST RESIDUAL FORCE            -4.021E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.250       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.733E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.365E-05   TIME AVG. MOMENT      7.955E-06
 LARGEST RESIDUAL MOMENT           -1.343E-13   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.663E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     8.827E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.150    
          TIME INCREMENT REDUCED TO   4.064E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.182E-02,  FRACTION OF STEP COMPLETED  0.441    
 STEP TIME COMPLETED       0.881    ,  TOTAL TIME COMPLETED        0.881    


  INCREMENT    17 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  4.064E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.377E-06   TIME AVG. FORCE       7.607E-07
 LARGEST RESIDUAL FORCE             1.320E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.249       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.478E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.435E-05   TIME AVG. MOMENT      8.332E-06
 LARGEST RESIDUAL MOMENT           -3.117E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.686E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.961E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.269E-06   TIME AVG. FORCE       7.524E-07
 LARGEST RESIDUAL FORCE            -5.464E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.242       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.597E-03   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.426E-05   TIME AVG. MOMENT      8.326E-06
 LARGEST RESIDUAL MOMENT            1.657E-09   AT NODE        111   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.362E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.249E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.727E-07
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


 AVERAGE FORCE                      1.185E-06   TIME AVG. FORCE       7.459E-07
 LARGEST RESIDUAL FORCE            -4.307E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.242       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.773E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.426E-05   TIME AVG. MOMENT      8.326E-06
 LARGEST RESIDUAL MOMENT            5.631E-11   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.352E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     9.096E-06   AT NODE        192   DOF  4
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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.185E-06   TIME AVG. FORCE       7.460E-07
 LARGEST RESIDUAL FORCE            -7.060E-12   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.242       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.131E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.426E-05   TIME AVG. MOMENT      8.326E-06
 LARGEST RESIDUAL MOMENT           -1.380E-13   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.352E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.652E-07   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.150    
          TIME INCREMENT REDUCED TO   3.953E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  4.064E-02,  FRACTION OF STEP COMPLETED  0.461    
 STEP TIME COMPLETED       0.922    ,  TOTAL TIME COMPLETED        0.922    


  INCREMENT    18 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.953E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.393E-06   TIME AVG. FORCE       7.922E-07
 LARGEST RESIDUAL FORCE             1.018E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.240       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.444E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.492E-05   TIME AVG. MOMENT      8.692E-06
 LARGEST RESIDUAL MOMENT           -2.425E-06   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.342E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.457E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.310E-06   TIME AVG. FORCE       7.863E-07
 LARGEST RESIDUAL FORCE            -3.988E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.234       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.617E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.485E-05   TIME AVG. MOMENT      8.689E-06
 LARGEST RESIDUAL MOMENT            1.117E-09   AT NODE        111   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.064E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.779E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -1.280E-07
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


 AVERAGE FORCE                      1.252E-06   TIME AVG. FORCE       7.821E-07
 LARGEST RESIDUAL FORCE            -1.376E-09   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.234       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.954E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.485E-05   TIME AVG. MOMENT      8.689E-06
 LARGEST RESIDUAL MOMENT            2.333E-11   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.059E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     5.140E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.150    
          TIME INCREMENT REDUCED TO   3.847E-02 TO MEET NEXT OUTPUT TIME POINT

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.953E-02,  FRACTION OF STEP COMPLETED  0.481    
 STEP TIME COMPLETED       0.962    ,  TOTAL TIME COMPLETED        0.962    


  INCREMENT    19 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.847E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.424E-06   TIME AVG. FORCE       8.249E-07
 LARGEST RESIDUAL FORCE             7.857E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.231       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.591E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.547E-05   TIME AVG. MOMENT      9.046E-06
 LARGEST RESIDUAL MOMENT           -1.889E-06   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -9.021E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.039E-04   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.358E-06   TIME AVG. FORCE       8.205E-07
 LARGEST RESIDUAL FORCE            -2.915E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.226       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.786E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.543E-05   TIME AVG. MOMENT      9.043E-06
 LARGEST RESIDUAL MOMENT            7.534E-10   AT NODE        111   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.783E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.378E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -9.486E-08
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


 AVERAGE FORCE                      1.318E-06   TIME AVG. FORCE       8.178E-07
 LARGEST RESIDUAL FORCE            -4.249E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -0.226       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.498E-05   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.543E-05   TIME AVG. MOMENT      9.043E-06
 LARGEST RESIDUAL MOMENT            9.502E-12   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.780E-03   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.858E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.150    

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.847E-02,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         1.00    


  INCREMENT    20 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.150    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.449E-03   TIME AVG. FORCE       1.538E-04
 LARGEST RESIDUAL FORCE            -0.115       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.73       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.844       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     9.208E-05   TIME AVG. MOMENT      1.320E-05
 LARGEST RESIDUAL MOMENT           -1.736E-03   AT NODE        114   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.208E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -4.785E-02   AT NODE        189   DOF  4
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


 AVERAGE FORCE                      1.051E-05   TIME AVG. FORCE       1.423E-06
 LARGEST RESIDUAL FORCE             1.887E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -1.72       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.423E-02   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.122E-05   TIME AVG. MOMENT      9.652E-06
 LARGEST RESIDUAL MOMENT           -2.642E-05   AT NODE         59   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -8.202E-02   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     2.181E-03   AT NODE         44   DOF  5
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
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.

     NEGATIVE EIGENVALUES MEAN THAT THE SYSTEM MATRIX IS NOT 
     POSITIVE DEFINITE:
     FOR EXAMPLE, A BIFURCATION (BUCKLING) LOAD MAY HAVE BEEN EXCEEDED.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.850E-02   TIME AVG. FORCE       5.532E-03
 LARGEST RESIDUAL FORCE             -6.06       AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         -7.25       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -5.53       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     3.498E-02   TIME AVG. MOMENT      1.757E-03
 LARGEST RESIDUAL MOMENT             1.10       AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.383       AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.301       AT NODE        189   DOF  4
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
 

 ***WARNING: THE SYSTEM MATRIX HAS 13 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.562E-02   TIME AVG. FORCE       1.602E-03
 LARGEST RESIDUAL FORCE             -2.95       AT NODE        165   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -7.47       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.921       AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     2.036E-02   TIME AVG. MOMENT      1.027E-03
 LARGEST RESIDUAL MOMENT            -1.59       AT NODE        170   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -0.611       AT NODE        165   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -0.273       AT NODE        188   DOF  5
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
 

 ***WARNING: THE SYSTEM MATRIX HAS 27 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                       1.09       TIME AVG. FORCE       6.831E-02
 LARGEST RESIDUAL FORCE             -79.1       AT NODE        128   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -12.5       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -6.95       AT NODE         73   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      1.08       TIME AVG. MOMENT      5.404E-02
 LARGEST RESIDUAL MOMENT            -87.1       AT NODE         90   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      -1.60       AT NODE        132   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION      1.49       AT NODE         89   DOF  4
   INSTANCE: PART-1-1                                                                        
          MOMENT    EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.
 

 ***NOTE: THE SOLUTION APPEARS TO BE DIVERGING. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    20 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  3.750E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.316E-05   TIME AVG. FORCE       4.714E-06
 LARGEST RESIDUAL FORCE            -2.084E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -5.573E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.165       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.686E-05   TIME AVG. MOMENT      9.434E-06
 LARGEST RESIDUAL MOMENT            4.531E-05   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.278E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     6.292E-03   AT NODE        189   DOF  4
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


 AVERAGE FORCE                      1.209E-06   TIME AVG. FORCE       8.423E-07
 LARGEST RESIDUAL FORCE            -7.908E-08   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -5.547E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.124E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.558E-05   TIME AVG. MOMENT      9.370E-06
 LARGEST RESIDUAL MOMENT            2.699E-08   AT NODE        163   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.263E-03   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.387E-05   AT NODE        132   DOF  5
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION    2.018E-08
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


 AVERAGE FORCE                      1.228E-06   TIME AVG. FORCE       8.434E-07
 LARGEST RESIDUAL FORCE            -1.855E-06   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -5.167E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.799E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.557E-05   TIME AVG. MOMENT      9.370E-06
 LARGEST RESIDUAL MOMENT           -6.246E-10   AT NODE         86   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.073E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.902E-04   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 ESTIMATE OF ROTATION CORRECTION   -4.402E-06
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


 AVERAGE FORCE                      1.207E-06   TIME AVG. FORCE       8.422E-07
 LARGEST RESIDUAL FORCE            -1.033E-12   AT NODE        179   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -5.167E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.622E-06   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.557E-05   TIME AVG. MOMENT      9.370E-06
 LARGEST RESIDUAL MOMENT           -3.674E-13   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.073E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -8.541E-08   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    5.625E-02

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.750E-02,  FRACTION OF STEP COMPLETED  0.519    
 STEP TIME COMPLETED        1.04    ,  TOTAL TIME COMPLETED         1.04    


  INCREMENT    21 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  5.625E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.297E-05   TIME AVG. FORCE       2.732E-06
 LARGEST RESIDUAL FORCE            -1.030E-03   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         3.372E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.111       AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.545E-05   TIME AVG. MOMENT      9.659E-06
 LARGEST RESIDUAL MOMENT            7.463E-07   AT NODE        187   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      1.354E-03   AT NODE        189   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.464E-03   AT NODE        190   DOF  4
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


 AVERAGE FORCE                      1.323E-06   TIME AVG. FORCE       8.705E-07
 LARGEST RESIDUAL FORCE             1.536E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         3.376E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.096E-04   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.548E-05   TIME AVG. MOMENT      9.660E-06
 LARGEST RESIDUAL MOMENT            1.664E-10   AT NODE         84   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      1.354E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     7.327E-07   AT NODE        133   DOF  5
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

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.324E-06   TIME AVG. FORCE       8.705E-07
 LARGEST RESIDUAL FORCE            -7.948E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         3.379E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.273E-05   AT NODE          1   DOF  3
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.548E-05   TIME AVG. MOMENT      9.660E-06
 LARGEST RESIDUAL MOMENT           -9.766E-13   AT NODE         59   DOF  6
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      1.355E-03   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     1.245E-06   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    8.438E-02

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  5.625E-02,  FRACTION OF STEP COMPLETED  0.547    
 STEP TIME COMPLETED        1.09    ,  TOTAL TIME COMPLETED         1.09    


  INCREMENT    22 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  8.438E-02

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  1.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.991E-06   TIME AVG. FORCE       1.044E-06
 LARGEST RESIDUAL FORCE            -9.729E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         1.635E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.433E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.544E-05   TIME AVG. MOMENT      9.923E-06
 LARGEST RESIDUAL MOMENT           -4.572E-07   AT NODE        183   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      6.569E-04   AT NODE        190   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.376E-03   AT NODE        192   DOF  4
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


 AVERAGE FORCE                      1.319E-06   TIME AVG. FORCE       8.954E-07
 LARGEST RESIDUAL FORCE             1.413E-10   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         1.635E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.021E-05   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.544E-05   TIME AVG. MOMENT      9.923E-06
 LARGEST RESIDUAL MOMENT            1.466E-11   AT NODE         80   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      6.569E-04   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.519E-08   AT NODE         64   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.127    

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  8.438E-02,  FRACTION OF STEP COMPLETED  0.589    
 STEP TIME COMPLETED        1.18    ,  TOTAL TIME COMPLETED         1.18    


  INCREMENT    23 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.127    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.728E-06   TIME AVG. FORCE       9.918E-07
 LARGEST RESIDUAL FORCE            -5.500E-05   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -1.176E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.570E-02   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.543E-05   TIME AVG. MOMENT      1.016E-05
 LARGEST RESIDUAL MOMENT           -5.563E-08   AT NODE        175   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.715E-05   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.032E-03   AT NODE        191   DOF  4
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


 AVERAGE FORCE                      1.319E-06   TIME AVG. FORCE       9.177E-07
 LARGEST RESIDUAL FORCE             4.527E-11   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -1.175E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.132E-05   AT NODE        189   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.544E-05   TIME AVG. MOMENT      1.016E-05
 LARGEST RESIDUAL MOMENT            7.348E-12   AT NODE         76   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -4.717E-05   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -2.637E-08   AT NODE         88   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.190    

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.127    ,  FRACTION OF STEP COMPLETED  0.652    
 STEP TIME COMPLETED        1.30    ,  TOTAL TIME COMPLETED         1.30    


  INCREMENT    24 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  1.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.321E-06   TIME AVG. FORCE       9.378E-07
 LARGEST RESIDUAL FORCE            -1.011E-07   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.        -6.607E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.101E-03   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.544E-05   TIME AVG. MOMENT      1.038E-05
 LARGEST RESIDUAL MOMENT            1.839E-11   AT NODE        186   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.653E-05   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.423E-05   AT NODE        191   DOF  4
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


 AVERAGE FORCE                      1.319E-06   TIME AVG. FORCE       9.378E-07
 LARGEST RESIDUAL FORCE             2.806E-12   AT NODE        165   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -6.607E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.058E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.544E-05   TIME AVG. MOMENT      1.038E-05
 LARGEST RESIDUAL MOMENT            1.241E-13   AT NODE        120   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION     -2.653E-05   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -9.330E-11   AT NODE         48   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED
          TIME INCREMENT MAY NOW INCREASE TO    0.285    

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.747    
 STEP TIME COMPLETED        1.49    ,  TOTAL TIME COMPLETED         1.49    


  INCREMENT    25 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.285    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.321E-06   TIME AVG. FORCE       9.560E-07
 LARGEST RESIDUAL FORCE            -8.266E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         5.238E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.963E-04   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                     1.544E-05   TIME AVG. MOMENT      1.059E-05
 LARGEST RESIDUAL MOMENT            3.677E-12   AT NODE        182   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.106E-07   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     4.000E-05   AT NODE        191   DOF  4
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


 AVERAGE FORCE                      1.319E-06   TIME AVG. FORCE       9.559E-07
 LARGEST RESIDUAL FORCE             2.349E-12   AT NODE          1   DOF  1
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         5.240E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.704E-08   AT NODE          1   DOF  2
   INSTANCE: REFPOINT-0                                                                      
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.544E-05   TIME AVG. MOMENT      1.059E-05
 LARGEST RESIDUAL MOMENT           -1.734E-13   AT NODE        181   DOF  5
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      2.105E-07   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -5.030E-11   AT NODE        101   DOF  4
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.285    ,  FRACTION OF STEP COMPLETED  0.890    
 STEP TIME COMPLETED        1.78    ,  TOTAL TIME COMPLETED         1.78    


  INCREMENT    26 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.221    

  NUMBER OF EQUATIONS =      1143     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  2.15E+06

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.319E-06   TIME AVG. FORCE       9.724E-07
 LARGEST RESIDUAL FORCE             1.991E-12   AT NODE          1   DOF  1
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         4.938E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.762E-07   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE MOMENT                     1.545E-05   TIME AVG. MOMENT      1.077E-05
 LARGEST RESIDUAL MOMENT           -1.673E-13   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      1.984E-07   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION     3.521E-08   AT NODE        191   DOF  4
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


 AVERAGE FORCE                      1.319E-06   TIME AVG. FORCE       9.724E-07
 LARGEST RESIDUAL FORCE             1.852E-12   AT NODE          1   DOF  1
   INSTANCE: REFPOINT-0                                                                      
 LARGEST INCREMENT OF DISP.         4.938E-06   AT NODE        189   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.667E-13   AT NODE        116   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                     1.545E-05   TIME AVG. MOMENT      1.077E-05
 LARGEST RESIDUAL MOMENT            1.464E-13   AT NODE        192   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF ROTATION      1.984E-07   AT NODE        191   DOF  4
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO ROTATION    -1.087E-12   AT NODE         32   DOF  5
   INSTANCE: PART-1-1                                                                        
          THE MOMENT    EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.221    ,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        2.00    ,  TOTAL TIME COMPLETED         2.00    


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF         26  INCREMENTS
                       1  CUTBACKS IN AUTOMATIC INCREMENTATION
                     115  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     115  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     115  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       2  WARNING MESSAGES DURING USER INPUT PROCESSING
                       7  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       3  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   3.2100    
       SYSTEM TIME (SEC)    =  0.21000    
       TOTAL CPU TIME (SEC) =   3.4200    
       WALLCLOCK TIME (SEC) =          4

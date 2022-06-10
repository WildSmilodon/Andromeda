subroutine solveStokes()
      
      USE mPar
      USE mEqns
      implicit none

!     Local vars
!
      REAL ts,te    

      if (parSLEtype .EQ. parLIS) then
            CALL WriteToLog("Preparing a SLE for all 3 directions!")
!
!           Set up X and RHS vectors (for 3 individual X,Y,Z equations)
!
            CALL sdSetUpXandB()
!
!           Make a big system of equation for all 3 directions
!                                     
            CALL StokesBigSystemXB()              
!
!           Set system matrix size and disribution across processors
!
            call plis_setup_A(stk%sle,stk%neq)
            call DivideRows()
!
!           Set solver options
!                              
            CALL WriteToLog("LIS solver options: "//TRIM(parLISslvSet))
            CALL plis_createSolver(stk%sle)            
!
!           Integrate & form system and r.h.s matrices
!
            CALL stokesFormLISsysMrhsMcsrDiv()
!
!           Use boundary conditions to set up left and right hand side vectors
!
            CALL formStokesLeftRightHandSideVectors()            
!            
!           Solve the system                  
!        
            CALL StokesBigSystemSOLVElis()
!
!           Use explicit calculation for pressure
!                                    
            CALL pressureStokesLIS()
      end if      


      if (parSLEtype .EQ. parLISfull) then
            CALL WriteToLog("Preparing a SLE for all 3 directions!")
!
!           Set up X and RHS vectors (for 3 individual X,Y,Z equations)
!
            CALL sdSetUpXandB()
!
!           Make a big system of equation for all 3 directions
!                                     
            CALL StokesBigSystemXB()              
!
!           Set system matrix size and disribution across processors
!
            call plis_setup_A(stk%sle,stk%neq)
!
!           Set solver options
!                              
            CALL plis_createSolver(stk%sle)            
!
!           Integrate
!
            CALL WriteToLog("Integration!")
            CALL CoR_Integrals_Stokes()                        
!
!           Form CRS system and rhs matrices
!                        
            CALL stokesFormLISsysMrhsMcsr()  ! LIS parallel VERSION
!            
!           Solve the system                     
!        
            CALL StokesBigSystemSOLVElis() ! LIS parallel VERSION

      end if     

      if (parSLEtype .EQ. parCRS .OR. parSLEtype .EQ. parFULL) then
!
!           Integrate
!
            CALL CPU_TIME(ts)
            CALL CoR_Integrals_Stokes()
            CALL CPU_TIME(te)
            WRITE (parLogTekst,'(A,F10.4)') "TIMER :: CoR_Integrals_Stokes [s] = ",te-ts
            CALL WriteToLog(parLogTekst)            
!
!           Set up X and RHS vectors (for 3 individual X,Y,Z equations)
!
            CALL sdSetUpXandB()
!
!           Make a big system of equation for all 3 directions
!                         
            CALL WriteToLog("Setting up System and RHS matrices from MEMORY!")
            CALL StokesBigSystemXB()                   
!
!           Solve
!                                    
            if (parSLEtype .EQ. parCRS) then
!
!                 Form CRS system and rhs matrices
!                        
                  CALL CPU_TIME(ts)
                  CALL stokesFormCRSsysMrhsM()  ! CRS VERSION
                  CALL CPU_TIME(te)
                  WRITE (parLogTekst,'(A,F10.4)') "TIMER :: stokesFormCRSsysMrhsM [s] = ",te-ts
                  CALL WriteToLog(parLogTekst)
               
                  IF (parTrii.EQ.parYes) THEN
                        CALL WriteToLog("Solving ... Stokes 3 sides inlet run!") 
                        CALL StokesInlet3sidesForceCRS() ! CRS VERSION
                        CALL StokesInlet3sidesTorqueOmegaCRS()
                        CALL StokesInlet3sidesTorquePiCRS()
                  END IF

                  IF (parSval.EQ.parYes) THEN
                        CALL WriteToLog("Solving ... Stokes from all sides run!") 
                        CALL StokesDragFromAllSidesCRS() 
                  END IF

                  IF (parFlop.EQ.parYes) THEN
                        CALL WriteToLog("Solving ... Stokes flow over rotated particle!") 
                        CALL StokesFLOPcrs()
                  END IF

                  IF (parSval.EQ.parNo.AND.parTrii.EQ.parNo.AND.parFlop.EQ.parNo) THEN
                        CALL WriteToLog("Solving ...") 
                        CALL CPU_TIME(ts)
                        CALL StokesBigSystemSOLVEcrs() ! CRS VERSION
                        CALL CPU_TIME(te)
                        WRITE (parLogTekst,'(A,F10.4)') "TIMER :: StokesBigSystemSOLVEcrs [s] = ",te-ts
                        CALL WriteToLog(parLogTekst)
         
                  END IF                  


            end if
            if (parSLEtype .EQ. parFULL) then
!
!                 Form FULL system and rhs matrices
!    
                  IF (parSval.EQ.parYes) THEN
                        CALL WriteToLog("Solving ... Stokes validation run!") 
                        CALL StokesEllipsoidValidationSOLVE() 
                  END IF
                  IF (parTrii.EQ.parYes) THEN
                        CALL WriteToLog("Solving ... Stokes 3 sides inlet run!") 
                        CALL StokesInlet3sidesForceTorque() 
                  END IF

                  IF (parSval.EQ.parNo.AND.parTrii.EQ.parNo) THEN
                        CALL WriteToLog("Solving ...") 
                        CALL StokesBigSystemSOLVE()
                  END IF

           endif

           CALL CPU_TIME(ts)
           call pressureStokes() 
           CALL CPU_TIME(te)
           WRITE (parLogTekst,'(A,F10.4)') "TIMER :: pressureStokes [s] = ",te-ts
           CALL WriteToLog(parLogTekst)

      END IF      

end subroutine      
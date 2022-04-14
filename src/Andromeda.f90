program Andromeda
!
!     Andromeda, a BEM code
!
!     author  - jure.ravnik@um.si
!
!
!     Include modules
!
      USE mMesh
      USE GaussIntegration
      USE Triangle
      USE mPar
      USE mEqns
      use linFlowFields
      use mCommon
      use parallel
      use plis 
#include "lisf.h" 

!     Local vars
!
      INTEGER i
      CHARACTER(123) SLEtype ! te
      REAL ts,te
!
!     Name and version of the code
!
      parIDname='Andromeda'
      parIDversion='1.6'
      parIDdate='April 2022'
!
!     Init parallel environment
!      
      call par_init()
!
!     Get start time and computer name
!
      CALL DatumInUra(parStartTime)
      CALL hostnm(parHostname)
!
!     Set up Log file
!
      call defineKeys()
      call readCLargs()
      CALL SetUpLogFile()
!
!     Read Parameters
!
      CALL ReadParameters()
!
!     Init integration routines
!
      CALL GaussPosWeights
      CALL Triangle_IntegrationInit(parTriInteg)
!
!     Read Mesh
!
      CALL ReadMesh(parMeshFullName)
      CALL CalMeshArea()
      CALL CalMeshNormals()
      !CALL VerifyMeshNormals()  ! daj pod IF !!! velja samo za primer delcev
      CALL CalQMesh()
!
!     Write mesh stats to log file
!

      CALL PrintMeshStats()
!
!     Read subdomain definitions from BIC file
!
      CALL WriteToLog("SetUpSubdomains!")
      CALL ReadsdBIC(parBiCFullName)
      CALL SetUpSubdomains()
!
!     Handle corners and edges
!
      CALL WriteToLog("HandleCornersEdges!")
      CALL HandleCornersEdges(parBiCFullName)
!
!     Read Boundary and Inital conditions file
!
      CALL WriteToLog("sdApplyInitalCond!")
      CALL ReadBIC(parBiCFullName)
      CALL sdApplyInitalCond()
!
!     Is this a postprocessing only run?
!
      IF (parPostOnly.EQ.parYes) THEN
        CALL PostOnly()
        GOTO 3141
      END IF
!
!     Output and.initial.vtu
!
      IF (parOutInit.EQ.parYes) CALL OutputInitialParaview()
!
!     Start
!
      SLEtype = "lis"
      !SLEtype = "lisFull"
      !SLEtype = "crs"
      !SLEtype = "full"
      !crs = .false.
      

      if (parPrType .EQ. parStokes .AND. SLEtype.EQ."lis") then
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


      if (parPrType .EQ. parStokes .AND. SLEtype.EQ."lisFull") then
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

      IF (parPrType .EQ. parStokes.AND..NOT.(SLEtype.EQ."lis").AND..NOT.(SLEtype.EQ."lisFull")) THEN
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
           if (SLEtype.EQ."crs") then
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
            if (SLEtype.EQ."full") then
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

      IF (parPrType .EQ. parLaplace) THEN
!
!           Integrate
!
            CALL CoR_Integrals()
!
!           Set up X and RHS vectors
!
            CALL sdSetUpXandB()
!
!           Set up System Matrix and RHS matrix
!
            IF (parfromDisk.EQ.parYes) THEN
!
!                 Set up matrices from disk
!
                  CALL WriteToLog("Setting up System and RHS matrices from HARD DISK!")
                  DO i=1,neq
                        CALL sdSetUpSysMrhsM_CRS_fromDisk(i)
                  END DO
            ELSE
!
!                 Set up matrices from memory
!
                  CALL WriteToLog("Setting up System and RHS matrices from MEMORY!")
                  DO i=1,neq
                        CALL sdSetUpSysMrhsM_CRS(i)
                  END DO
!
!                 Deallocate integral matrices
!
                  DO i=1,nosd
                        DEALLOCATE (subdomain(i)%Hmat)
                        DEALLOCATE (subdomain(i)%Gmat)
                  END DO
            END IF
!
!     Solve
!
            DO i=1,neq
                  CALL sdSolveCRS(i)
            END DO
      END IF

!
!     Write restart file
!
      IF (parOutRest.EQ.parYes.AND.amIroot) CALL WriteRestartFile()
!
!     Start postprocessing
!
3141  CONTINUE
!
!     Output results
!
      IF (parOutProfile.EQ.parYes.AND.amIroot) CALL linePostProcessing()
      IF (parOutPiec.EQ.parYes.AND.amIroot) CALL OutputPiecesParaview()
      IF (parOutMesh.EQ.parYes.AND.amIroot) CALL OutputMeshParaview()
      IF (parOutStoc.EQ.parYes.AND.amIroot) CALL OutputStochastic()
      IF (parOutDomainProfile.EQ.parYes.AND.amIroot) CALL OutputDomainProfile()

      IF (parQinw.EQ.parYes.AND.amIroot) CALL IntegrateFluxes()
      IF (parUinw.EQ.parYes.AND.amIroot) CALL IntegrateFunction()
      IF (parTinw.EQ.parYes.AND.amIroot) CALL IntegrateTorques()


123   continue      
!
!     Close log file
!
      CALL StopProgram()

end program Andromeda

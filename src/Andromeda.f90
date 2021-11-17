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
!
!     Local vars
!
      INTEGER i
      logical crs ! te
!
!     Name and version of the code
!
      parIDname='Andromeda'
      parIDversion='1.4'
      parIDdate='October 2021'
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
!     Integrate
!
      IF (parPrType .EQ. parStokes) THEN
!
!           Integrate
!
            CALL CoR_Integrals_Stokes()
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
           crs = .true.
           !crs = .false.

           if (crs) then
!
!           Form CRS system and rhs matrices
!                        
                  CALL stokesFormCRSsysMrhsM()  ! CRS VERSION

                  IF (parTrii.EQ.parYes) THEN
                        CALL WriteToLog("Solving ... Stokes 3 sides inlet run!") 
                        CALL StokesInlet3sidesForceCRS() ! CRS VERSION
                        CALL StokesInlet3sidesTorqueOmegaCRS()
                        CALL StokesInlet3sidesTorquePiCRS()
                  END IF

                  IF (parSval.EQ.parNo.AND.parTrii.EQ.parNo) THEN
                        CALL WriteToLog("Solving ...") 
                        CALL StokesBigSystemSOLVEcrs() ! CRS VERSION
                  END IF                  


           else

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

           call pressureStokes() 




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
      IF (parOutRest.EQ.parYes) CALL WriteRestartFile()
!
!     Start postprocessing
!
3141  CONTINUE
!
!     Output results
!
      IF (parOutProfile.EQ.parYes) CALL linePostProcessing()
      IF (parOutPiec.EQ.parYes) CALL OutputPiecesParaview()
      IF (parOutMesh.EQ.parYes) CALL OutputMeshParaview()
      IF (parOutStoc.EQ.parYes) CALL OutputStochastic()
      IF (parOutDomainProfile.EQ.parYes) CALL OutputDomainProfile()

      IF (parQinw.EQ.parYes) CALL IntegrateFluxes()
      IF (parUinw.EQ.parYes) CALL IntegrateFunction()
      IF (parTinw.EQ.parYes) CALL IntegrateTorques()

!
!     Close log file
!
      CALL StopProgram
END

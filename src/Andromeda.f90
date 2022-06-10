program Andromeda
!
!     Andromeda, a BEM code
!
!     author  - jure.ravnik@um.si
!
!
!     Include modules
!
      USE mPar
      USE parallel
      USE GaussIntegration
      USE Triangle     
!
!     Include LIS library header file
!                  
#include "lisf.h" 
!
!     Name and version of the code
!
      parIDname='Andromeda'
      parIDversion='1.6'
      parIDdate='June 2022'
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
!     Solve the Stokes equation in a single domain
!            
      IF (parPrType .EQ. parStokes) THEN
            call solveStokes()
      END IF
!
!     Solve the Laplace equation in several subdomains
!            
      IF (parPrType .EQ. parLaplace) THEN
            call solveLaplace()
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
!
!     Close log files & stop program
!
      CALL StopProgram()

end program Andromeda

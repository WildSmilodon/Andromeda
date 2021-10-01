!
!     parameters module
!     (c) Jure Ravnik, 2017
!
      MODULE mPar
      
      USE mCommon
      IMPLICIT NONE

!
!     Types
!
      TYPE DomainLineType

        CHARACTER(255) name ! profile name
        INTEGER nno  ! number of nodes
        REAL(rk) startPoint(3)
        REAL(rk) endPoint(3)

      END TYPE
!
! ----------------------------------------------------------------------
!

!
!     Variables
!
      CHARACTER(255) parMeshDir
      CHARACTER(255) parMeshFileName,parBiCFileName
      CHARACTER(255) parMeshFullName,parBiCFullName
      INTEGER parTriInteg             ! regular triangle integration (1..5)
      INTEGER parTriRecur             ! singular triangle integration (number of recursive steps)
      INTEGER parQuadIntegSing        ! singular quad integration (1..8)
      INTEGER parQuadIntegRegu        ! regualr quad integration (1..8)
      INTEGER parScreenOutput

      CHARACTER(39) parStartTime,parEndTime
      CHARACTER(255) parIDname
      CHARACTER(255) parIDversion
      CHARACTER(255) parIDdate
      CHARACTER(255) parHostname

!     Mesh transromation
      INTEGER pariMTrot
      INTEGER pariMTstr
      INTEGER pariMTtra
      INTEGER pariMTcyl,pariMTsel

      REAL(rk) parMTrot(3)
      REAL(rk) parMTstr(3)
      REAL(rk) parMTtra(3)
      REAL(rk) parMTcyl(6)
      REAL(rk) parMTsel(5)

      INTEGER parWriteIntegrals,parOutRest,parPostOnly,parConserveMemory
      INTEGER parfromDisk,parSval,parTrii

      INTEGER parPrType,parStokes,parLaplace

!     Internal parameters, fixed
      INTEGER parLogLun
      CHARACTER(255) parInputFileName
      CHARACTER(255) parResultsFileName,parStochasticFileName
      CHARACTER(255) parInitialFileName,parResultsWallName,parResultsForceErr
      CHARACTER(255) parLogFileName,parIntegralsFileName
      INTEGER parYes,parNo
      CHARACTER(255) parLogTekst
      CHARACTER(255) parPPFileName,parRstFileName

      INTEGER parOutInit
      INTEGER parOutMesh
      INTEGER parOutStoc,parOutPiec

      INTEGER parQinw
      integer parUinw,parTinw

      REAL(rk) parSuelLam1,parSuelLam2,parSuelE1,parSuelE2

      INTEGER parNOkeys
      CHARACTER(4), POINTER :: parKey(:)
      CHARACTER(255), POINTER :: parKeyDesc(:)

      PARAMETER (parInputFileName="and.inp")
      PARAMETER (parResultsFileName="and.results.vtu")
      PARAMETER (parResultsWallName="and.region.")
      PARAMETER (parResultsForceErr="and.forceErr.")
      PARAMETER (parIntegralsFileName="and.integrals.bin")
      PARAMETER (parStochasticFileName="and.stochastic.txt")
      PARAMETER (parInitialFileName="and.initial.vtu")
      PARAMETER (parLogFileName="and.log")
      PARAMETER (parRstFileName="and.rst")
      PARAMETER (parPPFileName="and.linePP.txt")

      PARAMETER (parLogLun=20)
      PARAMETER (parYes=1,parNo=0)
      PARAMETER (parStokes  = 1407)
      PARAMETER (parLaplace = 1973)

!     Post-processing (values along the line projected onto surface)
      INTEGER parOutProfile
      INTEGER parNumOfProf

!     Domain profile line export
      INTEGER parNoDPLE
      INTEGER parOutDomainProfile
      TYPE(DomainLineType), ALLOCATABLE :: parDomLine(:)

!
! ----------------------------------------------------------------------
!

!
!     Subroutines
!
      CONTAINS
!
! ----------------------------------------------------------------------
!
subroutine defineKeys()

!
!     Keywords
!
  parNOkeys=27
  allocate (parKey(parNOkeys))
  allocate (parKeyDesc(parNOkeys))
  parKey(1)="MDIR"
  parKeyDesc(1)="parMeshDir - directory where *.msh and *.bic are located."
  parKey(2)="LGEO"
  parKeyDesc(2)="parMeshFileName - name of the mesh file (*.msh)."
  parKey(3)="LBIC"
  parKeyDesc(3)="parBiCFileName - name of the boundary and initial conditions file (*.bic)."
  parKey(4)="ITRI"
  parKeyDesc(4)="parTriInteg, parTriRecur - triangle integration (1,3,5), no. of recursions."
  parKey(5)="IQUA"
  parKeyDesc(5)="parQuadIntegRegu,parQuadIntegSing - qudarilateral integration (reg: 1..8), (sing:1..8)."
  parKey(6)="MTRO"
  parKeyDesc(6)="parMTrot(1),parMTrot(2),parMTrot(3) - mesh rotation"
  parKey(7)="MTST"
  parKeyDesc(7)="parMTstr(1),parMTstr(2),parMTstr(3) - mesh stretch"
  parKey(8)="MTTR"
  parKeyDesc(8)="parMTtra(1),parMTtra(2),parMTtra(3) - mesh translate"
  parKey(9)="MTCY"
  parKeyDesc(9)="6 REAL(rk)s - cylinder mesh transform (TES project)."
  parKey(10)="SCRO"
  parKeyDesc(10)="Yes/No - Output info on screeen."
  parKey(11)="OINI"
  parKeyDesc(11)="Yes/No - output initial distribution of function and flux (and.initial.vtu)."
  parKey(12)="OMSH"
  parKeyDesc(12)="Yes/No - output results in Paraview format (and.results.vtu)."
  parKey(13)="OSTO"
  parKeyDesc(13)="Yes/No - output results for stochastic postprocessing (and.stochastic.txt)."
  parKey(14)="WINT"
  parKeyDesc(14)="Yes/No - Write integrals binary file (and.integrals.bin)."
  parKey(15)="OREG"
  parKeyDesc(15)="Yes/No - output results for each region in Paraview format (and.region.pvd)."
  parKey(16)="PPLI"
  parKeyDesc(16)=" - profile line export"
  parKey(17)="REST"
  parKeyDesc(17)="Yes/No - output binary restart file (and.rst)"
  parKey(18)="POST"
  parKeyDesc(18)="Yes/No - postprocessing only run."
  parKey(19)="DPLE"
  parKeyDesc(19)=" - domain profile line export"
  parKey(20)="COME"
  parKeyDesc(20)=" - conserve memory (Yes/No), only works with WINT Yes"
  parKey(21)="PRTY"
  parKeyDesc(21)=" - problem type : Laplace (default) or Stokes"
  parKey(22)="QINW"
  parKeyDesc(22)="Yes/No - integral of flux through the walls."
  parKey(23)="UINW"
  parKeyDesc(23)="Yes/No - integral of function through the walls."            
  parKey(24)="SUEL"
  parKeyDesc(23)="a,b,c,e1,e2 - mesh sphere -> superEllipsoid"
  parKey(25)="SVAL"
  parKeyDesc(25)="Yes/No - Run Stokes validation (particle in sphere)."
  parKey(26)="TRII"
  parKeyDesc(26)="lambda1 lambda2 e1 e2 - superellipsoid parameters"      
  parKey(27)="TINW"
  parKeyDesc(27)="Yes/No - integral of torque at the walls."      

end subroutine



!
! ----------------------------------------------------------------------
!
      SUBROUTINE ReadParameters()

      INTEGER lun,Cstring,i
      CHARACTER KeyWord*4,OneLine*255
      CHARACTER dummy*64

      
!
!     Default values
!
      parMeshDir="/home/ravnik/mesh/gmsh/"
      parTriInteg=5
      parTriRecur=20
      parQuadIntegSing=4
      parQuadIntegRegu=4
      parScreenOutput=parYes
      parSval = parNo
      parTrii = parNo

!     Mesh transromation
      pariMTrot=parNo
      pariMTstr=parNo
      pariMTtra=parNo
      pariMTcyl=parNo
      pariMTsel=parNo
      parMTrot=0.0_rk
      parMTstr=0.0_rk
      parMTtra=0.0_rk
      parMTcyl=0.0_rk
      parMTsel=0.0_rk

!     Output control
      parPostOnly = parNo
      parOutInit = parYes
      parOutMesh = parNo
      parOutStoc = parNo
      parOutPiec = parYes
      parWriteIntegrals = parYes
      parOutRest = parYes
      parConserveMemory = parNo

      parMeshFileName="cyl-1.msh"
      parBiCFileName="cyl-prevod.bic"

!     Post processing -line
      parOutProfile = parNo
      parNumOfProf = 0

!     Domain profile line export
      parNoDPLE = 0
      parOutDomainProfile = parNo

      parPrType = parLaplace

      parQinw = parYes
      parUinw = parYes
      parTinw = parNo


      parSuelLam1 = 1.0_rk
      parSuelLam2 = 1.0_rk
      parSuelE1 = 1.0_rk
      parSuelE2 = 1.0_rk


!
!     Open input file
!
      lun=11
      OPEN (lun,FILE=TRIM(parInputFileName),ERR=10,STATUS='OLD')


      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ (OneLine,*) KeyWord
!
!***    GEO AND BC FILE LOCATION :
!
        IF (KeyWord.EQ.parKey(1)) THEN ! MDIR
          READ(OneLine,'(A5,A)') dummy,parMeshDir

        ELSE IF (KeyWord.EQ.parKey(2)) THEN ! LGEO
          READ(OneLine,'(A5,A)') dummy,parMeshFileName

        ELSE IF (KeyWord.EQ.parKey(3)) THEN ! LBIC
          READ(OneLine,'(A5,A)') dummy,parBiCFileName

        ELSE IF (KeyWord.EQ.parKey(4)) THEN ! ITRI
          READ(OneLine,*) dummy,parTriInteg,parTriRecur

        ELSE IF (KeyWord.EQ.parKey(5)) THEN ! IQUA
          READ(OneLine,*) dummy,parQuadIntegRegu,parQuadIntegSing

        ELSE IF (KeyWord.EQ.parKey(6)) THEN ! MTRO
          READ(OneLine,*) dummy,parMTrot(1),parMTrot(2),parMTrot(3)
          pariMTrot=parYes
        ELSE IF (KeyWord.EQ.parKey(7)) THEN ! MTST
          READ(OneLine,*) dummy,parMTstr(1),parMTstr(2),parMTstr(3)
          pariMTstr=parYes
        ELSE IF (KeyWord.EQ.parKey(8)) THEN ! MTTR
          READ(OneLine,*) dummy,parMTtra(1),parMTtra(2),parMTtra(3)
          pariMTtra=parYes
        ELSE IF (KeyWord.EQ.parKey(9)) THEN ! MTCY
          READ(OneLine,*) dummy,parMTcyl(1),parMTcyl(2),parMTcyl(3),parMTcyl(4),parMTcyl(5),parMTcyl(6)
          parMTcyl(1)=parMTcyl(1)*0.5_rk ! transform to radius
          parMTcyl(2)=parMTcyl(2)*0.5_rk
          parMTcyl(4)=parMTcyl(4)*0.5_rk
          parMTcyl(5)=parMTcyl(5)*0.5_rk
          pariMTcyl=parYes
        ELSE IF (KeyWord.EQ.parKey(10)) THEN ! SCRO
          READ(OneLine,*) dummy,dummy
          parScreenOutput=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(11)) THEN ! OINI
          READ(OneLine,*) dummy,dummy
          parOutInit=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(12)) THEN ! OMSH
          READ(OneLine,*) dummy,dummy
          parOutMesh=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(13)) THEN ! OSTO
          READ(OneLine,*) dummy,dummy
          parOutStoc=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(14)) THEN ! WINT
          READ(OneLine,*) dummy,dummy
          parWriteIntegrals=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(15)) THEN ! OREG
          READ(OneLine,*) dummy,dummy
          parOutPiec=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(16)) THEN ! PPLI
          READ(OneLine,*) dummy, parNumOfProf
          IF (parNumOfProf.GT.0) parOutProfile = parYes
        ELSE IF (KeyWord.EQ.parKey(17)) THEN ! REST
          READ(OneLine,*) dummy,dummy
          parOutRest=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(18)) THEN ! POST
          READ(OneLine,*) dummy,dummy
          parPostOnly=Cstring(dummy,"YES")

        ELSE IF (KeyWord.EQ.parKey(19)) THEN ! DPLE
          READ(OneLine,*) dummy,parNoDPLE
          IF (parNoDPLE.GT.0) THEN
            parOutDomainProfile = parYes
            ALLOCATE (parDomLine(parNoDPLE))
            DO i=1,parNoDPLE
              CALL rOneTL(lun,OneLine)
              READ(OneLine,*)      parDomLine(i)%name &
                                  ,parDomLine(i)%startPoint(1),parDomLine(i)%startPoint(2),parDomLine(i)%startPoint(3) &
                                  ,parDomLine(i)%endPoint(1),parDomLine(i)%endPoint(2),parDomLine(i)%endPoint(3) &
                                  ,parDomLine(i)%nno
            END DO
          END IF

        ELSE IF (KeyWord.EQ.parKey(20)) THEN ! COME
          READ(OneLine,*) dummy,dummy
          parConserveMemory=Cstring(dummy,"YES")

        ELSE IF (KeyWord.EQ.parKey(21)) THEN ! COME
          READ(OneLine,*) dummy,dummy
          if (Cstring(dummy,"STOKES").EQ.parYes) parPrType = parStokes
          if (Cstring(dummy,"LAPLACE").EQ.parYes) parPrType = parLaplace
        ELSE IF (KeyWord.EQ.parKey(22)) THEN ! QINW
          READ(OneLine,*) dummy,dummy
          parQinw=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(23)) THEN ! UINTW
          READ(OneLine,*) dummy,dummy
          parUinw=Cstring(dummy,"YES")
        ELSE IF (KeyWord.EQ.parKey(24)) THEN ! SUEL a,b,c,e1,e2
          READ(OneLine,*) dummy,parMTsel(1),parMTsel(2),parMTsel(3),parMTsel(4),parMTsel(5)
          pariMTsel=parYes          
        ELSE IF (KeyWord.EQ.parKey(25)) THEN ! SVAL
          READ(OneLine,*) dummy,dummy
          parSval=Cstring(dummy,"YES")          
        ELSE IF (KeyWord.EQ.parKey(26)) THEN ! TRII
          READ(OneLine,*) dummy,parSuelLam1,parSuelLam2,parSuelE1,parSuelE2
          parTrii = parYes
        ELSE IF (KeyWord.EQ.parKey(27)) THEN ! TINW
          READ(OneLine,*) dummy,dummy
          parTinw=Cstring(dummy,"YES")          
        END IF

        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)
!
!     Some postprocessing
!
      WRITE (parMeshFullName,'(A,A)') TRIM(parMeshDir),TRIM(parMeshFileName)
      WRITE (parBiCFullName,'(A,A)') TRIM(parMeshDir),TRIM(parBiCFileName)

      RETURN

10    CONTINUE ! error when opening input file


      CALL WriteToLog("Could not open input file!")
      CALL ShowHelp()

      STOP


      END SUBROUTINE
!
! ----------------------------------------------------------------------
!
      SUBROUTINE ShowHelp()

      IMPLICIT NONE

      INTEGER i

      WRITE (*,'(/A)') "Structure of the 'and.inp' file: 4 letter keyword followed by parameters."
      WRITE (*,'(/A)') "Obligatory keywords:"

      DO i=1,3
        WRITE (*,'(A4,1X,A)') parKey(i),TRIM(parKeyDesc(i))
      END DO

      WRITE (*,'(/A)') "Optional keywords:"
      DO i=4,parNOkeys
        WRITE (*,'(A4,1X,A)') parKey(i),TRIM(parKeyDesc(i))
      END DO


      END SUBROUTINE

!
! ----------------------------------------------------------------------
!
subroutine readCLargs()

  if (COMMAND_ARGUMENT_COUNT().NE.0) THEN
    call ShowHelp()
    STOP
  end if

!
!     Try to open input file
!
  OPEN (11,FILE=TRIM(parInputFileName),ERR=10,STATUS='OLD')
  CLOSE (11)

  RETURN

10 CONTINUE
  write(*,*) "Could not open the input file!"
  call ShowHelp()
  STOP


end subroutine



!
! ----------------------------------------------------------------------
!
      SUBROUTINE WriteToLog(tekst)

      CHARACTER tekst*(*)
      CHARACTER*6 cas

      CALL logGetTime(cas)
      WRITE (parLogLun,'(A6,1X,A)') cas,TRIM(tekst)
      CALL FLUSH(parLogLun)

      IF (parScreenOutput.EQ.parYes) THEN
        WRITE (6,'(A)') TRIM(tekst)
      END IF

      END SUBROUTINE

!
! ----------------------------------------------------------------------
!
!
!     Set up log file
!

      SUBROUTINE SetUpLogFile()

      CHARACTER(255) tekst

      OPEN (parLogLun,FILE=TRIM(parLogFileName),ERR=10,STATUS='UNKNOWN')

      WRITE (tekst,'(6A)') trim(parIDname)," ",trim(parIDversion),", ",trim(parIDdate),"."
      CALL WriteToLog(tekst)
      WRITE (tekst,'(A)') parStartTime
      CALL WriteToLog(tekst)
      WRITE (tekst,'(A,A)') "Running on: ",TRIM(ParHostname)
      CALL WriteToLog(tekst)


      RETURN

10    CONTINUE ! error when opening input file

      WRITE (*,*) "Could not open LOG file!!"
      STOP


      END SUBROUTINE

! -----------------------------------------------------------------------------
      SUBROUTINE logGetTime(cas)
!
!     $: writes date and time to cas
!
! -----------------------------------------------------------------------------
      CHARACTER*50 D,T,Z
      INTEGER Value(8)
      CHARACTER*2 ura,min
      CHARACTER*6 cas

      CALL DATE_AND_TIME(D,T,Z,Value)

      IF (Value(5).LT.10) THEN
        WRITE (ura,'(A,I1)') "0",Value(5)
      else
        WRITE (ura,'(I2)') Value(5)
      end if

      IF (Value(6).LT.10) THEN
        WRITE (min,'(A,I1)') "0",Value(6)
      else
        WRITE (min,'(I2)') Value(6)
      end if

      WRITE (cas,'(A2,A1,A2,1X)') ura,':',min

      END SUBROUTINE

!
! ----------------------------------------------------------------------
!
!
!     Close log file
!

      SUBROUTINE StopProgram()

      CHARACTER(255) tekst

      CALL DatumInUra(parEndTime)
      WRITE (tekst,'(A,A)') "START:",TRIM(parStartTime)
      CALL WriteToLog(tekst)
      WRITE (tekst,'(A,A)') "END  :",TRIM(parEndTime)
      CALL WriteToLog(tekst)

      CLOSE (parLogLun)

      STOP

      END SUBROUTINE


      END MODULE

! -----------------------------------------------------------------------------

      SUBROUTINE OutputPiecesParaview()

!
!     $: Each wall in separate file
!
! -----------------------------------------------------------------------------
      USE mMesh
      USE mEqns
      USE mPar
      USE mCommon
      IMPLICIT NONE

      INTEGER itp,i,of,k,j
      CHARACTER*255 filename,tmp,vrstica

      INTEGER iWall,noe,non
      INTEGER, ALLOCATABLE :: nlist(:)


      itp=96

      WRITE (filename,'(2A)') TRIM(parResultsWallName),"pvd"
      OPEN (itp+1,FILE=TRIM(filename),STATUS='UNKNOWN')

      WRITE (itp+1,'(A)') '<?xml version="1.0"?>'
      WRITE (itp+1,'(A)') '<VTKFile type="Collection" version="0.1"'
      WRITE (itp+1,'(A)') '         byte_order="LittleEndian"'
      WRITE (itp+1,'(A)') '         compressor="vtkZLibDataCompressor">'
      WRITE (itp+1,'(A)') '  <Collection>'


      DO iWall = 1,nofw

      WRITE (filename,'(3A)') TRIM(parResultsWallName),TRIM(wall(iWall)%name),".vtu"
      WRITE (itp+1,'(3A)') '<DataSet timestep="0" group="" part="0" file="',TRIM(filename),'"/>'
      OPEN (itp,FILE=TRIM(filename),STATUS='UNKNOWN')

!     count number of nodes in wall

      noe = 0
      ALLOCATE (nlist(nnodes))
      nlist=0
      DO  i=1,nelem
        IF (element(i)%bcid.EQ.iWall) THEN
          noe = noe + 1
          DO j=1,element(i)%nno
            nlist(element(i)%con(j))=1
          END DO
        END IF
      END DO
      non=SUM(nlist)

      j=0
      DO  i=1,nnodes
        IF (nlist(i).EQ.1) THEN
          j=j+1
          nlist(i)=j
        END IF
      END DO





      WRITE (itp,'(A)') '<?xml version="1.0"?>'
      WRITE (itp,'(A)')&
     '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      WRITE (itp,'(A)') '<UnstructuredGrid>'
      WRITE (itp,'(A,I10,A,I10,A)') '<Piece NumberOfPoints="',non,'" NumberOfCells="',noe,'">'

      WRITE (tmp,'("(",I1,"A)")') neq+2
      WRITE (itp,tmp) '<PointData Scalars="node',(',u'//TRIM(eqn(i)%name),i=1,neq),'">'
!
!     PODATKI V VOZLISCIH
!
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="node" format="ascii">'
      DO  i=1,nnodes
        IF (nlist(i).NE.0) WRITE (itp,*) nlist(i)
      END DO
      WRITE (itp,'(A)') '</DataArray>'

      DO j=1,neq
        WRITE (itp,'(3A)') '<DataArray type="Float32" Name="',"u"//TRIM(eqn(j)%name),'" format="ascii">'
        DO  i=1,nnodes
          IF (nlist(i).NE.0) WRITE (itp,*) eqn(j)%u(i)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END DO

      WRITE (itp,'(A)') '</PointData>'
!
!     PODATKI NA ELEMENTE
      WRITE (tmp,'("(",I1,"A)")') neq+2

      IF (parPrType .EQ. parStokes) THEN
        WRITE (itp,tmp) '<CellData Vectors="normale" Scalars="pressure,elemBCid',(',q'//TRIM(eqn(i)%name),i=1,neq),'">'
      END IF
      IF (parPrType .EQ. parLaplace) THEN
        WRITE (itp,tmp) '<CellData Vectors="normale" Scalars="elemBCid',(',q'//TRIM(eqn(i)%name),i=1,neq),'">'
      END IF
!
            
!      NORMALE
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="normale" NumberOfComponents="3" format="ascii">'
      DO  i=1,nelem
        IF (element(i)%bcid.EQ.iWall)  THEN
          WRITE (vrstica,'(3G18.9)') element(i)%normal(1),element(i)%normal(2),element(i)%normal(3)
          CALL sqblnk(itp,vrstica)
        END IF
      END DO
      WRITE (itp,'(A)') '</DataArray>'

!     pressure
      IF (parPrType .EQ. parStokes) THEN
        WRITE (itp,'(A)') '<DataArray type="Float32" Name="pressure" format="ascii">'
        DO  i=1,nelem
          IF (element(i)%bcid.EQ.iWall)  THEN
            WRITE (itp,*) eqn(1)%p(i)
          END IF
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END IF

!     element BCid
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="elemBCid" format="ascii">'
      DO  i=1,nelem
        IF (element(i)%bcid.EQ.iWall)  THEN
          WRITE (itp,*) element(i)%bcid
        END IF
      END DO
      WRITE (itp,'(A)') '</DataArray>'

!     Fluksi
      DO j=1,neq
        WRITE (itp,'(3A)') '<DataArray type="Float32" Name="',"q"//TRIM(eqn(j)%name),'" format="ascii">'
        DO  i=1,nqnodes
          IF (element(i)%bcid.EQ.iWall) WRITE (itp,*) eqn(j)%q(i)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END DO


      WRITE (itp,'(A)') '</CellData>'


!      KOORDINATE VOZLISC
      WRITE (itp,'(A)') '<Points>'
      WRITE (itp,'(A)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'

      DO  i=1,nnodes
         IF (nlist(i).NE.0) THEN
           WRITE (vrstica,'(3G18.9)') node(i)%x(1),node(i)%x(2),node(i)%x(3)
           CALL sqblnk(itp,vrstica)
         END IF
      END DO

      WRITE (itp,'(A)') '</DataArray>'
      WRITE (itp,'(A)') '</Points>'

      WRITE (itp,'(A)') '<Cells>'
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="connectivity" format="ascii">'

        DO i=1,nelem
          IF (element(i)%bcid.EQ.iWall) THEN
            WRITE(vrstica,*) (nlist(element(i)%con(k))-1,k=1,element(i)%nno)
            CALL sqblnk(itp,vrstica)
          END IF
        END DO


        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="offsets" format="ascii">'
          of = 0
          DO i=1,nelem
            IF (element(i)%bcid.EQ.iWall) THEN
              of = of + element(i)%nno
              WRITE (itp,*) of
            END IF
          END DO
        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
          DO i=1,nelem
          IF (element(i)%bcid.EQ.iWall) THEN
            IF (element(i)%type.EQ.2) THEN ! three node triangle
              WRITE (itp,*) "5" ! 5 tikotnik, 9 quuad, 12 heksaeder
            ELSE IF (element(i)%type.EQ.3) THEN ! four node Quadrangle
             WRITE (itp,*) "9" ! 5 tikotnik, 9 quuad, 12 heksaeder
            ELSE
              CALL WriteToLog("Error :: Element type not supported!")
            END IF
          END IF
          END DO
        WRITE (itp,'(A)') '</DataArray>'
      WRITE (itp,'(A)') '</Cells>'


      WRITE (itp,'(A)') '</Piece>'
      WRITE (itp,'(A)') '</UnstructuredGrid>'
      WRITE (itp,'(A)') '</VTKFile>'

      CLOSE (itp)

      DEALLOCATE (nlist)

      END DO

      WRITE (itp+1,'(A)') '  </Collection>'
      WRITE (itp+1,'(A)') '</VTKFile>'
      CLOSE(itp+1)

      END


! -----------------------------------------------------------------------------

      SUBROUTINE OutputMeshParaview()

!
!     $: Outputs function
!
! -----------------------------------------------------------------------------
      USE mMesh
      USE mEqns
      USE mPar
      IMPLICIT NONE

      INTEGER itp,i,of,k,j
!      CHARACTER*6 cifra
      CHARACTER*255 vrstica,tmp


      itp=96

      OPEN (itp,FILE=TRIM(parResultsFileName),STATUS='UNKNOWN')


      WRITE (itp,'(A)') '<?xml version="1.0"?>'
      WRITE (itp,'(A)')&
     '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      WRITE (itp,'(A)') '<UnstructuredGrid>'
      WRITE (itp,'(A,I10,A,I10,A)') '<Piece NumberOfPoints="',nnodes,'" NumberOfCells="',nelem,'">'


      IF (parPrType .EQ. parStokes) THEN
        WRITE (tmp,'("(",I1,"A)")') neq+2
        WRITE (itp,tmp) '<PointData Vectors="u" Scalars="node',(",u" // TRIM(eqn(i)%name),i=1,neq),'">'
      END IF
      IF (parPrType .EQ. parLaplace) THEN
        WRITE (tmp,'("(",I1,"A)")') neq+2
        WRITE (itp,tmp) '<PointData Scalars="node',(",u" // TRIM(eqn(i)%name),i=1,neq),'">'
      END IF
!
!     PODATKI V VOZLISCIH
!
      IF (parPrType .EQ. parStokes) THEN
        WRITE (itp,'(A)') '<DataArray type="Float32" Name="u" NumberOfComponents="3" format="ascii">'
        DO  i=1,nnodes
          WRITE (vrstica,'(3G18.9)') eqn(1)%u(i),eqn(2)%u(i),eqn(3)%u(i)
          CALL sqblnk(itp,vrstica)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END IF

      WRITE (itp,'(A)') '<DataArray type="Float32" Name="node" format="ascii">'
      DO  i=1,nnodes
        WRITE (itp,*) i
      END DO
      WRITE (itp,'(A)') '</DataArray>'

      DO j=1,neq
        WRITE (itp,'(3A)') '<DataArray type="Float32" Name="',"u"//TRIM(eqn(j)%name),'" format="ascii">'
        DO  i=1,nnodes
          WRITE (itp,*) eqn(j)%u(i)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END DO

      WRITE (itp,'(A)') '</PointData>'
!
!     PODATKI NA ELEMENTE
!
      WRITE (tmp,'("(",I1,"A)")') neq+2


      IF (parPrType .EQ. parStokes) THEN
        WRITE (itp,tmp) '<CellData Vectors="normale" Scalars="pressure,elemBCid',(',q'//TRIM(eqn(i)%name),i=1,neq),'">'
      END IF
      IF (parPrType .EQ. parLaplace) THEN
        WRITE (itp,tmp) '<CellData Vectors="normale" Scalars="elemBCid',(',q'//TRIM(eqn(i)%name),i=1,neq),'">'
      END IF

!      NORMALE
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="normale" NumberOfComponents="3" format="ascii">'
      DO  i=1,nelem
        WRITE (vrstica,'(3G18.9)') element(i)%normal(1),element(i)%normal(2),element(i)%normal(3)
        CALL sqblnk(itp,vrstica)
      END DO
      WRITE (itp,'(A)') '</DataArray>'

!     pressure
      IF (parPrType .EQ. parStokes) THEN
        WRITE (itp,'(A)') '<DataArray type="Float32" Name="pressure" format="ascii">'
        DO  i=1,nelem
          WRITE (itp,*) eqn(1)%p(i)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END IF

!     element BCid
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="elemBCid" format="ascii">'
      DO  i=1,nelem
        WRITE (itp,*) element(i)%bcid
      END DO
      WRITE (itp,'(A)') '</DataArray>'


!     Fluksi
      DO j=1,neq
        WRITE (itp,'(3A)') '<DataArray type="Float32" Name="',"q"//TRIM(eqn(j)%name),'" format="ascii">'
        DO  i=1,nqnodes
          WRITE (itp,*) eqn(j)%q(i)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END DO





      WRITE (itp,'(A)') '</CellData>'


!      KOORDINATE VOZLISC
      WRITE (itp,'(A)') '<Points>'
      WRITE (itp,'(A)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'

      DO  i=1,nnodes
        WRITE (vrstica,'(3G18.9)') node(i)%x(1),node(i)%x(2),node(i)%x(3)
        CALL sqblnk(itp,vrstica)
      END DO

      WRITE (itp,'(A)') '</DataArray>'
      WRITE (itp,'(A)') '</Points>'

      WRITE (itp,'(A)') '<Cells>'
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="connectivity" format="ascii">'

        DO i=1,nelem
          WRITE(vrstica,*) (element(i)%con(k)-1,k=1,element(i)%nno)
          CALL sqblnk(itp,vrstica)
        END DO


        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="offsets" format="ascii">'
          of = 0
          DO i=1,nelem
            of = of + element(i)%nno
            WRITE (itp,*) of
          END DO
        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
          DO i=1,nelem
            IF (element(i)%type.EQ.2) THEN ! three node triangle
              WRITE (itp,*) "5" ! 5 tikotnik, 9 quuad, 12 heksaeder
            ELSE IF (element(i)%type.EQ.3) THEN ! four node Quadrangle
             WRITE (itp,*) "9" ! 5 tikotnik, 9 quuad, 12 heksaeder
            ELSE
              CALL WriteToLog("Error :: Element type not supported!")
            END IF
          END DO
        WRITE (itp,'(A)') '</DataArray>'
      WRITE (itp,'(A)') '</Cells>'


      WRITE (itp,'(A)') '</Piece>'
      WRITE (itp,'(A)') '</UnstructuredGrid>'
      WRITE (itp,'(A)') '</VTKFile>'

      CLOSE (itp)


      END




! -----------------------------------------------------------------------------

      SUBROUTINE OutputStochastic()

!
!     $: Outputs function
!
! -----------------------------------------------------------------------------
      USE mMesh
      USE mEqns
      USE mPar
      IMPLICIT NONE

      INTEGER itp,i,j

      itp=96

      OPEN (itp,FILE=TRIM(parStochasticFileName),STATUS='UNKNOWN')

      WRITE (itp,'(A)') '#'
      WRITE (itp,'(A)') '# Results for stochastic postprocessing'
      WRITE (itp,'(A)') '#'
      WRITE (itp,'(A,I0)') '# total number of datapoints, ',neq*(nnodes+nqnodes)
      WRITE (itp,'(A)') '#'
      WRITE (itp,'(A,I0)') '# potential, ',neq*nnodes
      WRITE (itp,'(A)') '#'
      DO j=1,neq
        DO  i=1,nnodes
          WRITE (itp,*) eqn(j)%u(i)
        END DO
      END DO
      WRITE (itp,'(A)') '#'
      WRITE (itp,'(A,I0)') '# flux, ',neq*nqnodes
      WRITE (itp,'(A)') '#'
      DO j=1,neq
        DO  i=1,nqnodes
          WRITE (itp,*) eqn(j)%q(i)
        END DO
      END DO

      CLOSE (itp)


      END


! -----------------------------------------------------------------------------

  SUBROUTINE OutputInitialParaview()

!
!     $: Outputs function
!
! -----------------------------------------------------------------------------
    USE mMesh
    USE mEqns
    USE mPar
    USE parallel
    IMPLICIT NONE
    INTEGER itp,i,of,k,j
    CHARACTER*255 vrstica,tmp

    IF (amIroot) THEN

      itp=96

      OPEN (itp,FILE=TRIM(parInitialFileName),STATUS='UNKNOWN')


      WRITE (itp,'(A)') '<?xml version="1.0"?>'
      WRITE (itp,'(A)')&
     '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      WRITE (itp,'(A)') '<UnstructuredGrid>'
      WRITE (itp,'(A,I10,A,I10,A)') '<Piece NumberOfPoints="',nnodes,'" NumberOfCells="',nelem,'">'

      WRITE (tmp,'("(",I1,"A)")') neq+1
      WRITE (itp,tmp) '<PointData Scalars="node',(",u" // TRIM(eqn(i)%name),i=1,neq),'">'

!
!     PODATKI V VOZLISCIH
!
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="node" format="ascii">'
      DO  i=1,nnodes
        WRITE (itp,*) i
      END DO
      WRITE (itp,'(A)') '</DataArray>'

      DO j=1,neq
        WRITE (itp,'(3A)') '<DataArray type="Float32" Name="',"u"//TRIM(eqn(j)%name),'" format="ascii">'
        DO  i=1,nnodes
          WRITE (itp,*) eqn(j)%u(i)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END DO

      WRITE (itp,'(A)') '</PointData>'
!
!     PODATKI NA ELEMENTE
!
      WRITE (tmp,'("(",I1,"A)")') neq+3
      WRITE (itp,'(A)') '<CellData Vectors="normale" Scalars="',(",q,elemBCid" // TRIM(eqn(i)%name),i=1,neq),'">'


!      NORMALE
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="normale" NumberOfComponents="3" format="ascii">'
      DO  i=1,nelem
        WRITE (vrstica,'(3G18.9)') element(i)%normal(1),element(i)%normal(2),element(i)%normal(3)
        CALL sqblnk(itp,vrstica)
      END DO
      WRITE (itp,'(A)') '</DataArray>'

!     Fluksi
      DO j=1,neq
        WRITE (itp,'(3A)') '<DataArray type="Float32" Name="',"q"//TRIM(eqn(j)%name),'" format="ascii">'
        DO  i=1,nqnodes
          WRITE (itp,*) eqn(1)%q(i)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
      END DO

!     element BCid
      WRITE (itp,'(A)') '<DataArray type="Float32" Name="elemBCid" format="ascii">'
      DO  i=1,nelem
        WRITE (itp,*) element(i)%bcid
      END DO
      WRITE (itp,'(A)') '</DataArray>'

      WRITE (itp,'(A)') '</CellData>'


!      KOORDINATE VOZLISC
      WRITE (itp,'(A)') '<Points>'
      WRITE (itp,'(A)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'

      DO  i=1,nnodes
        WRITE (vrstica,'(3G18.9)') node(i)%x(1),node(i)%x(2),node(i)%x(3)
        CALL sqblnk(itp,vrstica)
      END DO

      WRITE (itp,'(A)') '</DataArray>'
      WRITE (itp,'(A)') '</Points>'

      WRITE (itp,'(A)') '<Cells>'
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="connectivity" format="ascii">'

        DO i=1,nelem
          WRITE(vrstica,*) (element(i)%con(k)-1,k=1,element(i)%nno)
          CALL sqblnk(itp,vrstica)
        END DO


        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="offsets" format="ascii">'
          of = 0
          DO i=1,nelem
            of = of + element(i)%nno
            WRITE (itp,*) of
          END DO
        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
          DO i=1,nelem
            IF (element(i)%type.EQ.2) THEN ! three node triangle
              WRITE (itp,*) "5" ! 5 tikotnik, 9 quuad, 12 heksaeder
            ELSE IF (element(i)%type.EQ.3) THEN ! four node Quadrangle
             WRITE (itp,*) "9" ! 5 tikotnik, 9 quuad, 12 heksaeder
            ELSE
              CALL WriteToLog("Error :: Element type not supported!")
            END IF
          END DO
        WRITE (itp,'(A)') '</DataArray>'
      WRITE (itp,'(A)') '</Cells>'


      WRITE (itp,'(A)') '</Piece>'
      WRITE (itp,'(A)') '</UnstructuredGrid>'
      WRITE (itp,'(A)') '</VTKFile>'

      CLOSE (itp)

    END IF

  END


!
! ----------------------------------------------------------------------
!
      SUBROUTINE CalQMesh()
!
! ----------------------------------------------------------------------
      USE mMesh
      USE mPar
      USE mCommon
      IMPLICIT NONE

      INTEGER i,j,k
      REAL(rk) ETA1M,ETA1P,ETA2M,ETA2P,FIG4(4)

      nqnodes = nelem ! constant interpolation for flux

      ALLOCATE (qnode(nqnodes))
!
!     Calculate element center
!
      DO i=1,nelem

        IF (element(i)%type.EQ.2) THEN ! three node triangle

!         Barycentric coordinates = 1/3
          FIG4(1)=1.0_rk/3.0_rk
          FIG4(2)=1.0_rk/3.0_rk
          FIG4(3)=1.0_rk/3.0_rk

        ELSE IF (element(i)%type.EQ.3) THEN ! four node Quadrangle

          ETA1M=1.0_rk-0.0_rk
          ETA1P=1.0_rk+0.0_rk
          ETA2M=1.0_rk-0.0_rk
          ETA2P=1.0_rk+0.0_rk
!
          FIG4(1)=0.25_rk*ETA1M*ETA2M
          FIG4(2)=0.25_rk*ETA1P*ETA2M
          FIG4(3)=0.25_rk*ETA1P*ETA2P
          FIG4(4)=0.25_rk*ETA1M*ETA2P
!
        ELSE
           CALL WriteToLog("Error :: Element type not supported!")
        END IF

        DO j=1,3
          qnode(i)%x(j)=0.0_rk
          DO k=1,element(i)%nno
            qnode(i)%x(j)=qnode(i)%x(j)+FIG4(k)*node(element(i)%con(k))%x(j)
          END DO
        END DO

      END DO


      END SUBROUTINE


!
! ----------------------------------------------------------------------
!
SUBROUTINE InterpolateToQMesh(u,uq)
!
! ----------------------------------------------------------------------
      USE mMesh
      USE mPar
      USE mCommon
      USE mCommon
      IMPLICIT NONE
        
      INTEGER i,k
      REAL(rk) ETA1M,ETA1P,ETA2M,ETA2P,FIG4(4)
      REAL(rk) u(nnodes),uq(nnodes)

      nqnodes = nelem ! constant interpolation for flux
!
!     Interpolate to element center
!
      DO i=1,nelem

        IF (element(i)%type.EQ.2) THEN ! three node triangle

!         Barycentric coordinates = 1/3
          FIG4(1)=1.0_rk/3.0_rk
          FIG4(2)=1.0_rk/3.0_rk
          FIG4(3)=1.0_rk/3.0_rk

        ELSE IF (element(i)%type.EQ.3) THEN ! four node Quadrangle

          ETA1M=1.0_rk-0.0_rk
          ETA1P=1.0_rk+0.0_rk
          ETA2M=1.0_rk-0.0_rk
          ETA2P=1.0_rk+0.0_rk
!
          FIG4(1)=0.25_rk*ETA1M*ETA2M
          FIG4(2)=0.25_rk*ETA1P*ETA2M
          FIG4(3)=0.25_rk*ETA1P*ETA2P
          FIG4(4)=0.25_rk*ETA1M*ETA2P
!
        ELSE
           CALL WriteToLog("Error :: Element type not supported!")
        END IF

        
          uq(i)=0.0_rk
          DO k=1,element(i)%nno
            uq(i)=uq(i)+FIG4(k)*u(element(i)%con(k))
          END DO
        

      END DO


END SUBROUTINE
        


!
!
! ----------------------------------------------------------------------
!
      SUBROUTINE CalMeshNormals()
!
! ----------------------------------------------------------------------
      USE mMesh

      INTEGER i

!
!     Calculate element normal
!
      DO i=1,nelem
        CALL CalElementNormal(element(i))
      END DO

!      stop

      END SUBROUTINE

!
! -----------------------------------------------------------------------------------------
!
      SUBROUTINE CalElementNormal(e)

      USE mMesh
      USE mCommon
      IMPLICIT NONE

      TYPE(ElementType) e
      INTEGER i
      REAL(rk) t1x,t1y,t1z
      REAL(rk) t2x,t2y,t2z
      REAL(rk) t3x,t3y,t3z
      REAL(rk) v1x,v1y,v1z
      REAL(rk) v2x,v2y,v2z
      REAL(rk) d

        t1x=node(e%con(1))%x(1)   !p%x(p%ibc(i,1),1)
        t1y=node(e%con(1))%x(2)   !p%x(p%ibc(i,1),2)
        t1z=node(e%con(1))%x(3)   !p%x(p%ibc(i,1),3)

        t2x=node(e%con(2))%x(1)   !p%x(p%ibc(i,2),1)
        t2y=node(e%con(2))%x(2)   !p%x(p%ibc(i,2),2)
        t2z=node(e%con(2))%x(3)   !p%x(p%ibc(i,2),3)

        t3x=node(e%con(3))%x(1)   !p%x(p%ibc(i,3),1)
        t3y=node(e%con(3))%x(2)   !p%x(p%ibc(i,3),2)
        t3z=node(e%con(3))%x(3)   !p%x(p%ibc(i,3),3)

!       two vectors

        v1x=t2x-t1x
        v1y=t2y-t1y
        v1z=t2z-t1z

        v2x=t3x-t1x
        v2y=t3y-t1y
        v2z=t3z-t1z

!       cross product

        e%normal(1)= v1y*v2z-v1z*v2y
        e%normal(2)=-v1x*v2z+v1z*v2x
        e%normal(3)= v1x*v2y-v1y*v2x

!       normalization


        d=SQRT(e%normal(1)**2+e%normal(2)**2+e%normal(3)**2)
        DO i=1,3
          e%normal(i)=e%normal(i)/d
        END DO



      END


!
!
! ----------------------------------------------------------------------
!
SUBROUTINE VerifyMeshNormals()
!
! ----------------------------------------------------------------------
  USE mMesh
  INTEGER i
  
  !
  ! Verify element normal
  !
  DO i=1,nelem
    CALL VerifyElementNormal(element(i))
  END DO
  
END SUBROUTINE
        
!
! -----------------------------------------------------------------------------------------
!
SUBROUTINE VerifyElementNormal(e)

  USE mMesh
  USE mCommon
  IMPLICIT NONE
  TYPE(ElementType) e

  REAL(rk) t1x,t1y,t1z
  REAL(rk) t2x,t2y,t2z
  REAL(rk) t3x,t3y,t3z

  REAL(rk) cp(3), cpmag, dp

  t1x=node(e%con(1))%x(1)  
  t1y=node(e%con(1))%x(2)  
  t1z=node(e%con(1))%x(3)  
  t2x=node(e%con(2))%x(1)  
  t2y=node(e%con(2))%x(2)  
  t2z=node(e%con(2))%x(3)  
  t3x=node(e%con(3))%x(1)  
  t3y=node(e%con(3))%x(2)  
  t3z=node(e%con(3))%x(3)  
       

        
! check normal orientation
  
! facee center point
  cp(1) = (1.0D0/3.0D0) * ( t1x + t2x + t3x )
  cp(2) = (1.0D0/3.0D0) * ( t1y + t2y + t3y )
  cp(3) = (1.0D0/3.0D0) * ( t1z + t2z + t3z )
  !magnitude of cp vec
  cpmag=SQRT(cp(1)**2+cp(2)**2+cp(3)**2)
  !dot product of a vector from face center
  ! point to origin (-cp) and a normal vector (e%normal)
  CALL DotProduct(-cp,e%normal,dp)
  
  !cp mag check is so that the outer sphere boundary does not
  !interfere with particle boundary
  !logic is: particle size is << 10, boundary sphere size is >> 10)
  IF ( dp .LT. 0.0_rk .AND. cpmag .LT. 10.0_rk ) THEN
    PRINT *, "Normal wrong, turning around"
    stop
  !  e%normal(1) = -e%normal(1)
  !  e%normal(2) = -e%normal(2)
  !  e%normal(3) = -e%normal(3)
  END IF
  
END
        
        
        

!
!
! ----------------------------------------------------------------------
!
SUBROUTINE CalFluxIntegralWalls(q,name)
!
! ----------------------------------------------------------------------
  USE mMesh
  USE mPar


  implicit none

  real(rk) integral
  real(rk) q(nqnodes)
  integer i,j
  character(255) vrsta
  character(255) name
      
        
!
!     Calculate integral of flux
!
  DO j=1,nofw
    integral = 0.0_rk
    DO i=1,nelem
      IF (element(i)%bcid.EQ.wall(j)%id) THEN
        integral = integral + element(i)%area * q(i)
      END IF
    END DO
    write(vrsta,'(A,1X,A,G20.10)') trim(name),trim(wall(j)%name),integral
    call WriteToLog(vrsta)

  END DO

END SUBROUTINE



!
!
! ----------------------------------------------------------------------
!
SUBROUTINE CalFunctionIntegralWalls(u,name)
  !
  ! ----------------------------------------------------------------------
    USE mMesh
    USE mPar
    USE Triangle
  
    implicit none
  
    real(rk) integral,integ
    real(rk) u(nnodes)
    integer i,j,n1,n2,n3,n4
    character(255) vrsta
    character(255) name
        


    TYPE(ElementType) e    
          
  !
  !     Calculate integral of flux
  !
    DO j=1,nofw
      integral = 0.0_rk
      DO i=1,nelem
        IF (element(i)%bcid.EQ.wall(j)%id) THEN

          e = element(i)
          IF (e%type.EQ.2) THEN ! 3 node trangle
  !         list of nodes in triangle
            n1 = e%con(1)
            n2 = e%con(2)
            n3 = e%con(3)
            call Triangle_FunctionInt(node(n1)%x(1),node(n1)%x(2),node(n1)%x(3), &
                                      node(n2)%x(1),node(n2)%x(2),node(n2)%x(3), & 
                                      node(n3)%x(1),node(n3)%x(2),node(n3)%x(3), &
                                      u(n1),u(n2),u(n3), &
                                      integ)
          else
  !         list of nodes in quad
            n1 = e%con(1)
            n2 = e%con(2)
            n3 = e%con(3)            
            n4 = e%con(4)  
            call Quad_FunctionInt(node(n1)%x(1),node(n1)%x(2),node(n1)%x(3), &
                                  node(n2)%x(1),node(n2)%x(2),node(n2)%x(3), & 
                                  node(n3)%x(1),node(n3)%x(2),node(n3)%x(3), &
                                  node(n4)%x(1),node(n4)%x(2),node(n4)%x(3), &
                                  u(n1),u(n2),u(n3),u(n4), &
                                  integ)
          END IF
          integral = integral + integ
        END IF
      END DO
      write(vrsta,'(A,1X,A,G20.10)') trim(name),trim(wall(j)%name),integral
      call WriteToLog(vrsta)
  
    END DO
  
  END SUBROUTINE
          
  
!
!
! ----------------------------------------------------------------------
!
      SUBROUTINE CalMeshArea()
!
! ----------------------------------------------------------------------
      USE Triangle
      USE mMesh

      INTEGER i

!
!     Calculate element area
!
      DO i=1,nelem
        !print *, "Element ",i,": type=",element(i)%type," nno=",element(i)%nno," bcid=",element(i)%bcid, " con ",(element(i)%con(l),l=1,element(i)%nno)
        CALL CalElementArea(element(i))
      END DO
!
!     Calculate wall area
!
      DO j=1,nofw
        wall(j)%area = 0.0_rk
        DO i=1,nelem
          IF (element(i)%bcid.EQ.wall(j)%id) THEN
            wall(j)%area = wall(j)%area + element(i)%area
          END IF
        END DO
      END DO


      END SUBROUTINE



! ----------------------------------------------------------------------
!
      SUBROUTINE ReadsdBIC(fname)
!
! ----------------------------------------------------------------------
      USE mMesh
      USE mEqns
      USE mPar
      USE mCommon
      IMPLICIT NONE

      CHARACTER*(*) fname
      CHARACTER(255) OneLine,KeyWord

      INTEGER lun
      INTEGER i,j
      INTEGER Cstring

      INTEGER, ALLOCATABLE :: tmp(:)

      ALLOCATE (tmp(100))
!
!     Read BiC file
!
      lun=11

      OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)

      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')

        READ(Oneline,*) KeyWord

!
!       Subdomains
!
        IF (Cstring(KeyWord,"SUBDOMAINS").EQ.parYes) THEN
          READ(Oneline,*) KeyWord,nosd
          ALLOCATE (subdomain(nosd))
          DO i=1,nosd
            CALL rOneTL(lun,OneLine)
            READ(Oneline,*) subdomain(i)%name,subdomain(i)%diff ! first entry is the name
            subdomain(i)%nofw = 0
            DO j=1,nofw
              IF (INDEX(OneLine,TRIM(wall(j)%name)) .NE. 0) THEN
                subdomain(i)%nofw = subdomain(i)%nofw + 1
                tmp(subdomain(i)%nofw)=j
              END IF
            END DO

            ALLOCATE (subdomain(i)%loWalls(subdomain(i)%nofw))
            ALLOCATE (subdomain(i)%normMul(subdomain(i)%nofw))
            DO j=1,subdomain(i)%nofw
              subdomain(i)%loWalls(j)=tmp(j)
              subdomain(i)%normMul(j)=1.0
            END DO

          END DO

        END IF

        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)
      DEALLOCATE (tmp)
      RETURN

10    CALL WriteToLog("BiC error! - Can not open BiC file!")


      END SUBROUTINE




! ----------------------------------------------------------------------
!
      SUBROUTINE ReadBIC(fname)
!
! ----------------------------------------------------------------------
      USE mMesh
      USE mEqns
      USE mPar
      USE mCommon
      IMPLICIT NONE

      CHARACTER*(*) fname
      CHARACTER(255) OneLine,KeyWord,WallName,fq,bcform,EqName
      CHARACTER(255) sdName
      REAL(rk) multi
      INTEGER lun
      INTEGER i,j,wid,eid
      INTEGER Cstring

      INTEGER, ALLOCATABLE :: tmp(:)

      ALLOCATE (tmp(100))

!
!     Read BiC file for number of equations
!
      neq=-99
      lun=11
      OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)

      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')

        READ(Oneline,*) KeyWord
        IF (Cstring(KeyWord,"EQUATIONS").EQ.parYes) THEN
          READ(Oneline,*) KeyWord,neq
          CALL eqn_init(nofw,nnodes,nqnodes)
          DO i=1,neq
            CALL rOneTL(lun,OneLine)
            READ(Oneline,*) eqn(i)%name
          END DO
        END IF
        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)

      IF (neq.EQ.-99) THEN
        CALL WriteToLog("BiC error! - EQUATIONS missing!")
      END IF

!
!     Default - all walls zero Flux BC, initially all functions 0.0
!
      DO j=1,neq
        DO i=1,nofw
          eqn(j)%boundary(i)%known     = iFlux
          eqn(j)%boundary(i)%type      = iConst
          eqn(j)%boundary(i)%params(1) = 0.0_rk

          eqn(j)%initial(i)%known     = iFunction
          eqn(j)%initial(i)%Type      = iConst
          eqn(j)%initial(i)%params(1) = 0.0_rk
        END DO
        eqn(j)%slv%type=0
        eqn(j)%slv%pret=2
        eqn(j)%slv%prep=2
        eqn(j)%slv%maxit=500
        eqn(j)%slv%stopt=5
        eqn(j)%slv%eps=1.0E-15_rk
      END DO


!
!     Read BiC file
!
      lun=11

      OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)

      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')

        READ(Oneline,*) KeyWord

!
!       Flip normals
!
        IF (Cstring(KeyWord,"NORMMUL").EQ.parYes) THEN
          READ(Oneline,*) KeyWord,sdName,WallName,multi
          DO i=1,nosd
            IF (TRIM(sdName).EQ.TRIM(subdomain(i)%name)) THEN
              DO j=1,subdomain(i)%nofw
                IF (TRIM(WallName).EQ.TRIM(wall(subdomain(i)%loWalls(j))%name)) THEN
                   subdomain(i)%normMul(j)=multi
                END IF
              END DO
            END IF
          END DO
        END IF



!
!       Solver definitions
!
        IF (Cstring(KeyWord,"SOLVER").EQ.parYes) THEN

          READ(Oneline,*) KeyWord,EqName

!         which equation
          eid=-1
          DO i=1,neq
            IF (TRIM(EqName).EQ.TRIM(eqn(i)%name)) THEN
              eid=i
            END IF
          END DO
          IF (eid.EQ.-1) CALL WriteToLog("BiC error! - equation name")

          READ(Oneline,*) KeyWord,EqName, &
            eqn(eid)%slv%type, &
            eqn(eid)%slv%pret, &
            eqn(eid)%slv%prep, &
            eqn(eid)%slv%maxit, &
            eqn(eid)%slv%stopt, &
            eqn(eid)%slv%eps

        END IF
!
!       Boundary and initial conditions
!
        IF (Cstring(KeyWord,"BOUNDARY").EQ.parYes.OR.Cstring(KeyWord,"INITIAL").EQ.parYes) THEN

          READ(Oneline,*) KeyWord,EqName,WallName,fq,bcform

!         which equation
          eid=-1
          DO i=1,neq
            IF (TRIM(EqName).EQ.TRIM(eqn(i)%name)) THEN
              eid=i
            END IF
          END DO
          IF (eid.EQ.-1) CALL WriteToLog("BiC error! - equation name")

!         which wall
          wid=-1
          DO i=1,nofw
            IF (TRIM(WallName).EQ.TRIM(wall(i)%name)) THEN
              wid=i
            END IF
          END DO
          IF (wid.EQ.-1) CALL WriteToLog("BiC error! - wall name")

!         read condition for eqn/wall
          IF (Cstring(KeyWord,"BOUNDARY").EQ.parYes) THEN
            CALL RBiC(eqn(eid)%boundary(wid),lun,OneLine)
          ELSE
            CALL RBiC(eqn(eid)%initial(wid),lun,OneLine)
          END IF


        END IF

        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)
      DEALLOCATE (tmp)
      RETURN

10    CALL WriteToLog("BiC error! - Can not open BiC file!")


      END SUBROUTINE


! ----------------------------------------------------------------------
!
      SUBROUTINE OldReadBIC(fname)
!
! ----------------------------------------------------------------------
      USE mMesh
      USE mEqns
      USE mPar
      USE mCommon
      IMPLICIT NONE

      CHARACTER*(*) fname
      CHARACTER(255) OneLine,KeyWord,WallName,fq,bcform,EqName
      CHARACTER(255) WallName1,WallName2,WallName3,TargetWallName
      CHARACTER(255) sdName
      REAL(rk) multi
      INTEGER lun
      INTEGER i,j,k,l,jj,ll,wid,eid,wid1,wid2,twid,wid3
      INTEGER Cstring

      INTEGER, ALLOCATABLE :: tmp(:)

      ALLOCATE (tmp(100))

!
!     Read BiC file for number of equations
!
      neq=-99
      lun=11
      OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)

      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')

        READ(Oneline,*) KeyWord
        IF (Cstring(KeyWord,"EQUATIONS").EQ.parYes) THEN
          READ(Oneline,*) KeyWord,neq
          CALL eqn_init(nofw,nnodes,nqnodes)
          DO i=1,neq
            CALL rOneTL(lun,OneLine)
            READ(Oneline,*) eqn(i)%name
          END DO
        END IF
        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)

      IF (neq.EQ.-99) THEN
        CALL WriteToLog("BiC error! - EQUATIONS missing!")
      END IF

!
!     Default - all walls zero Flux BC, initially all functions 0.0
!
      DO j=1,neq
        DO i=1,nofw
          eqn(j)%boundary(i)%known     = iFlux
          eqn(j)%boundary(i)%type      = iConst
          eqn(j)%boundary(i)%params(1) = 0.0_rk

          eqn(j)%initial(i)%known     = iFunction
          eqn(j)%initial(i)%Type      = iConst
          eqn(j)%initial(i)%params(1) = 0.0_rk
        END DO
        eqn(j)%slv%type=0
        eqn(j)%slv%pret=2
        eqn(j)%slv%prep=2
        eqn(j)%slv%maxit=500
        eqn(j)%slv%stopt=5
        eqn(j)%slv%eps=1.0E-15_rk
      END DO
!
!     Read BiC file
!
      lun=11

      OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)

      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')

        READ(Oneline,*) KeyWord

!
!       Subdomains
!
        IF (Cstring(KeyWord,"SUBDOMAINS").EQ.parYes) THEN
          READ(Oneline,*) KeyWord,nosd
          ALLOCATE (subdomain(nosd))
          DO i=1,nosd
            CALL rOneTL(lun,OneLine)
            READ(Oneline,*) subdomain(i)%name,subdomain(i)%diff ! first entry is the name
            subdomain(i)%nofw = 0
            DO j=1,nofw
              IF (INDEX(OneLine,TRIM(wall(j)%name)) .NE. 0) THEN
                subdomain(i)%nofw = subdomain(i)%nofw + 1
                tmp(subdomain(i)%nofw)=j
              END IF
            END DO

            ALLOCATE (subdomain(i)%loWalls(subdomain(i)%nofw))
            ALLOCATE (subdomain(i)%normMul(subdomain(i)%nofw))
            DO j=1,subdomain(i)%nofw
              subdomain(i)%loWalls(j)=tmp(j)
              subdomain(i)%normMul(j)=1.0
            END DO

          END DO

        END IF

!
!       Flip normals
!
        IF (Cstring(KeyWord,"NORMMUL").EQ.parYes) THEN
          READ(Oneline,*) KeyWord,sdName,WallName,multi
          DO i=1,nosd
            IF (TRIM(sdName).EQ.TRIM(subdomain(i)%name)) THEN
              DO j=1,subdomain(i)%nofw
                IF (TRIM(WallName).EQ.TRIM(wall(subdomain(i)%loWalls(j))%name)) THEN
                   subdomain(i)%normMul(j)=multi
                END IF
              END DO
            END IF
          END DO
        END IF



!
!       Solver definitions
!
        IF (Cstring(KeyWord,"SOLVER").EQ.parYes) THEN

          READ(Oneline,*) KeyWord,EqName

!         which equation
          eid=-1
          DO i=1,neq
            IF (TRIM(EqName).EQ.TRIM(eqn(i)%name)) THEN
              eid=i
            END IF
          END DO
          IF (eid.EQ.-1) CALL WriteToLog("BiC error! - equation name")

          READ(Oneline,*) KeyWord,EqName, &
            eqn(eid)%slv%type, &
            eqn(eid)%slv%pret, &
            eqn(eid)%slv%prep, &
            eqn(eid)%slv%maxit, &
            eqn(eid)%slv%stopt, &
            eqn(eid)%slv%eps

        END IF
!
!       Boundary and initial conditions
!
        IF (Cstring(KeyWord,"BOUNDARY").EQ.parYes.OR.Cstring(KeyWord,"INITIAL").EQ.parYes) THEN

          READ(Oneline,*) KeyWord,EqName,WallName,fq,bcform

!         which equation
          eid=-1
          DO i=1,neq
            IF (TRIM(EqName).EQ.TRIM(eqn(i)%name)) THEN
              eid=i
            END IF
          END DO
          IF (eid.EQ.-1) CALL WriteToLog("BiC error! - equation name")

!         which wall
          wid=-1
          DO i=1,nofw
            IF (TRIM(WallName).EQ.TRIM(wall(i)%name)) THEN
              wid=i
            END IF
          END DO
          IF (wid.EQ.-1) CALL WriteToLog("BiC error! - wall name")

!         read condition for eqn/wall
          IF (Cstring(KeyWord,"BOUNDARY").EQ.parYes) THEN
            CALL RBiC(eqn(eid)%boundary(wid),lun,OneLine)
          ELSE
            CALL RBiC(eqn(eid)%initial(wid),lun,OneLine)
          END IF


        END IF
!
!       Node placement
!
        IF (Cstring(KeyWord,"EDGE").EQ.parYes) THEN


          READ(Oneline,*) KeyWord,WallName1,WallName2,TargetWallName


!         which wall 1
          wid1=-1
          DO i=1,nofw
            IF (TRIM(WallName1).EQ.TRIM(wall(i)%name)) THEN
              wid1=i
            END IF
          END DO
          IF (wid1.EQ.-1) CALL WriteToLog("BiC error! - EDGE wall name")

!         which wall 2
          wid2=-1
          DO i=1,nofw
            IF (TRIM(WallName2).EQ.TRIM(wall(i)%name)) THEN
              wid2=i
            END IF
          END DO
          IF (wid2.EQ.-1) CALL WriteToLog("BiC error! - EDGE wall name")

!         which wall TARGET
          twid=-1
          DO i=1,nofw
            IF (TRIM(TargetWallName).EQ.TRIM(wall(i)%name)) THEN
              twid=i
            END IF
          END DO
          IF (twid.EQ.-1) CALL WriteToLog("BiC error! - EDGE wall name")

!         apply edge rule to node%bcid

          DO i=1,nelem
            IF (element(i)%bcid.EQ.wid1) THEN
              DO j=1,nelem
                IF (element(j)%bcid.EQ.wid2) THEN
!                 Do elemnts share nodes?
                  DO k=1,element(i)%nno
                    DO l=1,element(j)%nno
                      IF (element(i)%con(k).EQ.element(j)%con(l)) THEN
                        node(element(i)%con(k))%bcid = twid
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END IF
          END DO

        END IF

        IF (Cstring(KeyWord,"CORNER").EQ.parYes) THEN


          READ(Oneline,*) KeyWord,WallName1,WallName2,WallName3,TargetWallName


!         which wall 1
          wid1=-1
          DO i=1,nofw
            IF (TRIM(WallName1).EQ.TRIM(wall(i)%name)) THEN
              wid1=i
            END IF
          END DO
          IF (wid1.EQ.-1) CALL WriteToLog("BiC error! - CORNER wall name")

!         which wall 2
          wid2=-1
          DO i=1,nofw
            IF (TRIM(WallName2).EQ.TRIM(wall(i)%name)) THEN
              wid2=i
            END IF
          END DO
          IF (wid2.EQ.-1) CALL WriteToLog("BiC error! - CORNER wall name")

!         which wall 3
          wid3=-1
          DO i=1,nofw
            IF (TRIM(WallName3).EQ.TRIM(wall(i)%name)) THEN
              wid3=i
            END IF
          END DO
          IF (wid3.EQ.-1) CALL WriteToLog("BiC error! - CORNER wall name")

!         which wall TARGET
          twid=-1
          DO i=1,nofw
            IF (TRIM(TargetWallName).EQ.TRIM(wall(i)%name)) THEN
              twid=i
            END IF
          END DO
          IF (twid.EQ.-1) CALL WriteToLog("BiC error! - CORNER wall name")

!         apply edge rule to node%bcid

          DO i=1,nelem
            IF (element(i)%bcid.EQ.wid1) THEN
              DO j=1,nelem
                IF (element(j)%bcid.EQ.wid2) THEN
                  DO jj=1,nelem
                    IF (element(jj)%bcid.EQ.wid3) THEN
!                     Do elemnts share nodes?
                      DO k=1,element(i)%nno
                        DO l=1,element(j)%nno
                          DO ll=1,element(jj)%nno
                            IF (element(i)%con(k).EQ.element(j)%con(l).AND.element(i)%con(k).EQ.element(jj)%con(ll)) THEN
                              node(element(i)%con(k))%bcid = twid
                            END IF
                          END DO
                        END DO
                      END DO
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END DO

        END IF








        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)
      DEALLOCATE (tmp)
      RETURN

10    CALL WriteToLog("BiC error! - Can not open BiC file!")


      END SUBROUTINE

!
!     ------------------------------------------------------------------
!
      SUBROUTINE RBiC(bic,lun,OneLine)

      USE mMesh
      USE mEqns
      IMPLICIT NONE

      INTEGER i,lun

      TYPE(ConditionType) bic

      CHARACTER(255) OneLine,KeyWord,WallName,fq,bcform,EqName


      READ(Oneline,*) KeyWord,EqName,WallName,fq,bcform

          IF (TRIM(fq).EQ."function") THEN
            bic%known = iFunction
          ELSE IF (TRIM(fq).EQ."flux") THEN
            bic%known = iFlux
          ELSE IF (TRIM(fq).EQ."contact") THEN
            bic%known = iContact
          ELSE
            Print *,"BC3 error!"
          END IF

          IF (TRIM(bcform).EQ."const") THEN
            bic%Type = iConst
          ELSE IF (TRIM(bcform).EQ."lin") THEN
            bic%Type = iLin
          ELSE IF (TRIM(bcform).EQ."quad") THEN
            bic%Type = iQuad
          ELSE
            Print *,"BC4 error!"
          END IF


          CALL rOneTL(lun,OneLine)

          IF ( bic%Type .EQ. iConst ) THEN
            READ(Oneline,*) bic%params(1)
          ELSE IF ( bic%Type .EQ. iLin ) THEN
            READ(Oneline,*) (bic%params(i),i=1,4)
          ELSE IF ( bic%Type .EQ. iQuad ) THEN
            READ(Oneline,*) (bic%params(i),i=1,10)
          END IF


      END SUBROUTINE
!
! ----------------------------------------------------------------------
!
      SUBROUTINE sdApplyInitalCond()

      USE mMesh
      USE mEqns

      IMPLICIT NONE

      INTEGER e,w,i,isd,inode
!
!     Set to zero
!
      DO e=1,neq
        DO i=1,nnodes
          eqn(e)%u(i)=0.0_rk
        END DO
        DO i=1,nqnodes
          eqn(e)%q(i)=0.0_rk
        END DO
      END DO


!
!     Initial conditions
!
      DO e=1,neq

        DO isd = 1,nosd
          DO i=1,subdomain(isd)%nnodes
            inode = subdomain(isd)%nodeList(i)
            w  = subdomain(isd)%BCidList(i)   ! wall to which this node belongs

!           znana funkcija
            IF ( eqn(e)%initial(w)%known .EQ. iFunction ) THEN
              CALL GetBiCValue(eqn(e)%initial(w),node(inode),eqn(e)%u(inode))
            END IF

            IF ( eqn(e)%boundary(w)%known .EQ. iFunction ) THEN
              CALL GetBiCValue(eqn(e)%boundary(w),node(inode),eqn(e)%u(inode))
            END IF

            IF ( eqn(e)%boundary(w)%known .EQ. iContact ) THEN
              CALL GetBiCValue(eqn(e)%boundary(w),node(inode),eqn(e)%u(inode))
            END IF

          END DO

          DO i=1,subdomain(isd)%nqnodes
            inode = subdomain(isd)%qnodeList(i)
            w  = subdomain(isd)%qBCidList(i)   ! wall to which this node belongs

!           znan fluks
            IF ( eqn(e)%initial(w)%known .EQ. iFlux ) THEN
              CALL GetBiCValue(eqn(e)%initial(w),qnode(inode),eqn(e)%q(inode))
            END IF

            IF ( eqn(e)%boundary(w)%known .EQ. iFlux ) THEN
              CALL GetBiCValue(eqn(e)%boundary(w),qnode(inode),eqn(e)%q(inode))
            END IF

            IF ( eqn(e)%boundary(w)%known .EQ. iContact ) THEN
              CALL GetBiCValue(eqn(e)%boundary(w),qnode(inode),eqn(e)%q(inode))
            END IF
          END DO

        END DO ! isd
      END DO ! e

      END SUBROUTINE



!
! ----------------------------------------------------------------------
!
      SUBROUTINE ApplyInitalCond()

      USE mMesh
      USE mEqns

      IMPLICIT NONE

      INTEGER e,w,i


!
!     Initial conditions
!
      DO e=1,neq
        DO w=1,nofw

!         znana funkcija
          IF ( eqn(e)%initial(w)%known .EQ. iFunction ) THEN
            DO i=1,nnodes
              IF (node(i)%bcid.EQ.w) THEN
                CALL GetBiCValue(eqn(e)%initial(w),node(i),eqn(e)%u(i))
              END IF
            END DO
          END IF
!         znan fluks
          IF ( eqn(e)%initial(w)%known .EQ. iFlux ) THEN
            DO i=1,nqnodes
              IF (qnode(i)%bcid.EQ.w) THEN
                CALL GetBiCValue(eqn(e)%initial(w),qnode(i),eqn(e)%q(i))
              END IF
            END DO
          END IF

        END DO
      END DO

!
!     Boundary conditions
!
      DO e=1,neq
        DO w=1,nofw

!         znana funkcija
          IF ( eqn(e)%boundary(w)%known .EQ. iFunction ) THEN
            DO i=1,nnodes
              IF (node(i)%bcid.EQ.w) THEN
                CALL GetBiCValue(eqn(e)%boundary(w),node(i),eqn(e)%u(i))
              END IF
            END DO
          END IF
!         znan fluks
          IF ( eqn(e)%boundary(w)%known .EQ. iFlux ) THEN

            DO i=1,nqnodes
              IF (qnode(i)%bcid.EQ.w) THEN
                CALL GetBiCValue(eqn(e)%boundary(w),qnode(i),eqn(e)%q(i))
              END IF
            END DO
          END IF
        END DO
      END DO



      END SUBROUTINE
!
! ----------------------------------------------------------------------
!
      SUBROUTINE GetBiCValue(bic,nod,val)

      USE mMesh
      USE mEqns
      USE mCommon
      TYPE(ConditionType) bic
      TYPE (NodeType) nod

      REAL(rk) val

      IF ( bic%type .EQ. iConst ) THEN
        val = bic%params(1)
      ELSE IF ( bic%type .EQ. iLin ) THEN
                  val =        bic%params(1)+ &
                               bic%params(2)*nod%x(1)+ &
                               bic%params(3)*nod%x(2)+ &
                               bic%params(4)*nod%x(3)
      ELSE IF ( bic%type .EQ. iQuad ) THEN
                  val =         bic%params(1)+ &
                               bic%params(2)*nod%x(1)+ &          ! x
                               bic%params(3)*nod%x(2)+ &          ! y
                               bic%params(4)*nod%x(3)+ &          ! z
                               bic%params(5)*nod%x(1)*nod%x(2)+ & ! xy
                               bic%params(6)*nod%x(1)*nod%x(3)+ & ! xz
                               bic%params(7)*nod%x(2)*nod%x(3)+ & ! yz
                               bic%params(8)*nod%x(1)*nod%x(1)+ & ! x^2
                               bic%params(9)*nod%x(2)*nod%x(2)+ & ! y^2
                              bic%params(10)*nod%x(3)*nod%x(3)    ! z^2
      END IF

      END SUBROUTINE


!______________________________________________________________________C
      SUBROUTINE rOneTL(lun,OneLine)
!     _    ___ _    _
!     Read One Text Line
!
!     Returns the first nonempty text line in file LUN, which does not
!     include the # character. If end of file is encoutered, it returns EOF
      CHARACTER*(*) OneLine
      INTEGER lun,i

10    READ(lun,'(A)',END=20) OneLine

!     Check if line is empty
      IF (len_trim(OneLine).EQ.0) GOTO 10

!     Check if line contains # character
      DO i=1,len_trim(OneLine)
        IF (OneLine(i:i).EQ.'#') GOTO 10
      ENDDO

      RETURN

20    OneLine='EOF'
      END

!
! ----------------------------------------------------------------------
!
      SUBROUTINE PrintMeshStats()
      USE mPar
      USE mMesh
      IMPLICIT NONE

      CHARACTER(255) tekst

      !INTEGER lun ! use lun=6 for screen
      INTEGER i

      WRITE (tekst,'(A)') "Mesh statistics"
      CALL WriteToLog(tekst)
      WRITE (tekst,'(A,A)') "Name   = ",TRIM(meshName)
      CALL WriteToLog(tekst)
      WRITE (tekst,'(A,I0)') "Number of nodes    = ",nnodes
      CALL WriteToLog(tekst)
      WRITE (tekst,'(A,I0)') "Number of elements = ",nelem
      CALL WriteToLog(tekst)
      WRITE (tekst,'(A,I0)') "Number of walls    = ",nofw
      CALL WriteToLog(tekst)

      DO i=1,nofw
        WRITE (tekst,*) wall(i)%id,wall(i)%area,TRIM(wall(i)%name)
        CALL WriteToLog(tekst)
      END DO

      WRITE (tekst,'(A,3(1X,F10.4))') "min:", mextMin
      CALL WriteToLog(tekst)
      WRITE (tekst,'(A,3(1X,F10.4))') "max:", mextMax
      CALL WriteToLog(tekst)

      END SUBROUTINE

!
!     ------------------------------------------------------------------
!
      SUBROUTINE WriteRestartFile()

      USE mMesh
      USE mPar
      USE mEqns

      IMPLICIT NONE

      INTEGER i,lun

      lun=12

      OPEN (lun,FILE=TRIM(parRstFileName),FORM='UNFORMATTED',STATUS='UNKNOWN')
!
!     Sizes
!
      DO i=1,nosd
        WRITE (lun) subdomain(i)%nsp,subdomain(i)%nnodes,subdomain(i)%nqnodes
      END DO
!
!     Solution fields
!

      DO i=1,neq
        CALL wrvec(lun,nnodes,eqn(i)%u)
        CALL wrvec(lun,nqnodes,eqn(i)%q)
      END DO

      CLOSE (lun)

      END


!
!     ------------------------------------------------------------------
!
      SUBROUTINE ReadRestartFile(ierr)


      USE mMesh
      USE mPar
      USE mEqns

      IMPLICIT NONE

      INTEGER i,lun,ierr,a,b,c

      ierr=1
      lun=12

      OPEN (lun,FILE=TRIM(parRstFileName),FORM='UNFORMATTED',STATUS='OLD',ERR=10)
!
!     Loop over subdomains
!
      DO i=1,nosd
        READ (lun,ERR=10) a,b,c
        IF ( (a.NE.subdomain(i)%nsp) .OR. (subdomain(i)%nnodes.NE.b) .OR. (subdomain(i)%nqnodes.NE.c) ) THEN
          GOTO 10
        END IF
      END DO
!
!     Solution fields
!
      DO i=1,neq
        CALL rdvec(lun,nnodes,eqn(i)%u)
        CALL rdvec(lun,nqnodes,eqn(i)%q)
      END DO


      ierr=0

10    CONTINUE
      CLOSE(lun)

      END SUBROUTINE



!
!     ------------------------------------------------------------------
!
      SUBROUTINE PostOnly()


      USE mMesh
      USE mPar
      USE mEqns

      IMPLICIT NONE

      INTEGER ierr

      CALL WriteToLog("Postprocessing run only!")


      CALL WriteToLog("Attempting to read restart file!")
      CALL ReadRestartFile(ierr)

      IF (ierr.NE.0) THEN
        CALL WriteToLog("Read failed - Exiting!")
        CALL StopProgram()
      ELSE
        CALL WriteToLog("Read sucessful!")
      END IF

      END SUBROUTINE






! -----------------------------------------------------------------------------

  SUBROUTINE OutputOpenFOAM()

!
!     $: Outputs data on domain mesh in OpenFOAM format
!
! -----------------------------------------------------------------------------
      USE mDomainMesh
      USE mMesh
      USE mDomainData
      USE mPar
      USE mEqns
      IMPLICIT NONE

      INTEGER itp,i,n1,n2,n3,iWall,ne
      CHARACTER*255 vrstica
      real(rk) Point(3)

      itp=96

      OPEN (itp,FILE=TRIM(parDomainOpenFOAMResultsFileName),STATUS='UNKNOWN')

!
!    HEADER
!      
      WRITE(itp,'(A)') '/*--------------------------------*- C++ -*----------------------------------*\'
      WRITE(itp,'(A)') '  =========                 |'
      WRITE(itp,'(A)') '  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox'
      WRITE(itp,'(A)') '   \\    /   O peration     | Website:  https://openfoam.org'
      WRITE(itp,'(A)') '    \\  /    A nd           | Version:  11'
      WRITE(itp,'(A)') '     \\/     M anipulation  |'
      WRITE(itp,'(A)') '\*---------------------------------------------------------------------------*/'
      WRITE(itp,'(A)') 'FoamFile'
      WRITE(itp,'(A)') '{'
      WRITE(itp,'(A)') '    format      ascii;'
      WRITE(itp,'(A)') '    class       volVectorField;'
      WRITE(itp,'(A)') '    location    "constant";'
      WRITE(itp,'(A)') '    object      U;'
      WRITE(itp,'(A)') '}'
      WRITE(itp,'(A)') '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
      WRITE(itp,'(A)') ''
      WRITE(itp,'(A)') 'dimensions      [0 1 -1 0 0 0 0];'
      WRITE(itp,'(A)') ''
      WRITE(itp,'(A)') 'internalField   nonuniform List<vector> '

      WRITE(itp,'(I0)') DMnelem
      WRITE(itp,'(A)') '('
!
!     ELEMENT BASED DATA
!              
!     Domain Velocity field      
      DO  i=1,DMnelem
          WRITE (vrstica,'(A,3G18.9,A)') "(",ddE(1)%val(i),ddE(2)%val(i),ddE(3)%val(i),")"
          CALL sqblnk(itp,vrstica)
      END DO
      WRITE(itp,'(A)') ')'
      WRITE(itp,'(A)') ';'
      WRITE(itp,'(A)') ''
!     Boundary velocity field
      WRITE(itp,'(A)') 'boundaryField'
      WRITE(itp,'(A)') '{'
      DO iWall = 1,nofw
        WRITE(itp,'(A,A)') '    ',TRIM(wall(iWall)%name)
        WRITE(itp,'(A)') '    {'
        WRITE(itp,'(A)') '        type            calculated;'
        WRITE(itp,'(A)') '        value           nonuniform List<vector> '
        ne = 0
        DO i = 1,nelem
          IF (element(i)%bcid.eq.iWall) THEN
            ne = ne + 1
          END IF
        END DO
        WRITE(itp,'(I0)') ne
        WRITE(itp,'(A)') '('
        
        do i = 1,nelem
          if (element(i)%bcid.eq.iWall) then
            ! get the coordinates of the element center
            n1 = element(i)%con(1)
            n2 = element(i)%con(2)
            n3 = element(i)%con(3)

            Point(1) = (eqn(1)%u(n1) + eqn(1)%u(n2) + eqn(1)%u(n3) ) / 3.0_rk
            Point(2) = (eqn(2)%u(n1) + eqn(2)%u(n2) + eqn(2)%u(n3) ) / 3.0_rk
            Point(3) = (eqn(3)%u(n1) + eqn(3)%u(n2) + eqn(3)%u(n3) ) / 3.0_rk
            WRITE (vrstica,'(A,3G18.9,A)') "(",Point(1),Point(2),Point(3),")"
            CALL sqblnk(itp,vrstica)
          end if
        end do

        WRITE(itp,'(A)') ')'
        WRITE(itp,'(A)') ';'
        WRITE(itp,'(A)') '    }'
      END DO
      WRITE(itp,'(A)') '}'
      WRITE(itp,'(A)') ''
      WRITE(itp,'(A)') ''
      WRITE(itp,'(A)') '// ************************************************************************* //'

    CLOSE (itp)
END
        
        


! -----------------------------------------------------------------------------

SUBROUTINE OutputOpenFOAMc()

  !
  !     $: Outputs element centers on domain mesh in OpenFOAM format
  !
  ! -----------------------------------------------------------------------------
        USE mDomainMesh
        USE mMesh
        USE mDomainData
        USE mPar
        IMPLICIT NONE
  
        INTEGER itp,i,n1,n2,n3,n4,iWall,ne
        CHARACTER*255 vrstica
        real(rk) Point(3)
  
        itp=96
  
        OPEN (itp,FILE="and.C",STATUS='UNKNOWN')
  
  !
  !    HEADER
  !      
        WRITE(itp,'(A)') '/*--------------------------------*- C++ -*----------------------------------*\'
        WRITE(itp,'(A)') '  =========                 |'
        WRITE(itp,'(A)') '  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox'
        WRITE(itp,'(A)') '   \\    /   O peration     | Website:  https://openfoam.org'
        WRITE(itp,'(A)') '    \\  /    A nd           | Version:  11'
        WRITE(itp,'(A)') '     \\/     M anipulation  |'
        WRITE(itp,'(A)') '\*---------------------------------------------------------------------------*/'
        WRITE(itp,'(A)') 'FoamFile'
        WRITE(itp,'(A)') '{'
        WRITE(itp,'(A)') '    format      ascii;'
        WRITE(itp,'(A)') '    class       volVectorField;'
        WRITE(itp,'(A)') '    location    "constant";'
        WRITE(itp,'(A)') '    object      C;'
        WRITE(itp,'(A)') '}'
        WRITE(itp,'(A)') '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'
        WRITE(itp,'(A)') ''
        WRITE(itp,'(A)') 'dimensions      [0 1 0 0 0 0 0];'
        WRITE(itp,'(A)') ''
        WRITE(itp,'(A)') 'internalField   nonuniform List<vector> '
  
        WRITE(itp,'(I0)') DMnelem
        WRITE(itp,'(A)') '('
  !
  !     ELEMENT BASED DATA
  !              
  !     Domain field      
        DO  i=1,DMnelem
  
            n1 = DMelement(i)%con(1)
            n2 = DMelement(i)%con(2)
            n3 = DMelement(i)%con(3)
            n4 = DMelement(i)%con(4)
          
            Point(1) = (DMnode(n1)%x(1) + DMnode(n2)%x(1) + DMnode(n3)%x(1) + DMnode(n4)%x(1)) / 4.0_rk
            Point(2) = (DMnode(n1)%x(2) + DMnode(n2)%x(2) + DMnode(n3)%x(2) + DMnode(n4)%x(2)) / 4.0_rk
            Point(3) = (DMnode(n1)%x(3) + DMnode(n2)%x(3) + DMnode(n3)%x(3) + DMnode(n4)%x(3)) / 4.0_rk        
            WRITE (vrstica,'(A,3G18.9,A)') "(",Point(1),Point(2),Point(3),")"
            CALL sqblnk(itp,vrstica)
        END DO
        WRITE(itp,'(A)') ')'
        WRITE(itp,'(A)') ';'
        WRITE(itp,'(A)') ''
  !     Boundary field
        WRITE(itp,'(A)') 'boundaryField'
        WRITE(itp,'(A)') '{'
        DO iWall = 1,nofw
          WRITE(itp,'(A,A)') '    ',TRIM(wall(iWall)%name)
          WRITE(itp,'(A)') '    {'
          WRITE(itp,'(A)') '        type            calculated;'
          WRITE(itp,'(A)') '        value           nonuniform List<vector> '
          ne = 0
          DO i = 1,nelem
            IF (element(i)%bcid.eq.iWall) THEN
              ne = ne + 1
            END IF
          END DO
          WRITE(itp,'(I0)') ne
          WRITE(itp,'(A)') '('
          
          do i = 1,nelem
            if (element(i)%bcid.eq.iWall) then
              ! get the coordinates of the element center
              n1 = element(i)%con(1)
              n2 = element(i)%con(2)
              n3 = element(i)%con(3)
              Point(1) = (node(n1)%x(1) + node(n2)%x(1) + node(n3)%x(1) ) / 3.0_rk
              Point(2) = (node(n1)%x(2) + node(n2)%x(2) + node(n3)%x(2) ) / 3.0_rk
              Point(3) = (node(n1)%x(3) + node(n2)%x(3) + node(n3)%x(3) ) / 3.0_rk
              WRITE (vrstica,'(A,3G18.9,A)') "(",Point(1),Point(2),Point(3),")"
              CALL sqblnk(itp,vrstica)
            end if
          end do
  
          WRITE(itp,'(A)') ')'
          WRITE(itp,'(A)') ';'
          WRITE(itp,'(A)') '    }'
        END DO
        WRITE(itp,'(A)') '}'
        WRITE(itp,'(A)') ''
        WRITE(itp,'(A)') ''
        WRITE(itp,'(A)') '// ************************************************************************* //'
  
      CLOSE (itp)
  END
  


! -----------------------------------------------------------------------------

SUBROUTINE OutputDomainList()

  !
  !     $: Outputs data on domain mesh
  !
  ! -----------------------------------------------------------------------------
        USE mDomainMesh
        USE mDomainData
        USE mPar
        IMPLICIT NONE
  
        INTEGER itp,i
        CHARACTER*255 vrstica
  
        itp=96
  
        OPEN (itp,FILE=TRIM(parDomainListResultsFileName),STATUS='UNKNOWN')
  
        WRITE (itp,'(A)') "#"
        WRITE (itp,'(A)') "# Domain values"
        WRITE (itp,'(A)') "#"
        WRITE (itp,'(A)') "# Number of values"
        WRITE (itp,'(I0)') DMnnodes
        WRITE (itp,'(A)') "#"
        WRITE (itp,'(8A)') "# x  y  z  ", trim(ddN(1)%name)," ",trim(ddN(2)%name)," ",trim(ddN(3)%name)," ",trim(ddN(4)%name)
        DO  i=1,DMnnodes
          WRITE (vrstica,'(7G18.9)') DMnode(i)%x(1),DMnode(i)%x(2),DMnode(i)%x(3),ddN(1)%val(i),ddN(2)%val(i),ddN(3)%val(i),ddN(4)%val(i)
          CALL sqblnk(itp,vrstica)
        END DO

        close(itp)   
  END



! -----------------------------------------------------------------------------

SUBROUTINE OutputDomainMeshParaview()

  !
  !     $: Outputs data on domain mesh
  !
  ! -----------------------------------------------------------------------------
        USE mDomainMesh
        USE mDomainData
        USE mPar
        IMPLICIT NONE
  
        INTEGER itp,i,of,k,j
        CHARACTER*255 vrstica,tmp
  
        itp=96
  
        OPEN (itp,FILE=TRIM(parDomainResultsFileName),STATUS='UNKNOWN')
  
  
        WRITE (itp,'(A)') '<?xml version="1.0"?>'
        WRITE (itp,'(A)')&
       '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
        WRITE (itp,'(A)') '<UnstructuredGrid>'
        WRITE (itp,'(A,I10,A,I10,A)') '<Piece NumberOfPoints="',DMnnodes,'" NumberOfCells="',DMnelem,'">'
  !
  !     PODATKI V VOZLISCIH
  !
        WRITE (itp,"(A)") '<PointData Scalars="node">'
        WRITE (itp,'(A,A,A)') '<DataArray type="Float32" Name="',trim(ddN(1)%name),'" format="ascii">'
        DO  i=1,DMnnodes
          WRITE (itp,*) ddN(1)%val(i)
        END DO
        WRITE (itp,'(A)') '</DataArray>'      
        WRITE (itp,'(A)') '</PointData>'
  
  !
  !     ELEMENT BASED DATA
  !              
        WRITE (itp,'(A)') '<CellData Scalars="pressure" Vectors="velocity">'
  !     Pressure field      
        WRITE (itp,'(A)') '<DataArray type="Float32"  Name="pressure" format="ascii">'
        DO  i=1,DMnelem
            WRITE (vrstica,'(3G18.9)') ddE(4)%val(i)
            CALL sqblnk(itp,vrstica)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
  !     Velocity field      
        WRITE (itp,'(A)') '<DataArray type="Float32"  Name="velocity" NumberOfComponents="3" format="ascii">'
        DO  i=1,DMnelem
            WRITE (vrstica,'(3G18.9)') ddE(1)%val(i),ddE(2)%val(i),ddE(3)%val(i)
            CALL sqblnk(itp,vrstica)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
  !      
        WRITE (itp,'(A)') '</CellData>'
  
  !
  !      KOORDINATE VOZLISC
  !      
        WRITE (itp,'(A)') '<Points>'
        WRITE (itp,'(A)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
  
        DO  i=1,DMnnodes
          WRITE (vrstica,'(3G18.9)') DMnode(i)%x(1),DMnode(i)%x(2),DMnode(i)%x(3)
          CALL sqblnk(itp,vrstica)
        END DO
  
        WRITE (itp,'(A)') '</DataArray>'
        WRITE (itp,'(A)') '</Points>'
  !
  !     ELEMENTS
  !                  
        WRITE (itp,'(A)') '<Cells>'
  !
  !     ELEMENT CONNECTIVITY
  !
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
  
        DO i=1,DMnelem
            WRITE(vrstica,*) (DMelement(i)%con(k)-1,k=1,DMelement(i)%nno)
            CALL sqblnk(itp,vrstica)
        END DO
        WRITE (itp,'(A)') '</DataArray>'
  !
  !     OFFSETS
  !                  
        WRITE (itp,'(A)') '<DataArray type="Int32" Name="offsets" format="ascii">'
        of = 0
        DO i=1,DMnelem
            of = of + DMelement(i)%nno
            WRITE (itp,*) of
        END DO
        WRITE (itp,'(A)') '</DataArray>'
  !
  !     ELEMENT TYPES
  !                  
        WRITE (itp,'(A)') '<DataArray type="UInt8" Name="types" format="ascii">'
        DO i=1,DMnelem
            !IF (element(i)%type.EQ.2) THEN ! three node triangle
            !  WRITE (itp,*) "5" ! 5 tikotnik, 9 quuad, 12 heksaeder, 10 tetraeder
            !ELSE IF (element(i)%type.EQ.3) THEN ! four node Quadrangle
            ! WRITE (itp,*) "9" ! 5 tikotnik, 9 quuad, 12 heksaeder
            !ELSE 
          IF (DMelement(i)%type.EQ.4) THEN ! four node tetraeder
              WRITE (itp,*) "10" ! 5 tikotnik, 9 quuad, 12 heksaeder, 10 tetraeder
          ELSE
            CALL WriteToLog("Error :: Element type not supported!")
          END IF
        END DO
        WRITE (itp,'(A)') '</DataArray>'
  
      WRITE (itp,'(A)') '</Cells>'
  
  
      WRITE (itp,'(A)') '</Piece>'
      WRITE (itp,'(A)') '</UnstructuredGrid>'
      WRITE (itp,'(A)') '</VTKFile>'
      CLOSE (itp)
  END
          
          
subroutine exportWallNodeList()
  use mMesh
  use mPar
  implicit none

  integer i, lun,isd,iW,ww,j
  character(255) vrstica

  lun = 12

  open(lun, file="and.wallNodeList", status='unknown')

  DO isd=1,nosd ! loop over subdomains
    WRITE(vrstica,*) "Subdomain name: ", subdomain(isd)%name
    CALL sqblnk(lun,vrstica)
    DO iW = 1,subdomain(isd)%nofw      
      ww = subdomain(isd)%loWalls(iW)
      WRITE(vrstica,*) "Wall name: ", wall(ww)%name
      CALL sqblnk(lun,vrstica)
      DO j=1,subdomain(isd)%nnodes        
        if (subdomain(isd)%BCidList(j).eq.ww) then
          i = subdomain(isd)%nodeList(j)
          WRITE(vrstica,*) node(i)%x(1), node(i)%x(2), node(i)%x(3)
          CALL sqblnk(lun,vrstica)
        end if
      END DO
    END DO                                      
  END DO

  close(lun)

end subroutine exportWallNodeList
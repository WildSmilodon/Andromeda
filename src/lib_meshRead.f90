
! ----------------------------------------------------------------------
!
SUBROUTINE ReadDomainMesh(fname)
!
! ----------------------------------------------------------------------
      USE mPar
      IMPLICIT NONE
            
      CHARACTER*(*) fname
      CHARACTER(255) OneLine
      INTEGER lun
      REAL(rk) version
      character(20) meshType

      lun=11
      version = 0.0_rk
      OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)

      CALL rOneTL(lun,OneLine)
      IF (trim(OneLine).EQ."$MeshFormat") THEN
        CALL rOneTL(lun,OneLine)
        READ(Oneline,*) version
      end if
      CLOSE(lun)

      IF (version .GE. 4.0_rk) then
        CALL getMeshTypev4(fname,meshType)
        if (meshType.EQ."volume") then
          CALL ReadGMSHv4Volume(fname)
        else
          CALL WriteToLog("Error :: could not read domain mesh file!")
          CALL StopProgram
        end if
      else
        goto 10
      end if
     
    
      RETURN
  
  
      10  CALL WriteToLog("Error :: could not open domain mesh file!")
        CALL StopProgram


END SUBROUTINE

! ----------------------------------------------------------------------
!
SUBROUTINE ReadMesh(fname)
!
! ----------------------------------------------------------------------
    USE mPar
    IMPLICIT NONE
          
    CHARACTER*(*) fname
    CHARACTER(255) OneLine
    INTEGER lun
    REAL(rk) version
    character(20) meshType
          
        
    lun=11
      
    OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)
    CALL rOneTL(lun,OneLine)
    IF (trim(OneLine).EQ."$MeshFormat") THEN
      CALL rOneTL(lun,OneLine)
      READ(Oneline,*) version
      CLOSE(lun)
      IF (version < 4.0_rk) then
        CALL ReadGMSHv2(fname)
      else
        CALL getMeshTypev4(fname,meshType)
        if (meshType.EQ."surface") then
          CALL ReadGMSHv4Surface(fname)
        else
          CALL WriteToLog("Error :: invalid surface mesh file!")
          CALL StopProgram
        end if
      end if
    ELSE
      CALL ReadVTKmesh(fname)
    END IF
   
  !
  !     Transform mesh (rotation, stretching, translation)
  !
    IF (pariMTstr.EQ.parYes) CALL MeshStretch()
    IF (pariMTrot.EQ.parYes) CALL MeshRotate()
    IF (pariMTtra.EQ.parYes) CALL MeshTranslate()
    IF (pariMTcyl.EQ.parYes) CALL MeshCylStretch()
    !IF (pariMTsel.EQ.parYes) call sphere2superE() ! JURE
    IF (pariMTsel.EQ.parYes) call projectPoints2superE() ! MITJA
  !
  !     Get mesh extents, output to log file.
  !
    CALL MeshGetExtents()
  
    RETURN
  
  
  10  CALL WriteToLog("Error :: could not open mesh file!")
    CALL StopProgram
  
  END SUBROUTINE
          

! ----------------------------------------------------------------------
!
  SUBROUTINE ReadVTKmesh(fname)
!
! ----------------------------------------------------------------------
  USE mMesh
  USE mPar
  USE mCommon
  IMPLICIT NONE

  CHARACTER*(*) fname

  CALL WriteToLog("Reading VTK mesh: "//trim(fname))
  CALL WriteToLog("Error :: VTK mesh reader not implemented yet!")
  CALL StopProgram()
END SUBROUTINE
                  

    

! ----------------------------------------------------------------------
!
  SUBROUTINE ReadGMSHv2(fname)  ! gmsh -2 -format msh2 ime.geo
    !
    ! ----------------------------------------------------------------------
          USE mMesh
          USE mPar
          IMPLICIT NONE
    
          CHARACTER*(*) fname
          CHARACTER(255) OneLine
          INTEGER lun,dummy
          INTEGER i,j,n,k,l,nd
    
          CALL WriteToLog("Reading v2 mesh: "//trim(fname))
    !
    !     remember mesh name
    !
          meshName = fname
    
          lun=11
    
          OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)
    
          CALL rOneTL(lun,OneLine)
          DO WHILE (OneLine(1:3).NE.'EOF')
    !
    !       Read nodes
    !
            IF (OneLine.EQ."$Nodes") THEN
              CALL rOneTL(lun,OneLine)
              READ(Oneline,*) n
              CALL InitNodes(n)
              DO i=1,nnodes
                CALL rOneTL(lun,OneLine)
                READ(Oneline,*) j,node(j)%x(1),node(j)%x(2),node(j)%x(3)
              END DO
            END IF
    !
    !       Read elements
    !
            IF (OneLine.EQ."$Elements") THEN
              CALL rOneTL(lun,OneLine)
              READ(Oneline,*) n
              CALL initElement(n)
              DO i=1,nelem
                CALL rOneTL(lun,OneLine)
                READ(Oneline,*) j,element(j)%type,nd,element(j)%bcid
    
                ! determine element type and number of nodes per element
                IF (element(j)%type.EQ.2) THEN ! three node triangle
                  element(j)%nno = 3 ! number of nodes in element
                ELSE IF (element(j)%type.EQ.3) THEN ! four node Quadrangle
                  element(j)%nno = 4 ! number of nodes in element
                ELSE
                  CALL WriteToLog("Error :: Element type not supported!")
                  CALL StopProgram
                END IF
                ! allocate connectivity
                ALLOCATE(element(j)%con(element(j)%nno))
                READ(Oneline,*) j,element(j)%type,nd,element(j)%bcid,(dummy,k=1,nd-1),(element(j)%con(l),l=1,element(j)%nno)
    
              END DO
            END IF
    !
    !       Read wall definitions
    !
            IF (OneLine.EQ."$PhysicalNames") THEN
              CALL rOneTL(lun,OneLine)
              READ(Oneline,*) n
              CALL InitWall(n)
              DO i=1,nofw
                CALL rOneTL(lun,OneLine)
                READ(Oneline,*) nd,wall(i)%id,wall(i)%name
              END DO
            END IF
    
    
            CALL rOneTL(lun,OneLine)
          END DO
    
          CLOSE (lun)
    
          RETURN
    
    
    10    CALL WriteToLog("Error :: could not open v2 mesh file!")
          CALL StopProgram
    
          END SUBROUTINE
    


! ----------------------------------------------------------------------
!
  SUBROUTINE getMeshTypev4(fname,meshType)
!
! ----------------------------------------------------------------------
!       
    USE mMesh
    USE mPar
    USE mCommon
    IMPLICIT NONE
              
    integer numPoints,numCurves,numSurfaces,numVolumes,lun
    CHARACTER*(*) fname  
    character(20) meshType
    CHARACTER(255) OneLine
    
    
    CALL WriteToLog("Determining type of v4 mesh: "//trim(fname))
          
    lun=11

    OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)
    
    CALL rOneTL(lun,OneLine)
    DO WHILE (OneLine(1:3).NE.'EOF')
      IF (OneLine.EQ."$Entities") THEN    
        CALL rOneTL(lun,OneLine)
        READ(Oneline,*) numPoints,numCurves,numSurfaces,numVolumes
      end if
      CALL rOneTL(lun,OneLine)
    end do
    close(lun)

    if (numVolumes.EQ.0) then
      meshType = "surface"
    else
      meshType = "volume"
    end if
    CALL WriteToLog("Found "//trim(meshType)//" mesh!")

    RETURN


10    CALL WriteToLog("Error :: could not open mesh file!")
    CALL StopProgram

  END SUBROUTINE

! ----------------------------------------------------------------------
!
SUBROUTINE ReadGMSHv4Surface(fname)
!
! ----------------------------------------------------------------------
  USE mMesh
  USE mPar
  USE mCommon
  IMPLICIT NONE
            
  CHARACTER*(*) fname
  CHARACTER(255) OneLine
  INTEGER lun,dummy
  REAL(rk) rdummy
  INTEGER i,j,n,k,l,nd
            
            
  INTEGER numPoints,numCurves,numSurfaces,numVolumes
  INTEGER surfaceTag,physicalTag
  INTEGER numEntityBlocks,numElements,minElementTag,maxElementTag
  INTEGER entityDim,entityTag,elType,numElementsInBlock    
  INTEGER numNodes, minNodeTag, maxNodeTag,parametric,numNodesInBlock       
          
  INTEGER, ALLOCATABLE :: ent2phyMap(:,:)
  INTEGER, ALLOCATABLE :: nodeTags(:)
          
  CALL WriteToLog("Reading v4 surface mesh: "//trim(fname))
  !
  !     remember mesh name
  !
  meshName = fname
          
  lun=11

  OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)
  
  CALL rOneTL(lun,OneLine)
  DO WHILE (OneLine(1:3).NE.'EOF')
!
!       Read entities
!

    IF (OneLine.EQ."$Entities") THEN    
      CALL rOneTL(lun,OneLine)
      READ(Oneline,*) numPoints,numCurves,numSurfaces,numVolumes

      DO i = 1,numPoints
!      pointTag(int) X(double) Y(double) Z(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!      ...

        CALL rOneTL(lun,OneLine)
      END DO

      DO i = 1,numCurves
!      curveTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundingPoints(size_t) pointTag(int) ...
!      ...      
        CALL rOneTL(lun,OneLine)
      END DO

      ALLOCATE (ent2phyMap(numSurfaces,2))
      DO i = 1,numSurfaces
!      surfaceTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundingCurves(size_t) curveTag(int) ...
!      ...        
        CALL rOneTL(lun,OneLine)
        READ(Oneline,*) surfaceTag,rdummy,rdummy,rdummy,rdummy,rdummy,rdummy,dummy,physicalTag
        ent2phyMap(i,1)=surfaceTag
        ent2phyMap(i,2)=abs(physicalTag) ! negativno pri zasukani normali
      END DO

      DO i = 1,numVolumes
!      volumeTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundngSurfaces(size_t) surfaceTag(int) ...
!      ...
        CALL rOneTL(lun,OneLine)
      END DO

    END IF


!
!       Read nodes
!
    IF (OneLine.EQ."$Nodes") THEN
      CALL rOneTL(lun,OneLine)
      READ(Oneline,*) numEntityBlocks, numNodes, minNodeTag, maxNodeTag
      CALL InitNodes(numNodes)

      DO i = 1,numEntityBlocks
        CALL rOneTL(lun,OneLine)
          READ(Oneline,*) entityDim,entityTag,parametric,numNodesInBlock       

        ALLOCATE (nodeTags(numNodesInBlock))
        DO k = 1,numNodesInBlock
          CALL rOneTL(lun,OneLine)
          READ(Oneline,*) nodeTags(k)
        END DO

        DO k = 1,numNodesInBlock
          j = nodeTags(k)
          CALL rOneTL(lun,OneLine) 
          READ(Oneline,*) node(j)%x(1),node(j)%x(2),node(j)%x(3)
        END DO
        DEALLOCATE (nodeTags)

      END DO
    END IF
!
!       Read elements
!
    IF (OneLine.EQ."$Elements") THEN
      CALL rOneTL(lun,OneLine)
      READ(Oneline,*) numEntityBlocks,numElements,minElementTag,maxElementTag
      CALL initElement(numElements)

      DO i = 1,numEntityBlocks
        CALL rOneTL(lun,OneLine)
        READ(Oneline,*) entityDim,entityTag,elType,numElementsInBlock       
        ! find physicalTab basedon entitiy tag
        DO k = 1,numSurfaces
          ! if (ent2phyMap(i,1).EQ.entityTag) physicalTag = ent2phyMap(i,2) 
          if (ent2phyMap(k,1).EQ.entityTag) physicalTag = ent2phyMap(k,2)
        END DO

        DO k = 1,numElementsInBlock
          CALL rOneTL(lun,OneLine)

          !elementTag(size_t) nodeTag(size_t) ...          
          READ(Oneline,*) j
          element(j)%type = elType
          element(j)%bcid = physicalTag

          ! determine element type and number of nodes per element
          IF (element(j)%type.EQ.2) THEN ! three node triangle
            element(j)%nno = 3 ! number of nodes in element
          ELSE IF (element(j)%type.EQ.3) THEN ! four node Quadrangle
            element(j)%nno = 4 ! number of nodes in element
          ELSE
            CALL WriteToLog("Error :: Element type not supported!")
            CALL StopProgram
          END IF
          ! allocate connectivity
          ALLOCATE(element(j)%con(element(j)%nno))
!         elementTag(size_t) nodeTag(size_t) ...          
          READ(Oneline,*) j,(element(j)%con(l),l=1,element(j)%nno)
!
        END DO
      END DO
    END IF
!
!       Read wall definitions
!
    IF (OneLine.EQ."$PhysicalNames") THEN
      CALL rOneTL(lun,OneLine)
      READ(Oneline,*) n
      CALL InitWall(n)
      DO i=1,nofw
        CALL rOneTL(lun,OneLine)
        READ(Oneline,*) nd,wall(i)%id,wall(i)%name
      END DO
    END IF


    CALL rOneTL(lun,OneLine)
  END DO

  CLOSE (lun)

  DEALLOCATE (ent2phyMap)

  RETURN


10    CALL WriteToLog("Error :: could not open mesh file!")
  CALL StopProgram

END SUBROUTINE
        

!      $MeshFormat // same as MSH version 2
!      version(ASCII double; currently 4.1)
!        file-type(ASCII int; 0 for ASCII mode, 1 for binary mode)
!        data-size(ASCII int; sizeof(size_t))
!      < int with value one; only in binary mode, to detect endianness >
!    $EndMeshFormat
!    
!    $PhysicalNames // same as MSH version 2
!      numPhysicalNames(ASCII int)
!      dimension(ASCII int) physicalTag(ASCII int) "name"(127 characters max)
!      ...
!    $EndPhysicalNames
!    
!    $Entities
!      numPoints(size_t) numCurves(size_t)
!        numSurfaces(size_t) numVolumes(size_t)
!      pointTag(int) X(double) Y(double) Z(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!      ...
!      curveTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundingPoints(size_t) pointTag(int) ...
!      ...
!      surfaceTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundingCurves(size_t) curveTag(int) ...
!      ...
!      volumeTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundngSurfaces(size_t) surfaceTag(int) ...
!      ...
!    $EndEntities
!    
!    
!    $Nodes
!      numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
!      entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
!        nodeTag(size_t)
!        ...
!        x(double) y(double) z(double)
!           < u(double; if parametric and entityDim >= 1) >
!           < v(double; if parametric and entityDim >= 2) >
!           < w(double; if parametric and entityDim == 3) >
!        ...
!      ...
!    $EndNodes
!    
!    $Elements
!      numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
!      entityDim(int) entityTag(int) elementType(int; see below) numElementsInBlock(size_t)
!        elementTag(size_t) nodeTag(size_t) ...
!        ...
!      ...
!    $EndElements
!    
            

! -----------------------------------------------------------------------------
SUBROUTINE MeshStretch()
!
!     $: stretches mesh
!
! -----------------------------------------------------------------------------
  USE mMesh
  USE mPar
  IMPLICIT NONE
  INTEGER i,j

  DO i=1,nnodes
    DO j=1,3
      node(i)%x(j)=node(i)%x(j)*parMTstr(j)
    END DO
  END DO

END
              
              
! -----------------------------------------------------------------------------
      SUBROUTINE MeshCylStretch()
!
!     $: stretches mesh
!
! -----------------------------------------------------------------------------
  USE mMesh
  USE mPar
  USE mCommon
  IMPLICIT NONE
              
  INTEGER i,j
  REAL(rk) rp,r
              
! Z coordinate
  DO i=1,nnodes
    node(i)%x(3)=node(i)%x(3)*parMTcyl(6)/parMTcyl(3)
  END DO

! X,Y
  DO i=1,nnodes

   r = SQRT(node(i)%x(1)**2+node(i)%x(2)**2)

   IF (r.LE.parMTcyl(1)) THEN  ! small cylinder
      DO j=1,2
        node(i)%x(j)=node(i)%x(j)*parMTcyl(4)/parMTcyl(1)
      END DO
    ELSE ! kolobar

     rp = parMTcyl(4) + (r-parMTcyl(1))*(parMTcyl(5)-parMTcyl(4))/(parMTcyl(2)-parMTcyl(1))

     DO j=1,2
        node(i)%x(j)=node(i)%x(j)*rp/r
      END DO
    END IF

 END DO


END

! -----------------------------------------------------------------------------
SUBROUTINE MeshGetExtents()
!
!     $: stretches mesh
!
! -----------------------------------------------------------------------------
  USE mMesh
  USE mPar
  IMPLICIT NONE
  INTEGER i,j

  mextMin= 1.0D10
  mextMax=-1.0D10

  DO i=1,nnodes
    DO j=1,3
      IF (node(i)%x(j) .GT. mextMax(j) ) mextMax(j) =  node(i)%x(j)
      IF (node(i)%x(j) .LT. mextMin(j) ) mextMin(j) =  node(i)%x(j)
    END DO
  END DO

END
                      
! -----------------------------------------------------------------------------
SUBROUTINE MeshTranslate()
!
!     $: translates mesh
!
! -----------------------------------------------------------------------------
  USE mMesh
  USE mPar
  IMPLICIT NONE

  INTEGER i,j

  DO i=1,nnodes
    DO j=1,3
      node(i)%x(j)=node(i)%x(j)+parMTtra(j)
    END DO
  END DO

END
                      
                      
                      
! -----------------------------------------------------------------------------
SUBROUTINE MeshRotate()
!
!     $: rotates mesh so that (1,0,0) points to inp%MTrot
!
! -----------------------------------------------------------------------------
  USE mMesh
  USE mPar
  USE mCommon
  IMPLICIT NONE

  INTEGER i

  REAL(rk) R(3,3),RT(3,3),vr(3)

  CALL NormVector(parMTrot)
  CALL GetRotationMatrix(R,RT,parMTrot)

  DO i=1,nnodes
    vr=MATMUL(RT,node(i)%x)
    node(i)%x=vr
  END DO

END
                      
                      
! -----------------------------------------------------------------------------
subroutine sphere2superE()
!
!   Transform sphere to superellipsoid
!     
! -----------------------------------------------------------------------------  
  use mPar
  use mMesh
  implicit none
  real(rk) selX(3),r
  integer i
  do i=1,nnodes
      r = node(i)%x(1)**2 + node(i)%x(2)**2 + node(i)%x(3)**2
      if (r.LT.100.0_rk) then
          call mapSphereToSuperE(parMTsel(1),parMTsel(2),parMTsel(3),parMTsel(4),parMTsel(5),selX,node(i)%x)
          node(i)%x = selX
      end if
  end do
end subroutine
                                                


! -----------------------------------------------------------------------------------------
subroutine projectPoints2superE()
  !
  !     Project mesh points to supere surface
  !
  ! -----------------------------------------------------------------------------------------  
    use mPar
    use mMesh
    implicit none
    real(rk) r
    integer i,k,maxk
    
    real(rk) u0,u,res,np(3)
    real(rk) a,b,c,e1,e2,f,df
    
    !PRINT *, "Entering Subroutine :: projectPoints2superE"
    
    DO i=1,nnodes
      r = node(i)%x(1)**2 + node(i)%x(2)**2 + node(i)%x(3)**2
      IF (r.LT.10.0_rk) THEN
    
          np(1) = node(i)%x(1)
          np(2) = node(i)%x(2)
          np(3) = node(i)%x(3)
    
          CALL NormVector(np)
    
          a = parMTsel(1)
          b = parMTsel(2)
          c = parMTsel(3)
          e1 = parMTsel(4)
          e2 = parMTsel(5)
          
          u0 = 1.0D0
          res = 1.0D0
        
          k = 0
          maxk = 1000
          DO WHILE ( res > 1.0e-3 )
        
            call calF(a,b,c,e1,e2,u0,np,f)
            call calDF(a,b,c,e1,e2,u0,np,df)
            u = u0 - f / df
            !PRINT *, "u=",u,f(u0,np),df(u0,np)
        
            call calF(a,b,c,e1,e2,u,np,f)
            res = abs( f )
            u0 = u
        
            IF (k > maxk) EXIT
            k = k+1
        
          END DO
        
          node(i)%x(1) = u * np(1)
          node(i)%x(2) = u * np(2)
          node(i)%x(3) = u * np(3)
        
          !PRINT *, "u=",u
          !PRINT *, "res=",res
          !PRINT *, "k=",k
          !PRINT *, "p=",node(i)%x 
    
      END IF
    END DO
    
  end subroutine
          
  subroutine calF(a,b,c,e1,e2,t,p,f)
    USE mPar
    IMPLICIT NONE
     
     real(rk) t,p(3),f
     real(rk) x,y,z,a,b,c,e1,e2
    !real(rk) sign_x, sign_y, sign_z
    !sign_x = SIGN( real(1,rk), x )
    !sign_y = SIGN( real(1,rk), y )
    !sign_z = SIGN( real(1,rk), z )
     x = abs( p(1) )
     y = abs( p(2) )
     z = abs( p(3) )
     f = ( (t*x/a)**(2/e2) + (t*y/b)**(2/e2) )**(e2/e1) + (t*z/c)**(2/e1) - 1
     !print *, "f+1=",f+1
     
     RETURN
   end subroutine
   
  subroutine  calDF(a,b,c,e1,e2,t,p,df)
   use mPar
   implicit none 
   real(rk) df,t,p(3)
   real(rk) x,y,z,a,b,c,e1,e2
  
    x = abs( p(1) )
    y = abs( p(2) )
    z = abs( p(3) )
    
    df = 2 * ( ((t*x/a)**(2/e2) + (t*y/b)**(2/e2))**(e2/e1) + (t*z/c)**(2/e1) )
    df = df / (e1*t)
  
    RETURN
  end subroutine   



! ----------------------------------------------------------------------
!
  SUBROUTINE ReadGMSHv4Volume(fname)
!
! ----------------------------------------------------------------------
  USE mDomainMesh
  USE mPar
  USE mCommon
  IMPLICIT NONE
                
  CHARACTER*(*) fname
  CHARACTER(255) OneLine
  INTEGER lun,dummy
  REAL(rk) rdummy
  INTEGER i,j,n,k,l,nd,ne
                
                
  INTEGER numPoints,numCurves,numSurfaces,numVolumes
  INTEGER surfaceTag,physicalTag,volumeTag
  INTEGER numEntityBlocks,numElements,minElementTag,maxElementTag
  INTEGER entityDim,entityTag,elType,numElementsInBlock    
  INTEGER numNodes, minNodeTag, maxNodeTag,parametric,numNodesInBlock       
              
  INTEGER, ALLOCATABLE :: ent2phyMap(:,:)
  INTEGER, ALLOCATABLE :: nodeTags(:)
          
  CALL WriteToLog("Reading v4 volume mesh: "//trim(fname))
!
! remember mesh name
!
  DMmeshName = fname
  lun=11


!
! Figure out the number of domain elements in the mesh
!
  OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)
  ne = 0
  CALL rOneTL(lun,OneLine)
  DO WHILE (OneLine(1:3).NE.'EOF')
    IF (OneLine.EQ."$Elements") THEN
      CALL rOneTL(lun,OneLine)
      READ(Oneline,*) numEntityBlocks,numElements,minElementTag,maxElementTag
      DO i = 1,numEntityBlocks
        CALL rOneTL(lun,OneLine)
        READ(Oneline,*) entityDim,entityTag,elType,numElementsInBlock       
        if (elType.eq.4) ne = ne + numElementsInBlock
        DO k = 1,numElementsInBlock
          CALL rOneTL(lun,OneLine)
        END DO
      END DO
    END IF
    CALL rOneTL(lun,OneLine)
  END DO
  CLOSE(lun)
  CALL DMinitElement(ne)
             
!
! Read mesh
!
    
  OPEN (unit=lun,file=TRIM(fname),status="OLD",ERR=10)
      
  CALL rOneTL(lun,OneLine)
  DO WHILE (OneLine(1:3).NE.'EOF')
!
!   Read entities
!
    IF (OneLine.EQ."$Entities") THEN    
      CALL rOneTL(lun,OneLine)
      READ(Oneline,*) numPoints,numCurves,numSurfaces,numVolumes
    
      DO i = 1,numPoints
!      pointTag(int) X(double) Y(double) Z(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!      ...
    
        CALL rOneTL(lun,OneLine)
      END DO
    
      DO i = 1,numCurves
!      curveTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundingPoints(size_t) pointTag(int) ...
!      ...      
        CALL rOneTL(lun,OneLine)
      END DO

      ALLOCATE (ent2phyMap(numSurfaces,2))
      DO i = 1,numSurfaces
!      surfaceTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundingCurves(size_t) curveTag(int) ...
!      ...        
        CALL rOneTL(lun,OneLine)
        READ(Oneline,*) surfaceTag,rdummy,rdummy,rdummy,rdummy,rdummy,rdummy,dummy,physicalTag
        ent2phyMap(i,1)=surfaceTag
        ent2phyMap(i,2)=abs(physicalTag) ! negativno pri zasukani normali
      END DO
    
      DO i = 1,numVolumes
!      volumeTag(int) minX(double) minY(double) minZ(double)
!        maxX(double) maxY(double) maxZ(double)
!        numPhysicalTags(size_t) physicalTag(int) ...
!        numBoundngSurfaces(size_t) surfaceTag(int) ...
!      ...
        CALL rOneTL(lun,OneLine)
        read(OneLine,*) volumeTag           
      END DO
    
    END IF
    
    
!
!   Read nodes
!
        IF (OneLine.EQ."$Nodes") THEN
          CALL rOneTL(lun,OneLine)
          READ(Oneline,*) numEntityBlocks, numNodes, minNodeTag, maxNodeTag
          CALL DMInitNodes(numNodes)
    
          DO i = 1,numEntityBlocks
            CALL rOneTL(lun,OneLine)
              READ(Oneline,*) entityDim,entityTag,parametric,numNodesInBlock       
    
            ALLOCATE (nodeTags(numNodesInBlock))
            DO k = 1,numNodesInBlock
              CALL rOneTL(lun,OneLine)
              READ(Oneline,*) nodeTags(k)
            END DO
    
            DO k = 1,numNodesInBlock
              j = nodeTags(k)
              CALL rOneTL(lun,OneLine) 
              READ(Oneline,*) DMnode(j)%x(1),DMnode(j)%x(2),DMnode(j)%x(3)
            END DO
            DEALLOCATE (nodeTags)
    
          END DO
        END IF
    !
    !       Read elements
    !
        IF (OneLine.EQ."$Elements") THEN
          CALL rOneTL(lun,OneLine)
          READ(Oneline,*) numEntityBlocks,numElements,minElementTag,maxElementTag
          ne = 0
          DO i = 1,numEntityBlocks
            CALL rOneTL(lun,OneLine)
            READ(Oneline,*) entityDim,entityTag,elType,numElementsInBlock       
            ! find physicalTab basedon entitiy tag
            DO k = 1,numSurfaces
              ! if (ent2phyMap(i,1).EQ.entityTag) physicalTag = ent2phyMap(i,2) 
              if (ent2phyMap(k,1).EQ.entityTag) physicalTag = ent2phyMap(k,2)
            END DO
    
            DO k = 1,numElementsInBlock
              !elementTag(size_t) nodeTag(size_t) ...          
              CALL rOneTL(lun,OneLine)              
              ! skip non domain elements                
              if (elType.EQ.4) then ! "Point","Line","Triangle","Quadrangle","Tetrahedron","Pyramid","Prism","Hexahedron"                
                ne = ne + 1     
                DMelement(ne)%type = elType
                DMelement(ne)%nno = 4 ! number of nodes in element 
                ! allocate connectivity
                ALLOCATE(DMelement(ne)%con(DMelement(ne)%nno))
                !         elementTag(size_t) nodeTag(size_t) ...          
                READ(Oneline,*) j,(DMelement(ne)%con(l),l=1,DMelement(ne)%nno)
              end if              

            END DO
          END DO
        END IF

        CALL rOneTL(lun,OneLine)
      END DO
    
      CLOSE (lun)
    
      DEALLOCATE (ent2phyMap)

      RETURN
    
    10    CALL WriteToLog("Error :: could not open domain mesh file!")
      CALL StopProgram
    
    END SUBROUTINE
            
    
    !      $MeshFormat // same as MSH version 2
    !      version(ASCII double; currently 4.1)
    !        file-type(ASCII int; 0 for ASCII mode, 1 for binary mode)
    !        data-size(ASCII int; sizeof(size_t))
    !      < int with value one; only in binary mode, to detect endianness >
    !    $EndMeshFormat
    !    
    !    $PhysicalNames // same as MSH version 2
    !      numPhysicalNames(ASCII int)
    !      dimension(ASCII int) physicalTag(ASCII int) "name"(127 characters max)
    !      ...
    !    $EndPhysicalNames
    !    
    !    $Entities
    !      numPoints(size_t) numCurves(size_t)
    !        numSurfaces(size_t) numVolumes(size_t)
    !      pointTag(int) X(double) Y(double) Z(double)
    !        numPhysicalTags(size_t) physicalTag(int) ...
    !      ...
    !      curveTag(int) minX(double) minY(double) minZ(double)
    !        maxX(double) maxY(double) maxZ(double)
    !        numPhysicalTags(size_t) physicalTag(int) ...
    !        numBoundingPoints(size_t) pointTag(int) ...
    !      ...
    !      surfaceTag(int) minX(double) minY(double) minZ(double)
    !        maxX(double) maxY(double) maxZ(double)
    !        numPhysicalTags(size_t) physicalTag(int) ...
    !        numBoundingCurves(size_t) curveTag(int) ...
    !      ...
    !      volumeTag(int) minX(double) minY(double) minZ(double)
    !        maxX(double) maxY(double) maxZ(double)
    !        numPhysicalTags(size_t) physicalTag(int) ...
    !        numBoundngSurfaces(size_t) surfaceTag(int) ...
    !      ...
    !    $EndEntities
    !    
    !    
    !    $Nodes
    !      numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)
    !      entityDim(int) entityTag(int) parametric(int; 0 or 1) numNodesInBlock(size_t)
    !        nodeTag(size_t)
    !        ...
    !        x(double) y(double) z(double)
    !           < u(double; if parametric and entityDim >= 1) >
    !           < v(double; if parametric and entityDim >= 2) >
    !           < w(double; if parametric and entityDim == 3) >
    !        ...
    !      ...
    !    $EndNodes
    !    
    !    $Elements
    !      numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)
    !      entityDim(int) entityTag(int) elementType(int; see below) numElementsInBlock(size_t)
    !        elementTag(size_t) nodeTag(size_t) ...
    !        ...
    !      ...
    !    $EndElements
    !    
                  
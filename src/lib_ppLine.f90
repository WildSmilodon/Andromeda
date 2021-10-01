! -----------------------------------------------------------------------------

       SUBROUTINE linePostProcessing()

!
!     $ computes points of the intersection of the line, starting at the 
!     $ line defined in the .inp, and the plane defined by elements (triangles)
!
! -----------------------------------------------------------------------------
      USE mMesh
      USE mPar
      USE mEqns
      USE mProfile
      USE mCommon
      IMPLICIT NONE


      REAL(rk)              :: t,d, denom
      REAL(rk)              :: interPoints(3)
      INTEGER           :: i, j, k,  m
      REAL(rk)           :: valAtInter
      REAL(rk)              :: primD
      LOGICAL           :: oneOrNil

!     READ IN THE POINTS OF PROFILES
      CALL ReadINPProf()
!     ALLOCATE THE VECTORS FOR LINES
      CALL lines_init(neq)

!     Defining points on the lines    
      DO k = 1, noLines
        DO i=1, prof(k)%numOfDiv
          t = 1.0_rk/(prof(k)%numOfDiv)*DBLE(i) ! parameter of the parametrised line
          DO j = 1, 3
              prof(k)%ptsOnLine(i+1,j) = prof(k)%firstLine(1,j)  * (1-t)+prof(k)%firstLine(2,j) * t
              prof(k)%ptsOnProj(i+1,j) = prof(k)%secndLine(1,j)  * (1-t)+prof(k)%secndLine(2,j) * t              
          END DO
        END DO
      END DO


      
!     Adding the first point into the vector
      DO i = 1, noLines
        prof(i)%ptsOnLine(1,:) = prof(i)%firstLine(1,:)
        prof(i)%ptsOnProj(1,:) = prof(i)%secndLine(1,:)
        prof(i)%results=-999999999.0_rk
      END DO

!       Compute the length ()
      DO k = 1, noLines
        DO i = 1, prof(k)%numOfDiv+1
          prof(k)%results(1,i) = NORM2(prof(k)%ptsOnLine(i,:) - prof(k)%ptsOnLine(1,:))
        END DO
      END DO
!     Compute intersections with each plane and find the one inside the triangle
      DO m = 1, noLines
        DO i = 1,prof(m)%numOfDiv+1
          primD = 1.0_rk
          DO j = 1,nelem
            IF (ANY(TRIM(wall(element(j)%bcid)%name).EQ.prof(m)%WallName(:))) THEN 
            d = DOT_PRODUCT((node(element(j)%con(2))%x(:)-prof(m)%ptsOnLine(i,:)), element(j)%normal(:))
            denom = DOT_PRODUCT((prof(m)%ptsOnProj(i,:)-prof(m)%ptsOnLine(i,:)), element(j)%normal(:))
              IF (denom.NE.0.0_rk) THEN
                d = d/denom
                interPoints = d*(prof(m)%ptsOnProj(i,:) - prof(m)%ptsOnLine(i,:)) + prof(m)%ptsOnLine(i,:)
                CALL isInTriangle(j,interPoints,oneOrNil)
                IF ( (oneOrNil.EQV..TRUE.) .AND. (d.LE.1.0_rk) .AND. (d.LT.primD)) THEN
                    primD = d
                    DO k = 1, neq
                      CALL compLambda(j, interPoints, valAtInter, k)
                      prof(m)%results( k+1, i ) = valAtInter
                    END DO
                    prof(m)%ptsOnBody(i,:) = interPoints
                END IF
              END IF
            END IF
          END DO
        END DO
      END DO


      CALL outputPostProc()
      DO i = 1, noLines
        DEALLOCATE (prof(i)%ptsOnLine)
        DEALLOCATE (prof(i)%ptsOnProj)
        DEALLOCATE (prof(i)%ptsOnBody)
        DEALLOCATE (prof(i)%results)
      END DO

      END

! -----------------------------------------------------------------------------
         SUBROUTINE isInTriangle(nTriangle, P, val)
!
!     $: Checks whether point is in triangle or not
!
! -----------------------------------------------------------------------------
    
      USE mPar
      USE mMesh
      USE mCommon
      implicit none
      REAL(rk), DIMENSION (3)    :: P0, P1, P2, P, new
      INTEGER                :: nTriangle
      INTEGER                :: i, samPlain
      REAL(rk)                   :: v0(3), v1(3), v2(3)
      REAL(rk)                   :: d00, d01, d11, d20, d21, denom
      REAL(rk)                   :: lambda1, lambda2, lambda3
      REAL(rk)                   :: pointsTog(3,3)
      REAL(rk)                   :: pointsTogNew(4,3)
      LOGICAL                :: val
          

      CALL threeDrotation( pointsTogNew, nTriangle, P )

      P0   = pointsTogNew(1,:)
      P1   = pointsTogNew(2,:)
      P2   = pointsTogNew(3,:)
      new  = pointsTogNew(4,:)

         lambda1 = ((P1(2)-P2(2))*(new(1)-P2(1))+(P2(1)-P1(1))*(new(2)-P2(2)))/ &
       ((P1(2)-P2(2))*(P0(1)-P2(1))+(P2(1)-P1(1))*(P0(2)-P2(2)))
         lambda2 = ((P2(2)-P0(2))*(new(1)-P2(1))+(P0(1)-P2(1))*(new(2)-P2(2)))/ &
       ((P1(2)-P2(2))*(P0(1)-P2(1))+(P2(1)-P1(1))*(P0(2)-P2(2)))

      IF (lambda1.LT.-1.0D-12.OR.lambda2.LT.-1.0D-12.OR.lambda1.GT.1.0_rk+ &
        1.0E-12.OR.lambda2.GT.1.0_rk+1.0E-12 &
                .OR.lambda1+lambda2.GT.1.0_rk) THEN
          val = .FALSE.
      ELSE
          val = .TRUE.
      END IF


      END
! -----------------------------------------------------------------------------
      SUBROUTINE compLambda(nTriangle, P, val, noEq)
!
!     $: Computes the values of the things on the sides
!
! -----------------------------------------------------------------------------
      USE mPar
      USE mMesh
      USE mEqns
      implicit none
      
      REAL(rk)  lambda1,lambda2,lambda3
      REAL(rk)  P1(3), P2(3), P3(3),P(3), v0(3), v1(3), v2(3)
      INTEGER nTriangle, noEq
      REAL(rk)  val
      REAL(rk)  valP1, valP2, valP3
      REAL(rk)  d00, d01, d11, d20, d21, denom

      P1 = node(element(nTriangle)%con(1))%x(:)
      P2 = node(element(nTriangle)%con(2))%x(:)
      P3 = node(element(nTriangle)%con(3))%x(:)

      valP1 = eqn(noEq)%u(element(nTriangle)%con(1))
      valP2 = eqn(noEq)%u(element(nTriangle)%con(2))
      valP3 = eqn(noEq)%u(element(nTriangle)%con(3))

      v0 = P2 - P1
      v1 = P3 - P1
      v2 = P  - P1

      d00 = DOT_PRODUCT(v0, v0)
      d01 = DOT_PRODUCT(v0, v1)
      d11 = DOT_PRODUCT(v1, v1)
      d20 = DOT_PRODUCT(v2, v0)
      d21 = DOT_PRODUCT(v2, v1)
      denom = d00 * d11 - d01 * d01

      lambda1 = (d11 * d20 - d01 * d21) / denom
      lambda2 = (d00 * d21 - d01 * d20) / denom
      lambda3 = 1 - lambda1 - lambda2
      !Compute the value on the very spot
      val = lambda3 * valP1 + lambda1 * valP2 + lambda2 * valP3        
      END

! -----------------------------------------------------------------------------

      SUBROUTINE outputPostProc()

!
!     $: PostProcessing output
!
! -----------------------------------------------------------------------------
        USE mMesh
        USE mPar
        USE mEqns
        USE mProfile

        IMPLICIT NONE

        INTEGER itp,i,of,k,j
        CHARACTER*255 vrstica,tmp



        itp=96
        WRITE(vrstica,*)
        DO k = 1,noLines
          OPEN (itp,FILE="and.prof."//TRIM(prof(k)%name)//".txt",STATUS='UNKNOWN')
          WRITE (itp,*) "x y z l",(" "//TRIM(eqn(i)%name),i=1,neq)

          DO j = 1,  prof(k)%numOfDiv + 1
            IF (prof(k)%results(2,j).NE.-999999999.0_rk) THEN
              WRITE (itp,*) prof(k)%ptsOnBody(j,1),prof(k)%ptsOnBody(j,2),prof(k)%ptsOnBody(j,3), &
                          (prof(k)%results(i,j), i = 1,neq + 1)
            END IF
          END DO
          CLOSE (itp)
        END DO

        END
! ----------------------------------------------------------------------
!
       SUBROUTINE ReadINPProf()
!
! ----------------------------------------------------------------------
      USE mProfile
      USE mPar
      IMPLICIT NONE

      CHARACTER(255) OneLine,KeyWord,line
      INTEGER lun, dummy,wallsNum
      INTEGER i,j

!     ALLOCATE THE MOD profiles
      CALL prof_init(parNumOfProf)

      lun=10
      OPEN (UNIT=lun ,FILE=parInputFileName,STATUS='OLD',ERR=10)
      
      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ(Oneline,*) KeyWord
        IF (KeyWord.EQ."PPLI") THEN          
          DO i = 1, parNumOfProf
            CALL rOneTL(lun,OneLine)
            READ(OneLine,*) prof(i)%name, wallsNum
            CALL walls_init(wallsNum, i)
            READ(OneLine,*)prof(i)%name, dummy, (prof(i)%WallName(j),j = 1,wallsNum), &
              prof(i)%firstLine(1,1), prof(i)%firstLine(1,2), prof(i)%firstLine(1,3), &
              prof(i)%firstLine(2,1), prof(i)%firstLine(2,2), prof(i)%firstLine(2,3), &
              prof(i)%secndLine(1,1), prof(i)%secndLine(1,2), prof(i)%secndLine(1,3), &
              prof(i)%secndLine(2,1), prof(i)%secndLine(2,2), prof(i)%secndLine(2,3), &
              prof(i)%numOfDiv
          END DO          
          EXIT
        END IF
        
        CALL rOneTL(lun,OneLine)
      END DO

      CLOSE (lun)
      
      RETURN
10    CALL WriteToLog("Error reading profile points ") 
      END

! ----------------------------------------------------------------------
!
      SUBROUTINE threeDrotation( newPoints, nTriangle, P )
!
! ----------------------------------------------------------------------
        USE mMesh
        USE mPar
        USE mEqns
        USE mProfile
        USE mCommon
        IMPLICIT NONE
        REAL(rk), DIMENSION(3,3)    ::  points
        REAL(rk), DIMENSION(4,3)    ::  newPoints
        REAL(rk), DIMENSION(3)      ::  angles
        REAL(rk), DIMENSION(3,3)    ::  normsMainPlain
        REAL(rk), DIMENSION(3)      ::  N
        REAL(rk), DIMENSION(3)      ::  plainNorm
        REAL(rk), DIMENSION(4)      ::  quat
        REAL(rk), DIMENSION(3,3)    ::  rotMat
        REAL(rk), DIMENSION(3)      ::  P, nov
        INTEGER                 ::  i, j, nTriangle
        REAL(rk)                    ::  angle2lines
!	First define points
        DO i = 1, 3
            points( i, : ) = node( element( nTriangle )%con(i))%x(:)
        END DO
!	Define normal, which is new z axis
        plainNorm = element(nTriangle)%normal(:)

        CALL    cartesianNormals( normsMainPlain )

!	Move coordinate system to the first point
        nov  = P
        nov  = nov - points( 1, : )
        DO i = 3, 1, -1
            points( i, : ) = points( i, : ) - points( 1, : )
        END DO


!       DEFINE ANGLES
!       FIRST COMPUTE N (z is in the direction of plain normal)
!       x axis direction \vec{13}

        CALL cross_product( normsMainPlain( 1, :), plainNorm ( : ),  N )


       	! Yaw, pitch and roll
        IF ( NORM2(N(:)).EQ.0.0_rk) THEN
            N ( : )   = normsMainPlain ( 3, : )
            angles(2) = 0.0_rk
        ELSE
        angles(2) = angle2lines(plainNorm(:),   normsMainPlain( 1,:) )
        END IF
        angles(1) = angle2lines( normsMainPlain(3,:),         N(:)   )
        angles(3) = angle2lines( N ( : ),     points( 3, : )         )
        


        ! COMPUTE QUATERNION
        CALL quaternion(angles, quat, rotMat)

        DO i = 1, 3
            newPoints( i, : ) = MATMUL( rotMat, points( i,: ) )
        END DO
            newPoints( 4, : ) = MATMUL( rotMat, nov )

      END SUBROUTINE
! ----------------------------------------------------------------------
!
      REAL(rk) FUNCTION angle2lines( line1, line2)
!
! ----------------------------------------------------------------------
          USE mCommon
          REAL(rk)               :: line1(3)
          REAL(rk)               :: line2(3)
          REAL(rk)               :: dummy

          IF (NORM2(line1)*NORM2(line2).EQ.0.0_rk) THEN
             angle2lines = 0.0_rk
          ELSE
              dummy = DOT_PRODUCT(line1,line2)/(NORM2(line1)*NORM2(line2))
              IF(dummy.LE.-1.0) THEN
                  angle2lines = 3.1415926535897931
              ELSE IF (dummy.GE.1.0) THEN
                  angle2lines = 0.0
              ELSE
                  angle2lines = ACOS(DOT_PRODUCT(line1,line2)/ &
                        (NORM2(line1)*NORM2(line2)))
              END IF
          END IF

      END FUNCTION

! ----------------------------------------------------------------------
!
      SUBROUTINE cross_product(a, b, newVector)
!
! ----------------------------------------------------------------------
            USE mCommon
            IMPLICIT NONE
          REAL(rk), DIMENSION(3), INTENT(IN)  :: a, b
          REAL(rk), DIMENSION(3)              ::   newVector
  
          newVector(1) = a(2)*b(3) - a(3)*b(2)
          newVector(2) = a(3)*b(1) - a(1)*b(3)
          newVector(3) = a(1)*b(2) - b(1)*a(2)
      END SUBROUTINE
! ----------------------------------------------------------------------
!
      SUBROUTINE quaternion( angles, qua, roMat )
!
! ----------------------------------------------------------------------
            USE mCommon
          IMPLICIT NONE
          REAL(rk)               :: angles(3)
          REAL(rk)               :: sx, sy, sz, cx, cy, cz, cp
          REAL(rk)               :: qua(4)
          REAL(rk)               :: roMat(3,3)

          sx = COS((angles(1) + angles(3))/2.0)
          sy = COS((angles(1) - angles(3))/2.0)
          cx = SIN((angles(1) + angles(3))/2.0)
          cy = SIN((angles(1) - angles(3))/2.0)
          cz = SIN(angles(2)/2.0)
          sz = COS(angles(2)/2.0)

          qua(1) = sx * sz
          qua(2) = sy * cz
          qua(3) = cy * cz
          qua(4) = cx * sz

          roMat(1,1) = qua(1)**2.0 + qua(2)**2.0-qua(3)**2.0-qua(4)**2.0
          roMat(2,2) = qua(1)**2.0 - qua(2)**2.0+qua(3)**2.0-qua(4)**2.0
          roMat(3,3) = qua(1)**2.0 - qua(2)**2.0-qua(3)**2.0+qua(4)**2.0
          roMat(1,2) = 2.0*(qua(2)*qua(3)+qua(1)*qua(4))
          roMat(1,3) = 2.0*(qua(2)*qua(4)-qua(1)*qua(3))
          roMat(2,1) = 2.0*(qua(2)*qua(3)-qua(1)*qua(4))
          roMat(2,3) = 2.0*(qua(3)*qua(4)+qua(1)*qua(2))
          roMat(3,1) = 2.0*(qua(2)*qua(4)+qua(1)*qua(3))
          roMat(3,2) = 2.0*(qua(3)*qua(4)-qua(1)*qua(2))


      END SUBROUTINE
! ----------------------------------------------------------------------
!
      SUBROUTINE cartesianNormals( normsMainPlain )
!
! ----------------------------------------------------------------------
            USE mCommon
        IMPLICIT NONE
        REAL(rk), DIMENSION(3,3)  ::  normsMainPlain
        INTEGER               ::  i, j

        DO i = 1, 3
          DO j = 1, 3
              normsMainPlain(i,j) = 0.0_rk
          END DO
        END DO

        normsMainPlain(1,3) = 1.0_rk
        normsMainPlain(2,2) = 1.0_rk
        normsMainPlain(3,1) = 1.0_rk
      END SUBROUTINE


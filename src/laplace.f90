subroutine solveLaplace()
      
      USE mMesh
      USE mPar
      USE mEqns
      implicit none
      INTEGER i

!
!     Integrate
!
      CALL CoR_Integrals()
!
!     Set up X and RHS vectors
!
      CALL sdSetUpXandB()
!
!     Set up System Matrix and RHS matrix
!
      IF (parfromDisk.EQ.parYes) THEN
!
!           Set up matrices from disk
!
            CALL WriteToLog("Setting up System and RHS matrices from HARD DISK!")
            DO i=1,neq
                  CALL sdSetUpSysMrhsM_CRS_fromDisk(i)
            END DO
      ELSE
!
!           Set up matrices from memory
!
            CALL WriteToLog("Setting up System and RHS matrices from MEMORY!")
            DO i=1,neq
                  CALL sdSetUpSysMrhsM_CRS(i)
            END DO
!
!           Deallocate integral matrices
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

end subroutine


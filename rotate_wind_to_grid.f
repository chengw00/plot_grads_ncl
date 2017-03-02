C NCLFORTSTART
      SUBROUTINE DCOMPUTEUVGRID(U,V,UVGRID,LONGCA,LONGCB,FLONG,FLAT,
     +                         CEN_LONG,CONE,RPD,NX,NY,NZ,NXP1,NYP1)
      IMPLICIT NONE
      INTEGER NX,NY,NZ,NXP1,NYP1,NL
      REAL U(NXP1,NY,NZ),V(NX,NYP1,NZ)
      REAL UVGRID(NX,NY,NZ,2)
      REAL FLONG(NX,NY),FLAT(NX,NY)
      REAL LONGCB(NX,NY),LONGCA(NX,NY)
      REAL CEN_LONG,CONE,RPD
C NCLEND

      INTEGER I,J,K
      DOUBLE PRECISION UK,VK


c      WRITE (6,FMT=*) ' in compute_uvmet ',NX,NY,NZ,NXP1,NYP1

      DO J = 1,NY
          DO I = 1,NX

              LONGCA(I,J) = FLONG(I,J) - CEN_LONG
              IF (LONGCA(I,J).GT.180.D0) THEN
                  LONGCA(I,J) = LONGCA(I,J) - 360.D0
              END IF
              IF (LONGCA(I,J).LT.-180.D0) THEN
                  LONGCA(I,J) = LONGCA(I,J) + 360.D0
              END IF
              IF (FLAT(I,J).LT.0.D0) THEN
                  LONGCB(I,J) = -LONGCA(I,J)*CONE*RPD
              ELSE
                  LONGCB(I,J) = LONGCA(I,J)*CONE*RPD
              END IF

              LONGCA(I,J) = COS(LONGCB(I,J))
              LONGCB(I,J) = SIN(LONGCB(I,J))

          END DO
      END DO

c      WRITE (6,FMT=*) ' computing velocities '

      DO K = 1,NZ
          DO J = 1,NY
              DO I = 1,NX
                  UK = U(I,J,K)
                  VK = V(I,J,K)
                  ! ====== original: rotate grid winds to Earth-relative winds
                  UVGRID(I,J,K,1) = VK*LONGCB(I,J) + UK*LONGCA(I,J)
                  UVGRID(I,J,K,2) = VK*LONGCA(I,J) - UK*LONGCB(I,J)

                  ! ====== rotate Earth-relative to grid-relative winds
                  !UVGRID(I,J,K,1) = UK*LONGCA(I,J) - VK*LONGCB(I,J)
                  !UVGRID(I,J,K,2) = UK*LONGCB(I,J) + VK*LONGCA(I,J) 
              END DO
          END DO
      END DO

      RETURN
      END



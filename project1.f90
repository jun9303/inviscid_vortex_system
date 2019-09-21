!====================================================================
    MODULE MOD_VARS
!====================================================================
    INTEGER *8, PARAMETER :: MAXP = 100
    REAL*8, PARAMETER :: PI = 4*ATAN(1.D0)
    INTEGER *8 :: NOV, NTST, NSUB        ! SETTINGS
    REAL *8    :: TIME, DT, TO           ! TIME VARIABLES
    REAL *8    :: RKC1, RKC2, RKC3, RKC4 ! RK4 COEFFICIENTS
    REAL *8    :: TM, XM, YM             ! RK4 INTERIM TIME & COORDINATES
    REAL *8    :: RKV1X, RKV2X, RKV3X, RKV4X ! RK4 VALUES FOR X
    REAL *8    :: RKV1Y, RKV2Y, RKV3Y, RKV4Y ! RK4 VALUES FOR Y

    REAL *8, DIMENSION(1:MAXP)    :: POSX, XO
    REAL *8, DIMENSION(1:MAXP)    :: POSY, YO
    REAL *8, DIMENSION(1:MAXP)    :: GAM

    END MODULE MOD_VARS
!====================================================================
    PROGRAM PROJECT_ONE
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE
    INTEGER *8 :: I, M, N

    !=== READ # OF VORTICES, INITIAL POSITION, GAMMA
    !=== SET # OF STEPS AND TIMESTEP
    CALL READ()

    TIME = 0.

    CALL RKINIT()

    !=== TIME-DEPENDENT CALCULATION
    DO M = 0, NTST

        !=== WRITE POSITION & VELOCITY
        IF (MOD(M,50)==0) THEN
            CALL WRITEHISTORY(M)
        ENDIF

        CALL STEPINIT()

        ! === RK4 CALCULATION FOR EACH VORTEX
        DO I = 1, NOV
            DO N = 1, NSUB
                CALL SUBSTEPINIT(I,N)
                CALL RK4_VAL(I,N)
            ENDDO
            !=== UPDATE POSITION
            CALL RK4_UPDATE(I)
        ENDDO
        ! =========================

        TIME = TIME + DT

    ENDDO

        CALL PRINT_CURRENT_TIME()

    STOP
    END PROGRAM PROJECT_ONE
!====================================================================
    SUBROUTINE READ
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE
    CHARACTER*10 :: DUMMY
    INTEGER*8 :: N

    OPEN(10,FILE='input.in')
    READ(10,*) DUMMY
    READ(10,*) NOV
    READ(10,*) DUMMY
    DO N = 1, NOV
        READ(10,*) POSX(N), POSY(N), GAM(N)
    ENDDO
    READ(10,*) DUMMY
    READ(10,*) NTST, DT
    CLOSE(10)

    RETURN
    END SUBROUTINE READ
!====================================================================
    SUBROUTINE RKINIT
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE

    ! USING CLASSIC RUNGE-KUTTA 4TH ORDER (RK4) METHOD
    RKC1 = 1./6.
    RKC2 = 1./3.
    RKC3 = 1./3.
    RKC4 = 1./6.

    NSUB = 4

    RETURN
    END SUBROUTINE RKINIT
!====================================================================
   SUBROUTINE WRITEHISTORY(ST)
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE
    INTEGER *8 :: ST, I

    WRITE(*,*) 'STEP: ', ST, '  TIME: ', TIME

    OPEN(2000,FILE='history.dat',POSITION='APPEND')
    WRITE(2000,100,ADVANCE='NO') TIME
    DO I = 1, NOV
        IF (I .NE. NOV) THEN
            WRITE(2000,101,ADVANCE='NO') POSX(I), POSY(I)
        ELSE
            WRITE(2000,102) POSX(I), POSY(I)
        ENDIF
    ENDDO

100 FORMAT(F12.4, ' ')
101 FORMAT(F12.4, F12.4, ' ')
102 FORMAT(F12.4, F12.4)

    CLOSE(2000)

    RETURN
    END SUBROUTINE WRITEHISTORY
!====================================================================
    SUBROUTINE STEPINIT
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE

    TO = TIME
    XO = POSX
    YO = POSY

    RETURN
    END SUBROUTINE STEPINIT
!====================================================================
    SUBROUTINE SUBSTEPINIT(II,SS)
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE
    INTEGER *8 :: II ! VORTEX INDEX
    INTEGER *8 :: SS ! RK4 SUBSTEP FROM 1 TO 4

    IF (SS .EQ. 1) THEN
        TM = TO
        XM = XO(II)
        YM = YO(II)
    ELSE IF (SS .EQ. 2) THEN
        TM = TO + .5*DT
        XM = XO(II) + .5*DT*RKV1X
        YM = YO(II) + .5*DT*RKV1Y
    ELSE IF (SS .EQ. 3) THEN
        TM = TO + .5*DT
        XM = XO(II) + .5*DT*RKV2X
        YM = YO(II) + .5*DT*RKV2Y
    ELSE IF (SS .EQ. 4) THEN
        TM = TO + DT
        XM = XO(II) + DT*RKV3X
        YM = YO(II) + DT*RKV3Y
    ELSE
        WRITE(*,*) 'ERROR: INVALID SUBSTEP (RK4)'
        STOP
    ENDIF

    RETURN
    END SUBROUTINE SUBSTEPINIT
!====================================================================
    SUBROUTINE RK4_VAL(II, SS)
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE
    INTEGER *8 :: II ! VORTEX INDEX
    INTEGER *8 :: SS ! RK4 SUBSTEP FROM 1 TO 4
    INTEGER *8 :: I
    REAL *8    :: RKVX, RKVY

    RKVX = 0.
    RKVY = 0.
    DO I = 1, NOV
        IF (I .NE. II) THEN
            RKVX = RKVX + GAM(I)/(((XM-XO(I))**2.+(YM-YO(I))**2.)) &
                                *(-(YM-YO(I)))
            RKVY = RKVY + GAM(I)/(((XM-XO(I))**2.+(YM-YO(I))**2.)) &
                                *(XM-XO(I))
        ENDIF
    ENDDO

    IF (SS .EQ. 1) THEN
        RKV1X = RKVX
        RKV1Y = RKVY
    ELSE IF (SS .EQ. 2) THEN
        RKV2X = RKVX
        RKV2Y = RKVY
    ELSE IF (SS .EQ. 3) THEN
        RKV3X = RKVX
        RKV3Y = RKVY
    ELSE IF (SS .EQ. 4) THEN
        RKV4X = RKVX
        RKV4Y = RKVY
    ELSE
        WRITE(*,*) 'ERROR: INVALID SUBSTEP (RK4)'
        STOP
    ENDIF

    RETURN
    END SUBROUTINE RK4_VAL
!====================================================================
    SUBROUTINE RK4_UPDATE(II)
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE
    INTEGER *8 :: II ! VORTEX INDEX

    POSX(II) = POSX(II) + DT*(RKC1*RKV1X+RKC2*RKV2X+RKC3*RKV3X+RKC4*RKV4X)
    POSY(II) = POSY(II) + DT*(RKC1*RKV1Y+RKC2*RKV2Y+RKC3*RKV3Y+RKC4*RKV4Y)

    RETURN
    END SUBROUTINE RK4_UPDATE
!====================================================================
    SUBROUTINE RK4_STABILITY
!====================================================================
    USE MOD_VARS
    IMPLICIT NONE

    RETURN
    END SUBROUTINE RK4_STABILITY
!=======================================================================
      SUBROUTINE PRINT_CURRENT_TIME
!=======================================================================
    IMPLICIT NONE
    CHARACTER*8 ::  DATE
    CHARACTER*10::  NOW
    CHARACTER*5 ::  ZONE
    INTEGER*8   ::  VALS(8)

    CALL DATE_AND_TIME(DATE,NOW,ZONE,VALS)

    WRITE(*,101) VALS(1),VALS(2),VALS(3),VALS(5),VALS(6),VALS(7)
    WRITE(*,*) ''
101 FORMAT(' @ 'I0.4,'-',I0.2,'-',I0.2,' ',I0.2,':',I0.2,':',I0.2)

    RETURN
    END SUBROUTINE PRINT_CURRENT_TIME
!=======================================================================






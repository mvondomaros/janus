PROGRAM ave_rho_interface
IMPLICIT NONE

!! *** Calculates the average oxygen within both interfacial regions ***

!!! MODIFIED TO READ FROM STDIN !

!!  *** Parameter declaration  ***

INTEGER, PARAMETER                      :: SP = SELECTED_REAL_KIND(p=6)
INTEGER, PARAMETER                      :: DP = SELECTED_REAL_KIND(p=12)

INTEGER, PARAMETER                      :: NT = 30000
INTEGER, PARAMETER                      :: NW = 839
INTEGER, PARAMETER                      :: NC = 767
INTEGER, PARAMETER                      :: MAXATOMS = 5000
REAL(KIND=DP), PARAMETER                :: ZMAX = 24.0_DP
REAL(KIND=DP), PARAMETER                :: ZMIN = 11.0_DP
REAL(KIND=DP), PARAMETER                :: XBOX = 31.92800_DP
REAL(KIND=DP), PARAMETER                :: YBOX = 31.19539_DP
REAL(KIND=DP), PARAMETER                :: ZBOX = 31.94800_DP
REAL(KIND=DP), PARAMETER                :: CENTER = 16.0_DP

!!  *** Variable Declaration  ***

INTEGER                                 :: i,j,k,l,m
INTEGER                                 :: o0,h1,h2
INTEGER                                 :: interfaceU, interfaceL
REAL(KIND=DP), DIMENSION(MAXATOMS)      :: x,y,z
CHARACTER(LEN=120)                      :: input, output

!!  *** Get Arguments  ***

CALL getarg(1, output)

!!  *** Execution Section  ***

!OPEN(UNIT=100, STATUS='OLD', ACTION='READ', FILE=input)
OPEN(UNIT=5)
OPEN(UNIT=200, STATUS='NEW', ACTION='WRITE', FILE=output)

DO i=1,NT
  DO j=1,9
    READ(5,*)
  END DO
  DO j=1,NW*3
    READ(5,*) k,l,x(k-NC),y(k-NC),z(k-NC)
    m=k-NC
    x(m)=x(m)*XBOX
    y(m)=y(m)*YBOX
    z(m)=z(m)*ZBOX
  END DO
  interfaceU = 0
  interfaceL = 0
  DO j=1,NW
  o0=3*j-2
  h1=3*j-1
  h2=3*j
  IF ( z(o0).GE.ZMIN .AND. z(o0).LE.ZMAX ) THEN
    IF ( z(o0).LT.CENTER ) THEN
      interfaceL = interfaceL + 1
    ELSE
      interfaceU = interfaceU + 1
    END IF
  END IF
  END DO
  WRITE(200,99) i, interfaceL, interfaceU
END DO

99 FORMAT (I6, 2X, 2I12)
CLOSE(5)
CLOSE(200)

END PROGRAM ave_rho_interface

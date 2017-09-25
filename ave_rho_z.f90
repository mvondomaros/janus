PROGRAM ave_rho_z
USE ISO_FORTRAN_ENV
IMPLICIT NONE

!! *** Calculates the average oxygen density within non-overlapping slabs in z-direction ***

!!! MODIFIED TO READ FROM STDIN !

!!  *** Parameter declaration  ***

INTEGER, PARAMETER                      :: SP = SELECTED_REAL_KIND(p=6)
INTEGER, PARAMETER                      :: DP = SELECTED_REAL_KIND(p=12)

INTEGER, PARAMETER                      :: NT = 30000
INTEGER, PARAMETER                      :: NW = 839
INTEGER, PARAMETER                      :: NC = 767
INTEGER, PARAMETER                      :: MAXATOMS = 5000
REAL(KIND=DP), PARAMETER                :: XBOX = 31.92800_DP
REAL(KIND=DP), PARAMETER                :: YBOX = 31.19539_DP
REAL(KIND=DP), PARAMETER                :: ZBOX = 31.94800_DP
INTEGER, PARAMETER                      :: NBINS = 100

!!  *** Variable Declaration  ***

INTEGER                                 :: i,j,k,l,m
INTEGER                                 :: ioxygen
INTEGER                                 :: o0,h1,h2
INTEGER, DIMENSION(0:NBINS-1)           :: noxygens
REAL(KIND=DP), DIMENSION(MAXATOMS)      :: x,y,z
CHARACTER(LEN=120)                      :: input, output

!!  *** Get Arguments  ***

CALL getarg(1, output)

!!  *** Execution Section  ***

OPEN(UNIT=INPUT_UNIT)
OPEN(UNIT=200, STATUS='NEW', ACTION='WRITE', FILE=output)

DO i=1,NT
  DO j=1,9
    READ(INPUT_UNIT,*)
  END DO
  DO j=1,NW*3
    READ(INPUT_UNIT,*) k,l,x(k-NC),y(k-NC),z(k-NC)
    m=k-NC
    x(m)=x(m)*XBOX
    y(m)=y(m)*YBOX
    z(m)=z(m)*ZBOX
  END DO
  noxygens = 0
  DO j=1,NW
    o0=3*j-2
    h1=3*j-1
    h2=3*j
    ioxygen = FLOOR(z(o0)*REAL(NBINS, DP)/ZBOX)
    IF (ioxygen < 0) THEN
        ioxygen = ioxygen + NBINS
    ELSE IF (ioxygen >= NBINS) THEN
        ioxygen = ioxygen - NBINS
    END IF
    noxygens(ioxygen) = noxygens(ioxygen) + 1
  END DO
  WRITE(200,99) noxygens
END DO

99 FORMAT (100(I12, 2X))
CLOSE(100)
CLOSE(200)

END PROGRAM ave_rho_z

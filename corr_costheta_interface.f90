PROGRAM corr_costheta_interface
USE shared
USE subroutines
IMPLICIT NONE

!! *** Calculates the time correlation function of cos(theta) within the interface ***

!!  *** Parameter declaration  ***

INTEGER, PARAMETER                      :: NT = 1660000
INTEGER, PARAMETER                      :: NCORR = 50000
INTEGER, PARAMETER                      :: EQUILIBRIUM = 15000

!!  *** Variable Declaration  ***

INTEGER                                 :: i,j
REAL(KIND=DP),DIMENSION(NT)             :: dipZU
REAL(KIND=DP),DIMENSION(NT)             :: dipZL
REAL(KIND=DP)                           :: dummy, averageL, averageU
REAL(KIND=DP),DIMENSION(0:NCORR)        :: corrU, corrL
CHARACTER(LEN=120)                      :: input, output

!!  *** Get Arguments  ***

CALL getarg(1, input)
CALL getarg(2, output)

!!  *** Execution Section  ***

OPEN(UNIT=100, STATUS='OLD', ACTION='READ', FILE=input)
OPEN(UNIT=200, STATUS='NEW', ACTION='WRITE', FILE=output)

DO i=1,EQUILIBRIUM
  READ(100,*)
END DO

DO i=1,NT
  READ(100,*) j, dummy, dummy, dipZL(i), dummy, dummy, dipZU(i) 
END DO

averageL = SUM(dipZL)/REAL(NT,dp)
averageU = SUM(dipZU)/REAL(NT,dp)
dipZL = dipZL - averageL
dipZU = dipZU - averageU

CALL autocorrelationF(NT,NCORR,dipZL,corrL)
CALL autocorrelationF(NT,NCORR,dipZU,corrU)

DO i=0,NCORR
  WRITE(200,*) REAL(i)/100., corrL(i), corrU(i)
END DO

CLOSE(100)
CLOSE(200)

END PROGRAM corr_costheta_interface

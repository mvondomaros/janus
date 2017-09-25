PROGRAM corr_rho_interface
USE shared
USE subroutines
IMPLICIT NONE

!! *** Calculates the time correlation function of the interfacial density ***
 
!!  *** Parameter declaration  ***

INTEGER, PARAMETER                      :: NT = 1000000
INTEGER, PARAMETER                      :: NCORR = 5000
INTEGER, PARAMETER                      :: EQUILIBRIUM = 0

!!  *** Variable Declaration  ***

INTEGER                                 :: i,j
REAL(KIND=DP),DIMENSION(NT)             :: dL, dU
REAL(KIND=DP)                           :: averageL, averageU
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
  READ(100,*) j, dL(i), dU(i)
END DO

averageL = SUM(dL)/NT
averageU = SUM(dU)/NT
dL = dL - averageL
dU = dU - averageU

CALL autocorrelationF(NT,NCORR,dL,corrL)
CALL autocorrelationF(NT,NCORR,dU,corrU)

DO i=0,NCORR
  WRITE(200,*) REAL(i)/100., corrL(i), corrU(i)
END DO

CLOSE(100)
CLOSE(200)

END PROGRAM corr_rho_interface

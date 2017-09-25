PROGRAM hist_costheta_interface
USE shared
IMPLICIT NONE

!! *** Calculates histograms of cos(theta) within the interfacial regions ***

!!  *** Parameter declaration  ***

INTEGER, PARAMETER                      :: NT = 1500000
INTEGER, PARAMETER                      :: EQUILIBRIUM = 30000
INTEGER, PARAMETER                      :: FINE = 50

!!  *** Variable Declaration  ***

INTEGER                                 :: i,j
INTEGER                                 :: temp 
REAL(KIND=DP),DIMENSION(-FINE:FINE)     :: pL,pU
REAL(KIND=SP),DIMENSION(NT)             :: dipXU, dipYU, dipZU
REAL(KIND=SP),DIMENSION(NT)             :: dipXL, dipYL, dipZL
CHARACTER(LEN=120)                      :: input, output

!!  *** Get Arguments  ***

CALL getarg(1, input)
CALL getarg(2, output)

!!  *** Execution Section  ***

OPEN(UNIT=100, STATUS='OLD', ACTION='READ', FILE=input)
OPEN(UNIT=200, STATUS='UNKNOWN', ACTION='WRITE', FILE=output)

DO i=1,EQUILIBRIUM
  READ(100,*)
END DO

pL=0.0_DP
pU=0.0_DP

DO i=1,NT
  READ(100,*) j, dipXL(i), dipYL(i), dipZL(i), dipXU(i), dipYU(i), dipZU(i) 
  temp=NINT(dipZL(i)*REAL(FINE))
  pL(temp)=pL(temp)+1.0_DP
  temp=NINT(dipZU(i)*REAL(FINE))
  pU(temp)=pU(temp)+1.0_DP
END DO

pL(-FINE)=REAL(2)*pL(-FINE)
pU(-FINE)=REAL(2)*pU(-FINE)
pL(FINE)=REAL(2)*pL(FINE)
pU(FINE)=REAL(2)*pU(FINE)

DO i=-FINE,FINE
  pL(i)=pL(i)/REAL(NT)
  pU(i)=pU(i)/REAL(NT)
  WRITE(200,*) REAL(i)/REAL(FINE), pL(i)*REAL(FINE), pU(i)*REAL(FINE)
END DO

CLOSE(100)
CLOSE(200)

END PROGRAM hist_costheta_interface

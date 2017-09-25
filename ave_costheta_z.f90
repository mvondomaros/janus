PROGRAM ave_costheta_z
IMPLICIT NONE

!! *** Calculates the average of cos(theta) within non-overlapping slabs in z-direction ***

!!  *** Parameter declaration  ***

INTEGER, PARAMETER                      :: SP = SELECTED_REAL_KIND(p=6)
INTEGER, PARAMETER                      :: DP = SELECTED_REAL_KIND(p=12)

INTEGER, PARAMETER                      :: NT = 30000
INTEGER, PARAMETER                      :: NW = 839
INTEGER, PARAMETER                      :: NC = 767
INTEGER, PARAMETER                      :: MAXATOMS = 5000
INTEGER, PARAMETER                      :: NRSLABS = 100
REAL(KIND=DP), PARAMETER                :: XBOX = 31.92800_DP
REAL(KIND=DP), PARAMETER                :: YBOX = 31.19539_DP
REAL(KIND=DP), PARAMETER                :: ZBOX = 31.94800_DP

!!  *** Variable Declaration  ***

INTEGER                                 :: i,j,k,l,m
INTEGER                                 :: o0,h1,h2
REAL(KIND=DP), DIMENSION(MAXATOMS)      :: x,y,z 
REAL(KIND=DP),DIMENSION(0:NRSLABS)      :: dX, dY, dZ, dipScale
REAL(KIND=DP),DIMENSION(0:NRSLABS)      :: aveDX, aveDY, aveDZ
CHARACTER(LEN=120)                      :: input, output

!!  *** Get Arguments  ***

CALL getarg(1, input)
CALL getarg(2, output)

!!  *** Execution Section  ***

OPEN(UNIT=100, STATUS='OLD', ACTION='READ', FILE=input)
OPEN(UNIT=200, STATUS='NEW', ACTION='WRITE', FILE=output)

dX = 0.0_DP
dY = 0.0_DP
dZ = 0.0_DP
aveDX = 0.0_DP
aveDY = 0.0_DP
aveDZ = 0.0_DP

DO i=1,NT
  DO j=1,9
    READ(100,*)
  END DO
  DO j=1,NW*3
    READ(100,*) k,l,x(k-NC),y(k-NC),z(k-NC) 
    m=k-NC
    x(m)=x(m)*XBOX
    y(m)=y(m)*YBOX
    z(m)=z(m)*ZBOX 
  END DO
  DO j=1,NW
    o0=3*j-2
    h1=3*j-1
    h2=3*j
    k=NINT(z(o0)/ZBOX*REAL(NRSLABS))
    dX(k)=x(h1)-x(o0)-NINT((x(h1)-x(o0))/XBOX)*XBOX &
         +x(h2)-x(o0)-NINT((x(h2)-x(o0))/XBOX)*XBOX + dX(k)
    dY(k)=y(h1)-y(o0)-NINT((y(h1)-y(o0))/YBOX)*YBOX &
         +y(h2)-y(o0)-NINT((y(h2)-y(o0))/YBOX)*YBOX + dY(k)
    dZ(k)=z(h1)-z(o0)-NINT((z(h1)-z(o0))/ZBOX)*ZBOX &
         +z(h2)-z(o0)-NINT((z(h2)-z(o0))/ZBOX)*ZBOX + dZ(k)
  END DO
  dX(0)=dX(0)*2.0_DP
  dY(0)=dY(0)*2.0_DP
  dZ(0)=dZ(0)*2.0_DP
  dX(NRSLABS)=dX(NRSLABS)*2.0_DP
  dY(NRSLABS)=dY(NRSLABS)*2.0_DP
  dZ(NRSLABS)=dZ(NRSLABS)*2.0_DP
  DO j=0,NRSLABS
    aveDX(j)=aveDX(j)+dX(j)
    aveDY(j)=aveDY(j)+dY(j)
    aveDZ(j)=aveDZ(j)+dZ(j)
  END DO   
END DO

DO j=0,NRSLABS
  dipScale(j)=SQRT(aveDX(j)**2+aveDY(j)**2+aveDZ(j)**2)
  IF (dipScale(j) .GT. 0.00001) THEN
  aveDX(j)=aveDX(j)/dipScale(j) 
  aveDY(j)=aveDY(j)/dipScale(j) 
  aveDZ(j)=aveDZ(j)/dipScale(j) 
  END IF
END DO   

DO i=0,NRSLABS
  WRITE(200,*) REAL(i)/REAL(NRSLABS)*ZBOX, aveDX(i), aveDY(i), aveDZ(i)
END DO

CLOSE(100)
CLOSE(200)

END PROGRAM ave_costheta_z

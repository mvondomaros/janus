PROGRAM ave_dipole_interface
IMPLICIT NONE

!! *** Calculates the average of cos(theta) within both interfacial regions ***

!!  *** Parameter declaration  ***

INTEGER, PARAMETER                      :: SP = SELECTED_REAL_KIND(p=6)
INTEGER, PARAMETER                      :: DP = SELECTED_REAL_KIND(p=12)

INTEGER, PARAMETER                      :: NT = 5000
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
INTEGER                                 :: sumNU, sumNL, interfaceU, interfaceL
REAL(KIND=DP), DIMENSION(MAXATOMS)      :: x,y,z 
REAL(KIND=DP)                           :: dipXU, dipYU, dipZU
REAL(KIND=DP)                           :: dipXL, dipYL, dipZL
REAL(KIND=DP)                           :: dX, dY, dZ, dipScale
CHARACTER(LEN=120)                      :: input, output

!!  *** Get Arguments  ***

CALL getarg(1, input)
CALL getarg(2, output)

!!  *** Execution Section  ***

OPEN(UNIT=100, STATUS='OLD', ACTION='READ', FILE=input)
OPEN(UNIT=200, STATUS='NEW', ACTION='WRITE', FILE=output)

sumNU = 0
sumNL = 0
dipXU = 0.0_DP
dipYU = 0.0_DP
dipZU = 0.0_DP
dipXL = 0.0_DP
dipYL = 0.0_DP
dipZL = 0.0_DP

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
  interfaceU = 0
  interfaceL = 0
  DO j=1,NW
  o0=3*j-2
  h1=3*j-1
  h2=3*j
  IF ( z(o0).GE.ZMIN .AND. z(o0).LE.ZMAX ) THEN
    dX=x(h1)-x(o0)-NINT((x(h1)-x(o0))/XBOX)*XBOX &
      +x(h2)-x(o0)-NINT((x(h2)-x(o0))/XBOX)*XBOX
    dY=y(h1)-y(o0)-NINT((y(h1)-y(o0))/YBOX)*YBOX &
      +y(h2)-y(o0)-NINT((y(h2)-y(o0))/YBOX)*YBOX
    dZ=z(h1)-z(o0)-NINT((z(h1)-z(o0))/ZBOX)*ZBOX &
      +z(h2)-z(o0)-NINT((z(h2)-z(o0))/ZBOX)*ZBOX
    IF ( z(o0).LT.CENTER ) THEN
      interfaceL = interfaceL + 1
      dipXL = dipXL + dX
      dipYL = dipYL + dY
      dipZL = dipZL + dZ
    ELSE
      interfaceU = interfaceU + 1
      dipXU = dipXU + dX
      dipYU = dipYU + dY
      dipZU = dipZU + dZ
    END IF
  END IF
  END DO
  dipXL = dipXL/REAL(interfaceL)
  dipYL = dipYL/REAL(interfaceL)
  dipZL = dipZL/REAL(interfaceL)
  dipXU = dipXU/REAL(interfaceU)
  dipYU = dipYU/REAL(interfaceU)
  dipZU = dipZU/REAL(interfaceU)
  dipScale = SQRT(dipXL**2 + dipYL**2 + dipZL**2)
  dipXL = dipXL/dipScale
  dipYL = dipYL/dipScale
  dipZL = dipZL/dipScale
  dipScale = SQRT(dipXU**2 + dipYU**2 + dipZU**2)
  dipXU = dipXU/dipScale
  dipYU = dipYU/dipScale
  dipZU = dipZU/dipScale
  sumNL = sumNL + interfaceL
  sumNU = sumNU + interfaceU
  WRITE(200,99) i, dipXL, dipYL, dipZL, dipXU, dipYU, dipZU
END DO

WRITE(*,*) '# of interfacial water', sumNL/NT, sumNU/NT

99 FORMAT (I6, 2X, 6F8.4)
CLOSE(100)
CLOSE(200)

END PROGRAM ave_dipole_interface

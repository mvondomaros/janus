PROGRAM hist_costheta_z
USE shared
IMPLICIT NONE

!! *** Calculates histograms of cos(theta) within non-overlapping slabs in z direction ***

    !!! Paramaters

    INTEGER, PARAMETER                      :: NT = 30000
    INTEGER, PARAMETER                      :: NRZBINS = 20
    INTEGER, PARAMETER                      :: NRCBINS = 10
    INTEGER, PARAMETER                      :: NW = 839
    INTEGER, PARAMETER                      :: NC = 767
    INTEGER, PARAMETER                      :: MAXATOMS = 5000
    REAL(KIND=DP), PARAMETER                :: XBOX = 31.92800_DP
    REAL(KIND=DP), PARAMETER                :: YBOX = 31.19539_DP
    REAL(KIND=DP), PARAMETER                :: ZBOX = 31.94800_DP

    !!! Variables
    CHARACTER(LEN=120)                          :: output, input
    INTEGER                                     :: i,j,k,l,m,f
    INTEGER                                     :: dipZ
    REAL(KIND=DP), DIMENSION(MAXATOMS)          :: x,y,z 
    REAL(KIND=DP), DIMENSION(NRZBINS)           :: overallX, overallY, overallZ
    REAL(KIND=DP), DIMENSION(NRZBINS)           :: cosine
    REAL(KIND=DP), DIMENSION(NRZBINS,-NRCBINS:NRCBINS)   :: propability
    INTEGER                         :: o0, h1, h2
    REAL(KIND=DP)                   :: dX, dY, dZ
    REAL(KIND=DP)                   :: norm
    INTEGER, DIMENSION(NW)          :: nroxygens
    
    !!! Get Arguments
    CALL getarg(2, output)
    OPEN(UNIT=200, STATUS='NEW', ACTION='WRITE', FILE=output)
    propability = 0.0_DP
    CALL getarg(1, input)
    !!! Execution section
    OPEN(UNIT=100, STATUS='OLD', ACTION='READ', FILE=input)
                                                                                                  
    DO i = 1, NT
        DO j = 1, 9
            READ(100,*)
        END DO
        DO j = 1, NW * 3
            READ(100,*) k,l,x(k-NC),y(k-NC),z(k-NC) 
            m = k - NC
            x(m) = x(m) * XBOX
            y(m) = y(m) * YBOX
            z(m) = z(m) * ZBOX 
        END DO

        cosine = 0.0_DP
        overallX = 0.0_DP
        overallY = 0.0_DP
        overallZ = 0.0_DP
        nroxygens = 0
        DO j=1,NW

            o0=3*j-2
            h1=3*j-1
            h2=3*j

            dX = x(h1) - x(o0) - NINT((x(h1) - x(o0)) / XBOX) * XBOX &
               + x(h2) - x(o0) - NINT((x(h2) - x(o0)) / XBOX) * XBOX
            dY = y(h1) - y(o0) - NINT((y(h1) - y(o0)) / YBOX) * YBOX &
               + y(h2) - y(o0) - NINT((y(h2) - y(o0)) / YBOX) * YBOX
            dZ = z(h1) - z(o0) - NINT((z(h1) - z(o0)) / ZBOX) * ZBOX &
               + z(h2) - z(o0) - NINT((z(h2) - z(o0)) / ZBOX) * ZBOX

            k = CEILING(z(o0) / ZBOX * REAL(NRZBINS))
            overallX(k) = overallX(k) + dX
            overallY(k) = overallY(k) + dY
            overallZ(k) = overallZ(k) + dZ
            nroxygens(k) = nroxygens(k) + 1
        END DO
        DO j = 1, NRZBINS
            IF (nroxygens(j) .GT. 0) THEN 
                norm = SQRT(overallX(j)**2 + overallY(j)**2 + overallZ(j)**2)
                cosine(j) = overallZ(j) /  norm
!                WRITE(*,*) j, cosine(j)
                dipZ = CEILING(cosine(j) * REAL(NRCBINS))
                propability(j, dipZ) = propability(j, dipZ) + 1.0_DP
            END IF
        END DO
    END DO
   !propability = propability / REAL(NT)
    DO i = 1, NRZBINS
        DO j = -NRCBINS + 1, NRCBINS
            WRITE (200,*) REAL(i)/REAL(NRZBINS)*ZBOX, REAL(j)/REAL(NRCBINS), propability(i,j)
        END DO
    END DO
    CLOSE(200)

END PROGRAM hist_costheta_z

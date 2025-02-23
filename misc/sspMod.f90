MODULE sspmod

  ! holds SSP input by user and associated variables

  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER     :: MaxSSP = 20001, MaxMedia = 501
  INTEGER, PRIVATE       :: N, iz, ILoc, Lay
  INTEGER                :: ISSP
  REAL (KIND=8)          :: alphaR = 1500, betaR = 0, alphaI = 0, betaI = 0, rhoR = 1
  REAL (KIND=8), PRIVATE :: h, z, R

  ! SSP
  TYPE SSPStructure
     INTEGER           :: Loc( MaxMedia ), NPts( MaxMedia ), NMedia
     REAL     (KIND=8) :: z( MaxSSP ), alphaR( MaxSSP ), alphaI( MaxSSP ), rho( MaxSSP ), betaR( MaxSSP ), betaI( MaxSSP )
     REAL     (KIND=8) :: Depth( MaxMedia ), sigma( MaxMedia ), beta( MaxMedia ), fT( MaxMedia )
     COMPLEX  (KIND=8) :: cp( MaxSSP ), cs( MaxSSP ), n2( MaxSSP ),  &
                          cpSpline( 4, MaxSSP ), csSpline( 4, MaxSSP ), rhoSpline( 4, MaxSSP )
     CHARACTER (LEN=1) :: Type
     CHARACTER (LEN=2) :: AttenUnit
     CHARACTER (LEN=8) :: Material( MaxMedia )
  END TYPE SSPStructure

  TYPE( SSPStructure ) :: SSP

  TYPE HSInfo
     CHARACTER (LEN=1) :: BC                            ! Boundary condition type
     REAL    (KIND=8)  :: alphaR, alphaI, betaR, betaI  ! P-wave, S-wave speeds (user units)
     REAL    (KIND=8)  :: beta, fT                      ! power law and transition frequency
     COMPLEX (KIND=8)  :: cP, cS                        ! P-wave, S-wave speeds (neper/m loss)
     REAL    (KIND=8)  :: rho, BumpDensity, eta, xi     ! density, boss parameters
  END TYPE

  TYPE( HSInfo )      :: HSTop, HSBot

CONTAINS

  SUBROUTINE EvaluateSSP( cP, cS, rho, Medium, N1, Freq, Task, ENVFile, PRTFile )

    ! Call the particular SSP routine specified by SSPType
    ! Performs two Tasks:
    !    Task = 'TAB'  then tabulate cP, cS, rho
    !    Task = 'INIT' then initialize
    ! Note that Freq is only needed if Task = 'INIT'

    INTEGER,           INTENT(IN)    :: ENVFile, PRTFile, Medium
    INTEGER,           INTENT(INOUT) :: N1
    REAL     (KIND=8), INTENT(OUT)   :: rho( * )
    REAL     (KIND=8), INTENT(IN)    :: Freq
    COMPLEX  (KIND=8), INTENT(OUT)   :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)    :: Task
    COMPLEX  (KIND=8)                :: cPT, cST

    SELECT CASE ( SSP%Type )
    CASE ( 'A' )  !  Analytic profile option 
       IF ( Task( 1 : 4 ) == 'INIT' ) THEN
          N1 = 21

          CALL ANALYT( cP, cS, rho, Medium, N1, Freq, Task )
          h = ( SSP%Depth( Medium + 1 ) - SSP%Depth( Medium ) ) / ( N1 - 1 )

          DO iz = 1, N1
             z   = SSP%Depth( Medium ) + ( iz - 1 ) * h
             cPT =  cP( iz )
             cST =  cS( iz )
             WRITE( PRTFile, FMT="( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )" ) &
                  z,  REAL( cPT ),  REAL( cST ), rho( iz ), AIMAG( cPT ), AIMAG( cST )
          END DO
       ELSE
          CALL ANALYT( cP, cS, rho, Medium, N1, Freq, Task )
       ENDIF
    CASE ( 'N' )  !  N2-linear profile option
       CALL n2Linear( cP, cS, rho, Medium, N1, Task, ENVFile, PRTFile )
    CASE ( 'C' )  !  C-linear profile option 
       CALL cLinear(  cP, cS, rho, Medium, N1, Task, ENVFile, PRTFile )
    CASE ( 'P' )  !  C-linear profile option 
       CALL PCHIP(    cP, cS, rho, Medium, N1, Task, ENVFile, PRTFile )
    CASE ( 'S' )  !  Cubic spline profile option 
       CALL CCubic(   cP, cS, rho, Medium, N1, Task, ENVFile, PRTFile )
    CASE DEFAULT  !  Non-existent profile option 
       WRITE( PRTFile, * ) 'Profile option: ', SSP%Type
       CALL ERROUT( PRTFile, 'F', 'EvaluateSSP', 'Unknown profile option' )
    END SELECT

    RETURN
  END SUBROUTINE EvaluateSSP

  !**********************************************************************!

  SUBROUTINE n2Linear( cP, cS, rho, Medium, N1, Task, ENVFile, PRTFile )

    ! Tabulate cP, cS, rho for specified Medium
    ! Uses N2-linear segments for P and S-wave speeds
    ! Uses rho-linear segments for density

    INTEGER,           INTENT(IN)  :: ENVFile, PRTFile
    INTEGER,           INTENT( INOUT ) :: N1
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    COMPLEX  (KIND=8)              :: N2Bot, N2Top

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1, ENVFile, PRTFile )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          z = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) z = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot

          DO WHILE ( z > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          iSSP = ILoc + Lay
          R = ( z - SSP%z( iSSP ) ) / ( SSP%z( iSSP + 1 ) - SSP%z( iSSP ) )

          ! P-wave
          N2Top    = 1.0 / SSP%cp( iSSP     )**2
          N2Bot    = 1.0 / SSP%cp( iSSP + 1 )**2
          cP( iz ) = 1.0 / SQRT( ( 1.0 - R ) * N2Top + R * N2Bot )

          ! S-wave
          IF ( SSP%cs( iSSP ) /= 0.0 ) THEN
             N2Top    = 1.0 / SSP%cs( iSSP     )**2
             N2Bot    = 1.0 / SSP%cs( iSSP + 1 )**2
             cS( iz ) = 1.0 / SQRT( ( 1.0 - R ) * N2Top + R * N2Bot )
          ELSE
             cS( iz ) = 0.0
          ENDIF

          rho( iz ) = ( 1.0 - R ) * SSP%rho( iSSP ) + R * SSP%rho( iSSP + 1 )
       END DO

    ENDIF

    RETURN
  END SUBROUTINE n2Linear
  
  !**********************************************************************!
  
  SUBROUTINE cLinear( cP, cS, rho, Medium, N1, Task, ENVFile, PRTFile  )

    ! Tabulate cP, cS, rho for specified Medium

    ! Uses c-linear segments for P and S-wave speeds
    ! Uses rho-linear segments for density
    INTEGER,           INTENT(IN)  :: ENVFile, PRTFile
    INTEGER,           INTENT( INOUT ) :: N1
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1, ENVFile, PRTFile )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          z = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) z = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot

          DO WHILE ( z > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          iSSP = ILoc + Lay
          R = ( z - SSP%z( iSSP ) ) / ( SSP%z( iSSP + 1 ) - SSP%z( iSSP ) )
          cP(  iz ) = ( 1.0 - R ) * SSP%cp(  iSSP ) + R * SSP%cp(  iSSP + 1 )
          cS(  iz ) = ( 1.0 - R ) * SSP%cs(  iSSP ) + R * SSP%cs(  iSSP + 1 )
          rho( iz ) = ( 1.0 - R ) * SSP%rho( iSSP ) + R * SSP%rho( iSSP + 1 )
       END DO
    ENDIF

    RETURN
  END SUBROUTINE cLinear

  !**********************************************************************!
  
  SUBROUTINE PCHIP( cP, cS, rho, Medium, N1, Task, ENVFile, PRTFile )

    ! Tabulate cP, cS, rho values for the specified Medium

    ! Monotone Piecewise Cubic Hermite Interpolating Polynomial of SSP data
    !  This is the variant of PCHIP described by Fritch-Butner in UCRL-87559
    !  Also see: F. N. Fritsch and R. E. Carlson. "Monotone Piecewise Cubic
    !  Interpolation. Siam Journal on Numerical Analysis, 17:238-246, Apr 1980.

    USE pchipMod
    INTEGER,           INTENT(IN)  :: ENVFile, PRTFile
    INTEGER,           INTENT( INOUT ) :: N1
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    REAL     (KIND=8) :: hSSP, hSSP1, hSSP2
    COMPLEX  (KIND=8) :: del1, del2, f1, f2, f1prime, f2prime

    ! If Task = 'INIT' then this is the first call and the SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization

       CALL ReadSSP( Medium, N1, ENVFile, PRTFile )

    ELSE   ! Task = 'TABULATE'

       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          z = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) z = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot

          DO WHILE ( z > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          iSSP = ILoc + Lay

          hSSP = SSP%z( iSSP + 1 ) - SSP%z( iSSP )

          ! *** cP ***

          f1 = SSP%cp( iSSP     )
          f2 = SSP%cp( iSSP + 1 )

          ! estimate the derivative at the left endpoint of this SSP segment

          IF ( iSSP == 1 ) THEN                        ! left (top)    endpoint (non-centered 3-point formula, see UCRL-85104)
             CALL h_del( SSP%z, SSP%cp, iSSP + 1, hSSP1, hSSP2, del1, del2 )
             f1prime = ( ( 2.0D0 * hSSP1 + hSSP2 ) * del1 - hSSP1 * del2 ) / ( hSSP1 + hSSP2 );
             f1prime = fprime_left_end_Cmplx( del1, del2, f1prime )
          ELSE   ! interior node (see UCRL-87559)
             CALL h_del( SSP%z, SSP%cp, iSSP, hSSP1, hSSP2, del1, del2 )
             f1prime = fprime_Cmplx( hSSP1, hSSP2, del1, del2 )
          END IF

          ! estimate the derivative at the right endpoint of this SSP segment

          IF ( iSSP == SSP%NPts( Medium ) - 1 ) THEN   ! right (bottom) endpoint (non-centered 3-point formula, see UCRL-85104)
             CALL h_del( SSP%z, SSP%cp, iSSP, hSSP1, hSSP2, del1, del2 )
             f2prime = ( -hSSP2 * del1 + ( hSSP1 + 2.0D0 * hSSP2 ) * del2 ) / ( hSSP1 + hSSP2 );
             f2prime = fprime_right_end_Cmplx( del1, del2, f2prime )
          ELSE   ! interior node (see UCRL-87559)
             CALL h_del( SSP%z, SSP%cp, iSSP + 1, hSSP1, hSSP2, del1, del2 )
             f2prime = fprime_Cmplx( hSSP1, hSSP2, del1, del2 )
          END IF

          cP( iz ) = Cubic( z - SSP%z( iSSP ), f1, f1prime, f2, f2prime, hSSP )   ! here's the Hermite cubic

          ! *** cS ***

          f1 = SSP%cs( iSSP     )
          f2 = SSP%cs( iSSP + 1 )

          ! estimate the derivative at the left endpoint of this SSP segment

          IF ( iSSP == 1 ) THEN   ! left (top) endpoint (non-centered 3-point formula, see UCRL-85104)
             CALL h_del( SSP%z, SSP%cs, iSSP + 1, hSSP1, hSSP2, del1, del2 )
             f1prime = ( ( 2.0D0 * hSSP1 + hSSP2 ) * del1 - hSSP1 * del2 ) / ( hSSP1 + hSSP2 );
             f1prime = fprime_left_end_Cmplx( del1, del2, f1prime )
          ELSE   ! interior node (see UCRL-87559)
             CALL h_del( SSP%z, SSP%cs, iSSP, hSSP1, hSSP2, del1, del2 )
             f1prime = fprime_Cmplx( hSSP1, hSSP2, del1, del2 )
          END IF

          ! estimate the derivative at the right endpoint of this SSP segment

          IF ( iSSP == SSP%NPts( Medium ) - 1 ) THEN
             ! handle right endpoint (non-centered 3-point formula, see UCRL-85104)
             CALL h_del( SSP%z, SSP%cs, iSSP, hSSP1, hSSP2, del1, del2 )
             f2prime = ( -hSSP2 * del1 + ( hSSP1 + 2.0D0 * hSSP2 ) * del2 ) / ( hSSP1 + hSSP2 );
             f2prime = fprime_right_end_Cmplx( del1, del2, f2prime )
          ELSE   ! interior node (see UCRL-87559)
             CALL h_del( SSP%z, SSP%cs, iSSP + 1, hSSP1, hSSP2, del1, del2 )
             f2prime = fprime_Cmplx( hSSP1, hSSP2, del1, del2 )
          END IF

          cS( iz ) = Cubic( z - SSP%z( iSSP ), f1, f1prime, f2, f2prime, hSSP )   ! here's the Hermite cubic

          ! density
          R = ( z - SSP%z( iSSP ) ) / hSSP
          rho( iz ) = ( 1.0 - R ) * SSP%rho( iSSP ) + R * SSP%rho( iSSP + 1 )   ! linear interpolation of density

       END DO
    ENDIF

    RETURN
  END SUBROUTINE PCHIP

  !**********************************************************************!
  
  SUBROUTINE CCubic( cP, cS, rho, Medium, N1, Task, ENVFile, PRTFile  )

    ! Tabulate cP, cS, rho for specified Medium
    ! using cubic spline interpolation

    INTEGER,           INTENT(IN)  :: ENVFile, PRTFile
    INTEGER,           INTENT( INOUT ) :: N1
    INTEGER,           INTENT(IN)  :: Medium
    REAL     (KIND=8), INTENT(OUT) :: rho( * )
    COMPLEX  (KIND=8), INTENT(OUT) :: cP( * ), cS( * )
    CHARACTER (LEN=8), INTENT(IN)  :: Task
    REAL     (KIND=8)              :: HSPLNE
    COMPLEX  (KIND=8)              :: SPLINE

    ! If Task = 'INIT' then this is the first call and SSP is read.
    ! Any other call is a request for SSP subtabulation.

    IF ( Task( 1 : 4 ) == 'INIT' ) THEN   ! Task 'INIT' for initialization
       CALL ReadSSP( Medium, N1, ENVFile, PRTFile )
    ELSE   ! Task = 'TABULATE'
       ILoc = SSP%Loc( Medium )
       N    = N1 - 1
       h    = ( SSP%z( ILoc + SSP%NPts( Medium ) ) - SSP%z( ILoc + 1 ) ) / N
       Lay  = 1

       DO iz = 1, N1
          z = SSP%z( ILoc + 1 ) + ( iz - 1 ) * h
          IF ( iz == N1 ) z = SSP%z( ILoc + SSP%NPts( Medium ) )   ! Make sure no overshoot
          DO WHILE ( z > SSP%z( ILoc + Lay + 1 ) )
             Lay = Lay + 1
          END DO

          iSSP = ILoc + Lay
          HSPLNE = z - SSP%z( iSSP )

          cP(  iz ) =       SPLINE( SSP%cpSpline(  1, iSSP ), HSPLNE )
          cS(  iz ) =       SPLINE( SSP%csSpline(  1, iSSP ), HSPLNE )
          rho( iz ) = DBLE( SPLINE( SSP%rhoSpline( 1, iSSP ), HSPLNE ) )

       END DO
    ENDIF

    RETURN
  END SUBROUTINE CCubic

!**********************************************************************!

  SUBROUTINE ReadSSP( Medium, N1, ENVFile, PRTFile )

    ! reads the SSP data from the environmental file for a given medium
    
    INTEGER, INTENT( IN    ) :: ENVFile, PRTFile
    INTEGER, INTENT( IN    ) :: Medium
    INTEGER, INTENT( INOUT ) :: N1
    INTEGER                  :: iSSP

    SSP%NPts( Medium ) = N1

    ! The variable SSP%Loc( Medium ) points to the starting point for the
    ! data in the arrays z, alpha, beta and rho
    IF ( Medium == 1 ) THEN
       SSP%Loc( Medium ) = 0
    ELSE
       SSP%Loc( Medium ) = SSP%Loc( Medium - 1 ) + SSP%NPts( Medium - 1 )
    ENDIF
    ILoc = SSP%Loc( Medium )

    !  Read in data and convert attenuation to Nepers/m 
    N1 = 1
    DO iSSP = 1, MaxSSP
       iz = SSP%Loc( Medium ) + iSSP

       READ(  ENVFile, *    ) SSP%z( iz ), alphaR, betaR, rhoR, alphaI, betaI
       WRITE( PRTFile, FMT="( F10.2,      3X, 2F10.2,       3X, F6.2, 3X, 2F10.4 )" ) &
                              SSP%z( iz ),    alphaR, betaR,    rhoR,     alphaI, betaI
       SSP%alphaR( iz ) = alphaR
       SSP%alphaI( iz ) = alphaI
       SSP%rho(    iz ) = rhoR
       SSP%betaR(  iz ) = betaR
       SSP%betaI(  iz ) = betaI

       ! Did we read the last point?
       IF ( ABS( SSP%z( iz ) - SSP%Depth( Medium + 1 ) ) < 100. * EPSILON( 1.0e0 ) ) THEN
          SSP%NPts( Medium ) = N1
          IF ( Medium == 1 ) SSP%Depth( 1 ) = SSP%z( 1 )
          IF ( SSP%NPts( Medium ) == 1 ) THEN
              WRITE( PRTFile, * ) '#SSP points: ', SSP%NPts( Medium )
              CALL ERROUT( PRTFile, 'F', 'ReadSSP', 'The SSP must have at least 2 points in each layer' )
          END IF

          RETURN
       ENDIF

       N1 = N1 + 1
    END DO

    ! Fall through means too many points in the profile
    WRITE( PRTFile, * ) 'Max. #SSP points: ', MaxSSP
    CALL ERROUT( PRTFile, 'F', 'ReadSSP', 'Number of SSP points exceeds limit' )

  END SUBROUTINE ReadSSP

  !**********************************************************************!

  SUBROUTINE UpdateSSPLoss( freq, freq0 )
    ! Updates the imaginary part of the sound speed based on the frequency
    ! The depth of the SSP point is also used if there is a bio layer
    USE AttenMod

    REAL     (KIND=8), INTENT(IN) :: freq, freq0   ! freq0 is the reference frequency where dB/m was specified
    INTEGER                       :: Medium
    INTEGER                       :: IBCBeg, IBCEnd

    DO Medium = 1, SSP%NMedia
       ILoc = SSP%Loc( Medium )

       DO iSSP = 1, SSP%NPts( Medium )
          iz = SSP%Loc( Medium ) + iSSP
          SSP%cp( iz ) = CRCI( SSP%z( iz ), SSP%alphaR( iz ),  SSP%alphaI( iz ), freq, freq0, &
             SSP%AttenUnit, SSP%beta( Medium), SSP%ft( Medium ) )
          SSP%cs( iz ) = CRCI( SSP%z( iz ), SSP%betaR(  iz ),  SSP%betaI(  iz ), freq, freq0, &
             SSP%AttenUnit, SSP%beta( Medium), SSP%ft( Medium ) )

          SSP%cpSpline(  1, iz ) = SSP%cp(  iz )
          SSP%csSpline(  1, iz ) = SSP%cs(  iz )
          SSP%rhoSpline( 1, iz ) = SSP%rho( iz )
       END DO

       ! Compute spline coefs if spline interpolation has been selected
       IF ( SSP%Type == 'S' ) THEN
          IBCBeg = 0
          IBCEnd = 0
          CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%cpSpline(  1, ILoc + 1 ), SSP%NPts( Medium ), IBCBeg, IBCEnd, SSP%NPts( Medium ) )
          CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%csSpline(  1, ILoc + 1 ), SSP%NPts( Medium ), IBCBeg, IBCEnd, SSP%NPts( Medium ) )
          CALL CSPLINE( SSP%z( ILoc + 1 ), SSP%rhoSpline( 1, ILoc + 1 ), SSP%NPts( Medium ), IBCBeg, IBCEnd, SSP%NPts( Medium ) )
       END IF
    END DO
  END SUBROUTINE UpdateSSPLoss

END MODULE sspmod

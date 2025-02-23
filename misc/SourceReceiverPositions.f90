MODULE SourceReceiverPositions
  ! Reads in source depths, receiver depths, receiver ranges, and receiver bearings

  USE monotonicMod
  USE SortMod
  USE SubTabulate

  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER         :: Number_to_Echo = 21
  INTEGER                    :: Nfreq          ! number of frequencies
  REAL (KIND=8), ALLOCATABLE :: freqVec( : )   ! frequency vector for braodband runs

  TYPE Position
     INTEGER              :: Nsx, Nsy, Nsz, Nrz, Nrr, Ntheta   ! number of x, y, z, r, theta coordinates
     REAL                 :: Delta_r, Delta_theta
     INTEGER, ALLOCATABLE :: iSz( : ), iRz( : )
     REAL,    ALLOCATABLE :: Sx( : ), Sy( : ), Sz( : )   ! source x, y, z coordinates
     REAL,    ALLOCATABLE :: Rz( : ), ws( : ), wr( : ), Rr( : )
     REAL,    ALLOCATABLE :: theta( : )         ! receiver bearings
  END TYPE Position

  TYPE (Position ) :: Pos   ! structure containing source and receiver positions

CONTAINS

  SUBROUTINE ReadfreqVec( ENVFile, PRTFile, freq, BroadbandOption )

    ! Optionally reads a vector of source frequencies for a broadband run
    ! allocating and creating a frequency vector
    !
    ! If the broadband option is not selected, then the input freq is stored in the frequency vector

    IMPLICIT NONE
    INTEGER,       INTENT( IN ) :: ENVFile, PRTFile
    REAL (KIND=8), INTENT( IN ) :: freq             ! default frequency
    CHARACTER,     INTENT( IN ) :: BroadbandOption*( 1 )
    INTEGER                     :: IAllocStat, ifreq

    Nfreq = 1

    ! Broadband run?
    IF ( BroadbandOption == 'B' ) THEN
       READ( ENVFile, * ) Nfreq
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of frequencies =', Nfreq
    END IF

    IF ( Nfreq <= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadEnvironment', 'Number of frequencies must be positive'  )

    IF ( ALLOCATED( FreqVec ) ) DEALLOCATE( FreqVec )
    ALLOCATE( FreqVec( MAX( 3, Nfreq ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadEnvironment', 'Too many frequencies'  )

    IF ( BroadbandOption == 'B' ) THEN
       WRITE( PRTFile, * ) 'Frequencies (Hz)'
       FreqVec( 3 ) = -999.9
       READ(  ENVFile, * ) FreqVec( 1 : Nfreq )
       CALL SubTab( FreqVec, Nfreq )

       WRITE( PRTFile, "( 5G14.6 )" ) ( freqVec( ifreq ), ifreq = 1, MIN( Nfreq, Number_to_Echo ) )
       IF ( Nfreq > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', freqVec( Nfreq )
    ELSE
       freqVec( 1 ) = freq
    END IF

    RETURN

  END SUBROUTINE ReadfreqVec

  !********************************************************************!

  SUBROUTINE ReadSxSy( ENVFile, PRTFile, ThreeD )

    ! Reads source x-y coordinates

    IMPLICIT NONE
    LOGICAL, INTENT( IN ) :: ThreeD   ! flag indicating whether this is a 3D run
    INTEGER, INTENT( IN ) :: ENVFile, PRTFile
    INTEGER               :: is, IAllocStat

    IF ( ThreeD ) THEN

       ! *** Read source x coordinates ***

       READ(  ENVFile, * ) Pos%NSx
       WRITE( PRTFile, * ) '__________________________________________________________________________'
       WRITE( PRTFile, * )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of source x coordinates = ', Pos%NSx

       IF ( Pos%NSx <= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSxSy', 'Number of source x coordinates must be positive' )

       IF ( ALLOCATED( Pos%Sx ) ) DEALLOCATE( Pos%Sx )
       ALLOCATE( Pos%Sx( MAX( 3, Pos%NSx ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'Readsxsy', 'Too many sources'  )

       WRITE( PRTFile, * ) 'Source x coordinate (km)'
       Pos%Sx( 3 ) = -999.9
       READ( ENVFile, * ) Pos%Sx( 1 : Pos%NSx )

       CALL SubTab( Pos%Sx, Pos%NSx )
       !CALL SORT(   Pos%Sx, Pos%NSx )

       WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%Sx( is ), is = 1, MIN( Pos%NSx, Number_to_Echo ) )
       IF ( Pos%NSx > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', Pos%Sx( Pos%NSx )

       ! *** Read source y coordinates ***

       READ(  ENVFile, * ) Pos%NSy
       WRITE( PRTFile, * )
       WRITE( PRTFile, * )
       WRITE( PRTFile, * ) 'Number of source y coordinates = ', Pos%NSy

       IF ( Pos%NSy <= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSxSy', 'Number of source y coordinates must be positive' )

       IF ( ALLOCATED( Pos%Sy ) ) DEALLOCATE( Pos%Sy )
       ALLOCATE( Pos%Sy( MAX( 3, Pos%NSy ) ), Stat = IAllocStat )
       IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSxSy', 'Too many sources'  )

       WRITE( PRTFile, * ) 'Source y coordinate (km)'
       Pos%sy( 3 ) = -999.9
       READ( ENVFile, * ) Pos%Sy( 1 : Pos%NSy )

       CALL SubTab( Pos%Sy, Pos%NSy )
       !CALL SORT(   Pos%Sy, Pos%NSy )

       WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%Sy( is ), is = 1, MIN( Pos%NSy, Number_to_Echo ) )
       IF ( Pos%NSy > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', Pos%Sy( Pos%NSy )

       Pos%Sx = 1000.0 * Pos%Sx   ! convert km to m
       Pos%Sy = 1000.0 * Pos%Sy

!!$       IF ( .NOT. monotonic( Pos%Sx, Pos%NSx ) ) THEN
!!$          CALL ERROUT( PRTFile, 'F', 'SzRzRMod', 'Source x-coordinates are not monotonically increasing' )
!!$       END IF 
!!$ 
!!$       IF ( .NOT. monotonic( Pos%Sy, Pos%NSy ) ) THEN
!!$          CALL ERROUT( PRTFile, 'F', 'SzRzRMod', 'Source y-coordinates are not monotonically increasing' )
!!$       END IF 
    ELSE
       Pos%NSx = 1
       Pos%NSy = 1
	   IF ( ALLOCATED( Pos%Sx ) ) DEALLOCATE( Pos%Sx )
	   IF ( ALLOCATED( Pos%Sy ) ) DEALLOCATE( Pos%Sy )
       ALLOCATE( Pos%Sx( 1 ), Pos%Sy( 1 ) )
       Pos%Sx( 1 ) = 0.
       Pos%Sy( 1 ) = 0.
    END IF

    RETURN
  END SUBROUTINE ReadSxSy

  !********************************************************************!

  SUBROUTINE ReadSzRz( ENVFile, PRTFile, zMin, zMax )

    ! Reads source and receiver depths
    ! zMin and zMax are limits for those depths; sources and receivers are shifted to be within those limits

    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: ENVFile, PRTFile
    REAL,    INTENT( IN ) :: zMin, zMax
    !LOGICAL               :: monotonic
    INTEGER               :: is, ir, IAllocStat

    ! *** Read source depths ***

    READ(  ENVFile, * ) Pos%NSz
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '__________________________________________________________________________'
    WRITE( PRTFile, * )
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of source   depths = ', Pos%NSz

    IF ( Pos%NSz <= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSzRz', 'Number of sources must be positive'  )

    IF ( ALLOCATED( Pos%Sz ) ) DEALLOCATE( Pos%Sz, Pos%ws, Pos%iSz )
    ALLOCATE( Pos%Sz( MAX( 3, Pos%NSz ) ), Pos%ws( Pos%NSz ), Pos%iSz( Pos%NSz ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSzRz', 'Too many sources'  )

    WRITE( PRTFile, * ) 'Source   depths (m)'
    Pos%sz( 3 ) = -999.9
    READ( ENVFile, * ) Pos%Sz( 1 : Pos%NSz )

    CALL SubTab( Pos%sz, Pos%NSz )
    !CALL SORT(   Pos%sz, Pos%NSz )

    WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%Sz( is ), is = 1, MIN( Pos%NSz, Number_to_Echo ) )
    IF ( Pos%NSz > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', Pos%sz( Pos%NSz )

    ! *** Read receiver depths ***

    READ(  ENVFile, * ) Pos%NRz
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of receiver depths = ', Pos%NRz

    IF ( Pos%NRz <= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSdRz', 'Number of receivers must be positive'  )

    IF ( ALLOCATED( Pos%rz ) ) DEALLOCATE( Pos%Rz, Pos%wr, Pos%iRz )
    ALLOCATE( Pos%Rz( MAX( 3, Pos%NRz ) ), Pos%wr( Pos%NRz ), Pos%iRz( Pos%NRz ), Stat = IAllocStat  )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadSzRz', 'Too many receivers'  )

    WRITE( PRTFile, * ) 'Receiver depths (m)'
    Pos%Rz( 3 ) = -999.9
    READ( ENVFile, * ) Pos%Rz( 1 : Pos%NRz )

    CALL SubTab( Pos%Rz, Pos%NRz )
    !CALL SORT(   Pos%Rz, Pos%NRz )

    WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%Rz( ir ), ir = 1, MIN( Pos%NRz, Number_to_Echo ) )
    IF ( Pos%NRz > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', Pos%Rz( Pos%NRz )

    ! *** Check for Sz/Rz in upper or lower halfspace ***

    IF ( ANY( Pos%Sz( 1 : Pos%NSz ) < zMin ) ) THEN
       WHERE ( Pos%Sz < zMin ) Pos%Sz = zMin
       CALL ERROUT( PRTFile, 'W', 'SzRzRMod', 'Source above or too near the top bdry has been moved down' )
    END IF

    IF ( ANY( Pos%Sz( 1 : Pos%NSz ) > zMax ) ) THEN
       WHERE( Pos%Sz > zMax ) Pos%Sz = zMax
       CALL ERROUT( PRTFile, 'W', 'SzRzRMod', 'Source below or too near the bottom bdry has been moved up' ) 
    END IF

    IF ( ANY( Pos%Rz( 1 : Pos%NRz ) < zMin ) ) THEN
       WHERE( Pos%Rz < zMin ) Pos%Rz = zMin
       CALL ERROUT( PRTFile, 'W', 'SzRzRMod', 'Receiver above or too near the top bdry has been moved down' ) 
    END IF

    IF ( ANY( Pos%Rz( 1 : Pos%NRz ) > zMax ) ) THEN
       WHERE( Pos%Rz > zMax ) Pos%Rz = zMax
       CALL ERROUT( PRTFile, 'W', 'SzRzRMod', 'Receiver below or too near the bottom bdry has been moved up' ) 
    END IF

!!$    IF ( .NOT. monotonic( Pos%sz, Pos%NSz ) ) THEN
!!$       CALL ERROUT( PRTFile, 'F', 'SzRzRMod', 'Source depths are not monotonically increasing' )
!!$    END IF 
!!$ 
!!$    IF ( .NOT. monotonic( Pos%rz, Pos%NRz ) ) THEN
!!$       CALL ERROUT( PRTFile, 'F', 'SzRzRMod', 'Receiver depths are not monotonically increasing' )
!!$    END IF 

    RETURN
  END SUBROUTINE ReadSzRz

  !********************************************************************!

  SUBROUTINE ReadRcvrRanges( ENVFile, PRTFile )

    ! Read receiver ranges

    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: ENVFile, PRTFile
    INTEGER               :: ir, IAllocStat

    READ(  ENVFile, * ) Pos%NRr
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of receiver ranges = ', Pos%NRr

    IF ( ALLOCATED( Pos%rr ) ) DEALLOCATE( Pos%rr )
    ALLOCATE( Pos%rr( MAX( 3, Pos%NRr ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadRcvrRanges', 'Too many range points' )

    WRITE( PRTFile, * ) 'Receiver ranges (km)'
    Pos%Rr( 3 ) = -999.9
    READ( ENVFile, * ) Pos%Rr( 1 : Pos%NRr )

    CALL SubTab( Pos%Rr, Pos%NRr )
    CALL Sort(   Pos%Rr, Pos%NRr )

    WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%Rr( ir ), ir = 1, MIN( Pos%NRr, Number_to_Echo ) )
    IF ( Pos%NRr > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', Pos%Rr( Pos%NRr )

    Pos%Rr( 1 : Pos%NRr ) = 1000.0 * Pos%Rr( 1 : Pos%NRr )   ! Convert ranges to meters

    ! calculate range spacing
    Pos%delta_r = 0.0
    IF ( Pos%NRr /= 1 ) Pos%delta_r = Pos%Rr( Pos%NRr ) - Pos%Rr( Pos%NRr - 1 )

    ! For a point source can't have receiver at origin
    ! IF ( OPT( 1 : 1 ) == 'R' .AND. Pos%rr( 1 ) <= 0.0 ) 
    ! IF ( Pos%rr( 1 ) <= 0.0 ) Pos%rr( 1 ) = MIN( 1.0, Pos%rr( 2 ) )

    IF ( .NOT. monotonic( Pos%rr, Pos%NRr ) ) THEN
       CALL ERROUT( PRTFile, 'F', 'SzRzRMod', 'Receiver ranges are not monotonically increasing' )
    END IF 
 
    RETURN
  END SUBROUTINE ReadRcvrRanges

  !********************************************************************!

  SUBROUTINE ReadRcvrBearings( ENVFile, PRTFile )

    ! Read receiver bearings

    IMPLICIT NONE
    INTEGER, INTENT( IN ) :: ENVFile, PRTFile
    INTEGER               :: itheta, IAllocStat

    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) '__________________________________________________________________________'
    WRITE( PRTFile, * )

    READ(  ENVFile, * ) Pos%Ntheta
    WRITE( PRTFile, * )
    WRITE( PRTFile, * ) 'Number of receiver bearings = ', Pos%Ntheta
    WRITE( PRTFile, * ) 'Receiver bearings (degrees)'

    IF ( ALLOCATED( Pos%theta ) ) DEALLOCATE( Pos%theta )
    ALLOCATE( Pos%theta( MAX( 3, Pos%Ntheta ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadBearings', 'Too many bearing angles' )

    Pos%theta( 3 ) = -999.9
    READ( ENVFile, * ) Pos%theta( 1 : Pos%Ntheta )

    CALL SubTab( Pos%theta, Pos%Ntheta )
    CALL Sort(   Pos%theta, Pos%Ntheta )

    WRITE( PRTFile, "( 5G14.6 )" ) ( Pos%theta( itheta ), itheta = 1, MIN( Pos%Ntheta, Number_to_Echo ) )
    IF ( Pos%Ntheta > Number_to_Echo ) WRITE( PRTFile,  "( G14.6 )" ) ' ... ', Pos%theta( Pos%Ntheta )

    ! full 360-degree sweep? remove duplicate angle
    IF ( Pos%theta( Pos%Ntheta ) == Pos%theta( 1 ) + 360.0D0 ) Pos%Ntheta = Pos%Ntheta - 1

    ! calculate angular spacing
    Pos%Delta_theta = 0.0
    IF ( Pos%Ntheta /= 1 ) Pos%Delta_theta = Pos%theta( Pos%Ntheta ) - Pos%theta( Pos%Ntheta - 1 )

    IF ( .NOT. monotonic( Pos%theta, Pos%Ntheta ) ) THEN
       CALL ERROUT( PRTFile, 'F', 'SzRzRMod', 'Receiver bearings are not monotonically increasing' )
    END IF 
 
    RETURN
  END SUBROUTINE ReadRcvrBearings

END MODULE SourceReceiverPositions

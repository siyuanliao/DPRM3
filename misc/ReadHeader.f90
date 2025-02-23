SUBROUTINE ReadHeader( SHDFile, FileName, Title, Atten, PlotType )

  ! Read header from disk file
  ! This routine is not currently used anywhere

  ! FileName is a SHDFIL for complex pressure or a GRNFIL for a Green's function
  ! Title   arbitrary title

  ! variables taken from SourceReceiverPositions:
  ! FreqVec vector of frequencies
  ! theta   vector of bearing lines,   theta( 1 : Ntheta )
  ! sz      vector of source   depths, sz(    1 : Nsz    )
  ! rz      vector of receiver depths, rz(    1 : Nrz    )
  ! r       vector of receiver ranges, r(     1 : Nr     )

  USE SourceReceiverPositions
  IMPLICIT NONE
  INTEGER, PARAMETER                :: PRTFile = 6
  REAL,               INTENT( OUT ) :: Atten           ! stabilizing attenuation for SCOOTER FFP runs
  INTEGER                           :: SHDFile         ! unit number of SHDFile
  CHARACTER (LEN=80), INTENT( OUT ) :: Title, FileName
  CHARACTER (LEN=10), INTENT( OUT ) :: PlotType
  INTEGER                           :: IAllocStat, IOStat, LRecL

  ! Open file, read header
  IF ( SHDFile == 0 ) SHDFile = 25
  IF ( FileName( 1 : 1 ) == ' ' ) FileName = 'SHDFIL'

  ! INQUIRE( FILE = FileName, RECL = IRECL )
  OPEN( UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4, &
        IOSTAT = IOStaT, ACTION = 'READ' )
  IF ( IOStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadHeader', 'Unable to open shade file' )
  
  READ( SHDFile, REC = 1 ) LRecl
  CLOSE( UNIT = SHDFile )
  OPEN(  UNIT = SHDFile,   FILE = FileName, STATUS = 'OLD', ACCESS = 'DIRECT', FORM = 'UNFORMATTED', RECL = 4 * LRecl )

  READ( SHDFile, REC = 1  ) LRecl, Title
  READ( SHDFile, REC = 2  ) PlotType
  READ( SHDFile, REC = 3  ) Nfreq, Pos%Ntheta, Pos%Nsx, Pos%Nsy, Pos%Nsz, Pos%Nrz, Pos%Nrr, atten

  ALLOCATE( FreqVec( Nfreq ), Pos%sz( Pos%Nsz ), Pos%rz( Pos%Nrz ), Pos%rr( Pos%Nrr ), Pos%theta( Pos%Ntheta ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'ReadHeader', 'Too many source/receiver combinations' )

  READ( SHDFile, REC = 4  ) FreqVec
  READ( SHDFile, REC = 5  ) Pos%theta
  READ( SHDFile, REC = 6  ) Pos%sx
  READ( SHDFile, REC = 7  ) Pos%sy
  READ( SHDFile, REC = 8  ) Pos%sz
  READ( SHDFile, REC = 9  ) Pos%rz
  READ( SHDFile, REC = 10 ) Pos%rr

  ! Pos%deltaR = Pos%r( Pos%Nrr ) - Pos%r( Pos%Nrr - 1 )

END SUBROUTINE ReadHeader

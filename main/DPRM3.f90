PROGRAM DPRM3
  
  USE MathConstants
  USE ReadEnvironmentBell
  USE BeamPattern
  USE bdryMod
  USE RefCoef
  USE SourceReceiverPositions
  USE angleMod
  USE sspMod
  USE influence
  USE ArrMod
  use SourceReceiverPositions
  use mpi
  USE omp_lib

  IMPLICIT NONE
  
  LOGICAL, PARAMETER   :: ThreeD = .FALSE., Inline = .FALSE.
  INTEGER              :: j, k, m, n, p1, p2, nt, maxNumarr, g, numnarr, y, maxindex, tt, lsy1, lsy2, lsy3, time1, time2
  INTEGER              :: time3, time4, time5, time6
  CHARACTER ( LEN=2 )  :: AttenUnit
  CHARACTER ( LEN=80 ) :: FileRoot
  CHARACTER*20 cTemp
  INTEGER, PARAMETER :: RLFile = 66
  INTEGER,allocatable :: timeindex(:)
  REAL,allocatable :: alldelay(:)
  REAL  :: delay1, delay2, del
  REAL,allocatable    :: allenergy(:), Arever(:), RL(:), tar1(:), tar2(:), myArever(:)
  REAL :: dr, maxdep, maxrange, timestep, A1, A2, Recang1, Recang2, theta, S, sigma, energy, tmax, factor, dz
  INTEGER :: mpiid, mpinps, ierr, worknum, ip, mywork, totaltime
  INTEGER :: stats(MPI_STATUS_SIZE)
  ! get the file root for naming all input and output files
  ! should add some checks here ...
  !mpi环境初始化
  call MPI_INIT(ierr)
  call system_clock(time1)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpiid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpinps, ierr)
  CALL GET_COMMAND_ARGUMENT( 1, FileRoot )
  read(FileRoot,*) nt !扇面循环
  g=1 !计数，用于存储能量-时延到达结构
  timestep=0.1 !混响计算时间步长
  maxrange=50000 !最大水平范围
  
  if (mpiid==0) then 
    write(0,*) "参与工作的进程数量为",mpinps
	do worknum= mpinps, nt+mpinps-1
		call MPI_RECV(ip,1,MPI_INT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,stats,ierr)
		call MPI_SEND(worknum,1,MPI_INT,ip,1,MPI_COMM_WORLD,ierr)
		!write(0,*),"进程",ip,"在完成第",worknum,"个任务"
	enddo
	!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	call MPI_RECV(maxindex,1,MPI_INT,1,2,MPI_COMM_WORLD,stats,ierr)
	allocate(Arever(maxindex))
	Arever(:)=0
	allocate(RL(maxindex))
  else
    totaltime=0
	mywork=mpiid
	do while(mywork<=nt)
		call system_clock(time3)
		write(cTemp,'(i3)')mywork
		FileRoot=AdjustL(TRIM("task"))//AdjustL(TRIM(cTemp))
		OPEN( UNIT = PRTFile, FILE = TRIM( FileRoot ) // '.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )
		! Read in or otherwise initialize inline all the variables used by BELLHOP
		CALL ReadEnvironment(  FileRoot, ThreeD )
		CALL ReadATI( FileRoot, Bdry%Top%HS%Opt( 5 : 5 ), Bdry%Top%HS%Depth, PRTFile )                          ! AlTImetry
		CALL ReadBTY( FileRoot, Bdry%Bot%HS%Opt( 2 : 2 ), Bdry%Bot%HS%Depth, PRTFile )                          ! BaThYmetry
		CALL ReadReflectionCoefficient( FileRoot, Bdry%Bot%HS%Opt( 1 : 1 ), Bdry%Top%HS%Opt( 2 : 2 ), PRTFile ) ! (top and bottom)
		SBPFlag = Beam%RunType( 3 : 3 )
		CALL ReadPat( FileRoot, PRTFile )   ! Source Beam Pattern
		! dummy bearing angles
		Pos%Ntheta = 1
		ALLOCATE( Pos%theta( Pos%Ntheta ), Stat = IAllocStat )
		Pos%theta( 1 ) = 0.
		CALL OpenOutputFiles( FileRoot, ThreeD )
		CALL BellhopCore
		maxNumarr=maxval(Narr)
		IF ( .not. ALLOCATED( allenergy )) then
			allocate(allenergy(size(NArr,dim=2)*maxNumarr**2))
			allenergy(:)=-1
		ENDIF
		IF ( .not. ALLOCATED( alldelay )) THEN
			allocate(alldelay(size(NArr,dim=2)*maxNumarr**2))
			alldelay(:)=-1
		ENDIF
		IF ( .not. ALLOCATED( timeindex )) allocate(timeindex(size(NArr,dim=2)*maxNumarr**2))
		IF ( .not. ALLOCATED( tar1 )) allocate(tar1(Pos%Nrz))
		IF ( .not. ALLOCATED( tar2 )) allocate(tar2(Pos%Nrz))
		
		!射线幅值修正
		do lsy1=1, Pos%Nrz
			do lsy2=1, Pos%Nrr
				IF ( Pos%Rr ( lsy2 ) == 0 ) THEN
					factor = 1e5   
				else
					factor = 1. / SQRT( Pos%Rr( lsy2 ) )
				endif
				do lsy3=1, Narr(lsy1,lsy2)
					Arr( lsy1, lsy2, lsy3)%A=factor*Arr( lsy1, lsy2, lsy3)%A
				enddo
			enddo
		enddo	
		
		dr=(Pos%Rr(2)-Pos%Rr(1))
		!maxdep=maxval(Bot( : )%x(2)) !最大海底深度
		!write(0,*) "NArr维度", size(NArr,dim=1), size(NArr,dim=2)
		
		theta=(360/nt)*4*ATAN2(1.0,1.0)/180 !单个散射单元所占的角度
		
		
		!海底散射单元循环
		do j=2, Pos%Nrr
			!Pos%Rr水平接收位置矩阵 单位是m，Pos%Nrr水平接收位置矩阵个数，Pos%Rz垂直接收位置矩阵 单位是m，Pos%Nrz垂直接收位置矩阵个数
			!Bot是地形，Bot( 1 )%x是水平位置和深度，可以直接打印出来  Bot( 1 )%x(1)是水平位置  Bot( 1 )%x(2)是深度位置
			!p是最接近海底的接收器位置
			tar1=abs(Pos%Rz-Bot(j-1)%x(2))!单位是m
			p1=maxval([minloc(tar1,1)-1, 1]) !减1是为了得到海底边界以上的接收器位置
			tar2=abs(Pos%Rz-Bot(j)%x(2))!单位是m
			p2=maxval([minloc(tar2,1)-1, 1]) !减1是为了得到海底边界以上的接收器位置
			numnarr=Narr(p2,j)
			S=theta*dr*Pos%Rr(j) !散射单元的面积
			!根据坡度修正散射单元的面积
			dz=Pos%Rz(p2)-Pos%Rz(p1)
			S=(sqrt(dz**2+dr**2)/dr)*S 
			
			do m=1, numnarr
				do n=1, numnarr
					!入射声压/时延/角度
					A1=Arr(p2,j,m)%A !类型是REAL
					delay1=abs(Arr(p2,j,m)%delay) !类型是COMPLEX
					Recang1=Arr(p2,j,m)%RcvrDeclAngle
					!出射声压/时延/角度
					A2=Arr(p2,j,n)%A
					delay2=abs(Arr(p2,j,n)%delay)
					Recang2=Arr(p2,j,n)%RcvrDeclAngle
					!剔除假声线
					if (A1==0 .or. A2==0 .or. delay1==0 .or. delay2==0) then
						continue
					endif
					
					!求解在第j个单元，每组组合射线的三维散射系数
					call Scatter3D(Recang1,Recang2,sigma)
					!射线组合最终到达接收器的能量和时延
					energy=A1**2*A2**2*S*sigma**2
					del=delay1+delay2
					allenergy(g)=energy
					alldelay(g)=del
					g=g+1
				enddo
			enddo
			tmax=2*sqrt(4000**2+maxrange**2)/1500
			maxindex=nint( tmax/timestep )
			do tt=1,g-1
				timeindex(tt)=nint(abs(alldelay(tt)/timestep))
			enddo
			IF ( .not. ALLOCATED( Arever )) then 
				allocate(Arever(maxindex))
				Arever(:)=0
			endif
			!$OMP PARALLEL DEFAULT(shared)  &
			!$OMP PRIVATE(y,tt)
			!$OMP DO
			do y=1, maxindex
				do tt=1,g-1
					if (timeindex(tt)==y) then
						Arever(y)=Arever(y)+allenergy(tt)
					endif
				enddo
			enddo
			!$OMP END DO
			!$OMP END PARALLEL
			g=1
			
		enddo
		call system_clock(time4)
		totaltime=totaltime+time4-time3
		call MPI_SEND(mpiid,1,MPI_INT,0,1,MPI_COMM_WORLD,ierr)
		call MPI_RECV(mywork,1,MPI_INT,0,1,MPI_COMM_WORLD,stats,ierr)
	enddo
	!call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	WRITE( 0, * ) '进程', mpiid, '的计算时间为', totaltime, 'ms'
	
	if (mpiid==1) then
		call MPI_SEND(maxindex,1,MPI_INT,0,2,MPI_COMM_WORLD,ierr)
	endif	
  endif
  
  allocate(myArever(maxindex))
  call MPI_REDUCE(Arever(1),myArever(1),maxindex,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  
  if (mpiid==0) then 
	do y=1, maxindex
		RL(y)=10*log10(abs(myArever(y)))
	enddo
	OPEN ( FILE = TRIM( 'myreverberation' ) // '.rl', UNIT = RLFile, FORM = 'FORMATTED' )
	WRITE( RLFile, * ) timestep, maxindex
	WRITE( RLFile, * ) RL(:)
	close(RLFile)
	call system_clock(time2)
	WRITE( 0, * ) '总的计算时间为:', time2-time1, 'ms'
  endif
  call MPI_FINALIZE(ierr)
  

  CONTAINS

! **********************************************************************!

SUBROUTINE Scatter3D(Recang1,Recang2,sigma)

  REAL :: Recang1, Recang2
  REAL :: sigma
  
  sigma=0.05*sin(Recang1*3.14159/180)*sin(Recang2*3.14159/180)

END SUBROUTINE Scatter3D


SUBROUTINE BellhopCore

  USE ArrMod
  USE AttenMod
  USE omp_lib

  INTEGER, PARAMETER   :: SHDFile = 25, RAYFile = 21, ArrivalsStorage = 20000000
  INTEGER, PARAMETER   :: RAYTRACEFILE = 50
  INTEGER              :: IBPvec( 1 ), ibp, is, iBeamWindow2, Irz1, Irec, NalphaOpt, iSeg, &
			Myid, threadNUM,is0_red, is1_red, is2_red
  REAL                 :: Tstart, Tstop, TStart_openmp, TStop_openmp, Time_openmp
  REAL        (KIND=8) :: Amp0, alpha0, DalphaOpt, xs( 2 ), RadMax, s, &
                          c, cimag, gradc( 2 ), crr, crz, czz, rho
  COMPLEX, ALLOCATABLE :: U( :, : )
  COMPLEX, ALLOCATABLE :: U_shared( :, :, : )
  COMPLEX     (KIND=8) :: epsilon

  CHARACTER (LEN=80)   :: filename
  INTEGER              :: BeamNsteps
  TYPE( ray2DPt )      :: ray2D( MaxN )
  INTEGER              :: iSegz=1
  INTEGER              :: iSegr=1
  INTEGER              :: iSegx=1
  INTEGER              :: iSegy=1
  

  filename = 'ray2dstep'

  open(Unit=RAYTRACEFILE, FILE=trim(filename)//'.prt', status='unknown', IOSTAT=iostat)

  CALL CPU_TIME( Tstart )

  omega = 2.0 * pi * freq

  IF ( Beam%deltas == 0.0 ) THEN
     Beam%deltas = ( Bdry%Bot%HS%Depth - Bdry%Top%HS%Depth ) / 10.0   ! Automatic step size selection
     WRITE( PRTFile, * )
     WRITE( PRTFile, fmt = '(  '' Step length,       deltas = '', G11.4, '' m (automatically selected)'' )' ) Beam%deltas
  END IF

  Angles%alpha  = DegRad * Angles%alpha   ! convert to radians
  Angles%Dalpha = 0.0
  IF ( Angles%Nalpha /= 1 ) &
       Angles%Dalpha = ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) / ( Angles%Nalpha - 1 )  ! angular spacing between beams

  ! convert range-dependent geoacoustic parameters from user to program units
  IF ( atiType( 2 : 2 ) == 'L' ) THEN
     DO iSeg = 1, NatiPts
        Top( iSeg )%HS%cp = CRCI( 1D20, Top( iSeg )%HS%alphaR, Top( iSeg )%HS%alphaI, freq, freq, 'W ', &
             betaPowerLaw, ft )   ! compressional wave speed
        Top( iSeg )%HS%cs = CRCI( 1D20, Top( iSeg )%HS%betaR,  Top( iSeg )%HS%betaI,  freq, freq, 'W ', &
             betaPowerLaw, ft )   ! shear         wave speed
     END DO
  END IF
   
  IF ( btyType( 2 : 2 ) == 'L' ) THEN
     DO iSeg = 1, NbtyPts
        Bot( iSeg )%HS%cp = CRCI( 1D20, Bot( iSeg )%HS%alphaR, Bot( iSeg )%HS%alphaI, freq, freq, 'W ', &
             betaPowerLaw, ft )   ! compressional wave speed 
        Bot( iSeg )%HS%cs = CRCI( 1D20, Bot( iSeg )%HS%betaR,  Bot( iSeg )%HS%betaI,  freq, freq, 'W ', &
             betaPowerLaw, ft )   ! shear         wave speed
     END DO
  END IF

  SELECT CASE ( Beam%RunType( 5 : 5 ) )
  CASE ( 'I' )
     Nrz_per_range = 1         ! irregular grid
  CASE DEFAULT
     Nrz_per_range = Pos%Nrz   ! rectilinear grid
  END SELECT

  ! for a TL calculation, allocate space for the pressure matrix
  IF ( ALLOCATED( U ) ) DEALLOCATE( U )
  ALLOCATE ( U( 1, 1 ), Stat = IAllocStat )   ! open a dummy variable
  IF ( ALLOCATED( U_shared ) ) DEALLOCATE( U )
  ALLOCATE ( U_shared( 1, 1, Maxthreads ), Stat = IAllocStat )   ! open a dummy variable

  ! for an arrivals run, allocate space for arrivals matrices
  MaxNArr = MAX( ArrivalsStorage / ( Nrz_per_range * Pos%Nrr ), 10 )   ! allow space for at least 10 arrivals
  WRITE( PRTFile, * )
  WRITE( PRTFile, * ) '( Maximum # of arrivals = ', MaxNArr, ')'
  IF ( ALLOCATED( Arr ) ) DEALLOCATE( Arr )
  IF ( ALLOCATED( NArr ) ) DEALLOCATE( NArr )
  ALLOCATE ( Arr( Nrz_per_range, Pos%Nrr, MaxNArr ), NArr( Nrz_per_range, Pos%Nrr ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( PRTFile, 'F', 'BELLHOP', &
          'Insufficient memory to allocate arrivals matrix; reduce parameter ArrivalsStorage' )

  NArr( 1:Nrz_per_range, 1:Pos%Nrr ) = 0

  WRITE( PRTFile, * )

  SourceDepth: DO is = 1, Pos%Nsz
     xs = [ 0.0, Pos%sz( is ) ]   ! source coordinate

     SELECT CASE ( Beam%RunType( 1 : 1 ) )
     CASE ( 'C', 'S', 'I' ) ! TL calculation, zero out pressure matrix
        U = 0.0
        U_shared = 0.0
     CASE ( 'A', 'a' )      ! Arrivals calculation, zero out arrival matrix
        NArr = 0
     END SELECT

     CALL EvaluateSSP( xs, c, cimag, gradc, crr, crz, czz, rho, Freq, 'TAB', iSegz, iSegr, iSegx, iSegy )
     RadMax = 5 * c / freq  ! 5 wavelength max radius

     ! Are there enough beams?
     DalphaOpt = SQRT( c / ( 6.0 * freq * Pos%rr( Pos%Nrr ) ) )
     NalphaOpt = 2 + INT( ( Angles%alpha( Angles%Nalpha ) - Angles%alpha( 1 ) ) / DalphaOpt )

     IF ( Beam%RunType( 1 : 1 ) == 'C' .AND. Angles%Nalpha < NalphaOpt ) THEN
        CALL ERROUT( PRTFile, 'W', 'BELLHOP', 'Too few beams' )
        WRITE( PRTFile, * ) 'Nalpha should be at least = ', NalphaOpt
     ENDIF

     !openmp并行开始
     !TStart_openmp = OMP_get_wtime()
     ! Trace successive beams
     !$OMP PARALLEL DEFAULT(shared)  &
     !$OMP PRIVATE(ialpha, alpha0,IBPvec, IBP, s, Amp0, ibeamwindow2, RadMax, BeamNsteps, epsilon, Myid, ray2D) &
     !$OMP FIRSTPRIVATE (iSegz, iSegr, iSegx, iSegy)
         Myid = OMP_GET_THREAD_NUM()
	 threadNUM = OMP_GET_NUM_THREADS()
         !write(0, *), 'thread id = ',  Myid
         !OPEN( UNIT = (RAYTRACEFile+Myid), FILE = trim(filename) // '.prt', STATUS = 'UNKNOWN', IOSTAT = iostat )
     !$OMP DO
     ElevationAngle: DO ialpha = 1, Angles%Nalpha

        IF ( Angles%iSingle == 0 .OR. ialpha == Angles%iSingle ) THEN    ! Single beam run?

           alpha0 = Angles%alpha( ialpha ) * RadDeg   ! take-off angle in degrees

           IBPvec = maxloc( SrcBmPat( :, 1 ), mask = SrcBmPat( :, 1 ) < alpha0 )  ! index of ray angle in beam pattern
           IBP    = IBPvec( 1 )
           IBP    = MAX( IBP, 1 )               ! don't go before beginning of table
           IBP    = MIN( IBP, NSBPPts - 1 )     ! don't go past end of table

           ! linear interpolation to get amplitude
           s    = ( alpha0  - SrcBmPat( IBP, 1 ) ) / ( SrcBmPat( IBP + 1, 1 ) - SrcBmPat( IBP, 1 ) )
           Amp0 = ( 1 - s ) * SrcBmPat( IBP, 2 ) + s * SrcBmPat( IBP + 1, 2 )

           ! Lloyd mirror pattern for semi-coherent option
           IF ( Beam%RunType( 1 : 1 ) == 'S' ) &
              Amp0 = Amp0 * SQRT( 2.0 ) * ABS( SIN( omega / c * xs( 2 ) * SIN( Angles%alpha( ialpha ) ) ) )

           ! show progress ...
           IF ( MOD( ialpha - 1, max( Angles%Nalpha / 50, 1 ) ) == 0 ) THEN
              !WRITE( PRTFile, FMT = "( 'Tracing beam ', I7, F10.2 )" ) ialpha, alpha0
               !WRITE( 0, FMT = "( 'Tracing beam ', I7, F10.2 )" ) ialpha, alpha0
           END IF

           CALL TraceRay2D( xs, Angles%alpha( ialpha ), Amp0, BeamNsteps, ray2D, iSegz, iSegr, iSegx, iSegy )   ! Trace a ray

           IF (ialpha == 32) THEN
		write(RAYTRACEFILE, *) ray2D(512)
           END IF


           IF ( Beam%RunType( 1 : 1 ) == 'R' ) THEN   ! Write the ray trajectory to RAYFile
              !CALL WriteRay2D( alpha0, Beam%Nsteps )
           ELSE                                       ! Compute the contribution to the field

              epsilon = PickEpsilon( Beam%Type( 1 : 2 ), omega, c, gradc, Angles%alpha( ialpha ), &
                   Angles%Dalpha, Beam%rLoop, Beam%epsMultiplier ) ! 'optimal' beam constant

              SELECT CASE ( Beam%Type( 1 : 1 ) )
              CASE ( 'R' )
                 iBeamWindow2 = Beam%iBeamWindow **2
                 RadMax       = 50 * c / freq  ! 50 wavelength max radius
                 CALL InfluenceCervenyRayCen(   U_shared(:,:,(Myid+1)), epsilon, Angles%alpha( ialpha ), &
				iBeamWindow2, RadMax, BeamNsteps, ray2D )
              CASE ( 'C' )
                 iBeamWindow2 = Beam%iBeamWindow **2
                 RadMax       = 50 * c / freq  ! 50 wavelength max radius
                 CALL InfluenceCervenyCart(     U_shared(:,:,(Myid+1)), epsilon, Angles%alpha( ialpha ), &
				iBeamWindow2, RadMax, BeamNsteps, ray2D, iSegz, iSegr, iSegx, iSegy )
              CASE ( 'g' )
                 CALL InfluenceGeoHatRayCen(    U_shared(:,:,(Myid+1)), Angles%alpha( ialpha ), Angles%Dalpha, &
				 BeamNsteps, ray2D )              
              CASE ( 'S' )
                 CALL InfluenceSGB(             U_shared(:,:,(Myid+1)), Angles%alpha( ialpha ), Angles%Dalpha, &
				BeamNsteps, ray2D )
              CASE ( 'B' )
                 CALL InfluenceGeoGaussianCart( U_shared(:,:,(Myid+1)), Angles%alpha( ialpha ), Angles%Dalpha, &
				BeamNsteps, ray2D )
              CASE DEFAULT
                 CALL InfluenceGeoHatCart(      U_shared(:,:,(Myid+1)), Angles%alpha( ialpha ), Angles%Dalpha, &
				BeamNsteps, ray2D )
              END SELECT

           END IF
        END IF
     END DO ElevationAngle
     !$OMP END DO
     !$OMP END PARALLEL

     !TStop_openmp = OMP_get_wtime()
     !Time_openmp = TStop_openmp - TStart_openmp
     !WRITE( 0, "( /, ' 射线追踪  = ', G15.6, 's' )" ) Time_openmp

  END DO SourceDepth

  CALL CPU_TIME( Tstop )
  WRITE( PRTFile, "( /, ' CPU Time = ', G15.3, 's' )" ) Tstop - Tstart
END SUBROUTINE BellhopCore

! **********************************************************************!

COMPLEX (KIND=8 ) FUNCTION PickEpsilon( BeamType, omega, c, gradc, alpha, Dalpha, rLoop, EpsMultiplier )

  ! Picks the optimum value for epsilon

  INTEGER,            PARAMETER     :: PRTFile = 6
  COMPLEX,            PARAMETER     :: i = ( 0.0, 1.0 )
  REAL      (KIND=8), INTENT( IN  ) :: omega, c, gradc( 2 ) ! angular frequency, sound speed and gradient
  REAL      (KIND=8), INTENT( IN  ) :: alpha, Dalpha        ! angular spacing for ray fan
  REAL      (KIND=8), INTENT( IN  ) :: epsMultiplier, Rloop ! multiplier, loop range
  CHARACTER (LEN= 2), INTENT( IN  ) :: BeamType
  LOGICAL, SAVE      :: INIFlag = .TRUE.
  REAL      (KIND=8) :: HalfWidth
  REAL      (KIND=8) :: cz
  COMPLEX   (KIND=8) :: epsilonOpt
  CHARACTER (LEN=40) :: TAG

  SELECT CASE ( BeamType( 1 : 1 ) )
  CASE ( 'C', 'R' )   ! Cerveny beams
     TAG    = 'Cerveny style beam'
     SELECT CASE ( BeamType( 2 : 2 ) )
     CASE ( 'F' )
        TAG       = 'Space filling beams'
        halfwidth = 2.0 / ( ( omega / c ) * Dalpha )
        epsilonOpt    = i * 0.5 * omega * halfwidth ** 2
     CASE ( 'M' )
        TAG       = 'Minimum width beams'
        halfwidth = SQRT( 2.0 * c * 1000.0 * rLoop / omega )
        epsilonOpt    = i * 0.5 * omega * halfwidth ** 2
     CASE ( 'W' )
        TAG       = 'WKB beams'
        halfwidth = HUGE( halfwidth )
        cz        = gradc( 2 )
        IF ( cz == 0.0 ) THEN
           epsilonOpt = 1.0D10
        ELSE
           epsilonOpt = ( -SIN( alpha ) / COS( alpha ** 2 ) ) * c * c / cz
        ENDIF
     END SELECT

  CASE ( 'G', 'g' )   ! geometric hat beams
     TAG        = 'Geometic hat beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2

  CASE ( 'B' )   ! geometric Gaussian beams
     TAG        = 'Geometric Gaussian beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2

  CASE ( 'S' )   ! simple Gaussian beams
     TAG        = 'Simple Gaussian beams'
     halfwidth  = 2.0 / ( ( omega / c ) * Dalpha )
     epsilonOpt = i * 0.5 * omega * halfwidth ** 2
  END SELECT

  PickEpsilon = EpsMultiplier * epsilonOpt

  ! On first call write info to prt file
  IF ( INIFlag ) THEN
     WRITE( PRTFile, * )
     WRITE( PRTFile, * ) TAG
     WRITE( PRTFile, * ) 'halfwidth  = ', halfwidth
     WRITE( PRTFile, * ) 'epsilonOpt = ', epsilonOpt
     WRITE( PRTFile, * ) 'EpsMult    = ', EpsMultiplier
     WRITE( PRTFile, * )
     INIFlag = .FALSE.
  END IF

END FUNCTION PickEpsilon

! **********************************************************************!

SUBROUTINE TraceRay2D( xs, alpha, Amp0, BeamNsteps, ray2D, iSegz, iSegr, iSegx, iSegy)

  ! Traces the beam corresponding to a particular take-off angle

  USE WriteRay

  REAL     (KIND=8), INTENT( IN ) :: xs( 2 )      ! x-y coordinate of the source
  REAL     (KIND=8), INTENT( IN ) :: alpha, Amp0  ! initial angle, amplitude

  INTEGER  , INTENT(OUT)	  :: BeamNsteps
  TYPE( ray2DPt ), INTENT(OUT)    :: ray2D( MaxN )
  INTEGER           :: istep
  INTEGER            :: IsegTop  , IsegBot
  REAL (KIND=8)      :: rTopseg( 2 ), rBotseg( 2 )

  INTEGER           :: is, is1                    ! index for a step along the ray
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho
  REAL     (KIND=8) :: dEndTop( 2 ), dEndBot( 2 ), TopnInt( 2 ), BotnInt( 2 ), ToptInt( 2 ), BottInt( 2 )
  REAL     (KIND=8) :: DistBegTop, DistEndTop, DistBegBot, DistEndBot ! Distances from ray beginning, end to top and bottom
  REAL     (KIND=8) :: sss

  INTEGER,           INTENT(INOUT) :: iSegz
  INTEGER,           INTENT(INOUT) :: iSegr
  INTEGER,           INTENT(INOUT) :: iSegx
  INTEGER,           INTENT(INOUT) :: iSegy

  ! Initial conditions

  iSmallStepCtr = 0
  CALL EvaluateSSP( xs, c, cimag, gradc, crr, crz, czz, rho, Freq, 'TAB', iSegz, iSegr, iSegx, iSegy )
  ray2D( 1 )%c         = c
  ray2D( 1 )%x         = xs
  ray2D( 1 )%t         = [ COS( alpha ), SIN( alpha ) ] / c
  ray2D( 1 )%p         = [ 1.0, 0.0 ]
  ray2D( 1 )%q         = [ 0.0, 1.0 ]
  ray2D( 1 )%tau       = 0.0
  ray2D( 1 )%Amp       = Amp0
  ray2D( 1 )%Phase     = 0.0
  ray2D( 1 )%NumTopBnc = 0
  ray2D( 1 )%NumBotBnc = 0

  ! second component of qv is not used in geometric beam tracing
  ! set I.C. to 0 in hopes of saving run time
  IF ( Beam%RunType( 2 : 2 ) == 'G' ) ray2D( 1 )%q = [ 0.0, 0.0 ]

  CALL GetTopSeg( xs( 1 ), IsegTop, rTopseg )   ! identify the top    segment above the source
  CALL GetBotSeg( xs( 1 ), IsegBot, rBotseg )   ! identify the bottom segment below the source

  ! convert range-dependent geoacoustic parameters from user to program units
  ! compiler is not accepting the copy of the whole structure at once ...
  IF ( atiType( 2 : 2 ) == 'L' ) THEN
     Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp   ! grab the geoacoustic info for the new segment
     Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
     Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
  END IF

  IF ( btyType( 2 : 2 ) == 'L' ) THEN
     Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp
     Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
     Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
  END IF

  ! Trace the beam (note that Reflect alters the step index is)
  is = 0
  CALL Distances2D( ray2D( 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
       Top( IsegTop )%n, Bot( IsegBot )%n, DistBegTop, DistBegBot )

  IF ( DistBegTop <= 0 .OR. DistBegBot <= 0 ) THEN
     BeamNsteps = 1
     WRITE( PRTFile, * ) 'Terminating the ray trace because the source is on or outside the boundaries'
     RETURN       ! source must be within the medium
  END IF

  Stepping: DO istep = 1, MaxN - 1
     is  = is + 1
     is1 = is + 1

     CALL Step2D( ray2D( is ), ray2D( is1 ),  &
          Top( IsegTop )%x, Top( IsegTop )%n, &
          Bot( IsegBot )%x, Bot( IsegBot )%n, rTopSeg, rBotSeg, iSegz, iSegr, iSegx, iSegy )

     ! New altimetry segment?
     IF ( ray2D( is1 )%x( 1 ) < rTopSeg( 1 ) .OR. &
          ray2D( is1 )%x( 1 ) > rTopSeg( 2 ) ) THEN
        CALL GetTopSeg( ray2D( is1 )%x( 1 ), IsegTop, rTopseg )
        IF ( atiType( 2 : 2 ) == 'L' ) THEN
           Bdry%Top%HS%cp  = Top( IsegTop )%HS%cp   ! grab the geoacoustic info for the new segment
           Bdry%Top%HS%cs  = Top( IsegTop )%HS%cs
           Bdry%Top%HS%rho = Top( IsegTop )%HS%rho
        END IF
     END IF

     ! New bathymetry segment?
     IF ( ray2D( is1 )%x( 1 ) < rBotSeg( 1 ) .OR. &
          ray2D( is1 )%x( 1 ) > rBotSeg( 2 ) ) THEN
        CALL GetBotSeg( ray2D( is1 )%x( 1 ), IsegBot, rBotseg )
        IF ( btyType( 2 : 2 ) == 'L' ) THEN
           Bdry%Bot%HS%cp  = Bot( IsegBot )%HS%cp   ! grab the geoacoustic info for the new segment
           Bdry%Bot%HS%cs  = Bot( IsegBot )%HS%cs
           Bdry%Bot%HS%rho = Bot( IsegBot )%HS%rho
        END IF
     END IF

     ! Reflections?
     ! Tests that ray at step is is inside, and ray at step is+1 is outside
     ! to detect only a crossing from inside to outside
     ! DistBeg is the distance at step is,   which is saved
     ! DistEnd is the distance at step is+1, which needs to be calculated

     CALL Distances2D( ray2D( is1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
          Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     IF      ( DistBegTop > 0.0d0 .AND. DistEndTop <= 0.0d0 ) THEN  ! test top reflection

        IF ( atiType == 'C' ) THEN
           sss     = DOT_PRODUCT( dEndTop, Top( IsegTop )%t ) / Top( IsegTop )%Len   ! proportional distance along segment
           TopnInt = ( 1 - sss ) * Top( IsegTop )%Noden + sss * Top( 1 + IsegTop )%Noden
           ToptInt = ( 1 - sss ) * Top( IsegTop )%Nodet + sss * Top( 1 + IsegTop )%Nodet
        ELSE
           TopnInt = Top( IsegTop )%n   ! normal is constant in a segment
           ToptInt = Top( IsegTop )%t
        END IF

        CALL Reflect2D( is, Bdry%Top%HS, 'TOP', ToptInt, TopnInt, Top( IsegTop )%kappa, RTop, NTopPTS &
			, ray2D, iSegz, iSegr, iSegx, iSegy )
        ray2D( is + 1 )%NumTopBnc = ray2D( is )%NumTopBnc + 1

        CALL Distances2D( ray2D( is + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot,  &
             Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     ELSE IF ( DistBegBot > 0.0d0 .AND. DistEndBot <= 0.0d0 ) THEN  ! test bottom reflection

        IF ( btyType == 'C' ) THEN
           sss     = DOT_PRODUCT( dEndBot, Bot( IsegBot )%t ) / Bot( IsegBot )%Len   ! proportional distance along segment
           BotnInt = ( 1 - sss ) * Bot( IsegBot )%Noden + sss * Bot( 1 + IsegBot )%Noden
           BottInt = ( 1 - sss ) * Bot( IsegBot )%Nodet + sss * Bot( 1 + IsegBot )%Nodet
        ELSE
           BotnInt = Bot( IsegBot )%n   ! normal is constant in a segment
           BottInt = Bot( IsegBot )%t
        END IF

        CALL Reflect2D( is, Bdry%Bot%HS, 'BOT', BottInt, BotnInt, Bot( IsegBot )%kappa, RBot, NBotPTS &
			, ray2D, iSegz, iSegr, iSegx, iSegy )
        ray2D( is + 1 )%NumBotBnc = ray2D( is )%NumBotBnc + 1
        CALL Distances2D( ray2D( is + 1 )%x, Top( IsegTop )%x, Bot( IsegBot )%x, dEndTop,    dEndBot, &
             Top( IsegTop )%n, Bot( IsegBot )%n, DistEndTop, DistEndBot )

     END IF

     ! Has the ray left the box, lost its energy, escaped the boundaries, or exceeded storage limit?
     IF ( ABS( ray2D( is + 1 )%x( 1 ) ) > Beam%Box%r .OR. &
          ABS( ray2D( is + 1 )%x( 2 ) ) > Beam%Box%z .OR. ray2D( is + 1 )%Amp < 0.005 .OR. &
          ( DistBegTop < 0.0 .AND. DistEndTop < 0.0 ) .OR. &
          ( DistBegBot < 0.0 .AND. DistEndBot < 0.0 ) ) THEN
          ! ray2D( is + 1 )%t( 1 ) < 0 ) THEN ! this last test kills off a backward traveling ray
        BeamNsteps = is + 1
        EXIT Stepping
     ELSE IF ( is >= MaxN - 3 ) THEN
        CALL ERROUT( PRTFile, 'W', 'TraceRay2D', 'Insufficient storage for ray trajectory' )
        BeamNsteps = is
        EXIT Stepping
     END IF

     DistBegTop = DistEndTop
     DistBegBot = DistEndBot

  END DO Stepping
END SUBROUTINE TraceRay2D

! **********************************************************************!

SUBROUTINE Distances2D( rayx, Topx, Botx, dTop, dBot, Topn, Botn, DistTop, DistBot )

  ! Calculates the distances to the boundaries
  ! Formula differs from JKPS because code uses outward pointing normals

  REAL (KIND=8), INTENT( IN  ) :: rayx( 2 )              ! ray coordinate
  REAL (KIND=8), INTENT( IN  ) :: Topx( 2 ), Botx( 2 )   ! top, bottom coordinate
  REAL (KIND=8), INTENT( IN  ) :: Topn( 2 ), Botn( 2 )   ! top, bottom normal vector (outward)
  REAL (KIND=8), INTENT( OUT ) :: dTop( 2 ), dBot( 2 )   ! vector pointing from top, bottom bdry to ray
  REAL (KIND=8), INTENT( OUT ) :: DistTop, DistBot       ! distance (normal to bdry) from the ray to top, bottom boundary

  dTop    = rayx - Topx  ! vector pointing from top    to ray
  dBot    = rayx - Botx  ! vector pointing from bottom to ray
  DistTop = -DOT_PRODUCT( Topn, dTop )
  DistBot = -DOT_PRODUCT( Botn, dBot )

END SUBROUTINE Distances2D

! **********************************************************************!

SUBROUTINE Step2D( ray0, ray2, Topx, Topn, Botx, Botn, rTopSeg, rBotSeg, iSegz, iSegr, iSegx, iSegy )

  ! Does a single step along the ray
  ! x denotes the ray coordinate, (r,z)
  ! t denotes the scaled tangent to the ray (previously (rho, zeta))
  ! c * t would be the unit tangent

  TYPE( ray2DPt )    :: ray0, ray1, ray2
  REAL (KIND=8 ), INTENT( IN ) :: Topx( 2 ), Topn( 2 ), Botx( 2 ), Botn( 2 )

  REAL (KIND=8), INTENT( IN  )      :: rTopseg( 2 ), rBotseg( 2 )
  INTEGER,           INTENT(INOUT) :: iSegz
  INTEGER,           INTENT(INOUT) :: iSegr
  INTEGER,           INTENT(INOUT) :: iSegx
  INTEGER,           INTENT(INOUT) :: iSegy

  INTEGER            :: iSegz0, iSegr0
  REAL     (KIND=8 ) :: gradc0( 2 ), gradc1( 2 ), gradc2( 2 ), &
       c0, cimag0, crr0, crz0, czz0, csq0, cnn0_csq0, &
       c1, cimag1, crr1, crz1, czz1, csq1, cnn1_csq1, &
       c2, cimag2, crr2, crz2, czz2, urayt0( 2 ), urayt1( 2 ), &
       h, halfh, hw0, hw1, ray2n( 2 ), RM, RN, gradcjump( 2 ), cnjump, csjump, w0, w1, rho 



  ! The numerical integrator used here is a version of the polygon (a.k.a. midpoint, leapfrog, or Box method), and similar
  ! to the Heun (second order Runge-Kutta method).
  ! However, it's modified to allow for a dynamic step change, while preserving the second-order accuracy).

  ! *** Phase 1 (an Euler step)

  CALL EvaluateSSP( ray0%x, c0, cimag0, gradc0, crr0, crz0, czz0, rho, Freq, 'TAB', iSegz, iSegr, iSegx, iSegy )

  csq0      = c0 * c0
  cnn0_csq0 = crr0 * ray0%t( 2 )**2 - 2.0 * crz0 * ray0%t( 1 ) * ray0%t( 2 ) + czz0 * ray0%t( 1 )**2
  iSegz0    = iSegz     ! make note of current layer
  iSegr0    = iSegr

  h = Beam%deltas       ! initially set the step h, to the basic one, deltas
  urayt0 = c0 * ray0%t  ! unit tangent

  CALL ReduceStep2D( ray0%x, urayt0, iSegz0, iSegr0, Topx, Topn, Botx, Botn, Beam%deltas, h &
			, rTopSeg, rBotSeg ) ! reduce h to land on boundary
  halfh = 0.5 * h   ! first step of the modified polygon method is a half step

  ray1%x = ray0%x + halfh * urayt0
  ray1%t = ray0%t - halfh * gradc0 / csq0
  ray1%p = ray0%p - halfh * cnn0_csq0 * ray0%q
  ray1%q = ray0%q + halfh * c0        * ray0%p

  ! *** Phase 2

  CALL EvaluateSSP( ray1%x, c1, cimag1, gradc1, crr1, crz1, czz1, rho, Freq, 'TAB', iSegz, iSegr, iSegx, iSegy )
  csq1      = c1 * c1
  cnn1_csq1 = crr1 * ray1%t( 2 )**2 - 2.0 * crz1 * ray1%t( 1 ) * ray1%t( 2 ) + czz1 * ray1%t( 1 )**2

  ! The Munk test case with a horizontally launched ray caused problems.
  ! The ray vertexes on an interface and can ping-pong around that interface.
  ! Have to be careful in that case about big changes to the stepsize (that invalidate the leap-frog scheme) in phase II.
  ! A modified Heun or Box method could also work.

  urayt1 = c1 * ray1%t   ! unit tangent

  CALL ReduceStep2D( ray0%x, urayt1, iSegz0, iSegr0, Topx, Topn, Botx, Botn, Beam%deltas, h &
			, rTopSeg, rBotSeg ) ! reduce h to land on boundary

  ! use blend of f' based on proportion of a full step used.
  w1  = h / ( 2.0d0 * halfh )
  w0  = 1.0d0 - w1
  hw0 = h * w0
  hw1 = h * w1

  ray2%x   = ray0%x   + hw0 * urayt0              + hw1 * urayt1
  ray2%t   = ray0%t   - hw0 * gradc0 / csq0       - hw1 * gradc1 / csq1
  ray2%p   = ray0%p   - hw0 * cnn0_csq0 * ray0%q  - hw1 * cnn1_csq1 * ray1%q
  ray2%q   = ray0%q   + hw0 * c0        * ray0%p  + hw1 * c1        * ray1%p
  ray2%tau = ray0%tau + hw0 / CMPLX( c0, cimag0, KIND=8 ) + hw1 / CMPLX( c1, cimag1, KIND=8 )

  ray2%Amp       = ray0%Amp
  ray2%Phase     = ray0%Phase
  ray2%NumTopBnc = ray0%NumTopBnc
  ray2%NumBotBnc = ray0%NumBotBnc

  ! If we crossed an interface, apply jump condition

  CALL EvaluateSSP( ray2%x, c2, cimag2, gradc2, crr2, crz2, czz2, rho, Freq, 'TAB', iSegz, iSegr, iSegx, iSegy )
  ray2%c = c2

  IF ( iSegz /= iSegz0 .OR. iSegr /= iSegr0 ) THEN
     gradcjump =  gradc2 - gradc0
     ray2n     = [ -ray2%t( 2 ), ray2%t( 1 ) ]   ! ray normal

     cnjump    = DOT_PRODUCT( gradcjump, ray2n  )
     csjump    = DOT_PRODUCT( gradcjump, ray2%t )

     IF ( iSegz /= iSegz0 ) THEN         ! crossing in depth
        RM = +ray2%t( 1 ) / ray2%t( 2 )  ! this is tan( alpha ) where alpha is the angle of incidence
     ELSE                                ! crossing in range
        RM = -ray2%t( 2 ) / ray2%t( 1 )  ! this is tan( alpha ) where alpha is the angle of incidence
     END IF

     RN     = RM * ( 2 * cnjump - RM * csjump ) / c2
     ray2%p = ray2%p - ray2%q * RN

  END IF

END SUBROUTINE Step2D

! **********************************************************************!

SUBROUTINE ReduceStep2D( x0, urayt, iSegz0, iSegr0, Topx, Topn, Botx, Botn, deltas, h &
			, rTopSeg, rBotSeg )

  ! calculate a reduced step size, h, that lands on any points where the environment changes

  INTEGER,       INTENT( IN  ) :: iSegz0, iSegr0                             ! SSP layer the ray is in
  REAL (KIND=8), INTENT( IN  ) :: x0( 2 ), urayt( 2 )                        ! ray coordinate and tangent
  REAL (KIND=8), INTENT( IN  ) :: Topx( 2 ), Topn( 2 ), Botx( 2 ), Botn( 2 ) ! Top, bottom coordinate and normal
  REAL (KIND=8), INTENT( IN  ) :: deltas                                     ! default step size

  REAL (KIND=8), INTENT( IN  )      :: rTopseg( 2 ), rBotseg( 2 )

  REAL (KIND=8), INTENT( OUT ) :: h                                          ! reduced step size 
  REAL (KIND=8)                :: x( 2 ), d( 2 ), d0( 2 ), h1, h2, h3, h4, h5, h10

  ! Detect interface or boundary crossing and reduce step, if necessary, to land on that crossing.
  ! Keep in mind possibility that user put source right on an interface
  ! and that multiple events can occur (crossing interface, top, and bottom in a single step).

  x = x0 + h * urayt ! make a trial step

  ! interface crossing in depth
  h1 = huge( h1 )
  IF ( ABS( urayt( 2 ) ) > EPSILON( h1 ) ) THEN
     IF      ( SSP%z( iSegz0     ) > x(  2 ) ) THEN
        h1 = ( SSP%z( iSegz0     ) - x0( 2 ) ) / urayt( 2 )
     ELSE IF ( SSP%z( iSegz0 + 1 ) < x(  2 ) ) THEN
        h1 = ( SSP%z( iSegz0 + 1 ) - x0( 2 ) ) / urayt( 2 )
     END IF
  END IF

  ! top crossing
  h2 = huge( h2 )
  d  = x - Topx              ! vector from top to ray
  IF ( DOT_PRODUCT( Topn, d ) > EPSILON( h2 ) ) THEN
     d0  = x0 - Topx         ! vector from top    node to ray origin
     h2 = -DOT_PRODUCT( d0, Topn ) / DOT_PRODUCT( urayt, Topn )
  END IF

  ! bottom crossing
  h3 = huge( h3 )
  d  = x - Botx              ! vector from bottom to ray
  IF ( DOT_PRODUCT( Botn, d ) > EPSILON( h2 ) ) THEN
     d0  = x0 - Botx         ! vector from bottom node to ray origin
     h3 = -DOT_PRODUCT( d0, Botn ) / DOT_PRODUCT( urayt, Botn )
  END IF

  !!!! following tests could all be merged based on a single min or max range
  
  ! top segment crossing in range
  h4 = huge( h4 )
  IF ( ABS( urayt( 1 ) )  > EPSILON( h4 ) ) THEN
     IF       ( x(  1 ) < rTopSeg( 1 ) ) THEN
        h4 = -( x0( 1 ) - rTopSeg( 1 ) ) / urayt( 1 )
     ELSE IF  ( x(  1 ) > rTopSeg( 2 ) ) THEN
        h4 = -( x0( 1 ) - rTopSeg( 2 ) ) / urayt( 1 )
     END IF
  END IF

  ! bottom segment crossing in range
  h5 = huge( h5 )
  IF ( ABS( urayt( 1 ) )  > EPSILON( h5 ) ) THEN
     IF       ( x(  1 ) < rBotSeg( 1 ) ) THEN
        h5 = -( x0( 1 ) - rBotSeg( 1 ) ) / urayt( 1 )
     ELSE IF  ( x(  1 ) > rBotSeg( 2 ) ) THEN
        h5 = -( x0( 1 ) - rBotSeg( 2 ) ) / urayt( 1 )
     END IF
  END IF

  ! ocean segment crossing in r
  h10 = huge( h10 )

  IF ( SSP%Type == 'Q' ) THEN
     IF ( ABS( urayt( 1 ) ) > EPSILON( h1 ) ) THEN
        IF        ( x(  1 ) < SSP%Seg%r( iSegr0     ) ) THEN
           h10 = -( x0( 1 ) - SSP%Seg%r( iSegr0     ) ) / urayt( 1 )
           ! write( *, * ) 'ocean segment crossing in r'

        ELSE IF   ( x(  1 ) > SSP%Seg%r( iSegr0 + 1 ) ) THEN
           h10 = -( x0( 1 ) - SSP%Seg%r( iSegr0 + 1 ) ) / urayt( 1 )
           ! write( *, * ) 'ocean segment crossing in r'

        END IF
     END IF

  END IF

  h = MIN( h, h1, h2, h3, h4, h5, h10 )  ! take limit set by shortest distance to a crossing

  IF ( h < 1.0d-4 * deltas ) THEN   ! is it taking an infinitesimal step?
     h = 1.0d-5 * deltas            ! make sure we make some motion
     iSmallStepCtr = iSmallStepCtr + 1   ! keep a count of the number of sequential small steps
  ELSE
     iSmallStepCtr = 0   ! didn't do a small step so reset the counter
  END IF
  
END SUBROUTINE ReduceStep2D

! **********************************************************************!

SUBROUTINE Reflect2D( is, HS, BotTop, tBdry, nBdry, kappa, RefC, Npts, ray2D, iSegz, iSegr, iSegx, iSegy )

  INTEGER,              INTENT( IN ) :: Npts
  REAL     (KIND=8),    INTENT( IN ) :: tBdry( 2 ), nBdry( 2 )  ! Tangent and normal to the boundary
  REAL     (KIND=8),    INTENT( IN ) :: kappa                   ! Boundary curvature
  CHARACTER (LEN=3),    INTENT( IN ) :: BotTop                  ! Flag indicating bottom or top reflection
  TYPE( HSInfo ),       INTENT( IN ) :: HS                      ! half-space properties
  TYPE(ReflectionCoef), INTENT( IN ) :: RefC( NPts )            ! reflection coefficient
  INTEGER,              INTENT( INOUT ) :: is

  TYPE( ray2DPt ), INTENT(INOUT)    :: ray2D( MaxN )
  INTEGER,           INTENT(INOUT) :: iSegz
  INTEGER,           INTENT(INOUT) :: iSegr
  INTEGER,           INTENT(INOUT) :: iSegx
  INTEGER,           INTENT(INOUT) :: iSegy

  INTEGER           :: is1
  REAL     (KIND=8) :: c, cimag, gradc( 2 ), crr, crz, czz, rho       ! derivatives of sound speed
  REAL     (KIND=8) :: RM, RN, Tg, Th, rayt( 2 ), rayn( 2 ), rayt_tilde( 2 ), rayn_tilde( 2 ), cnjump, csjump  ! for curvature change
  REAL     (KIND=8) :: ck, co, si, cco, ssi, pdelta, rddelta, sddelta, theta_bot ! for beam shift
  COMPLEX  (KIND=8) :: gamma1, gamma2, gamma1Sq, gamma2Sq, GK, Refl   ! for tabulated reflection coef.
  COMPLEX  (KIND=8) :: ch, a, b, d, sb, delta, ddelta                 ! for beam shift
  TYPE(ReflectionCoef) :: RInt

  is  = is + 1
  is1 = is + 1

  Tg = DOT_PRODUCT( ray2D( is )%t, tBdry )  ! component of ray tangent, along boundary
  Th = DOT_PRODUCT( ray2D( is )%t, nBdry )  ! component of ray tangent, normal to boundary

  ray2D( is1 )%NumTopBnc = ray2D( is )%NumTopBnc
  ray2D( is1 )%NumBotBnc = ray2D( is )%NumBotBnc
  ray2D( is1 )%x         = ray2D( is )%x
  ray2D( is1 )%t         = ray2D( is )%t - 2.0 * Th * nBdry  ! changing the ray direction

  ! Calculate the change in curvature
  ! Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

  CALL EvaluateSSP( ray2D( is )%x, c, cimag, gradc, crr, crz, czz, rho, Freq, 'TAB', iSegz, iSegr, iSegx, iSegy )   ! just to get c

  ! incident unit ray tangent and normal
  rayt = c * ray2D( is )%t                              ! unit tangent to ray
  rayn = [ -rayt( 2 ), rayt( 1 ) ]                      ! unit normal  to ray

  ! reflected unit ray tangent and normal (the reflected tangent, normal system has a different orientation)
  rayt_tilde = c * ray2D( is1 )%t                       ! unit tangent to ray
  rayn_tilde = -[ -rayt_tilde( 2 ), rayt_tilde( 1 ) ]   ! unit normal  to ray

  RN = 2 * kappa / c ** 2 / Th    ! boundary curvature correction

  ! get the jumps (this could be simplified, e.g. jump in rayt is roughly 2 * Th * nbdry
  cnjump = -DOT_PRODUCT( gradc, rayn_tilde - rayn  )
  csjump = -DOT_PRODUCT( gradc, rayt_tilde - rayt )

  IF ( BotTop == 'TOP' ) THEN
     cnjump = -cnjump   ! this is because the (t,n) system of the top boundary has a different sense to the bottom boundary
     RN     = -RN
  END IF

  RM = Tg / Th   ! this is tan( alpha ) where alpha is the angle of incidence
  RN = RN + RM * ( 2 * cnjump - RM * csjump ) / c ** 2

  SELECT CASE ( Beam%Type( 3 : 3 ) )
  CASE ( 'D' )
     RN = 2.0 * RN
  CASE ( 'Z' )
     RN = 0.0
  END SELECT

  ray2D( is1 )%c   = c
  ray2D( is1 )%tau = ray2D( is )%tau
  ray2D( is1 )%p   = ray2D( is )%p + ray2D( is )%q * RN
  ray2D( is1 )%q   = ray2D( is )%q

  ! account for phase change

  SELECT CASE ( HS%BC )
  CASE ( 'R' )                 ! rigid
     ray2D( is1 )%Amp   = ray2D( is )%Amp
     ray2D( is1 )%Phase = ray2D( is )%Phase
  CASE ( 'V' )                 ! vacuum
     ray2D( is1 )%Amp   = ray2D( is )%Amp
     ray2D( is1 )%Phase = ray2D( is )%Phase + pi
  CASE ( 'F' )                 ! file
     RInt%theta = RadDeg * ABS( ATAN2( Th, Tg ) )           ! angle of incidence (relative to normal to bathymetry)
     IF ( RInt%theta > 90 ) RInt%theta = 180. - RInt%theta  ! reflection coefficient is symmetric about 90 degrees
     CALL InterpolateReflectionCoefficient( RInt, RefC, Npts, PRTFile )
     ray2D( is1 )%Amp   = ray2D( is )%Amp * RInt%R
     ray2D( is1 )%Phase = ray2D( is )%Phase + RInt%phi
  CASE ( 'A', 'G' )            ! half-space
     GK       = omega * Tg     ! wavenumber in direction parallel to bathymetry
     gamma1Sq = ( omega / c     ) ** 2 - GK ** 2 - i * tiny( omega )   ! tiny prevents g95 giving -zero, and wrong branch cut
     gamma2Sq = ( omega / HS%cP ) ** 2 - GK ** 2 - i * tiny( omega )
     gamma1   = SQRT( -gamma1Sq )
     gamma2   = SQRT( -gamma2Sq )

     Refl = ( HS%rho * gamma1 - rho * gamma2 ) / ( HS%rho * gamma1 + rho * gamma2 )

     IF ( ABS( Refl ) < 1.0E-5 ) THEN   ! kill a ray that has lost its energy in reflection
        ray2D( is1 )%Amp   = 0.0
        ray2D( is1 )%Phase = ray2D( is )%Phase
     ELSE
        ray2D( is1 )%Amp   = ABS( Refl ) * ray2D(  is )%Amp
        ray2D( is1 )%Phase = ray2D( is )%Phase + ATAN2( AIMAG( Refl ), REAL( Refl ) )

        ! compute beam-displacement Tindle, Eq. (14)
        ! needs a correction to beam-width as well ...
        !  IF ( REAL( gamma2Sq ) < 0.0 ) THEN
        !     rhoW   = 1.0   ! density of water
        !     rhoWSq  = rhoW  * rhoW
        !     rhoHSSq = rhoHS * rhoHS
        !     DELTA = 2 * GK * rhoW * rhoHS * ( gamma1Sq - gamma2Sq ) /
        ! &( gamma1 * i * gamma2 *
        ! &( -rhoWSq * gamma2Sq + rhoHSSq * gamma1Sq ) )
        !     RV( is + 1 ) = RV( is + 1 ) + DELTA
        !  END IF

        if ( Beam%Type( 4 : 4 ) == 'S' ) then   ! beam displacement & width change (Seongil's version)

           ch = ray2D( is )%c / conjg( HS%cP )
           co = ray2D( is )%t( 1 ) * ray2D( is )%c
           si = ray2D( is )%t( 2 ) * ray2D( is )%c
           ck = omega / ray2D( is )%c

           a   = 2 * HS%rho * ( 1 - ch * ch )
           b   = co * co - ch * ch
           d   = HS%rho * HS%rho * si * si + b
           sb  = sqrt( b )
           cco = co * co
           ssi = si * si

           delta   = a * co / si / ( ck * sb * d )    
           pdelta  = real( delta ) / ( ray2D( is )%c / co)
           ddelta  = -a / ( ck*sb*d ) - a*cco / ssi / (ck*sb*d) + a*cco / (ck*b*sb*d) &
                     -a*co / si / (ck*sb*d*d) * (2* HS%rho * HS%rho *si*co-2*co*si)
           rddelta = -real( ddelta )
           sddelta = rddelta / abs( rddelta )        

           ! next 3 lines have an update by Diana McCammon to allow a sloping bottom
           theta_bot = datan( tBdry( 2 ) / tBdry( 1 ))  ! bottom angle
           ray2D( is1 )%x( 1 ) = ray2D( is1 )%x( 1 ) + real( delta ) * dcos( theta_bot )   ! range displacement
           ray2D( is1 )%x( 2 ) = ray2D( is1 )%x( 2 ) + real( delta ) * dsin( theta_bot )   ! depth displacement
           ray2D( is1 )%tau    = ray2D( is1 )%tau + pdelta             ! phase change
           ray2D( is1 )%q      = ray2D( is1 )%q + sddelta * rddelta * si * c * ray2D( is )%p   ! beam-width change
        endif

     ENDIF
  CASE DEFAULT
     WRITE( PRTFile, * ) 'HS%BC = ', HS%BC
     CALL ERROUT( PRTFile, 'F', 'Reflect2D', 'Unknown boundary condition type' )
  END SELECT

END SUBROUTINE Reflect2D

END PROGRAM DPRM3



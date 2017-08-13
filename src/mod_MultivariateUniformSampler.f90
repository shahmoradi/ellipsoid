module MultivariateUniformSampler
  
  use Parameters, only: RK

  implicit none

  !***************************************************************************
  !***************************************************************************
  ! These NLOPT variables are taken from nlopt.f
  integer, parameter :: NLOPT_GN_DIRECT = 0
  integer, parameter :: NLOPT_GN_DIRECT_L=1
  integer, parameter :: NLOPT_GN_DIRECT_L_RAND=2
  integer, parameter :: NLOPT_GN_DIRECT_NOSCAL=3
  integer, parameter :: NLOPT_GN_DIRECT_L_NOSCAL=4
  integer, parameter :: NLOPT_GN_DIRECT_L_RAND_NOSCAL=5
  integer, parameter :: NLOPT_GN_ORIG_DIRECT=6
  integer, parameter :: NLOPT_GN_ORIG_DIRECT_L=7
  integer, parameter :: NLOPT_GD_STOGO=8
  integer, parameter :: NLOPT_GD_STOGO_RAND=9
  integer, parameter :: NLOPT_LD_LBFGS_NOCEDAL=10
  integer, parameter :: NLOPT_LD_LBFGS=11
  integer, parameter :: NLOPT_LN_PRAXIS=12
  integer, parameter :: NLOPT_LD_VAR1=13
  integer, parameter :: NLOPT_LD_VAR2=14
  integer, parameter :: NLOPT_LD_TNEWTON=15
  integer, parameter :: NLOPT_LD_TNEWTON_RESTART=16
  integer, parameter :: NLOPT_LD_TNEWTON_PRECOND=17
  integer, parameter :: NLOPT_LD_TNEWTON_PRECOND_RESTART=18
  integer, parameter :: NLOPT_GN_CRS2_LM=19
  integer, parameter :: NLOPT_GN_MLSL=20
  integer, parameter :: NLOPT_GD_MLSL=21
  integer, parameter :: NLOPT_GN_MLSL_LDS=22
  integer, parameter :: NLOPT_GD_MLSL_LDS=23
  integer, parameter :: NLOPT_LD_MMA=24
  integer, parameter :: NLOPT_LN_COBYLA=25
  integer, parameter :: NLOPT_LN_NEWUOA=26
  integer, parameter :: NLOPT_LN_NEWUOA_BOUND=27
  integer, parameter :: NLOPT_LN_NELDERMEAD=28
  integer, parameter :: NLOPT_LN_SBPLX=29
  integer, parameter :: NLOPT_LN_AUGLAG=30
  integer, parameter :: NLOPT_LD_AUGLAG=31
  integer, parameter :: NLOPT_LN_AUGLAG_EQ=32
  integer, parameter :: NLOPT_LD_AUGLAG_EQ=33
  integer, parameter :: NLOPT_LN_BOBYQA=34
  integer, parameter :: NLOPT_GN_ISRES=35
  integer, parameter :: NLOPT_AUGLAG=36
  integer, parameter :: NLOPT_AUGLAG_EQ=37
  integer, parameter :: NLOPT_G_MLSL=38
  integer, parameter :: NLOPT_G_MLSL_LDS=39
  integer, parameter :: NLOPT_LD_SLSQP=40
  integer, parameter :: NLOPT_LD_CCSAQ=41
  integer, parameter :: NLOPT_GN_ESCH=42
  integer, parameter :: NLOPT_FAILURE=-1
  integer, parameter :: NLOPT_INVALID_ARGS=-2
  integer, parameter :: NLOPT_OUT_OF_MEMORY=-3
  integer, parameter :: NLOPT_ROUNDOFF_LIMITED=-4
  integer, parameter :: NLOPT_FORCED_STOP=-5
  integer, parameter :: NLOPT_SUCCESS=1
  integer, parameter :: NLOPT_STOPVAL_REACHED=2
  integer, parameter :: NLOPT_FTOL_REACHED=3
  integer, parameter :: NLOPT_XTOL_REACHED=4
  integer, parameter :: NLOPT_MAXEVAL_REACHED=5
  integer, parameter :: NLOPT_MAXTIME_REACHED=6
  !***************************************************************************
  !***************************************************************************
  
  character(len=300)      :: directory, os, fileName
  integer , parameter     :: nd_glob=3
  real(RK), allocatable   :: CovOG(:,:), ScaleMat(:,:)
  real(RK)                :: radOG, volOG

contains


!****************************************************************************
!****************************************************************************

  subroutine test_getMVUfromEllipsoids(nd, np, ntest)
    use Statistics     , only: getRandInt,getRandIntLecuyer, getRandPointOnEllipsoid
    use Misc           , only: num2str, getSlash
    implicit none
    integer , parameter   :: ncMinGen=2,ncMaxGen=4,ndMin=2,ndMax=2,npMax=10,npb=560,npog=1000000
    integer , intent(in)  :: nd, np, ntest
    integer               :: i,ic,ip,ipb,itest,nc,npTot,ncMax, npbTot, istart, nRuns
    character(len=63)     :: ndFormat,ndMemberFormat,fileBase,ndMemberFormatOG,scaleFormat
    character(len=1)      :: slash
    character(len=127)    :: fileNameMembershipTrue, fileNameBoundTrue, fileNameCenterTrue, fileNameAMVETrue, fileNameMembershipOG
    character(len=127)    :: fileNameScale

    integer , allocatable :: Membership(:),MembershipOG(:)            ! true cluster memberships
    integer , allocatable :: CMembership(:)                           ! cluster memberships according to kvolume clustering
    real(RK), allocatable :: Point(:,:)                               ! the input data (points) to be clustered. on output, it is ordered.
    real(RK), allocatable :: Bound(:,:),CenterCrd(:,:)                ! real bourndaries from which points were sampled
    real(RK), allocatable :: EllVolTrue(:)                            ! volumes of the true clusters (sqrt(determinant))
    real(RK), allocatable :: BoundAMVE(:,:),CenterAMVE(:,:)
    real(RK), allocatable :: InvCovMat(:,:,:),InvCovMatOG(:,:,:)
    real(RK), allocatable :: PointOG(:,:),CenterOG(:,:),BoundOG(:,:),EllVolOG(:),ScaleMat(:,:)
    character(8)          :: date
    character(10)         :: time
    character(18)         :: dateNtime
    integer               :: idum      ! debugging
    integer               :: seedSize
    integer, allocatable  :: seed(:)

    call random_seed(size=seedSize)
    allocate(seed(seedSize))
    allocate(CovOG(nd,nd))
    idum = -777

    slash = getSlash()
    directory = '..'//slash//'out'//slash
    !call execute_command_line( 'mkdir ' // trim(adjustl(directory)) )
    call date_and_time(date,time)
    dateNtime = trim(adjustl(date)) // trim(adjustl(time))
    nRuns = npMax
    fileBase = 'allUniform'

    ! get one million points from unity sphere
    allocate(MembershipOG(npog),PointOG(nd,npog),BoundOG(nd,npb),CenterOG(nd,1),EllVolOG(1),ScaleMat(6,ntest),InvCovMatOG(nd,nd,1))
    call getMVUfromEllipsoids(1,nd,npog,npb,trim(adjustl(fileBase)),EllVolOG,CenterOG,BoundOG,PointOG,MembershipOG,InvCovMatOG)
    fileNameMembershipOG = '_'//'_'//trim(adjustl(fileBase))//'_'// &
                               'MembershipOG'//'.txt'
    !open(unit=25,file=trim(adjustl(directory))//trim(adjustl(fileNameMembershipOG)),status='replace')
    !ndMemberFormatOG = '(' // num2str(nd) // 'F30.15,1I30)'
    !do ip=1,npog
    !  write(25,ndMemberFormatOG) PointOG(1:nd,ip), MembershipOG(ip)
    !end do


    do itest = 1,ntest                                 ! do [ntest] times

      seed = -getRandIntLecuyer(1,10000,idum)
      call random_seed(put=seed)

      !write(*,*) 'itest: ', itest

      npTot = np                                    ! total number of points
      ncMax = npTot/(nd+1)
      nc = 1                                        ! number of elipsoids generated
      npbTot = npb*ncMax                            ! total points in each bounding ellipse

      allocate( Membership(npTot) , Point(nd,npTot) , Bound(nd,nc*npb) , CenterCrd(nd,nc) , EllVolTrue(nc) )
      allocate( BoundAMVE(nd,nc*npb) , CenterAMVE(nd,nc) , InvCovMat(nd,nd,nc) )

      fileBase = 'allUniform'
      ! Generate some random clusters
      call getMVUfromEllipsoids(nc,nd,npTot,npb,trim(adjustl(fileBase)),EllVolTrue,CenterCrd,Bound,Point,Membership,InvCovMat)
      call getSamAMVE(nc,nd,npTot,npb,fileBase,CenterAMVE,BoundAMVE,Point,PointOG,InvCovMatOG,CenterOG,ScaleMat(:,itest))! generate AMVE from point sample

      ! file naming
      ndFormat = '(' // num2str(nd) // 'F30.15)'
      ndMemberFormat = '(' // num2str(nd) // 'F30.15,1I30)'
      !fileNameMembershipTrue = '_'//num2str(itest)//'_'//trim(adjustl(fileBase))//'_'// &
      !                         'MembershipTrue'//'_'//num2str(nc)//'_'//num2str(nd)//'_'//num2str(npTot)//'.txt'
      !fileNameBoundTrue = '_'//num2str(itest)//'_'//trim(adjustl(fileBase))//'_'// &
      !                    'BoundTrue'//'_'//num2str(nc)//'_'//num2str(nd)//'_'//num2str(npTot)//'.txt'
      !fileNameCenterTrue = '_'//num2str(itest)//'_'//trim(adjustl(fileBase))//'_'// &
      !                     'CenterTrue'//'_'//num2str(nc)//'_'//num2str(nd)//'_'//num2str(npTot)//'.txt'
      !fileNameAMVETrue = '_'//num2str(itest)//'_'//trim(adjustl(fileBase))//'_'// &
      !                     'AMVETrue'//'_'//num2str(nc)//'_'//num2str(nd)//'_'//num2str(npTot)//'.txt'
      !! write data to file
      !open(unit=21,file=trim(adjustl(directory))//trim(adjustl(fileNameMembershipTrue)),status='replace')
      !open(unit=22,file=trim(adjustl(directory))//trim(adjustl(fileNameBoundTrue)),status='replace')
      !open(unit=23,file=trim(adjustl(directory))//trim(adjustl(fileNameCenterTrue)),status='replace')
      !open(unit=24,file=trim(adjustl(directory))//trim(adjustl(fileNameAMVETrue)),status='replace')
      !! write point data
      !do ip=1,npTot
      !  write(21,ndMemberFormat) Point(1:nd,ip), Membership(ip)
      !end do
      !! write bound data
      !do ip=1,nc*npb
      !  write(22,ndFormat) Bound(1:nd,ip)
      !end do
      !! write center data
      !do ic=1,nc
      !  write(23,ndFormat) CenterCrd(1:nd,ic)
      !end do
      !! write AMVE bound data
      !do ip=1,nc*npb
      !  write(24,ndFormat) BoundAMVE(1:nd,ip)
      !end do
      !close(21); close(22); close(23); close(24);

      deallocate(Point,Membership,Bound,CenterCrd,EllVolTrue,BoundAMVE,CenterAMVE,InvCovMat)
      !deallocate(EllVolOG,CenterOG,BoundOG,PointOG,MembershipOG)
      !nRuns = nRuns + 10

    end do

    fileNameScale = trim(adjustl(fileBase))//'_'// 'ScaleData'//'_nd_'//num2str(nd)//'_np_'//num2str(np)//'.txt'
    open(unit=26,file=trim(adjustl(directory))//trim(adjustl(fileNameScale)),status='replace')
    scaleFormat = '(' // num2str(6) // 'F30.15)'
    do itest=1,ntest
      write(26,scaleFormat) ScaleMat(:,itest)
    end do

    deallocate( CovOG, ScaleMat )

  end subroutine test_getMVUfromEllipsoids


!****************************************************************************
!****************************************************************************

  subroutine getMVUfromEllipsoids(nc,nd,npTot,npb,fileBase,EllVol,CenterCrd,Bound,Point,Membership,InvCovMat)
    use Matrix    , only: getCholeskyFactor, getInvMatFromCholFac
    use Statistics, only: getRandMVU, getRandInt, getCovMatFromCorMat &
                        , getRandPointOnEllipsoid, getRandCorMatRejection, getMahalSq
    use SpecialFunctions, only: getEllVolCoef
    implicit none
    real(RK), parameter          :: maxRange=1.d1,PI=4*atan(1.0)
    integer , intent(in)         :: nc,nd,npb            ! number of points for each bounding ellipsoid
    integer , intent(in)         :: npTot                ! number of points in all clusters together
    integer , intent(out)        :: Membership(npTot)
    real(RK), intent(out)        :: CenterCrd(nd,nc), Bound(nd,npb*nc), Point(nd,npTot), EllVol(nc),InvCovMat(nd,nd,nc)
    real(RK)                     :: SigmaVec(nd,nc),Diagonal(nd,nc),EllVolCDF(nc),CSize(nc),CorMat(nd,nd)
    real(RK)                     :: CovMat(nd,nd,nc)
    real(RK)                     :: dummy, mahalSq(1), rad
    integer                      :: ic, id, ip, ipb, istart, icc, jointCounter, k
    character(len=*), intent(in) :: fileBase
    real(RK)                     :: myStd

    ! assign sphere center at origin
    do ic = 1,nc
      do id = 1,nd
        CenterCrd(id,ic) = 0.0_RK
      end do
    end do
    
    ! get radius of unity sphere in nd dimensions
    if (mod(nd, 2) == 0) then
      k = nd / 2
      rad = (factorial(k))**(1/(nd)) / sqrt(PI)
    else
      k = nd / 2
      rad = (factorial(factorial(2*k+1)) / (2**(k+1) * PI**k))**(1/(2*k+1))
    end if


    ! assign cluster standard deviations
    SigmaVec = rad**(1/nd)
    radOG = rad**(1/nd)
    !SigmaVec = getEllVolCoef(nd)

    ! assign cluster covmat, and choleskky factors
    do ic = 1,nc
      CorMat = getRandCorMatRejection(nd,-0.0d0,0.0d0)
      CovMat(:,:,ic) = getCovMatFromCorMat(nd,SigmaVec(:,ic),CorMat)
      call getCholeskyFactor(nd,CovMat(:,:,ic),Diagonal(:,ic))
      InvCovMat(:,:,ic) = getInvMatFromCholFac(nd,CovMat(:,:,ic),Diagonal(:,ic))
      EllVol(ic) = product(Diagonal(:,ic))
    end do
    
    CovOG(:,:) = InvCovMat(:,:,1)

    do ic = 1,nc
      EllVolCDF(ic) = sum(EllVol(1:ic)) / sum(EllVol)
    end do
    volOG = EllVol(1)

    ! assign points "randomly equally" to each cluster
    !if (trim(adjustl(fileBase))=='mvuEll') then
    !  do ip = 1,npTot
    !    !ic = getRandIntLecuyer(1,nc,idum)
    !    ic = getRandInt(1,nc)
    !    Membership(ip) = ic
    !    Point(:,ip) = getRandMVU(nd,CenterCrd(:,ic),CovMat(:,:,ic),Diagonal(:,ic))
    !  end do
    !elseif (trim(adjustl(fileBase))=='mvuEllW') then
    !  CSize = 0
    !  do ip = 1,npTot
    !    call random_number(dummy)
    !    ic = minloc(EllVolCDF,1,mask=EllVolCDF>dummy)
    !    !ic = minloc(EllVolCDF,1,mask=EllVolCDF>getRandRealLecuyer(idum))
    !    CSize(ic) = CSize(ic) + 1
    !    Membership(ip) = ic
    !    Point(:,ip) = getRandMVU(nd,CenterCrd(:,ic),CovMat(:,:,ic),Diagonal(:,ic))
    !  end do
    if (trim(adjustl(fileBase))=='allUniform') then
      CSize = 0
      do ip = 1,npTot
        do  ! check if the point is not common between multiple clusters, is so, reduce the probability of accepting the point correctly.
          call random_number(dummy)
          ic = minloc(EllVolCDF,1,mask=EllVolCDF>dummy)
          !ic = minloc(EllVolCDF,1,mask=EllVolCDF>getRandRealLecuyer(idum))
          Point(:,ip) = getRandMVU(nd,CenterCrd(:,ic),CovMat(:,:,ic),Diagonal(:,ic))
          jointCounter = 0
          do icc = 1, nc
            mahalSq = getMahalSq(nd,1,CenterCrd(:,icc),InvCovMat(:,:,icc),Point(:,ip))
            if ( mahalSq(1) <= 1.d0 ) jointCounter = jointCounter + 1
          end do
          call random_number(dummy)
          if ( dummy > 1.d0/jointCounter ) cycle
          !if ( getRandRealLecuyer(idum) > 1.d0/jointCounter ) cycle
          CSize(ic) = CSize(ic) + 1
          Membership(ip) = ic
          exit
        end do
      end do
    else
      write(*,*) 'Requested sampling style not supported.'
      stop
    end if

    do ic = 1,nc
      istart = (ic-1)*npb
      !istart = (ic-1)*npb + ic
      !Bound(:,istart) = CenterCrd(:,ic)
      do ipb = istart+1,istart+npb
        Bound(:,ipb) = getRandPointOnEllipsoid(nd,CenterCrd(:,ic),CovMat(:,:,ic),Diagonal(:,ic))
      end do
    end do

  end subroutine getMVUfromEllipsoids

!****************************************************************************
!****************************************************************************

  subroutine getSamAMVE(nc,nd,npTot,npb,fileBase,CenterAMVE,BoundAMVE,Point,PointOG,InvCovMatOG,CenterOG,ScaleVec)
    use Matrix    , only: getCholeskyFactor, getInvMatFromCholFac, getInvPosDefMat 
    use Statistics, only: getRandMVU, getSamCholFac, getCovMatFromCorMat, getSamCovMean &
                        , getRandPointOnEllipsoid, getRandCorMatRejection, getMahalSq, getMean

    implicit none
    real(RK), parameter          :: maxRange=1.d1
    integer , intent(in)         :: nc,nd,npb            ! number of points for each bounding ellipsoid
    integer , intent(in)         :: npTot                ! number of points in all clusters together
    real(RK), intent(in)         :: Point(nd,npTot)      ! points sampled from original ellipse
    real(RK), intent(in)         :: PointOG(2,1000000)
    real(RK), intent(out)        :: BoundAMVE(nd,npb*nc),CenterAMVE(nd,nc)
    real(RK)                     :: SigmaVec(nd,nc),Diagonal(nd,nc),EllVolCDF(nc),CSize(nc),CorMat(nd,nd)
    real(RK)                     :: CovMat(nd,nd,nc),InvCovMat(nd,nd,nc),Mean(nd),counter
    real(RK)                     :: CovMat0(nd,nd,nc),Diagonal0(nd,nc),PointAMVE(nd,1000000)
    real(RK), intent(in)         :: InvCovMatOG(nd,nd,nc),CenterOG(nd,1)
    real(RK)                     :: dummy, MahalSqMax(1), sqrtDetInvCovMat, mahalSq(1), maxMahal
    real(RK)                     :: scaleFac, Xnf(nd), F_data(nd+1,nd), unscaledVol, scaledVol
    real(RK), intent(inout)      :: ScaleVec(6,1)
    integer                      :: ic, id, ip, ipb, istart, icc
    integer                      :: ires, maxeval
    real(RK)                     :: x(nd), maxf, lb(nd), ub(nd)
    integer*8                    :: opt
    character(len=*), intent(in) :: fileBase
    
    ic = 1
    counter = 0
    scaleFac = 0
    Mean = getMean(nd,npTot,Point)
    do id = 1,nd
      CenterAMVE(id,1) = Mean(id)
    end do
    call getSamCholFac(nd,npTot,Mean,Point,CovMat,Diagonal)
    InvCovMat(:,:,ic) = getInvMatFromCholFac(nd,CovMat(:,:,ic),Diagonal(:,ic))
    MahalSqMax(1) = maxval(getMahalSq(nd,npTot,Mean,InvCovMat,Point))
    CovMat = sqrt(MahalSqMax(1)) * CovMat
    Diagonal = sqrt(MahalSqMax(1))* Diagonal
    InvCovMat(:,:,ic) = getInvMatFromCholFac(nd,CovMat(:,:,ic),Diagonal(:,ic))
    unscaledVol = product(Diagonal(:,ic))

    
    !do ip = 1,1000000
    !  ic = 1
    !  PointAMVE(:,ip) = getRandMVU(nd,CenterAMVE(:,ic),CovMat(:,:,ic),Diagonal(:,ic))
    !end do
    !counter = 0
    !do ip = 1,1000000
    !  ic = 1
    !  mahalSq = getMahalSq(nd,1,CenterOG(:,ic),InvCovMatOG(:,:,ic),PointAMVE(:,ip))
    !  if ( sqrt(mahalSq(1)) <= 1.d0 ) counter = counter + 1
    !end do
    !ScaleVec(3,1) = 1 - (counter / 1000000)
    !write(*,*) 'Blank Space: ', ScaleVec(3,1)
    
    F_data(1:nd,1:nd) = InvCovMat(:,:,1)
    F_data(nd+1,:) = Mean
    
    ! global optimization with NLOPT's DIRECT algoritm
    call nlo_create(opt, NLOPT_GN_ORIG_DIRECT, nd)
    call nlo_set_max_objective(ires, opt, myfunc, F_data)
    call nlo_get_lower_bounds(ires, opt, lb)
    call nlo_get_upper_bounds(ires,opt, ub)
    lb = -1.0_RK*radOG
    ub = radOG
    call nlo_set_lower_bounds(ires, opt, lb)
    call nlo_set_upper_bounds(ires, opt, ub)
    call nlo_add_inequality_constraint(ires, opt, myconstraint1, CovOG, 1.D-6)
    call nlo_set_xtol_rel(ires, opt, 1.D-4)
    maxeval = 50000
    call nlo_set_maxeval(ires, opt, maxeval)
    call nlo_get_maxeval(maxeval, opt)
    x = 0.0_RK
    call nlo_optimize(ires, opt, x, maxf)
    call nlo_destroy(opt)
    
    ! local optimization with NLOPT's COBYLA algorithm
    call nlo_create(opt, NLOPT_LN_COBYLA, nd)
    call nlo_set_max_objective(ires, opt, myfunc, F_data)
    call nlo_get_lower_bounds(ires, opt, lb)
    call nlo_get_upper_bounds(ires,opt,ub)
    lb = -1.0_RK*radOG
    ub = radOG
    call nlo_set_lower_bounds(ires, opt, lb)
    call nlo_set_upper_bounds(ires, opt, ub)
    call nlo_add_inequality_constraint(ires, opt, myconstraint1, CovOG, 1.D-14)
    call nlo_add_inequality_constraint(ires, opt, myconstraint2, CovOG, 1.D-14)
    call nlo_set_xtol_rel(ires, opt, 1.D-8)
    call nlo_set_maxeval(ires, opt, maxeval)
    call nlo_get_maxeval(maxeval, opt)
    call nlo_optimize(ires, opt, x, maxf)
    call nlo_destroy(opt)
    
    ! record results
    scaleFac = sqrt(maxf)
    !write(*,*) 'Scale Factor = ', scaleFac

    ! scale AMVE with value found from optimization
    CovMat = scaleFac * CovMat
    Diagonal = scaleFac * Diagonal
    scaledVol = product(Diagonal(:,ic))

    ScaleVec(1,1) = nd
    ScaleVec(2,1) = npTot
    ScaleVec(3,1) = volOG
    ScaleVec(4,1) = unscaledVol
    ScaleVec(5,1) = scaledVol
    ScaleVec(6,1) = scaleFac
    
    ! scale AMVE with value found from optimization
    CovMat = scaleFac * CovMat
    Diagonal = scaleFac * Diagonal

    do ic = 1,nc
      istart = (ic-1)*npb
      do ipb = istart+1,istart+npb
        BoundAMVE(:,ipb) = getRandPointOnEllipsoid(nd,Mean,CovMat(:,:,ic),Diagonal(:,ic))
      end do
    end do

  end subroutine getSamAMVE

!****************************************************************************
!****************************************************************************

  integer function factorial(n)
    implicit none
    integer , intent(in)         :: n
    integer                      :: i, res

    res = 1
    do i = 1,n
      res = res * i
    end do
    factorial = res
  end function factorial

!****************************************************************************
!****************************************************************************

  subroutine myfunc(val, n, x, grad, need_gradient, f_data)
    double precision val, x(n), grad(n), f_data(n+1,n), xc(1,n), cen(n), S(n,n)
    integer n, need_gradient, i
    
    cen = f_data(n+1,:)
    S = f_data(1:n,:)

    val = dot_product(x - cen, matmul(S, x - cen))

  end subroutine myfunc

!****************************************************************************
!****************************************************************************

  subroutine myconstraint1(val, n, x, grad, need_gradient, d)
    integer n, need_gradient
    double precision val, x(n), grad(n), d(n,n), xc(1,n)

    xc(1,:) = x
    val = norm2(matmul(matmul(xc,d),transpose(xc))) - 1.0

  end subroutine myconstraint1

!****************************************************************************
!****************************************************************************
  
  subroutine myconstraint2(val, n, x, grad, need_gradient, d)
    integer n, need_gradient
    double precision val, x(n), grad(n), d(n,n), xc(1,n)

    xc(1,:) = x
    val = -1.0*norm2(matmul(matmul(xc,d),transpose(xc))) + 1.0

  end subroutine myconstraint2



end module MultivariateUniformSampler

module Statistics

  use Parameters, only: RK
  implicit none

contains

!****************************************************************************
!****************************************************************************
  
  function getMean(nd,np,Point)
    implicit none
    integer, intent(in)  :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
    real(RK), intent(in) :: Point(nd,np)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
    real(RK)             :: getMean(nd)            ! output mean vector
    integer              :: ip
    getMean = 0._RK
    do ip = 1,np
      getMean = getMean + Point(1:nd,ip)
    end do
    getMean = getMean / np
  end function getMean
    
!****************************************************************************
!****************************************************************************
  
  function getNormData(nd,np,Mean,Point)
    implicit none
    integer , intent(in)  :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
    real(RK), intent(in)  :: Mean(nd)               ! Mean vector
    real(RK), intent(in)  :: Point(nd,np)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
    real(RK)              :: getNormData(np,nd)
    integer               :: i
    do i = 1,np
      getNormData(i,1:nd) = Point(1:nd,i) - Mean
    end do
  end function getNormData
    
!****************************************************************************
!****************************************************************************

  ! Returns the lower triangle Cholesky Factor of the covariance matrix of a set of points in the lower part of CholeskyLower.
  ! The upper part of CholeskyLower, including the diagonal elements of it, contains the Covariance matrix of the sample.
  ! Diagonal, contains on output, the diagonal terms of Cholesky Factor.
  subroutine getSamCholFac(nd,np,Mean,Point,CholeskyLower,Diagonal)
    use matrix, only: getCholeskyFactor
    implicit none
    integer , intent(in)  :: nd,np                  ! np is the number of observations, nd is the number of parameters for each observation
    real(RK), intent(in)  :: Mean(nd)               ! Mean vector
    real(RK), intent(in)  :: Point(nd,np)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
    real(RK), intent(out) :: CholeskyLower(nd,nd)   ! Lower Cholesky Factor of the covariance matrix
    real(RK), intent(out) :: Diagonal(nd)           ! Diagonal elements of the Cholesky factorization
    real(RK)              :: NormedData(np,nd)
    integer               :: i,j
    
    do i = 1,np
      NormedData(i,1:nd) = Point(1:nd,i) - Mean
    end do
    
    ! Only upper half of CovMat is needed
    do j = 1,nd
      do i = 1,j
        CholeskyLower(i,j) = dot_product( NormedData(1:np,i) , NormedData(1:np,j) ) / real(np-1,kind=RK) ! Only the upper half of CovMat is needed
      end do
    end do

    call getCholeskyFactor(nd,CholeskyLower,Diagonal)
    
  end subroutine getSamCholFac

!****************************************************************************
!****************************************************************************
  
  ! This subroutine finds the elements of the symmetric nd*nd sample covariance matrix of a given set of
  ! nd observations, each one with nd parameters in the format of a "" np*nd "" matrix.
  ! For a review refer to Geisser & Cornfield (1963) "Posterior distributions for multivariate normal parameters".
  ! Also refer to Box and Tiao (1973), "Bayesian Inference in Statistical Analysis" Page 421.
  ! Amir Shahmoradi, Oct 16, 2009, 11:14 AM, MTU
  subroutine getSamCovMean(np,nd,Point,CovMat,Mean,MahalSq,InvCovMat,sqrtDetInvCovMat)
    
    use matrix, only: getInvPosDefMatSqrtDet
    implicit none
    integer , intent(in)            :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
    real(RK), intent(in)            :: Point(np,nd)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
    real(RK), intent(out)           :: CovMat(nd,nd)          ! Covariance matrix of the input data
    real(RK), intent(out)           :: Mean(nd)               ! Mean vector
    real(RK), intent(out), optional :: MahalSq(np)            ! Vector of Mahalanobis Distances Squared, with respect to the sample's mean position
    real(RK), intent(out), optional :: InvCovMat(nd,nd)       ! Inverse Covariance matrix of the input data
    real(RK), intent(out), optional :: sqrtDetInvCovMat       ! sqrt determinant of the inverse covariance matrix
    real(RK), dimension(nd)         :: DummyVec
    real(RK), dimension(np,nd)      :: NormedData
    integer                         :: i,j,k
    
    do j = 1,nd
      Mean(j) = sum(Point(1:np,j)) / real(np,kind=RK)
      NormedData(1:np,j) = Point(1:np,j) - Mean(j)
    end do
    do i = 1,nd
        do j = 1,nd
            CovMat(i,j) = dot_product(NormedData(1:np,i),NormedData(1:np,j))/real(np-1,kind=RK)
        end do
    end do
    
    if (present(sqrtDetInvCovMat)) then   ! Calculate InCovMat, MahalSq, and sqrt Determinant of InCovMat
      do j = 1,nd
        do i = 1,j
          InvCovMat(i,j) = CovMat(i,j)    ! Only the upper half of CovMat is needed
        end do
      end do
      call getInvPosDefMatSqrtDet(nd,InvCovMat,sqrtDetInvCovMat)
      do i = 1,np
        do j = 1,nd
          DummyVec(j) = dot_product(InvCovMat(1:nd,j),NormedData(i,1:nd))
        end do
        MahalSq(i) = dot_product(NormedData(i,1:nd),DummyVec)
      end do
    end if
    
  end subroutine getSamCovMean

!****************************************************************************
!****************************************************************************
  
  ! This subroutine is exact same thing as getSamCovMeanLU, with the only difference that input data is transposed here on input.
  subroutine getSamCovMeanTrans(np,nd,Point,CovMat,Mean,MahalSq,InvCovMat,sqrtDetInvCovMat)

    use matrix, only: getInvPosDefMatSqrtDet
    implicit none
    integer, intent(in)             :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
    real(RK), intent(in)            :: Point(nd,np)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
    real(RK), intent(out)           :: CovMat(nd,nd)          ! Covariance matrix of the input data
    real(RK), intent(out)           :: Mean(nd)               ! Mean vector
    real(RK), intent(out), optional :: InvCovMat(nd,nd)       ! Inverse Covariance matrix of the input data
    real(RK), intent(out), optional :: MahalSq(np)            ! Vector of Mahalanobis Distances Squared, with respect to the sample's mean position
    real(RK), intent(out), optional :: sqrtDetInvCovMat       ! sqrt determinant of the inverse covariance matrix
    real(RK), dimension(nd)         :: DummyVec
    real(RK), dimension(nd,np)      :: NormedData
    integer                         :: i,j,k
    
    Mean = 0._RK
    do i = 1,np
      do j = 1,nd
        Mean(j) = Mean(j) + Point(j,i)
      end do
    end do
    Mean = Mean / real(np,kind=RK)
    
    do i = 1,np
      NormedData(1:nd,i) = Point(1:nd,i) - Mean
    end do
    
    do i = 1,nd
        do j = 1,nd
            CovMat(i,j) = dot_product(NormedData(i,1:np),NormedData(j,1:np)) / real(np-1,kind=RK)
        end do
    end do

    if (present(sqrtDetInvCovMat)) then   ! Calculate InCovMat, MahalSq, and sqrt Determinant of InCovMat
      do j = 1,nd
        do i = 1,j
          InvCovMat(i,j) = CovMat(i,j)    ! Only the upper half of CovMat is needed
        end do
      end do
      call getInvPosDefMatSqrtDet(nd,InvCovMat,sqrtDetInvCovMat)
      do i = 1,np
        do j = 1,nd
          DummyVec(j) = dot_product(InvCovMat(1:nd,j),NormedData(1:nd,i))
        end do
        MahalSq(i) = dot_product(NormedData(1:nd,i),DummyVec)
        !MahalSq = dot_product(NormedData(1:nd,i),DummyVec)
        !if (maxMahal<MahalSq) maxMahal = MahalSq
      end do
      !maxMahal = maxval(MahalSq)
    end if
    
  end subroutine getSamCovMeanTrans

!****************************************************************************
!****************************************************************************
  
  ! This subroutine combines the mean and covariance of two samples given in the input
  ! It uses the recursion equation from http://stats.stackexchange.com/questions/97642/how-to-combine-sample-means-and-sample-variances
  ! See also https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Covariance to update the covariance:
  subroutine combineCovMean(nd,npA,MeanA,CovA,npB,MeanB,CovB,Mean,CovMat)
    implicit none
    integer, intent(in)       :: nd
    integer, intent(in)       :: npA,npB
    real(RK), intent(in)      :: MeanA(nd),CovA(nd,nd)
    real(RK), intent(in)      :: MeanB(nd),CovB(nd,nd)
    real(RK), intent(out)     :: Mean(nd),CovMat(nd,nd)
    real(RK), dimension(nd,1) :: MatMeanA,MatMeanB,MatMean
    real(RK)                  :: npAB
    
    npAB = real(npA + npB,kind=RK)
    MatMeanA(1:nd,1) = MeanA
    MatMeanB(1:nd,1) = MeanB
    
    ! First find Mean:
    Mean = ( real(npA,kind=RK)*MeanA + real(npB,kind=RK)*MeanB ) / npAB
    MatMean(1:nd,1) = Mean
     
    ! Now find new Covariance matrix:
    CovMat = real(npA,kind=RK) * ( CovA + matmul(MatMeanA,transpose(MatMeanA)) ) &
           + real(npB,kind=RK) * ( CovB + matmul(MatMeanB,transpose(MatMeanB)) )
    CovMat = CovMat/npAB - matmul(MatMean,transpose(MatMean))
    
  end subroutine combineCovMean
  
!****************************************************************************
!****************************************************************************
  
  function getRandGaus()
    ! This code returns a normally distributed deviate with zero mean and unit variance.
    implicit none
    integer , save :: iset=0
    real(RK), save :: gset
    real(RK)       :: fac,rsq,vec(2)  
    real(RK)       :: getRandGaus
  
    if (iset == 0) then
      do
        call random_number(vec)
        vec = 2._RK*vec - 1._RK
        rsq = vec(1)**2 + vec(2)**2
        if ( rsq > 0._RK .and. rsq < 1._RK ) exit
      end do
      fac = sqrt(-2._RK*log(rsq)/rsq)
      gset = vec(1)*fac
      getRandGaus = vec(2)*fac
      iset = 1
    else
      getRandGaus = gset
      iset = 0
    endif
  
  end function getRandGaus
  
!****************************************************************************
!****************************************************************************
  
  ! This subroutine is legacy. use getRandMVN()
  ! Given the mean vector MeanVec and the covariance matrix CovMat, this subroutine generates a random vector x (of length nd>=2)
  ! from an nd-dimensional multivariate normal distribution.
  ! First a vector of nd standard normal random deviates is generated,
  ! then this vector is multiplied by the Cholesky decomposition of the covariance matrix,
  ! then the vector MeanVec is added to the resulting vector, and is stored in the output vector x.
  ! ATTENTION: Only the upper half of the covariance matrix (plus the diagonal terms) need to be given in the input.
  ! in the ouput, the upper half and diagonal part will still be the covariance matrix, while the lower half will be
  ! the Cholesky decomposition elements (excluding its diagonal terms that are provided only in the vector Diagonal).
  ! USES choldc.f90, getRandGaus.f90
  ! Amir Shahmoradi, March 22, 2012, 2:21 PM, IFS, UTEXAS
  subroutine getMVNDev(nd,MeanVec,CovMatIn,X)

    use matrix, only: getCholeskyFactor
    
    implicit none
    integer , intent(in)  :: nd
    real(RK), intent(in)  :: MeanVec(nd), CovMatIn(nd,nd)
    real(RK), intent(out) :: X(nd)
    real(RK)              :: CovMat(nd,nd), Diagonal(nd), DummyVec(nd)
    integer               :: i,j
  
    CovMat = CovMatIn
    call getCholeskyFactor(nd,CovMat,Diagonal)
    if (Diagonal(1)<0._RK) then
      write(*,*) 'getCholeskyFactor() failed in getMVNDev()'
      stop
    end if
    do i=1,nd
      DummyVec(i) = getRandGaus()
      x(i) = DummyVec(i) * Diagonal(i)
    end do
    do i=2,nd
      x(i) = x(i) + dot_product(CovMat(i,1:i-1),DummyVec(1:i-1))
    end do
    x = x + MeanVec
  
  end subroutine getMVNDev

!****************************************************************************
!****************************************************************************

  ! This is legacy. use getRandMVU()
  subroutine getMVUDev(nd,MeanVec,CovMatIn,X)

    ! Given the mean vector MeanVec and the covariance matrix CovMat, this subroutine generates a random vector X (of length nd>=2)
    ! from an nd-dimensional multivariate normal distribution.
    ! First a vector of nd standard normal random deviates is generated,
    ! then this vector is multiplied by the Cholesky decomposition of the covariance matrix,
    ! then the vector MeanVec is added to the resulting vector, and is stored in the output vector X.
    ! ATTENTION: Only the upper half of the covariance matrix (plus the diagonal terms) need to be given in the input.
    ! in the ouput, the upper half and diagonal part will still be the covariance matrix, while the lower half will be
    ! the Cholesky decomposition elements (excluding its diagonal terms that are provided only in the vector Diagonal).

    ! USES getCholeskyFactor.f90, getRandGaus.f90
    ! Amir Shahmoradi, April 25, 2016, 2:21 PM, IFS, UTEXAS
    use matrix, only: getCholeskyFactor
    
    implicit none
    integer , intent(in)  :: nd
    real(RK), intent(in)  :: MeanVec(nd)
    real(RK), intent(in)  :: CovMatIn(nd,nd)
    real(RK), intent(out) :: X(nd)
    real(RK)              :: Diagonal(nd), DummyVec(nd), CovMat(nd,nd), dummy
    integer               :: i,j
    
    CovMat = CovMatIn
    call getCholeskyFactor(nd,CovMat,Diagonal)
    if (Diagonal(1)<0._RK) then
      write(*,*) 'getCholeskyFactor() failed in getMVUDev()'
      stop
    end if
    do i=1,nd
      DummyVec(i) = getRandGaus()
    end do
    call random_number(dummy)
    dummy = (dummy**(1._RK/real(nd,kind=RK)))/norm2(DummyVec)  ! Now DummyVec is a uniformly drawn point from inside of nd-D sphere.
    DummyVec = dummy*DummyVec
    
    ! Now transform this point to a point inside the ellipsoid.
    do i=1,nd
      X(i) = DummyVec(i)*Diagonal(i)
    end do
    
    do i=2,nd
      X(i) = X(i) + dot_product(CovMat(i,1:i-1),DummyVec(1:i-1))
    end do
    
    X = X + MeanVec
  
  end subroutine getMVUDev
  
!****************************************************************************
!****************************************************************************
  
  ! Amir Shahmoradi, April 23, 2017, 12:36 AM, ICES, UTEXAS
  function getRandMVN(nd,MeanVec,CholeskyLower,Diagonal)
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: MeanVec(nd)
    real(RK), intent(in) :: CholeskyLower(nd,nd), Diagonal(nd)   ! Cholesky lower triangle and its diagonal terms, calculated from the MVN CovMat.
    real(RK)             :: getRandMVN(nd), dummy
    integer              :: j,i
    getRandMVN = 0._RK
    do j = 1,nd
      dummy = getRandGaus()
      getRandMVN(j) = getRandMVN(j) + Diagonal(j) * dummy
      do i = j+1,nd
        getRandMVN(i) = getRandMVN(i) + CholeskyLower(i,j) * dummy
      end do
    end do
    getRandMVN = getRandMVN + MeanVec
  end function getRandMVN

!****************************************************************************
!****************************************************************************
  
  ! Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
  ! Given the Cholesky Lower triangle and diagonals of a given covariance matrix, this function return one point
  ! uniformly randomly drawn from an nd-ellipsoid, whose nd elements getRandMVU(i), i=1,nd are guaranteed to be in range
  ! MeanVec(i) - sqrt(CovMat(i,i)) < getRandMVU(i) < MeanVec(i) + sqrt(CovMat(i,i))
  function getRandMVU(nd,MeanVec,CholeskyLower,Diagonal)
    implicit none
    integer, intent(in)  :: nd
    real(RK), intent(in) :: MeanVec(nd)
    real(RK), intent(in) :: CholeskyLower(nd,nd), Diagonal(nd)   ! Cholesky lower triangle and its diagonal terms, calculated from the MVN CovMat.
    real(RK)             :: getRandMVU(nd), dummy, DummyVec(nd), sumSqDummyVec
    integer              :: i,j
    sumSqDummyVec = 0._RK
    do j=1,nd
      DummyVec(j) = getRandGaus()
      sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
    end do
    call random_number(dummy)
    dummy = (dummy**(1._RK/dble(nd)))/sqrt(sumSqDummyVec)    ! DummyVec(j) * dummy is a uniform random point from inside of nd-sphere
    DummyVec = DummyVec * dummy
    getRandMVU = 0._RK
    do j = 1,nd
      getRandMVU(j) = getRandMVU(j) + Diagonal(j) * DummyVec(j)
      do i = j+1,nd
        getRandMVU(i) = getRandMVU(i) + CholeskyLower(i,j) * DummyVec(j)
      end do
    end do
    getRandMVU = getRandMVU + MeanVec
  end function getRandMVU

!****************************************************************************
!****************************************************************************
  
  ! Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
  ! This is algorithm is similar to getRandMVU, with the only difference that points are drawn randomly from the surface of the ellipsoid instead of inside of its interior.
  function getRandPointOnEllipsoid(nd,MeanVec,CholeskyLower,Diagonal)
    implicit none
    integer, intent(in)  :: nd
    real(RK), intent(in) :: MeanVec(nd)
    real(RK), intent(in) :: CholeskyLower(nd,nd), Diagonal(nd)   ! Cholesky lower triangle and its diagonal terms, calculated from the MVN CovMat.
    real(RK)             :: getRandPointOnEllipsoid(nd), dummy, DummyVec(nd), sumSqDummyVec
    integer              :: i,j
    sumSqDummyVec = 0._RK
    do j=1,nd
      DummyVec(j) = getRandGaus()
      sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
    end do
    DummyVec = DummyVec / sqrt(sumSqDummyVec)    ! DummyVec is a uniform random point on the surface of nd-sphere.
    getRandPointOnEllipsoid = 0._RK
    do j = 1,nd
      getRandPointOnEllipsoid(j) = getRandPointOnEllipsoid(j) + Diagonal(j) * DummyVec(j)
      do i = j+1,nd
        getRandPointOnEllipsoid(i) = getRandPointOnEllipsoid(i) + CholeskyLower(i,j) * DummyVec(j)
      end do
    end do
    getRandPointOnEllipsoid = getRandPointOnEllipsoid + MeanVec
  end function getRandPointOnEllipsoid

!****************************************************************************
!****************************************************************************

  function erfcc(x)
    implicit none
    ! returns the complementary error function with fractional error everywhere less than 1.2*10^-7.
    ! amir shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    real(RK), intent(in) :: x
    real(RK)             :: erfcc
    real(RK)             :: t,z
    z = abs(x)
    t = 1._RK/(1._RK+0.5_RK*z)
    erfcc = t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+t*&
    (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*&
    (1.48851587+t*(-.82215223+t*.17087277)))))))))
    if (x < 0._RK) erfcc = 2._RK - erfcc
  end function erfcc

!****************************************************************************
!****************************************************************************
    
  function getMahalSq(nd,np,MeanVec,InvCovMat,Point)
    ! returns sqruare of Mahalanobis distance
    implicit none
    integer, intent(in)                  :: nd,np
   !real(RK), contiguous, intent(in) :: MeanVec(:)
   !real(RK), contiguous, intent(in) :: InvCovMat(:,:)     ! Inverse of the covariance matrix
   !real(RK), contiguous, intent(in) :: Point(:,:)         ! input data points
    real(RK), intent(in)             :: MeanVec(nd)
    real(RK), intent(in)             :: InvCovMat(nd,nd)     ! Inverse of the covariance matrix
    real(RK), intent(in)             :: Point(nd,np)         ! input data points
    real(RK)                         :: getMahalSq(np)       ! function return
    integer                          :: ip
    
    do ip = 1,np
      getMahalSq(ip) = dot_product( Point(1:nd,ip)-MeanVec , matmul(InvCovMat,Point(1:nd,ip)-MeanVec) )
      if (getMahalSq(ip)<0._RK) then
        write(*,*) 'in getMahalSq(): ', getMahalSq(ip)
        stop
      end if
    end do

  end function getMahalSq

!****************************************************************************
!****************************************************************************
    
  function getMahalSqS(nd,MeanVec,InvCovMat,Point)
    ! returns sqruare of Mahalanobis distance for a single point. output is Scalar variable
    implicit none
    integer, intent(in)  :: nd
    real(RK), intent(in) :: MeanVec(nd)
    real(RK), intent(in) :: InvCovMat(nd,nd)  ! Inverse of the covariance matrix
    real(RK), intent(in) :: Point(nd)         ! input data points
    real(RK)             :: getMahalSqS       ! function return
    getMahalSqS = dot_product( Point-MeanVec , matmul(InvCovMat,Point-MeanVec) )
    if (getMahalSqS<0._RK) then
      write(*,*) 'in getMahalSqS(): ', getMahalSqS
      stop
    end if
  end function getMahalSqS

!****************************************************************************
!****************************************************************************
    
  function getProbMVN(nd,np,MeanVec,InvCovMat,sqrtDetInvCovMat,Point)
  
    !use Matrix, only: getSqrtDetPosDefMat
    use Parameters, only: INVSQRT2PI
    implicit none
    integer , intent(in) :: nd,np
    real(RK), intent(in) :: MeanVec(nd)
    real(RK), intent(in) :: InvCovMat(nd,nd)               ! inverse of covariance matrix
    real(RK), intent(in) :: sqrtDetInvCovMat               ! 
    real(RK), intent(in) :: Point(nd,np)
    real(RK)             :: dummy
    real(RK)             :: getProbMVN(np)                 ! function return
    integer              :: ip
    getProbMVN = INVSQRT2PI**nd * sqrtDetInvCovMat * exp( -0.5_RK * getMahalSq(nd,np,MeanVec,InvCovMat,Point) )
  end function getProbMVN

!****************************************************************************
!****************************************************************************

    !	Press states that this function gives good random numbers up to ~2*10**18 deviates.
  function getRandRealLecuyer(idum)     ! do not change idum between calls
    implicit none
    integer , intent(inout) :: idum
    integer , parameter     :: im1=2147483563, im2=2147483399, imm1=im1-1, ia1=40014, ia2=40692
    integer , parameter     :: iq1=53668, iq2=52774, ir1=12211, ir2=3791, ntab=32, ndiv=1+imm1/ntab
    real(RK), parameter     :: am=1._RK/dble(im1), eps=1.2d-7, rnmx=1._RK-eps
    real(RK)                :: getRandRealLecuyer
    integer                 :: idum2,j,k,iv(ntab),iy
    save                    :: iv, iy, idum2

    data idum2/123456789/, iv/ntab*0/, iy/0/
    if (idum <= 0) then
      idum = max(-idum,1)
      idum2 = idum
      do j = ntab+8,1,-1
        k = idum/iq1
        idum = ia1*(idum-k*iq1)-k*ir1
        if (idum < 0) idum = idum+im1
        if (j <= ntab) iv(j) = idum
      end do
      iy = iv(1)
    endif
    k = idum/iq1
    idum = ia1*(idum-k*iq1)-k*ir1
    if (idum < 0) idum=idum+im1
    k = idum2/iq2
    idum2 = ia2*(idum2-k*iq2)-k*ir2
    if (idum2 < 0) idum2=idum2+im2
    j = 1+iy/ndiv
    iy = iv(j)-idum2
    iv(j) = idum
    if(iy < 1)iy = iy+imm1
    getRandRealLecuyer = min(am*iy,rnmx)
  end function getRandRealLecuyer

!****************************************************************************
!****************************************************************************

  integer function getRandIntLecuyer(lowerBound,upperBound,idum)
    implicit none
    integer, intent(in)    :: lowerBound,upperBound
    integer, intent(inout) :: idum
    getRandIntLecuyer = lowerBound + nint( getRandRealLecuyer(idum)*dble(upperBound-lowerBound) )
  end function getRandIntLecuyer

!****************************************************************************
!****************************************************************************

  integer function getRandInt(lowerBound,upperBound)
    implicit none
    integer, intent(in) :: lowerBound,upperBound
    real(RK)            :: dummy
    
    call random_number(dummy)
    getRandInt = lowerBound + nint( dummy*dble(upperBound-lowerBound) )
  end function getRandInt

!****************************************************************************
!****************************************************************************

  ! alpha must be > 0, else a negative random Gamma variable will be output.
  function getRandGamma(alpha)  ! as used in GSL library
    !use Parameters, only: NAPIER
    implicit none
    real(RK), intent(in) :: alpha
    real(RK)             :: getRandGamma
    real(RK)             :: c,u,v,z    !,residual,residualInv,frac,expo,eta,Vector(3)
    if (alpha<=0._RK) then  ! illegal value of alpha
      getRandGamma = -1._RK
      return
    else
      getRandGamma = alpha
      if (getRandGamma<1._RK) getRandGamma = getRandGamma + 1._RK
      getRandGamma = getRandGamma - 0.3333333333333333
      c = 3._RK*sqrt(getRandGamma)
      c = 1._RK / c
      do
        do
          z = getRandGaus()
          v = 1._RK + c*z
          if (v<=0._RK) cycle
          exit
        end do
        v = v**3
        call random_number(u)
        if ( log(u)>=0.5_RK*z**2+getRandGamma*(1._RK-v+log(v)) ) cycle
        getRandGamma = getRandGamma*v
        exit
      end do
      if (alpha<1._RK) then
        call random_number(u)
        getRandGamma = getRandGamma * u**(1._RK/alpha)
      end if
    end if
  end function getRandGamma

!****************************************************************************
!****************************************************************************

  ! alpha must be > 1, else a negative random Gamma variable will be output.
  function getRandGammaIntShape(alpha)
    implicit none
    integer , intent(in) :: alpha
    real(RK)             :: getRandGammaIntShape
    real(RK)             :: am,e,h,s,x,y,Vector(2),Array(5)
    if (alpha < 1) then  ! illegal value of alpha
      getRandGammaIntShape = -1._RK
      return
    elseif (alpha < 6) then
      call random_number(Array(1:alpha))
      x = -log(product(Array(1:alpha)))
    else    ! use rejection sampling
      do
        call random_number(Vector)
        Vector(2) = 2._RK*Vector(2)-1._RK
        if (dot_product(Vector,Vector) > 1._RK) cycle
        y = Vector(2) / Vector(1)
        am = alpha - 1
        s = sqrt(2._RK*am + 1._RK)
        x = s*y + am
        if (x <= 0.0) cycle
        e = (1._RK+y**2) * exp(am*log(x/am)-s*y)
        call random_number(h)
        if (h <= e) exit
      end do
    end if
    getRandGammaIntShape = x
  end function getRandGammaIntShape

!****************************************************************************
!****************************************************************************

  ! both alpha, beta must be > 1, else a negative random Beta variable will be output.
  function getRandBeta(alpha,beta)
    implicit none
    real(RK), intent(in) :: alpha,beta
    real(RK)             :: getRandBeta
    real(RK)             :: x
    if ( alpha>0._RK .and. beta>0._RK ) then  ! illegal value of alpha
      x = getRandGamma(alpha)
      getRandBeta = x / ( x + getRandGamma(beta) )
    else
      getRandBeta = -1._RK
    end if
  end function getRandBeta

!****************************************************************************
!****************************************************************************

  function getRandExp()
    implicit none
    real(RK) :: getRandExp
    call random_number(getRandExp)
    getRandExp = -log(getRandExp)
  end function getRandExp

!****************************************************************************
!****************************************************************************

  function getRandCorMat(nd,eta)
    use Matrix, only: getCholeskyFactor
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: eta
    real(RK)             :: getRandCorMat(nd,nd), dummy
    real(RK)             :: beta,sumSqDummyVec,DummyVec(nd-1),W(nd-1),Diagonal(nd-1)
    integer              :: m,j,i
    
    if (nd<2 .or. eta<=0._RK) then  ! illegal value for eta. set getRandCorMat=0, return
      getRandCorMat = -1._RK
      return
    end if
    
    do m = 1,nd
      getRandCorMat(m,m) = 1._RK
    end do
    beta = eta + 0.5_RK*(nd-2._RK)
    dummy = getRandBeta(beta,beta)
    if (dummy<=0._RK .or. dummy>=1._RK) then
      write(*,*) "random Beta variable out of bound.", dummy
      stop
    end if
    getRandCorMat(1,2) = 2._RK * dummy - 1._RK    ! for the moment, only the upper half of getRandCorMat is needed, the lower half will contain cholesky lower triangle.
    
    do m = 2,nd-1
      beta = beta - 0.5_RK    
      sumSqDummyVec = 0._RK
      do j=1,m
        DummyVec(j) = getRandGaus()
        sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
      end do
      DummyVec(1:m) = DummyVec(1:m) / sqrt(sumSqDummyVec)   ! DummyVec is now a uniform random point from inside of m-sphere
      dummy = getRandBeta(0.5e0_RK*m,beta)
      W(1:m) = sqrt(dummy) * DummyVec(1:m)
      call getCholeskyFactor(m,getRandCorMat(1:m,1:m),Diagonal(1:m))
      if (Diagonal(1)<0._RK) then
        write(*,*) "Cholesky factorization failed in getRandCorMat()."
        write(*,*) m
        write(*,*) getRandCorMat(1:m,1:m)
        stop
      end if
      DummyVec(1:m) = 0._RK
      do j = 1,m
        DummyVec(j) = DummyVec(j) + Diagonal(j) * W(j)
        do i = j+1,m
          DummyVec(i) = DummyVec(i) + getRandCorMat(i,j) * DummyVec(j)
        end do
      end do
      getRandCorMat(1:m,m+1) = DummyVec(1:m)
    end do
    do i=1,nd-1
      getRandCorMat(i+1:nd,i) = getRandCorMat(i,i+1:nd)
    end do
  end function getRandCorMat

!****************************************************************************
!****************************************************************************

!  function getRandCorMat(nd,eta)    ! based on the idea of LKJ (2007). But there is something wrong with this routine
!    use Matrix, only: getCholeskyFactor
!    implicit none
!    !integer, intent(in) :: nd,eta
!    integer, intent(in) :: nd
!    real(RK), intent(in) :: eta
!    integer :: m,mNew,j,i
!    real(RK) :: getRandCorMat(nd,nd), dummy, failureCounter
!    real(RK) :: beta,sumSqDummyVec,DummyVec(nd-1),W(nd-1),Diagonal(nd-1)
!    
!    if (nd<2 .or. eta<=0._RK) then  ! illegal value for eta. set getRandCorMat=0, return
!      getRandCorMat = -1._RK
!      return
!    end if
!    
!    do m = 1,nd
!      getRandCorMat(m,m) = 1._RK
!    end do
!    beta = eta + 0.5_RK*(nd-2._RK)
!    
!    do
!      dummy = getRandBeta(beta,beta)
!      if (dummy>0._RK .and. dummy<1._RK) exit
!      write(*,*) "**Warning** random Beta variable out of bound.", dummy
!      write(*,*) "Something is wrong with getRandBeta()."
!      cycle
!    end do
!    getRandCorMat(1,2) = 2._RK * dummy - 1._RK    ! for the moment, only the upper half of getRandCorMat is needed, the lower half will contain cholesky lower triangle.
!    
!    m = 2
!    call getCholeskyFactor(m,getRandCorMat(1:m,1:m),Diagonal(1:m))
!    
!    failureCounter = 0
!    onionLayer: do
!      
!      beta = beta - 0.5_RK
!      
!      sumSqDummyVec = 0._RK
!      do j=1,m
!        DummyVec(j) = getRandGaus()
!        sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
!      end do
!      DummyVec(1:m) = DummyVec(1:m) / sqrt(sumSqDummyVec)   ! DummyVec is now a uniform random point from inside of m-sphere
!      
!      mNew = m + 1
!      posDefCheck: do
!      
!        do
!          dummy = getRandBeta(0.5_RK*m,beta)
!          if (dummy>0._RK .and. dummy<1._RK) exit
!          write(*,*) "**Warning** random Beta variable out of bound.", dummy
!          write(*,*) "Something is wrong with getRandBeta()."
!          read(*,*)
!          cycle
!        end do
!        W(1:m) = sqrt(dummy) * DummyVec(1:m)
!        
!        getRandCorMat(1:m,mNew) = 0._RK
!        do j = 1,m
!          getRandCorMat(j,mNew) = getRandCorMat(j,mNew) + Diagonal(j) * W(j)
!          do i = j+1,m
!            getRandCorMat(i,mNew) = getRandCorMat(i,mNew) + getRandCorMat(i,j) * getRandCorMat(j,mNew)
!          end do
!        end do
!        
!        
!        call getCholeskyFactor(mNew,getRandCorMat(1:mNew,1:mNew),Diagonal(1:mNew))  ! Now check if the new matrix is positive-definite, then proceed with the next layer
!        if (Diagonal(1)<0._RK) then
!          failureCounter = failureCounter + 1
!          cycle posDefCheck
!          !write(*,*) "Cholesky factorization failed in getRandCorMat()."
!          !write(*,*) m
!          !write(*,*) getRandCorMat(1:m,1:m)
!          !stop
!        end if
!        exit posDefCheck
!
!      end do posDefCheck
!
!      if (mNew==nd) exit onionLayer
!      m = mNew
!        
!    end do onionLayer
!    
!    if (failureCounter>0) write(*,*) 'failureRatio: ', dble(failureCounter)/dble(nd-2)
!    do i=1,nd-1
!      getRandCorMat(i+1:nd,i) = getRandCorMat(i,i+1:nd)
!    end do
!    
!  end function getRandCorMat

!****************************************************************************
!****************************************************************************
  
  ! Only the upper half is the correlatrion matrix, lower half is giberish
  function getRandCorMatRejection(nd,minRho,maxRho)
    use Matrix, only: isPosDef
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: minRho,maxRho
    real(RK)             :: getRandCorMatRejection(nd,nd),RhoVec(nd*(nd-1))
    integer              :: i,j,irho
    if (maxRho<minRho .or. nd<1) then
      write(*,*) 'invalid input to getRandCorMatRejection()'
      stop
    end if
    if (nd==1) then
      getRandCorMatRejection = 1._RK
    else
      rejection: do
        call random_number(RhoVec)
        RhoVec = minRho + RhoVec*(maxRho-minRho)
        irho = 0
        do j=1,nd
          getRandCorMatRejection(j,j) = 1._RK
          do i=1,j-1
            irho = irho + 1
            getRandCorMatRejection(i,j) = RhoVec(irho)
          end do
        end do
        if (isPosDef(nd,getRandCorMatRejection)) exit rejection
        cycle rejection
      end do rejection
    end if
    do j=1,nd-1
      getRandCorMatRejection(j+1:nd,j) = getRandCorMatRejection(j,j+1:nd)
    end do
  end function getRandCorMatRejection
  
!****************************************************************************
!****************************************************************************
  
  ! Only the upper half is the correlatrion matrix, lower half is giberish
  function getCovMatFromCorMat(nd,SigmaVec,CorMat)
    use Matrix, only: isPosDef
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: SigmaVec(nd), CorMat(nd,nd)   ! only upper half needed
    real(RK)             :: getCovMatFromCorMat(nd,nd)
    integer              :: i,j
    do j=1,nd
      getCovMatFromCorMat(j,j) = SigmaVec(j)**2
      do i=1,j-1
        getCovMatFromCorMat(i,j) = CorMat(i,j) * SigmaVec(j) * SigmaVec(i)
      end do
    end do
  end function getCovMatFromCorMat
  
!****************************************************************************
!****************************************************************************

  function getNormCDF(mean,stdev,x)
    use Parameters, only: SP,SQRT2
    implicit none
    real(RK), intent(in) :: mean,stdev,x
    real(RK)             :: getNormCDF
    getNormCDF = 0.5_RK * ( 1._RK + erf( real((x-mean)/(SQRT2*stdev),kind=SP) ) )
  end function getNormCDF

!****************************************************************************
!****************************************************************************

  function getSNormCDF(x)
    use Parameters, only: SP,SQRT2
    implicit none
    real(RK), intent(in) :: x
    real(RK)             :: getSNormCDF
    getSNormCDF = 0.5_RK * ( 1._RK + erf( real(x/SQRT2,kind=SP) ) )
  end function getSNormCDF

!****************************************************************************
!****************************************************************************

  ! if x is not in [0,1], a negative value for getBetaCDF will be return as a sign of error
  ! The efficiency of this code can be improved by making x a vector on input.
  function getBetaCDF(a,b,x)
    use Parameters, only : SP
    implicit none
    real(RK), intent(in) :: a,b,x
    real(RK)             :: bt
    real(RK)             :: getBetaCDF
    if (x < 0._RK .or. x > 1._RK) then
      getBetaCDF = -1._RK
    end if
    if (x == 0._RK .or. x == 1._RK) then
      bt = 0._RK
    else
      bt = exp( log_gamma(real(a+b,kind=SP)) - log_gamma(real(a,kind=SP)) - log_gamma(real(b,kind=SP)) &
              + a*log(x) + b*log(1._RK-x) )
    end if
    if ( x < (a+1._RK) / (a+b+2._RK) ) then
      getBetaCDF = bt * getBetaContinuedFraction(a,b,x) / a
    else
      getBetaCDF = 1._RK - bt * getBetaContinuedFraction(b,a,1._RK-x) / b
    end if
  end function getBetaCDF

!****************************************************************************
!****************************************************************************

  function getBetaContinuedFraction(a,b,x)
    implicit none
    real(RK), intent(in) :: a,b,x
    real(RK), parameter  :: eps = epsilon(x), fpmin = tiny(x)/eps
    integer , parameter  :: maxit = 100
    real(RK)             :: aa,c,d,del,qab,qam,qap
    real(RK)             :: getBetaContinuedFraction
    integer              :: m,m2
    qab = a+b
    qap = a+1._RK
    qam = a-1._RK
    c = 1._RK
    d = 1._RK-qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1._RK/d
    getBetaContinuedFraction = d
    do m = 1,maxit
      m2 = 2*m
      aa = m*(b-m)*x/((qam+m2)*(a+m2))
      d = 1._RK+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1._RK+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1._RK/d
      getBetaContinuedFraction = getBetaContinuedFraction*d*c
      aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
      d = 1._RK+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1._RK+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1._RK/d
      del = d*c
      getBetaContinuedFraction = getBetaContinuedFraction*del
      if (abs(del-1._RK) <=  eps) exit
    end do
    if (m > maxit) then
      write(*,"('In getBetaContinuedFraction() in Statistics: a or b too big, or maxit too small.')")
      stop
    end if
  end function getBetaContinuedFraction

!****************************************************************************
!****************************************************************************

  subroutine doKS1(np,Point,getCDF,statKS,probKS)
    use SpecialFunctions, only : sortAscending
    implicit none
    integer , intent(in)    :: np
    real(RK), intent(out)   :: statKS,probKS
    real(RK), intent(inout) :: Point(np)
    real(RK)                :: npSqrt
    real(RK)                :: cdf,cdfObserved,dt,frac
    integer                 :: j
    interface
      function getCDF(x)
        use Parameters, only: RK
        real(RK), intent(in) :: x
        real(RK)             :: getCDF
      end function getCDF
    end interface
    
    call sortAscending(np,Point)
    statKS = 0._RK
    cdfObserved = 0._RK
    npSqrt = np
    do j = 1,np
      frac = j/npSqrt
      cdf = getCDF(Point(j))
      dt = max( abs(cdfObserved-cdf) , abs(frac-cdf) )
      if( dt > statKS ) statKS = dt
      cdfObserved = frac
    end do
    npSqrt = sqrt(npSqrt)
    probKS = getProbKS( (npSqrt+0.12_RK+0.11_RK/npSqrt)*statKS )
  end subroutine doKS1

!****************************************************************************
!****************************************************************************

  subroutine doUniformKS1(np,Point,statKS,probKS)   ! This assumes that points are coming from a uniform distribution in [0,1]. So, all Point must be in [0,1] on input.
    use SpecialFunctions, only : sortAscending
    implicit none
    integer, intent(in)     :: np
    real(RK), intent(out)   :: statKS,probKS
    real(RK), intent(inout) :: Point(np)
    real(RK)                :: npSqrt
    real(RK)                :: cdf,cdfObserved,dt,frac
    integer                 :: j    
    call sortAscending(np,Point)
    statKS = 0._RK
    cdfObserved = 0._RK
    npSqrt = np
    do j = 1,np
      frac = j/npSqrt
      cdf = Point(j)
      dt = max( abs(cdfObserved-cdf) , abs(frac-cdf) )
      if( dt > statKS ) statKS = dt
      cdfObserved = frac
    end do
    npSqrt = sqrt(npSqrt)
    probKS = getProbKS( (npSqrt+0.12_RK+0.11_RK/npSqrt)*statKS )
  end subroutine doUniformKS1

!****************************************************************************
!****************************************************************************

  function getProbKS(lambda)
    implicit none
    real(RK), intent(in) :: lambda
    real(RK), parameter  :: EPS1 = 0.001_RK, EPS2 = 1.e-8_RK
    integer , parameter  :: NITER=100
    integer              :: j
    real                 :: a2,fac,term,termbf
    real(RK)             :: getProbKS
    a2 = -2._RK*lambda**2
    fac = 2._RK
    getProbKS = 0._RK
    termbf = 0._RK
    do j = 1, NITER
      term = fac*exp(a2*j**2)
      getProbKS = getProbKS+term
      if (abs(term) <= EPS1*termbf .or. abs(term) <= EPS2*getProbKS) return
      fac = -fac
      termbf = abs(term)
    end do
    getProbKS = 1.0
  end function getProbKS

!****************************************************************************
!****************************************************************************

  function getUniformCDF(x)     ! returns the uniform CDF on support [0,1)
    implicit none
    real(RK), intent(in) :: x
    real(RK)             :: getUniformCDF
    getUniformCDF = x
  end function getUniformCDF

!****************************************************************************
!****************************************************************************

end module Statistics

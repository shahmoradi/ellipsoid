module Matrix
  
  use Parameters, only: RK
  
  implicit none

contains

!****************************************************************************
!****************************************************************************

  ! Returns on the lower triangle of PosDefMat, the Cholesky factorization L of PosDefMat=L.L^T.
  ! On input, the upper upper triangle of PosDefMat should be given, which remains intact on output.
  subroutine getCholeskyFactor(nd,PosDefMat,Diagonal)
    implicit none
    integer , intent(in)    :: nd
    real(RK), intent(inout) :: PosDefMat(nd,nd) ! Upper triangle + diagonal is input matrix, lower is output.
    real(RK), intent(out)   :: Diagonal(nd)
    real(RK)                :: summ
    integer             :: i
    do i=1,nd
      summ = PosDefMat(i,i) - dot_product(PosDefMat(i,1:i-1),PosDefMat(i,1:i-1))
      if (summ <= 0._RK) then
        Diagonal(1) = -1._RK
        return
        !write (*,*) 'getCholeskyFactor() failed'
        !stop
      end if
      Diagonal(i) = sqrt(summ)
      PosDefMat(i+1:nd,i) = ( PosDefMat(i,i+1:nd) - matmul(PosDefMat(i+1:nd,1:i-1),PosDefMat(i,1:i-1)) ) / Diagonal(i)
    end do
  end subroutine getCholeskyFactor
    
!****************************************************************************
!****************************************************************************

  ! Solve the linear equation system of the form: PosDefMat.InputSolution = Intercept
  ! PosDefMat and Diagonal are the output of subroutine getCholeskyFactor (i.g., only the lower triangle of PosDefMat is used).
  subroutine solveLinearPosDefSystem(nd,PosDefMat,Diagonal,Intercept,InputSolution)
    implicit none
    integer , intent(in)    :: nd
    real(RK), intent(in)    :: PosDefMat(nd,nd),Diagonal(nd),Intercept(nd)
    real(RK), intent(inout) :: InputSolution(nd)
    integer                 :: i
    do i=1,nd
      InputSolution(i) = ( Intercept(i) - dot_product(PosDefMat(i,1:i-1),InputSolution(1:i-1)) ) / Diagonal(i)
    end do
    do i = nd,1,-1
      InputSolution(i) = ( InputSolution(i) - dot_product(PosDefMat(i+1:nd,i),InputSolution(i+1:nd)) ) / Diagonal(i)
    end do
  end subroutine solveLinearPosDefSystem

!****************************************************************************
!****************************************************************************

  ! This code returns the inverse matrix of a symmetric-positive-definite input matrix which is given in the upper triangle of MatInvMat.
  ! On output MatInvMat is completely overwritten by in the inverse of the matrix.
  ! Also returns: InverseMatrix of input matrix
  ! Determinant of 
  ! Amir Shahmoradi, Apr 21, 2017, 1:54 AM, ICES, UT
  ! USES getCholeskyFactor.f90 & getCholeskyFactor.f90
  subroutine getInvPosDefMatSqrtDet(nd,MatInvMat,sqrtDetInvPosDefMat)
    implicit none
    integer , intent(in)    :: nd
    real(RK), intent(inout) :: MatInvMat(nd,nd)           ! input: upper half is covariance matrix, output: inverse matrix
    real(RK), intent(out)   :: sqrtDetInvPosDefMat        ! determinant of the inverse matrix
    real(RK)                :: CholeskyLower(nd,nd)       ! Cholesky factor
    real(RK)                :: Diagonal(nd)               ! diagonal terms of the inverse matrix
    real(RK)                :: summ
    integer                 :: i,j,k
    do j=1,nd
      do i=1,j
        CholeskyLower(i,j) = MatInvMat(i,j)
      end do
    end do
    call getCholeskyFactor(nd,CholeskyLower,Diagonal)
    if (Diagonal(1)<0._RK) then
      sqrtDetInvPosDefMat = -1._RK
      return
    end if
    sqrtDetInvPosDefMat = 1._RK / product(Diagonal)
    do i = 1,nd
      CholeskyLower(i,i) = 1._RK / Diagonal(i)
      do j = i+1,nd
          summ = 0._RK
          do k = i,j-1
            summ = summ - CholeskyLower(j,k) * CholeskyLower(k,i)
          end do
        CholeskyLower(j,i) = summ / Diagonal(j)
      end do
    end do
    do i = 1,nd
      do j = i,nd
        MatInvMat(j,i) = dot_product(CholeskyLower(j:nd,j),CholeskyLower(j:nd,i))
      end do
      MatInvMat(i,i:nd) = MatInvMat(i:nd,i)
    end do
  end subroutine getInvPosDefMatSqrtDet

!****************************************************************************
!****************************************************************************
  
  ! Returns the inverse of a matrix whose Cholesky Lower traingle is given in the lower part of CholeskyLower, and its diagonals in Diagonal
  function getInvMatFromCholFac(nd,CholeskyLower,Diagonal)
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: CholeskyLower(nd,nd),Diagonal(nd)
    real(RK)             :: getInvMatFromCholFac(nd,nd)
    real(RK)             :: summ
    integer              :: i,j,k
    getInvMatFromCholFac = 0._RK
    do j=1,nd-1
      do i=j+1,nd
        getInvMatFromCholFac(i,j) = CholeskyLower(i,j)
      end do
    end do
    do i = 1,nd
      getInvMatFromCholFac(i,i) = 1._RK / Diagonal(i)
      do j = i+1,nd
          summ = 0._RK
          do k = i,j-1
            summ = summ - getInvMatFromCholFac(j,k) * getInvMatFromCholFac(k,i)
            !write(*,*) 'Used:'
            !write(*,*) j,k
            !write(*,*) k,i
          end do
        getInvMatFromCholFac(j,i) = summ / Diagonal(j)
        !write(*,*) 'Assigned:'
        !write(*,*) j,i
      end do
    end do
    do i = 1,nd
      do j = i,nd
        getInvMatFromCholFac(j,i) = dot_product(getInvMatFromCholFac(j:nd,j),getInvMatFromCholFac(j:nd,i))
      end do
      getInvMatFromCholFac(i,i:nd) = getInvMatFromCholFac(i:nd,i)
    end do
  end function getInvMatFromCholFac

!****************************************************************************
!****************************************************************************

  ! exact same thing as subroutine getInvPosDefMatSqrtDet, but returns the result as a function.
  ! NOTE: according to my timing tests, compare_InvMatRoutines_1(), the function version seems to be 15-30% faster than the subroutine version above.
  ! Amir Shahmoradi, Apr 8, 2017, 1:54 PM, ICES, UT
  ! USES getLU.f90 & solveLinearSystem.f90
  function getInvPosDefMat(nd,PosDefMat)
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: PosDefMat(nd,nd)
    real(RK)             :: getInvPosDefMat(nd,nd),CholeskyLower(nd,nd)
    real(RK)             :: Diagonal(nd)
    real(RK)             :: summ
    integer              :: i,j,k
    do j=1,nd
      do i=1,j
        CholeskyLower(i,j) = PosDefMat(i,j)
      end do
    end do
    call getCholeskyFactor(nd,CholeskyLower,Diagonal)
    if (Diagonal(1)<0._RK) then
      write(*,*) 'getCholeskyFactor() failed in getInvPosDefMat()'
      stop
    end if
    do i = 1,nd
      CholeskyLower(i,i) = 1._RK / Diagonal(i)
      do j = i+1,nd
          summ = 0._RK
          do k = i,j-1
            summ = summ - CholeskyLower(j,k) * CholeskyLower(k,i)
          end do
        CholeskyLower(j,i) = summ / Diagonal(j)
      end do
    end do
    do i = 1,nd
      do j = i,nd
        getInvPosDefMat(j,i) = dot_product(CholeskyLower(j:nd,j),CholeskyLower(j:nd,i))
      end do
      getInvPosDefMat(i,i:nd) = getInvPosDefMat(i:nd,i)
    end do
  end function getInvPosDefMat

!****************************************************************************
!****************************************************************************

  ! This code replaces returns the inverse matrix Y of a n*n matrix S, and its determinant.
  ! NOTE: according to my timing tests, compare_InvMatRoutines_1(), the function version of this code below seems to be 15-30% faster than this subroutine version.
  ! Amir Shahmoradi, Oct 18, 2009, 1:54 AM, MTU
  ! USES getLU.f90 & solveLinearSystem.f90
  subroutine getInvMatDet(n,LU,InverseMatrix,detInvMat)
    implicit none
    integer , intent(in)    :: n
    real(RK), intent(inout) :: LU(n,n)     ! on input it is the matrix, on output it is LU decomposition
    real(RK), intent(out)   :: InverseMatrix(n,n)
    real(RK), intent(out)   :: detInvMat   ! determinant of the inverse matrix
    integer                 :: i,j,Permutation(n)
    do i = 1,n
      do j = 1,n
        InverseMatrix(i,j) = 0._RK
      end do
      InverseMatrix(i,i) = 1._RK
    end do
    call getLU(n,LU,Permutation,detInvMat)
    do j = 1,n
      detInvMat = detInvMat * LU(j,j)
      call solveLinearSystem(n,LU,Permutation,InverseMatrix(1:n,j))
    end do
    detInvMat = 1._RK/detInvMat
  end subroutine getInvMatDet

!****************************************************************************
!****************************************************************************
    
  ! exact same thing as subroutine getInvMatDet, but returns the result as a function.
  ! NOTE: according to my timing tests, compare_InvMatRoutines_1(), the function version seems to be 15-30% faster than the subroutine version above.
  ! This code replaces returns the inverse matrix Y of a n*n matrix S of physical dimension np*np.
  ! Amir Shahmoradi, Apr 8, 2017, 1:54 PM, MTU
  ! USES getLU.f90 & solveLinearSystem.f90
  function getInvMat(n,S)
  implicit none
  integer , intent(in) :: n
  real(RK), intent(in) :: S(n,n)
  real(RK)             :: getInvMat(n,n),DummyMat(n,n)
  integer              :: i,j,Permutation(n)
  real(RK)             :: parity
  do i = 1,n
      do j = 1,n
          getInvMat(i,j) = 0._RK
      end do
      getInvMat(i,i) = 1._RK
  end do
  DummyMat = S
  call getLU(n,DummyMat,Permutation,parity)
  do j = 1,n
      call solveLinearSystem(n,DummyMat,Permutation,getInvMat(1:n,j))
  end do
  end function getInvMat

!****************************************************************************
!****************************************************************************

  ! Solves the set of n linear equations A.X = B. Here a is input, not as the matrix A but
  ! rather as its LU decomposition, determined by the routine getLU. indx is input as the
  ! permutation vector returned by getLU. b(1:n) is input as the right-hand side vector B,
  ! and returns with the solution vector X. a, n, np, and indx are not modiffied by this routine
  ! and can be left in place for successive calls with different right-hand sides b. This routine
  ! takes into account the possibility that b will begin with many zero elements, so it is efficient
  ! for use in matrix inversion.
  ! Amir Shahmoradi, Oct 18, 2009, 1:54 PM, MTU
  subroutine solveLinearSystem(nd,Array,Permutation,InputSolution)
    integer , intent(in)    :: nd
    integer , intent(in)    :: Permutation(nd)
    real(RK), intent(in)    :: Array(nd,nd)
    real(RK), intent(inout) :: InputSolution(nd)
    integer                 :: i,ii
    real(RK)                :: summ
    ii = 0
    do i=1,nd
      summ = InputSolution(Permutation(i))
      InputSolution(Permutation(i)) = InputSolution(i)
      if (ii /= 0) then
        summ = summ - dot_product( Array(i,ii:i-1) , InputSolution(ii:i-1) )
      else if (summ /= 0._RK) then
        ii = i
      end if
      InputSolution(i) = summ
    end do
    do i = nd,1,-1
      InputSolution(i) = ( InputSolution(i) - dot_product( Array(i,i+1:nd) , InputSolution(i+1:nd)) ) / Array(i,i)
    end do
  end subroutine solveLinearSystem

!****************************************************************************
!****************************************************************************

  ! Returns the LU decomposition of the input matrix Array(nd,nd).
  ! Permutation(1:n) is an output vector that records the row permutation effected by the partial pivoting
  ! Parity is output as +-1 depending on whether the number of row interchanges was even or odd, respectively.
  ! This routine is used in combination with solveLinearSystem to solve linear equations or invert a matrix.
  ! Amir Shahmoradi, Apr 21, 2017, 1:43 PM, ICES, UT
  subroutine getLU(n,Array,Permutation,parity) !,errorMessage)
    integer , intent(in)    :: n
    real(RK), intent(inout) :: Array(n,n)
    integer , intent(out)   :: Permutation(n)
    real(RK), intent(out)   :: parity
    real(RK), parameter     :: TINY = 1.e-20_RK
    integer                 :: i,imax,j,k
    real(RK)                :: aamax,dum,summ,vv(n)
    !character(len=63)   ``     :: errorMessage
    parity = 1._RK
    do i = 1,n
      aamax = 0._RK
      do j=1,n
        if ( abs(Array(i,j)) > aamax ) aamax = abs( Array(i,j) )
      end do
      if (aamax == 0._RK) then
          write(*,*) 'singular matrix in getLU()'
          STOP
      end if
      vv(i) = 1._RK/aamax
    end do
    
    do j=1,n
      do i=1,j-1
        summ = Array(i,j)
        do k=1,i-1
          summ = summ - Array(i,k)*Array(k,j)
        end do
        Array(i,j) = summ
      end do
      aamax = 0._RK
      do i = j, n
        summ = Array(i,j)
        do k = 1, j-1
          summ = summ - Array(i,k)*Array(k,j)
        end do
        Array(i,j) = summ
        dum = vv(i) * abs(summ)
        if (dum >= aamax) then
          imax = i
          aamax = dum
        endif
      end do
      if (j /= imax)then
        do k=1,n
          dum = Array(imax,k)
          Array(imax,k) = Array(j,k)
          Array(j,k) = dum
        end do
        parity = - parity
        vv(imax) = vv(j)
      endif
      Permutation(j) = imax
      if (Array(j,j) == 0._RK) Array(j,j) = TINY
      if (j /= n) then
        dum = 1._RK / Array(j,j)
        do i = j+1,n
          Array(i,j) = Array(i,j) * dum
        end do
      endif
    end do
  end subroutine getLU

!****************************************************************************
!****************************************************************************

!    subroutine to find product of two matrices
!    Amir Shahmoradi, Oct 20, 2009, 10:56 PM, MTU
!    Product of two matrices is defined by
!    c(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j)+........+a(i,n)*b(n,j)
    subroutine multiplyMatrix(A,rowsA,colsA,B,rowsB,colsB,C)
    implicit none
    integer , intent(in)  :: rowsA, colsA, rowsB, colsB !Matrix Dimensions
    real(RK), intent(in)  :: A(rowsA,colsA) !Matrix A
    real(RK), intent(in)  :: B(rowsB,colsB) !Matrix B
    real(RK), intent(out) :: C(rowsA,colsB) !Matrix C
    integer  :: i,j,k,rowsC,colsC !Counters
    IF (colsA /= rowsB) THEN
        write(*,*) 'Error! Order of matrices donot match'
!        Two matrices can be multiplied if and only if the number of columns
!        of the first matrix equals the number of rows of the second
        STOP
    ELSE
        rowsC = rowsA
        colsC = colsB
    ENDIF
!    Initialize product matrix to 0
    do i = 1, rowsC
        do j = 1, colsC
            C(i,j) = 0._RK
        end do
    end do
!    Find product as per above formula
    do i = 1, rowsA
        do j = 1, colsB
            do k = 1, colsA
                C(i,j) = C(i,j) + A(i,k)*B(k,j)
            end do
        end do
    end do
!    ErrCode = 1 => subroutine returned successfully
    end subroutine multiplyMatrix 

!****************************************************************************
!****************************************************************************

  ! This subroutine finds the determinant of a given nd*nd matrix InputMat.
  ! Amir Shahmoradi, Oct 18, 2009, 4:10 PM, MTU
  ! USES getLU.f90
  real(RK)function getDeterminant(nd,InputMat)
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: InputMat(nd,nd)
    integer              :: Permutation(nd)
    real(RK)             :: DummyMat(nd,nd)
    integer              :: i,j
    DummyMat = InputMat
    call getLU(nd,DummyMat,Permutation,getDeterminant) !    This returns getDeterminant as +-1.
    do j=1,nd
      getDeterminant = getDeterminant*DummyMat(j,j)
    end do
  end function getDeterminant

!****************************************************************************
!****************************************************************************

  ! This subroutine finds the determinant of a given positive-definite nd*nd matrix PosDefMat.
  ! Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
  ! USES getCholeskyFactor.f90
  real(RK) function getSqrtDetPosDefMat(nd,PosDefMat)
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: PosDefMat(nd,nd)
    real(RK)             :: Diagonal(nd),DummyMat(nd,nd)
    integer              :: i,j
    do j=1,nd
      do i=1,j
        DummyMat(i,j) = PosDefMat(i,j)
      end do
    end do
    call getCholeskyFactor(nd,DummyMat,Diagonal)
    if (Diagonal(1)<0._RK) then
      getSqrtDetPosDefMat = -1._RK
      return
    end if
    getSqrtDetPosDefMat = product(Diagonal)
  end function getSqrtDetPosDefMat

!****************************************************************************
!****************************************************************************

  pure logical function isPosDef(nd,ArrayIn)
    ! this subroutines returns False value for isPosDef, if the cholesky decomposition fails (i.e. matrix is not positive definite), otherwise isPosDef=true)
    implicit none
    integer , intent(in) :: nd
    real(RK), intent(in) :: ArrayIn(nd,nd)
    real(RK)             :: Array(nd,nd),p(nd)
    real(RK)             :: dummySum
    integer              :: i,j,k
    isPosDef = .true.
    Array = ArrayIn
    do i = 1,nd
      do j = i,nd
        dummySum = Array(i,j)
        do k = i-1,1,-1
          dummySum = dummySum - Array(i,k) * Array(j,k)
        end do
        if(i==j)then
          if(dummySum<=0._RK) then 
        isPosDef = .false.
            return
          end if
          p(i) = sqrt(dummySum)
        else
          Array(j,i)=dummySum/p(i)
        endif
      end do
    end do
  end function isPosDef

!****************************************************************************
!****************************************************************************
  
  function getOuterProd(Array1,Array2)
    real(RK), intent(in) :: Array1(:),Array2(:)
    real(RK)             :: getOuterProd(size(Array1),size(Array2))
    getOuterProd = spread(Array1,dim=2,ncopies=size(Array2)) * spread(Array2,dim=1,ncopies=size(Array1))
  end function getOuterProd

!****************************************************************************
!****************************************************************************

end module Matrix
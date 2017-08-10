module SpecialFunctions

  use Parameters, only: RK
  implicit none

contains

!****************************************************************************
!****************************************************************************

  function getFactorial(intNum) ! calculates the Gamma function for whole integer input. This is basically factorial(intNum-1)
    implicit none
    integer, intent(in) :: intNum
    integer             :: i
    real(RK)            :: getFactorial
    getFactorial = 1._RK
    do i = 2,intNum
      getFactorial = getFactorial * i
    end do
  end function getFactorial

!****************************************************************************
!****************************************************************************

  function getGammaHalfInt(halfIntNum) ! calculates the Gamma function for a half integer input.
    use Parameters, only: SQRTPI
    implicit none
    real(RK), intent(in) :: halfIntNum
    real(RK)             :: getGammaHalfInt
    integer              :: i,k
    getGammaHalfInt = SQRTPI
    k = nint(halfIntNum-0.5_RK)    ! halfIntNum = k + 1/2
    do i = k+1, 2*k
      getGammaHalfInt = getGammaHalfInt * i / 4._RK
    end do
  end function getGammaHalfInt

!****************************************************************************
!****************************************************************************

  function getEllVolCoef(nd)
    ! calculates the coefficient for the calculation of an nd-dimensional ellipsoid.
    ! reference: http://math.stackexchange.com/questions/606184/volume-of-n-dimensional-ellipsoid
    ! reference: https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    ! reference: https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function
    ! getEllVolCoef = PI^(nd/2) / gamma(nd/2+1) where n is just a positive integer
    use Parameters, only: PI
    implicit none
    integer, intent(in) :: nd
    real(RK)            :: getEllVolCoef
    integer             :: i,k
    if (mod(nd,2)==0) then  ! nd is even
      getEllVolCoef = PI
      do i = 2, nd/2
        getEllVolCoef = getEllVolCoef * PI / i    ! nd = 2k ; getEllVolCoef = PI^(k) / Factorial(k)
      end do
    else    ! nd is an odd integer
      k = (nd+1)/2   ! nd = 2k-1 ; gamma(nd/2 + 1) = gamma(k + 1/2) ; gamma(k+1/2) = sqrt(PI) * (2k)! / (4^k * k!)
      getEllVolCoef = 4._RK / (k+1)     ! This is to avoid an extra unnecessary division of getEllVolCoef by PI
      do i = k+2, 2*k
        getEllVolCoef = getEllVolCoef * PI * 4._RK / i   ! getEllVolCoef = PI^(k-1/2) / gamma(k+1/2) = PI^(k+1) * 4^k * k! / (2k)!
      end do
    end if
  end function getEllVolCoef

!****************************************************************************
!****************************************************************************
  
  ! This function has not been tested yet
  ! There is really no reason to use this function for HPC calculations, as it can be likely done more efficiently.
  function getEllipsoidVolume(nEllipsoid,ndim,sqrtDeterminant)
    use Parameters, only: SQRTPI
    implicit none
    integer         :: i
    integer, intent(in) :: nEllipsoid                   ! number of ellipsoids
    integer, intent(in) :: ndim(nEllipsoid)             ! dimensions of the ellipsoids
    integer, intent(in) :: sqrtDeterminant(nEllipsoid)  ! sqrtDeterminant of the covariance matrices of the ellipsoids
    real(RK)            :: getEllipsoidVolume(nEllipsoid)
    real(RK)            :: factor
    do i=1,nEllipsoid
      getEllipsoidVolume(i) = getEllVolCoef(ndim(i)) / sqrtDeterminant(i)
    end do 
  end function getEllipsoidVolume

!****************************************************************************
!****************************************************************************

  recursive subroutine sortArray(array)
    real(RK), intent(inout) :: array(:)
    integer             :: iq
  
    if(size(array) > 1) then
       call partition(array, iq)
       call sortArray(array(:iq-1))
       call sortArray(array(iq:))
    endif
  end subroutine sortArray

  subroutine partition(array, marker)
    implicit none
    real(RK), intent(inout) :: array(:)
    integer , intent(out)   :: marker
    integer                 :: i, j
    real(RK)                :: temp
    real(RK)                :: x      ! pivot point
    x = array(1)
    i= 0
    j= size(array) + 1
    do
       j = j-1
       do
          if (array(j) <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (array(i) >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange array(i) and array(j)
          temp = array(i)
          array(i) = array(j)
          array(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do
  end subroutine partition
  
!****************************************************************************
!****************************************************************************

  subroutine sortAscending(np,Point)
    implicit none
    integer , parameter     :: nn = 15, nstack = 100
    integer , intent(in)    :: np
    real(RK), intent(inout) :: Point(np)
    real(RK)                :: dummy
    integer                 :: k,i,j,jstack,m,r,istack(nstack)
    jstack=0
    m = 1
    r = np
    do
      if (r-m < nn) then
        do j = m+1,r
          dummy = Point(j)
          do i = j-1,m,-1
            if (Point(i) <=  dummy) exit
            Point(i+1) = Point(i)
          end do
          Point(i+1) = dummy
        end do
        if (jstack == 0) return
        r = istack(jstack)
        m = istack(jstack-1)
        jstack = jstack-2
      else
        k = (m+r)/2
        call swapReal(Point(k),Point(m+1))
        if (Point(m)>Point(r)) call swapReal(Point(m),Point(r))
        if (Point(m+1)>Point(r)) call swapReal(Point(m+1),Point(r))
        if (Point(m)>Point(m+1)) call swapReal(Point(m),Point(m+1))
        i = m+1
        j = r
        dummy = Point(m+1)
        do
          do
            i = i+1
            if (Point(i) >= dummy) exit
          end do
          do
            j = j-1
            if (Point(j) <= dummy) exit
          end do
          if (j < i) exit
          call swapReal(Point(i),Point(j))
        end do
        Point(m+1) = Point(j)
        Point(j) = dummy
        jstack = jstack+2
        if (jstack > nstack) then
          write(*,*) 'sortAscending() failed: nstack too small'
          stop
        end if
        if (r-i+1 >= j-m) then
          istack(jstack) = r
          istack(jstack-1) = i
          r = j-1
        else
          istack(jstack) = j-1
          istack(jstack-1) = m
          m = i
        end if
      end if
    end do
  end subroutine sortAscending

!****************************************************************************
!****************************************************************************

  subroutine swapReal(a,b)
    real(RK), intent(inout) :: a,b
    real(RK)                :: dummy
    dummy = a
    a = b
    b = dummy
  end subroutine swapReal

!****************************************************************************
!****************************************************************************

  subroutine swapInt(a,b)
    integer, intent(inout) :: a,b
    integer                :: dummy
    dummy = a
    a = b
    b = dummy
  end subroutine swapInt

!****************************************************************************
!****************************************************************************

  subroutine indexArrayReal(n,Array,index)
    implicit none
    real(RK), intent(in)  :: Array(n)
    integer , intent(out) :: index(n)
    integer , parameter   :: nn=15, nstack=50
    integer               :: n,k,i,j,indext,jstack,l,r
    integer               :: istack(nstack)
    real(RK)              :: a
    do j = 1,n
      index(j) = j
    end do
    jstack=0
    l=1
    r=n
    do
      if (r-l < nn) then
        do j=l+1,r
          indext=index(j)
          a=Array(indext)
          do i=j-1,l,-1
            if (Array(index(i)) <= a) exit
            index(i+1)=index(i)
          end do
          index(i+1)=indext
        end do
        if (jstack == 0) return
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+r)/2
        call swapInt(index(k),index(l+1))
        call exchangeIndex(index(l),index(r))
        call exchangeIndex(index(l+1),index(r))
        call exchangeIndex(index(l),index(l+1))
        i=l+1
        j=r
        indext=index(l+1)
        a=Array(indext)
        do
          do
            i=i+1
            if (Array(index(i)) >= a) exit
          end do
          do
            j=j-1
            if (Array(index(j)) <= a) exit
          end do
          if (j < i) exit
          call swapInt(index(i),index(j))
        end do
        index(l+1)=index(j)
        index(j)=indext
        jstack=jstack+2
        if (jstack > nstack) then
          write(*,*) 'NSTACK too small in indexArrayReal'
          stop
        end if
        if (r-i+1 >= j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
        else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
        end if
      end if
    end do
  
  contains

    subroutine exchangeIndex(i,j)
      integer, intent(inout) :: i,j
      integer            :: swp
      if (Array(j) < Array(i)) then
        swp=i
        i=j
        j=swp
      end if
    end subroutine exchangeIndex

  end subroutine indexArrayReal

!****************************************************************************
!****************************************************************************

  subroutine indexArrayInt(n,Array,index)
    implicit none
    integer, intent(in)  :: Array(n)
    integer, intent(out) :: index(n)
    integer, parameter   :: nn=15, nstack=50
    integer              :: n,k,i,j,indext,jstack,l,r
    integer              :: istack(nstack)
    integer              :: a
    do j = 1,n
      index(j) = j
    end do
    jstack=0
    l=1
    r=n
    do
      if (r-l < nn) then
        do j=l+1,r
          indext=index(j)
          a=Array(indext)
          do i=j-1,l,-1
            if (Array(index(i)) <= a) exit
            index(i+1)=index(i)
          end do
          index(i+1)=indext
        end do
        if (jstack == 0) return
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+r)/2
        call swapInt(index(k),index(l+1))
        call exchangeIndex(index(l),index(r))
        call exchangeIndex(index(l+1),index(r))
        call exchangeIndex(index(l),index(l+1))
        i=l+1
        j=r
        indext=index(l+1)
        a=Array(indext)
        do
          do
            i=i+1
            if (Array(index(i)) >= a) exit
          end do
          do
            j=j-1
            if (Array(index(j)) <= a) exit
          end do
          if (j < i) exit
          call swapInt(index(i),index(j))
        end do
        index(l+1)=index(j)
        index(j)=indext
        jstack=jstack+2
        if (jstack > nstack) then
          write(*,*) 'NSTACK too small in indexArrayInt'
          stop
        end if
        if (r-i+1 >= j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
        else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
        end if
      end if
    end do
  
  contains

    subroutine exchangeIndex(i,j)
      integer, intent(inout) :: i,j
      integer                :: swp
      if (Array(j) < Array(i)) then
        swp=i
        i=j
        j=swp
      end if
    end subroutine exchangeIndex

  end subroutine indexArrayInt

!****************************************************************************
!****************************************************************************

end module SpecialFunctions
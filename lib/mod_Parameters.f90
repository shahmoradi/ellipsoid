module Parameters
  use, intrinsic :: iso_fortran_env, only: RequestedRealKind => REAL64
  implicit none
  
  ! Constants from comutational accuracy
  integer , parameter :: SP = selected_real_kind(6), DP = selected_real_kind(15)     ! single precision, double precision
  integer , parameter :: RK = RequestedRealKind
  integer , parameter :: RKP = precision(1._RK)

  ! Mathematical constants
  real(RK), parameter :: SQRT2 = sqrt(2._RK)                                  ! Square root of 2.
  real(RK), parameter :: SQRTPI = sqrt(acos(-1._RK))                          ! Square root of Pi.
  real(RK), parameter :: SQRT2PI = sqrt(2._RK*acos(-1._RK))                   ! Square root of 2Pi.
  real(RK), parameter :: NAPIER = exp(1._RK)                                  ! Napier number e.
  real(RK), parameter :: LOG10NAPIER = 0.434294481903259_RK                   ! Log10 of Napier constant.
  real(RK), parameter :: LN10 = 2.302585092994046_RK                          ! Log10 of Napier constant.
  real(RK), parameter :: SQRT_HALF_PI = 1.2533141373155_RK                    ! Square root of PI/2. This is needed in the integration of Epk and ELR density functions.
  real(RK), parameter :: HUGE_NUMBER = huge(LN10)                             ! Square root of PI/2. This is needed in the integration of Epk and ELR density functions.
  real(RK), parameter :: PI = acos(-1._RK)                                    ! The irrational number PI.
  real(RK), parameter :: INVSQRT2PI = 0.398942280401432_RK                    ! 1/sqrt(2*Pi)
  real(RK), parameter :: LN2PI = 1.837877066409345_RK                         ! ln(2pi)

  ! Physical constants
  real(RK), parameter :: LIGHT_SPEED = 3.e5_RK,HUBBLE_CONST = 7.1e1_RK        ! LIGHT_SPEED is the speed of light (Km/s), HUBBLE_CONST is the Hubble constant in units of km/s/MPc.
  real(RK), parameter :: LS2HC = LIGHT_SPEED / HUBBLE_CONST                   ! the speed of light in units of km/s divided by the Hubble constant.
  real(RK), parameter :: MPC2CM = 3.09e24_RK                                  ! 1 Mega Parsec = MPC2CM centimeters.
  real(RK), parameter :: LOG10MPC2CMSQ4PI = log10(4._RK*PI) + 2*log10(MPC2CM) ! Log10(MPC2CM centimeters.
  real(RK), parameter :: OMEGA_DE = 0.7_RK                                    ! Dark Energy density.
  real(RK), parameter :: OMEGA_DM = 0.3_RK                                    ! Dark Matter density.

  ! Parameters required for the calculation of Luminosity distance, that are private to this module
  real(RK), parameter, private :: alpha0 = 1._RK + 2._RK * OMEGA_DE / OMEGA_DM
  real(RK), parameter, private :: x0 = log( alpha0 + sqrt( alpha0**2 - 1._RK ) )
  real(RK), parameter, private :: psi0=x0**0.33333333*(1.5874010519682-6.2992105236833e-3_RK*x0**2+7.5375168659459e-5_RK*x0**4)
  real(RK), parameter, private :: lumDisFactor = LS2HC / ( OMEGA_DE**0.16666666 * OMEGA_DM**0.33333333 )
  
contains

  !*********************************************************************
  !*********************************************************************
  
  !    This function calculates the luminosity distance in units of Mpc.
  !    The input parameters is the redshift (z), while the required parameters are taken from the two modules constants:
  !        LS2HC: which is the speed of light in units of km/s divided by the Hubble Constant in units of km/s/MPc.
  !    Also from module COSMOparameters:
  !        OMEGA_DE: the dark energy density,
  !        OMEGA_DM: the dark matter density.
  !    For input redshifts z>zlim=0.1, the luminosity distance will be calculated according to approximation algorithm of Wickramasinghe & Okwatta (2010).
  !    Note that for redshifts less than 0.1, the error in the calculated luminosity distance grows to more than 0.001.
  !    This algorithm should therefore not be used for z<0.1, and only when the desired accuracy in Lum. dis. is 0.001.
  !    Amir Shahmoradi, Wednesday July 14, 2016, 9:59 PM, ICES, The University of Texas at Austin.
  real(RK) pure function getLumDisWickram(z)
    implicit none
    real(RK), intent(in) :: z
    real(RK)             :: alpha,x,psi
    alpha = 1._RK + 2._RK * OMEGA_DE / ( OMEGA_DM * (1._RK+z)**3 )
    x     = log( alpha  + sqrt( alpha**2 - 1._RK ) )
    psi   = x**0.33333333  * ( 1.5874010519682 - 6.2992105236833e-3_RK * x**2  + 7.5375168659459e-5_RK * x**4  )
    getLumDisWickram = lumDisFactor * ( 1._RK + z ) * ( psi0 - psi )
  end function getLumDisWickram

  !*********************************************************************
  !*********************************************************************  
  
end module Parameters

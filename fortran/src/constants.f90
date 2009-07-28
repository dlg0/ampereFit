module constants
  implicit none

  integer, parameter :: SGL = selected_real_kind ( p = 6, r = 37 )
  integer, parameter :: DBL = selected_real_kind ( p = 13, r = 200 )

  real(kind=DBL), parameter :: rE = 6356.75e3
  real(kind=DBL), parameter :: rI = 1.02
  real(kind=DBL), parameter :: pi = 3.141593
  real(kind=DBL), parameter :: degToRad = pi / 180.0
  real(kind=DBL), parameter :: radToDeg = 180.0 / pi
  real(kind=DBL), parameter :: u0_ = 1.2566e-6;
  real(kind=DBL), parameter :: e0_ = 8.8541878e-12;
  real(kind=DBL), parameter :: c_  = 3e8;
  real(kind=DBL), parameter :: u0  = u0_ / rE;
  real(kind=DBL), parameter :: e0  = e0_ * rE ** 3;
  real(kind=DBL), parameter :: c = c_ / rE;

end module constants



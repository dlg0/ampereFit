module constants

    implicit none

    integer, parameter :: SGL = selected_real_kind ( p = 6, r = 37 )
    integer, parameter :: DBL = selected_real_kind ( p = 13, r = 200 )

    real(kind=DBL), parameter :: rE = 6356.75e3
    real(kind=DBL), parameter :: rSat = 780.0e3
    real(kind=DBL), parameter :: rI = 1.02
    real(kind=DBL), parameter :: pi = 3.1415926535897932384626
    real(kind=DBL), parameter :: degToRad = pi / 180d0
    real(kind=DBL), parameter :: radToDeg = 180d0 / pi
    real(kind=DBL), parameter :: u0_ = 1.2566e-6;
    real(kind=DBL), parameter :: e0_ = 8.8541878e-12;
    real(kind=DBL), parameter :: c_  = 3e8;
    real(kind=DBL), parameter :: u0  = u0_ / rE;
    real(kind=DBL), parameter :: e0  = e0_ * rE ** 3;
    real(kind=DBL), parameter :: c = c_ / rE;

    type :: ampBasis

        real(DBL) :: PLM, dPLM
        real(DBL) :: Y, br, bTh, bPh

    end type ampBasis

    type :: ampData
        real(kind=DBL) :: utc           ! UT time in dec hours
        integer :: isat                 ! coded SV number (for Haje)
        integer :: iPln                 ! orbit track number (0->5)
        real :: qual                    ! data quality from Lars
        integer :: splice               ! flag for spliced data - where missing data estimated
        integer :: typ                  ! data type (for CLW) - This is NOT a helpful description !!
        real(kind=DBL) :: x, y, z
        real(kind=DBL) :: bX, bY, bZ
        real(kind=DBL) :: &
            R, & ! [km]
            T, & ! [coLat, rad]
            P, & ! [lon, rad]
            bR, bT, bP
        real(kind=DBL) :: jPar ! [uAm^{-2}]

    end type ampData

end module constants



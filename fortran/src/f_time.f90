module f_time

use constants, only : SGL, DBL
use ISO_C_BINDING, only : C_INT, C_DOUBLE

interface

        integer(C_INT) function f_TimeJulianToYMDHMS(jd,yr,mo,dy,hr,mt,sc) bind(C,name='TimeJulianToYMDHMS')
            use ISO_C_BINDING, only: C_INT, C_DOUBLE, C_PTR
            implicit none
            type(C_PTR), value, intent(in) :: yr,mo,dy,hr,mt,sc
            real(C_DOUBLE), value, intent(in) :: jd
        end function f_TimeJulianToYMDHMS

        subroutine f_TimeEpochToYMDHMS(tme,yr,mo,dy,hr,mn,sc) bind(C,name='TimeEpochToYMDHMS')
            use ISO_C_BINDING, only: C_INT, C_DOUBLE, C_PTR
            implicit none
            real(C_DOUBLE), value, intent(in) :: tme
            type(C_PTR), value, intent(in) :: yr,mo,dy,hr,mn,sc
        end subroutine f_TimeEpochToYMDHMS

        real(C_DOUBLE) function f_TimeYMDHMSToEpoch(yr,mo,dy,hr,mt,sc) bind(C,name='TimeYMDHMSToEpoch')
            use ISO_C_BINDING, only: C_INT, C_DOUBLE
            implicit none
            integer(C_INT), value, intent(in) :: yr,mo,dy,hr,mt
            real(C_DOUBLE), value, intent(in) :: sc
        end function f_TimeYMDHMSToEpoch

        real(C_DOUBLE) function f_TimeYMDHMSToJulian(yr,mo,dy,hr,mt,sc) bind(C,name='TimeYMDHMSToJulian')
            use ISO_C_BINDING, only: C_INT, C_DOUBLE
            implicit none
            integer(C_INT), value, intent(in) :: yr,mo,dy,hr,mt
            real(C_DOUBLE), value, intent(in) :: sc
        end function f_TimeYMDHMSToJulian

end interface

end module f_time

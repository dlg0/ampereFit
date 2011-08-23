module f_aacgm

use ISO_C_BINDING, only : C_INT, C_DOUBLE

interface

        real(C_DOUBLE) function f_mltconvertyrsec (yr, yrsec, mlon) bind(C, name='MLTConvertYrsec')
                use ISO_C_BINDING, only : C_INT, C_DOUBLE, C_PTR
                implicit none
                integer(C_INT), value, intent(IN) :: yr, yrsec
                real(C_DOUBLE), value :: mlon
        end function f_mltconvertyrsec


        integer(C_INT) function f_aacgmconvert (iLat,iLon,hgt,oLat,oLon,r,flg) bind(C, name='AACGMConvert')
                use ISO_C_BINDING, only : C_INT, C_DOUBLE, C_PTR
                implicit none
                real(C_DOUBLE), value, intent(IN) :: iLat, iLon, hgt
                type(C_PTR), value :: oLat, oLon, r
                integer(C_INT), value :: flg
        end function f_aacgmconvert

        integer(C_INT) function f_aacgminit (yr) bind(C, name='AACGMInit')
                use ISO_C_BINDING, only : C_INT, C_DOUBLE, C_PTR
                implicit none
                integer(C_INT), value :: yr 
        end function f_aacgminit


end interface

end module f_aacgm

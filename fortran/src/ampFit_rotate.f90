module ampFit_rotate

use constants
use dlg

implicit none

contains

subroutine XYZ_from_SPH ( dataIn )

    implicit none
    type(ampData), intent(inout) :: dataIn(:)
    
    integer :: j, i

    write(*,*) 'Generating XYZ coords from SPH coords ...'

    j = 1

    coords_XYZ_from_SPH: &
    do i=1,size(dataIn)

        ! Double precision GEOPACK only
        call sphcar_08 ( &
                    dataIn(i)%R, &
                    dataIn(i)%T, &
                    dataIn(i)%P, &
                    dataIn(i)%x, &
                    dataIn(i)%y, &
                    dataIn(i)%z, &
                    j )

    enddo coords_XYZ_from_SPH

    write(*,*) 'DONE'

end subroutine XYZ_from_SPH

subroutine XYZ_to_SPH ( dataIn )

    implicit none

    type(ampData), intent(inout) :: dataIn(:)
    integer :: i,j

    ! Rotate XYZ GEI db to spherical GEI
    ! ----------------------------------

    write(*,*) '    Rotating GEI b vectors from XYZ to SPH ...'

    vectors_GEI_XYZ_to_SPH: &
    do i=1,size(dataIn)
        
        ! Double precision GEOPACK only
        call bcarsp_08 ( &
                    dataIn(i)%x, &
                    dataIn(i)%y, &
                    dataIn(i)%z, &
                    dataIn(i)%bx, &
                    dataIn(i)%by, &
                    dataIn(i)%bz, &
                    dataIn(i)%bR, &
                    dataIn(i)%bT, &
                    dataIn(i)%bP )

    enddo vectors_GEI_XYZ_to_SPH

    write(*,*) '    DONE'


    ! Get spherical coords of the GEI XYZ locations
    ! -----------------------------------------------

    write(*,*) '    Calculating GEI SPH coords from XYZ ...'

    j = -1

    coords_GEI_XYZ_to_SPH: &
    do i=1,size(dataIn)

        ! Double precision GEOPACK only
        call sphcar_08 ( &
                    dataIn(i)%R, &
                    dataIn(i)%T, &
                    dataIn(i)%P, &
                    dataIn(i)%x, &
                    dataIn(i)%y, &
                    dataIn(i)%z, &
                    j )

    end do coords_GEI_XYZ_to_SPH

    write(*,*) '    DONE'

    write(*,*) 'DONE'

end subroutine XYZ_to_SPH

end module ampFit_rotate

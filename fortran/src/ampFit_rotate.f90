module ampFit_rotate

use constants
use dlg

implicit none

contains

subroutine XYZ_from_SPH ( dataIn )

    implicit none
    type(ampData), intent(inout) :: dataIn(:)
    
    integer :: j, i

    write(*,*) '    Updating XYZ coords from SPH coords ...'

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

    write(*,*) '    DONE'

end subroutine XYZ_from_SPH

subroutine XYZ_to_SPH ( dataIn )

    implicit none

    type(ampData), intent(inout) :: dataIn(:)
    integer :: i,j

    ! Rotate XYZ GEI db to spherical GEI
    ! ----------------------------------

    write(*,*) '    Updating SPH vecs and coords from XYZ...'

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


    ! Get spherical coords of the GEI XYZ locations
    ! -----------------------------------------------

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

end subroutine XYZ_to_SPH

subroutine SPH_to_XYZ ( dataIn )

    implicit none

    type(ampData), intent(inout) :: dataIn(:)
    integer :: i,j

    ! Rotate XYZ GEI db to spherical GEI
    ! ----------------------------------

    write(*,*) '    Updating XYZ vecs and coords from SPH ...'

    vectors_GEI_XYZ_to_SPH: &
    do i=1,size(dataIn)
        
        ! Double precision GEOPACK only
        call bspcar_08 ( &
                    dataIn(i)%T, &
                    dataIn(i)%P, &
                    dataIn(i)%bR, &
                    dataIn(i)%bT, &
                    dataIn(i)%bP, &
                    dataIn(i)%bX, &
                    dataIn(i)%bY, &
                    dataIn(i)%bZ )

    enddo vectors_GEI_XYZ_to_SPH


    ! Get spherical coords of the GEI XYZ locations
    ! -----------------------------------------------

    j = 1

    coords_GEI_SPH_to_XYZ: &
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

    end do coords_GEI_SPH_to_XYZ

    write(*,*) '    DONE'

end subroutine SPH_to_XYZ


end module ampFit_rotate

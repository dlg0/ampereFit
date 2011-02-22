module ampFit_rotate

use constants

implicit none

contains

subroutine XYZ_to_SPH ( dataIn )

    implicit none

    type(ampData), intent(inout) :: dataIn(:)
    integer :: i,j, nObs

    ! Get spherical coords of the GEI XYZ locations
    ! ---------------------------------------------

    write(*,*) '    Calculating GEI SPH coords from XYZ ...'

    nObs    = size ( dataIn )
    
    write(*,*) '        for this many pts: ', nObs

    j = -1

    coords_GEI_XYZ_to_SPH: &
    do i=1,nObs

        ! Double precision GEOPACK only
        call sphcar_08 ( &
                    dataIn(i)%GEI_R_km, &
                    dataIn(i)%GEI_coLat_rad, &
                    dataIn(i)%GEI_lon_rad, &
                    dataIn(i)%px, &
                    dataIn(i)%py, &
                    dataIn(i)%pz, &
                    j )

    end do coords_GEI_XYZ_to_SPH

    dataIn%GEI_coLat_deg = dataIn%GEI_coLat_rad * radToDeg
    dataIn%GEI_lon_deg = dataIn%GEI_lon_rad * radToDeg 

    write(*,*) '    DONE'

    ! Rotate XYZ GEI db to spherical GEI
    ! ----------------------------------

    write(*,*) '    Rotating GEI b vectors from XYZ to SPH ...'

    vectors_GEI_XYZ_to_SPH: &
    do i=1,nObs
        
        ! Double precision GEOPACK only
        call bcarsp_08 ( &
                    dataIn(i)%px, &
                    dataIn(i)%py, &
                    dataIn(i)%pz, &
                    dataIn(i)%dbx, &
                    dataIn(i)%dby, &
                    dataIn(i)%dbz, &
                    dataIn(i)%br_GEI, &
                    dataIn(i)%bTheta_GEI, &
                    dataIn(i)%bPhi_GEI )

    enddo vectors_GEI_XYZ_to_SPH

    write(*,*) '    DONE'

    write(*,*) 'DONE'

end subroutine XYZ_to_SPH

end module ampFit_rotate

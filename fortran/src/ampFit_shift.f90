module ampFit_shift

use constants

real(kind=DBL) :: R_intersection_GEI_km, &
        T_intersection_GEI_coLat_rad, &
        P_intersection_GEI_rad
real(kind=DBL) :: x_intersection, y_intersection, z_intersection

real, target :: shift_rot_mat(3,3)

contains

subroutine calculate_intersection_point ( dataIn )

    use la_precision, only: WP=>DP
    use f95_lapack, only: la_gelss, la_lange

    implicit none

    type(ampData), intent(in) :: dataIn(:)

    integer :: nTrackPts, p, i, M, N, pp, cnt, j
    real(kind=WP), allocatable :: la_A(:,:), la_B(:), la_S(:)
    integer, allocatable :: iiSubSet(:)
    logical, allocatable :: mask(:)
    real(kind=DBL), allocatable :: parabola_coeffs(:,:)
    integer :: la_rank, la_info
    real(kind=DBL) :: AA, BB, CC, x1, x2, y1, y2
!    real(kind=DBL) :: x_intersection_points(0:5,5), y_intersection_points(0:5,5), &
    real(kind=DBL) :: x_intersection_points(15), y_intersection_points(15), &
        distance1, distance2, r_intersect
    integer :: nCloseDataPts

    write(*,*) 'Calculating intersection point ...'

    allocate ( mask ( size(dataIn) ) )
    allocate ( parabola_coeffs(0:5,3) )

    fit_2d_poly_to_track: &
    do p=0,5

        mask    = dataIn%iPln==p .and. dataIn%T*radToDeg < 35
        nTrackPts    = count ( mask )
        M = nTrackPts 
        N = 3

        !write(*,*) nTrackPts

        allocate ( iiSubSet(nTrackPts) )
        allocate ( la_A(M,N), la_B(M), la_S(N) )

        iiSubSet    = pack ( (/ (i, i=1, size(dataIn)) /), mask = mask )  

        ! y = Ax^2 + Bx + C

        la_A(:,1)   = dataIn(iiSubSet)%x**2
        la_A(:,2)   = dataIn(iiSubSet)%x
        la_A(:,3)   = 1

        la_B    = dataIn(iiSubSet)%y

        !write(*,*) '    Calling la_gelss ...'
        call la_gelss ( la_A, la_B, la_rank, la_S, info = la_info ) 
    
        !write(*,*) '    Rank: ', la_rank
        !write(*,*) '    Info: ', la_info
        !write(*,*) '    A: ', la_B(1)
        !write(*,*) '    B: ', la_B(2)
        !write(*,*) '    C: ', la_B(3)

        parabola_coeffs(p,:) = la_B(1:N)

        deallocate (iiSubSet, la_A, la_B, la_S )

    enddo fit_2d_poly_to_track

    cnt=1
    find_intersection_points: &
    do p=0,4

!        cnt = 1
        all_other_tracks: &
        do pp=p+1,5

!             if ( p /= pp ) then
!             if (pp .gt. p) then

!    For each track pair (15 of them) solve the track parabola eqns as simul.
!    eqns

                AA  = parabola_coeffs(p,1) - parabola_coeffs(pp,1)
                BB  = parabola_coeffs(p,2) - parabola_coeffs(pp,2)
                CC  = parabola_coeffs(p,3) - parabola_coeffs(pp,3)

                x1  = ( -BB + sqrt ( BB**2 - 4 * AA * CC ) ) / ( 2 * AA )
                x2  = ( -BB - sqrt ( BB**2 - 4 * AA * CC ) ) / ( 2 * AA )

                y1  = parabola_coeffs(p,1)*x1**2 + parabola_coeffs(p,2)*x1 + parabola_coeffs(p,3)
                y2  = parabola_coeffs(p,1)*x2**2 + parabola_coeffs(p,2)*x2 + parabola_coeffs(p,3)

                !write(*,*) '        (x1,y1): ', x1,y1
                !write(*,*) '        (x2,y2): ', x2,y2

                ! Select out solution closest to pole

                distance1   = sqrt ( x1**2 + y1**2 )
                distance2   = sqrt ( x2**2 + y2**2 )
 
                if (distance1<distance2) then 

!                    x_intersection_points(p,cnt) = x1
!                    y_intersection_points(p,cnt) = y1
                    x_intersection_points(cnt) = x1
                    y_intersection_points(cnt) = y1

                else

!                    x_intersection_points(p,cnt) = x2
!                    y_intersection_points(p,cnt) = y2
                    x_intersection_points(cnt) = x2
                    y_intersection_points(cnt) = y2

                endif

                cnt = cnt + 1
 
!            endif

        enddo all_other_tracks

    enddo find_intersection_points

!    x_intersection  = sum ( x_intersection_points ) / 30d0
!    y_intersection  = sum ( y_intersection_points ) / 30d0
    x_intersection  = sum ( x_intersection_points ) / 15d0
    y_intersection  = sum ( y_intersection_points ) / 15d0
!    print*,'cnt,xint,yint=',cnt,x_intersection,y_intersection

    nCloseDataPts = 0
    cnt = 0
    do while ( nCloseDataPts < 20 ) 

        mask    = sqrt((x_intersection-dataIn%x)**2+(y_intersection-dataIn%x)**2) < (10 + cnt * 10) &
                    .and. dataIn%T*radToDeg < 35
        nCloseDataPts   = count ( mask )

        cnt = cnt + 1 

    enddo

    allocate ( iiSubSet(nCloseDataPts) )
    iiSubSet    = pack ( (/ (i, i=1, size(dataIn)) /), mask = mask )  

    r_intersect=sum (dataIn(iiSubSet)%R ) / nCloseDataPts
    z_intersection = sqrt(r_intersect**2 - x_intersection**2-y_intersection**2)
! Set for Nth hemis (+z) at the moment)
!    z_intersection  = sum ( dataIn(iiSubSet)%pz ) / nCloseDataPts

    j = -1

    ! Double precision GEOPACK only
    call sphcar_08 ( &
                R_intersection_GEI_km, &
                T_intersection_GEI_coLat_rad, &
                P_intersection_GEI_rad, &
                x_intersection, &
                y_intersection, &
                z_intersection, &
                j )


    !write(*,*) dataIn(iiSubSet)%pz

    write(*,*) '    X intersection: ', x_intersection
    write(*,*) '    Y intersection: ', y_intersection
    write(*,*) '    Z intersection: ', z_intersection

    write(*,*) '    R intersection: ', R_intersection_GEI_km
    write(*,*) '    T intersection: ', T_intersection_GEI_coLat_rad * radToDeg
    write(*,*) '    P intersection: ', P_intersection_GEI_rad * radToDeg

    write(*,*) 'DONE' 

end subroutine calculate_intersection_point 


subroutine calculate_shift_rotation_matrix ()

    implicit none

    real :: u_x, u_y, u_z, uMag, angle
    real, pointer :: R(:,:)
    real :: detR

    write(*,*) 'Calculating rotation matrix ...'

    ! Calculate a rotation matrix that shifts the z-axis to the
    ! vector [x,y,z] where x,y,z are the x,y,z coords of the intersection 
    ! point calculated by the calculate_intersection_point subroutine.

    ! The rotation will be around an axis perpendicular to the
    ! plane defined by the z-axis and intersection point vector.
    ! i.e., the cross product of those two vectors gives the 
    ! axis of rotation (u_x,u_y,u_z) ...

! cross product (0,0,r) x (x, y, z)
    u_x = 0 * z_intersection - 1 * y_intersection
    u_y = -1 * ( 0 * z_intersection - 1 * x_intersection )
    u_z = 0 * y_intersection - 0 * x_intersection

    uMag    = sqrt(u_x**2+u_y**2+u_z**2)

    u_x = u_x / uMag
    u_y = u_y / uMag
    u_z = u_z / uMag

    angle   = -atan2 ( sqrt ( x_intersection**2 + y_intersection**2 ), z_intersection )
    
    write(*,*) '    Rotation angle: ', angle * radToDeg

    ! See http://en.wikipedia.org/wiki/Rotation_matrix

    shift_rot_mat(1,1)  = cos ( angle ) + u_x**2 * ( 1 - cos (angle) )
    shift_rot_mat(2,1)  = u_y * u_x * (1-cos(angle)) + u_z*sin(angle)
    shift_rot_mat(3,1)  = u_z*u_x*(1-cos(angle)) - u_y*sin(angle)

    shift_rot_mat(1,2)  = u_x * u_y * (1-cos(angle)) - u_z * sin(angle)
    shift_rot_mat(2,2)  = cos(angle) + u_y**2*(1-cos(angle))
    shift_rot_mat(3,2)  = u_z*u_y*(1-cos(angle)) + u_x*sin(angle)

    shift_rot_mat(1,3)  = u_x*u_z*(1-cos(angle))+u_y*sin(angle)
    shift_rot_mat(2,3)  = u_y*u_z*(1-cos(angle))-u_x*sin(angle)
    shift_rot_mat(3,3)  = cos(angle)+u_z**2*(1-cos(angle))

    write(*,*) '    Rotation matrix: '
    write(*,*) '        ',  shift_rot_mat(1,:)
    write(*,*) '        ',  shift_rot_mat(2,:)
    write(*,*) '        ',  shift_rot_mat(3,:)

    !shift_rot_mat   = 0
    !shift_rot_mat(1,1)  = 1
    !shift_rot_mat(2,2)  = 1
    !shift_rot_mat(3,3)  = 1

    R => shift_rot_mat

    ! Check determinant == 1

!    detR    = R(1,1)*R(2,2)*R(3,3) &
!        + R(1,2)*R(2,3)*R(3,1) &
!        + R(1,3)*R(2,1)*R(3,2) &
!        - R(1,1)*R(2,3)*R(3,2) &
!        - R(1,2)*R(2,1)*R(3,3) &
!        - R(1,3)*R(2,2)*R(3,1)

!    write(*,*) '    Det of rotation matrix: ', detR

    write(*,*) 'DONE'

end subroutine calculate_shift_rotation_matrix


subroutine create_dataShifted ( dataIn, dataOut )

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    type(ampData), intent(inout) :: dataOut(:)
    integer :: i
    real :: pos_vec(3), shifted_pos_vec(3)
    real :: b_vec(3), shifted_b_vec(3)

    write(*,*) 'Shifting data ...'

    write(*,*) '    This many pts: ', size(dataIn)

    dataOut = dataIn

    shift_coordinates_and_vectors: &
    do i=1,size(dataIn)

        pos_vec(1)  = dataIn(i)%x
        pos_vec(2)  = dataIn(i)%y
        pos_vec(3)  = dataIn(i)%z

        b_vec(1)  = dataIn(i)%bx
        b_vec(2)  = dataIn(i)%by
        b_vec(3)  = dataIn(i)%bz

        shifted_pos_vec = matMul ( shift_rot_mat, pos_vec )
        shifted_b_vec = matMul ( shift_rot_mat, b_vec )

        dataOut(i)%x  = shifted_pos_vec(1)
        dataOut(i)%y  = shifted_pos_vec(2)
        dataOut(i)%z  = shifted_pos_vec(3)

        dataOut(i)%bx  = shifted_b_vec(1)
        dataOut(i)%by  = shifted_b_vec(2)
        dataOut(i)%bz  = shifted_b_vec(3)

    enddo shift_coordinates_and_vectors

    write(*,*) 'DONE'

end subroutine create_dataShifted

end module ampFit_shift

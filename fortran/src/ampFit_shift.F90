module ampFit_shift
!
! Routines for shifting the dB data
! C.L. Waters and D.L. Green
! University of Newcastle
! Australia
!
! Nov 2011
!
! SUBROUTINES:
! calculate_intersection_point
! calculate_shift_rotation_matrix
! ShiftData
! UnShiftData
!
! Ver : 201202
!

  use constants
  use ampFit_rotate

  real(kind=DBL) :: R_intersection_km, &
                    T_intersection_coLat_rad, &
                    P_intersection_rad
  real(kind=DBL) :: x_intersection, y_intersection, z_intersection

  real(kind=DBL), target :: shift_rot_mat(3,3)

  contains
! calculate_intersection_point
! calculate_shift_rotation_matrix
! ShiftData
! UnShiftData

  subroutine calculate_intersection_point ( dataIn, south )

#ifdef _intel_
    use mkl95_precision, only: WP=>DP
    use mkl95_lapack, only: la_gelss=>gelss
#else
    use la_precision, only: WP=>DP
    use f95_lapack, only: la_gelss
#endif

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    integer, intent(in) :: south

    integer :: nTrackPts, p, i, M, N, pp, cnt, j
    real(kind=WP), allocatable :: la_A(:,:), la_B(:)
    real(kind=WP), allocatable :: la_S(:)
    integer, allocatable :: iiSubSet(:)
    logical, allocatable :: mask(:)
    real(kind=DBL), allocatable :: parabola_coeffs(:,:)
    integer :: la_rank, la_info
    real(kind=WP) :: la_rcond
    real(kind=DBL) :: AA, BB, CC, x1, x2, y1, y2
    real(kind=DBL) :: x_intersection_points(15), y_intersection_points(15), &
                      distance1, distance2, r_intersect
    integer :: nCloseDataPts

    allocate ( mask ( size(dataIn) ) )
    allocate ( parabola_coeffs(0:5,3) )

    fit_2d_poly_to_track: &
    do p=0,5                   ! for each Iridium track

      if (south .eq. 1) then
!        mask = dataIn%iPln==p .and. dataIn%T*radToDeg > 145
        mask = dataIn%iPln==p .and. dataIn%T*radToDeg > 95.0   ! same as IDL
      else
!        mask = dataIn%iPln==p .and. dataIn%T*radToDeg < 35
        mask = dataIn%iPln==p .and. dataIn%T*radToDeg < 85.0
      end if
      nTrackPts    = count ( mask )
      M = nTrackPts 
      N = 3

      if (nTrackPts .lt. 20) then
        err_stat = 8
        err_msg = 'Insufficient data in Calc_Intersection'
        write(*,*) err_msg
        print*,'Track Num:',p
        return
      endif

      allocate ( iiSubSet(nTrackPts) )
      allocate ( la_A(M,N), la_B(M), la_S(N) )

      iiSubSet    = pack ( (/ (i, i=1, size(dataIn)) /), mask = mask )  

! y = Ax^2 + Bx + C
      la_A(:,1)   = dataIn(iiSubSet)%x**2
      la_A(:,2)   = dataIn(iiSubSet)%x
      la_A(:,3)   = 1

      la_B    = dataIn(iiSubSet)%y

!write(*,*) '    Calling la_gelss ...'
      call la_gelss ( la_A, la_B, la_rank, la_S, 0.00001_WP, la_info ) 
    
      if (la_info .ne. 0) then
        err_stat = 9
        err_msg = 'Error in la_gelss in Calc_Intersection'
        write(*,*) err_msg
        print*,'la_info=',la_info
        return
      endif

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

      all_other_tracks: &
      do pp=p+1,5

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

          x_intersection_points(cnt) = x1
          y_intersection_points(cnt) = y1

        else

          x_intersection_points(cnt) = x2
          y_intersection_points(cnt) = y2

        endif

        cnt = cnt + 1
 
      enddo all_other_tracks

    enddo find_intersection_points

    x_intersection  = sum ( x_intersection_points ) / 15d0
    y_intersection  = sum ( y_intersection_points ) / 15d0

!    r_intersect = sum (dataIn(iiSubSet)%R ) / nTrackPts
    r_intersect = (rE + rSat)/1000.0d0

    z_intersection = sqrt(r_intersect**2 - &
       x_intersection**2 - y_intersection**2)
    if (south .eq. 1) z_intersection = -z_intersection
!
! calc x,y,z of intersection point
!    j = -1
!    z_intersection = (rE + rSat)/1.0d3  ! km
!    if (south .eq. 1) z_intersection = -z_intersection

! Double precision GEOPACK, j=-1 -> x,y,z => r,t,p
!    call sphcar_08 ( &
!              R_intersection_km, &
!              T_intersection_coLat_rad, &
!              P_intersection_rad, &
!              x_intersection, &
!              y_intersection, &
!              z_intersection, &
!              j )


!write(*,*) dataIn(iiSubSet)%pz

    write(*,*) '    X intersection: ', x_intersection
    write(*,*) '    Y intersection: ', y_intersection
    write(*,*) '    Z intersection: ', z_intersection
!
!    write(*,*) '    R intersection: ', R_intersection_km
!    write(*,*) '    T intersection: ', T_intersection_coLat_rad * radToDeg
!    write(*,*) '    P intersection: ', P_intersection_rad * radToDeg

  end subroutine calculate_intersection_point 
!
! ------------------------------------------------------------------
!
  subroutine calculate_shift_rotation_matrix (south)

    implicit none

    integer, intent(in) :: south
    real(kind=DBL), dimension(3,3) :: p_a, i_a, q_a
    real(kind=DBL) :: u_x, u_y, u_z, uMag, angle, r0, sh_mag
    
!    real, pointer :: R(:,:)
!    real :: detR

    ! Calculate a rotation matrix that shifts the z-axis to the
    ! vector [x,y,z] where x,y,z are the x,y,z coords of the intersection 
    ! point calculated by the calculate_intersection_point subroutine.

    ! The rotation will be around an axis perpendicular to the
    ! plane defined by the z-axis and intersection point vector.
    ! i.e., the cross product of those two vectors gives the 
    ! axis of rotation (u_x,u_y,u_z) ...a
!
! C.L. Waters and D.L. Green
!
! Last modified:
! 21 Feb 2012 - changed to matrix mult to do both Nth and Sth hemis data -CLW
!
! Ver : 201202
!

! cross product (x,y,z) x (0,0,r)
    r0 = (rE+rSat)/1000.0d0
    if (south .eq. 1) r0 = -r0

    u_x =  y_intersection*r0
    u_y = -x_intersection*r0
    u_z =  0.0d0

    uMag = sqrt(u_x**2+u_y**2+u_z**2)

    u_x = u_x / uMag
    u_y = u_y / uMag
    u_z = u_z / uMag

    sh_mag = sqrt(x_intersection**2+y_intersection**2+z_intersection**2)
    angle   = acos ( r0*z_intersection/(sh_mag*abs(r0)))
    
!    write(*,*) '    Rotation angle: ', angle * radToDeg

! See http://en.wikipedia.org/wiki/Rotation_matrix

! Identity matrix
    i_a = reshape( (/1,0,0,0,1,0,0,0,1 /), (/3,3/) )

! cos(angle) matrix
    p_a(1,1) = u_x*u_x    ! row, column
    p_a(1,2) = u_y*u_x
    p_a(1,3) = u_z*u_x

    p_a(2,1) = u_x*u_y    ! row, column
    p_a(2,2) = u_y*u_y
    p_a(2,3) = u_z*u_y

    p_a(3,1) = u_x*u_z    ! row, column
    p_a(3,2) = u_y*u_z
    p_a(3,3) = u_z*u_z

! sin(angle) matrix
    q_a(1,1) = 0.0d0
    q_a(1,2) = -u_z
    q_a(1,3) = u_y

    q_a(2,1) = u_z
    q_a(2,2) = 0.0d0
    q_a(2,3) = -u_x

    q_a(3,1) = -u_y
    q_a(3,2) = u_x
    q_a(3,3) = 0.0d0

    shift_rot_mat = p_a + (i_a-p_a)*cos(angle) + q_a * sin(angle)

!    write(*,*) '    Rotation matrix: '
!    write(*,*) '        ',  shift_rot_mat(1,:)
!    write(*,*) '        ',  shift_rot_mat(2,:)
!    write(*,*) '        ',  shift_rot_mat(3,:)

    !shift_rot_mat   = 0
    !shift_rot_mat(1,1)  = 1
    !shift_rot_mat(2,2)  = 1
    !shift_rot_mat(3,3)  = 1

!    R => shift_rot_mat

    ! Check determinant == 1

!    detR    = R(1,1)*R(2,2)*R(3,3) &
!        + R(1,2)*R(2,3)*R(3,1) &
!        + R(1,3)*R(2,1)*R(3,2) &
!        - R(1,1)*R(2,3)*R(3,2) &
!        - R(1,2)*R(2,1)*R(3,3) &
!        - R(1,3)*R(2,2)*R(3,1)

!    write(*,*) '    Det of rotation matrix: ', detR
  end subroutine calculate_shift_rotation_matrix
!
! ---------------------------------------------------------------------------
!
  subroutine ShiftData ( dataIn, dataOut )

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    type(ampData), intent(out), allocatable :: dataOut(:)
    integer :: i, j
    real(kind=DBL) :: pos_vec(3), shifted_pos_vec(3)
    real(kind=DBL) :: b_vec(3), shifted_b_vec(3)

    allocate (dataOut(size(dataIn)) )
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

! print*,'x,y,z:',dataIn(i)%x,dataIn(i)%y,dataIn(i)%z
! print*,'xs,ys,zs:',dataOut(i)%x,dataOut(i)%y,dataOut(i)%z
!
! if (i.eq.10) stop

    enddo shift_coordinates_and_vectors

    j  = -1   ! j<0 is xyz->rtp
    call rtp_xyz_coord ( dataOut, j )
    call xyz_to_rtp_vec ( dataOut )

  end subroutine ShiftData
!
! --------------------------------------------------------------------
!
  subroutine UnShiftData ( dataIn, dataOut )

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    type(ampData), intent(inout) :: dataOut(:)
    integer :: i, j
    real :: pos_vec(3), shifted_pos_vec(3)
    real :: b_vec(3), shifted_b_vec(3)

    dataOut = dataIn

    shift_coordinates_and_vectors: &
    do i=1,size(dataIn)

        pos_vec(1)  = dataIn(i)%x
        pos_vec(2)  = dataIn(i)%y
        pos_vec(3)  = dataIn(i)%z

        b_vec(1)  = dataIn(i)%bx
        b_vec(2)  = dataIn(i)%by
        b_vec(3)  = dataIn(i)%bz

        shifted_pos_vec = matMul ( transpose(shift_rot_mat), pos_vec )
        shifted_b_vec = matMul ( transpose(shift_rot_mat), b_vec )

        dataOut(i)%x  = shifted_pos_vec(1)
        dataOut(i)%y  = shifted_pos_vec(2)
        dataOut(i)%z  = shifted_pos_vec(3)

        dataOut(i)%bx  = shifted_b_vec(1)
        dataOut(i)%by  = shifted_b_vec(2)
        dataOut(i)%bz  = shifted_b_vec(3)

    enddo shift_coordinates_and_vectors
    
    j  = -1
    call rtp_xyz_coord ( dataOut, j )
    call xyz_to_rtp_vec ( dataOut )

  end subroutine UnShiftData 

end module ampFit_shift
!
! ==============================================================================

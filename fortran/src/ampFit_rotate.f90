module ampFit_rotate
!
! Coordinate and Vector conversion: XYZ <-> RTP
!
! rtp_xyz_coord
! rtp_to_xyz_vec
! xyz_to_rtp_vec
!
! C.L. Waters, D.L. Green
! University of Newcastle, Australia
! Nov, 2011
!
! Ver : 201202
!
  use constants

  implicit none

  contains

subroutine rtp_xyz_coord ( dataIn, j )

  implicit none

  type(ampData), intent(inout) :: dataIn(:)

! Switch, j>0 is rtp->xyz
!         j<0 is xyz->rtp
  integer, intent(in) :: j
    
  integer :: i

  coords_rtp_xyz: &
  do i=1,size(dataIn)

! Double precision GEOPACK only
    call sphcar_08 ( &
            dataIn(i)%R, &
            dataIn(i)%T, &
            dataIn(i)%P, &
            dataIn(i)%x, &
            dataIn(i)%y, &
            dataIn(i)%z, j ) 

  enddo coords_rtp_xyz

end subroutine rtp_xyz_coord
!
! ----------------------------------------------------------------------------
!
subroutine rtp_to_xyz_vec ( dataIn )
! Convert spherical vector components to Cartesian (xyz)

  implicit none

  type(ampData), intent(inout) :: dataIn(:)
  integer :: i

  vectors_rtp_to_xyz: &
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

  enddo vectors_rtp_to_xyz

end subroutine rtp_to_xyz_vec
!
! -----------------------------------------------------------------------------
!
subroutine xyz_to_rtp_vec (dataIn )
! Convert XYZ vector components to spherical ones
  
  implicit none

  type(ampData), intent(inout) :: dataIn(:)
  integer :: i

  vectors_xyz_to_rtp: &
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

  enddo vectors_xyz_to_rtp

end subroutine xyz_to_rtp_vec

end module ampFit_rotate
!
! ============================================================================

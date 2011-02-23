module dlg

contains
    subroutine dlg_check ( status )
        use netcdf
        integer, intent ( in) :: status
      
        if(status /= nf90_noerr) then 
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if

    end subroutine dlg_check 

    subroutine dlg_bcarsp ( x, y, z, bx, by, bz, bR, bT, bP )

        use constants

        implicit none

        real(kind=DBL), intent(in) :: x, y, z
        real(kind=DBL), intent(in) :: bx, by, bz
        real(kind=DBL), intent(out) :: bR, bT, bP

        real(kind=DBL) :: r, th, ph
        real(kind=DBL) :: bXYZ(3), bSPH(3)

        real(kind=DBL) :: rot(3,3)
   
        r   = sqrt ( x**2 + y**2 + z**2 ) 
        th  = acos ( z / r )
        ph  = atan2 ( y, x )

        rot(1,1) = sin(th)*cos(ph)
        rot(2,1) = cos(th)*cos(ph)
        rot(3,1) = -sin(ph)

        rot(1,2) = sin(th)*sin(ph)
        rot(2,2) = cos(th)*sin(ph)
        rot(3,2) = cos(ph)

        rot(1,3) = cos(th)
        rot(2,3) = -sin(th)
        rot(3,3) = 0  

        bXYZ(1) = bx
        bXYZ(2) = by
        bXYZ(3) = bz

        bSPH = matMul ( rot, bXYZ )
    
        bR = bSPH(1)
        bT = bSPH(2)
        bP = bSPH(3)

    end subroutine dlg_bcarsp
 
    function dlg_pDeriv ( array, dir, dS ) result ( dArray )
        implicit none
        
        real, dimension (:,:), allocatable :: dArray
        integer :: nX, nY, i, j
        integer, intent(IN) :: dir
        real, intent(IN) :: array (:,:), dS

        nX  = size ( array, 1 )
        nY  = size ( array, 2 )

        allocate ( dArray ( nX, nY ) )
      
        if ( dir == 2 ) then
            
            do i=1,nX
                do j=1,nY

                    if ( j > 1 .AND. j < nY ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( array(i,j+1) - array(i,j-1) )

                    if ( j == 1 ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            -3.0 * array(i,j) + 4.0 * array(i,j+1) &
                            - array(i,j+2) )

                    if ( j == nY ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            3.0 * array(i,j) - 4.0 * array(i,j-1) &
                            + array(i,j-2) )

                end do
            end do
       
        else if ( dir == 1 ) then
            
            do i=1,nX
                do j=1,nY

                    if ( i > 1 .AND. i < nX ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( array(i+1,j) - array(i-1,j) )

                    if ( i == 1 ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            -3.0 * array(i,j) + 4.0 * array(i+1,j) &
                            - array(i+2,j) )

                    if ( i == nX ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            3.0 * array(i,j) - 4.0 * array(i-1,j) &
                            + array(i-1,j) )

                end do
            end do
       
        end if 
        
    end function dlg_pDeriv

end module dlg

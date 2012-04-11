module ampFit_solve
!
! subroutine ampFit_solve_svd ( basis, dataIn, coeffs ) 
! subroutine ampFit_sumBasis ( basis, dataIn, coeffs )
!
! D.L. Green and C.L. Waters
! Centre for Space Physics
! University of Newcastle
! Australia
!
! Sep 2011
!
! ver : 201202
!
#ifdef _intel_
    use mkl95_precision, only: WP=>DP
    use mkl95_lapack, only: la_gelss=>gelss
#else
    use la_precision, only: WP=>DP
    use f95_lapack, only: la_gelss, la_lange, la_gesdd
#endif

use constants
use basisFns, only: nBFns, nArr
use ampFit_rotate

contains

subroutine ampFit_solve_svd ( basis, dataIn, coeffs ) 

#ifdef _parallel_
    use parallel
#endif

    implicit none

    type(ampData), intent(in) :: dataIn(:)
    type(ampBasis), intent(in) :: basis(:,:)
    real(DBL), allocatable, intent(inout) ::  coeffs(:)

    real(kind=WP), allocatable :: U(:,:), VT(:,:)
    real(kind=WP), allocatable :: V(:,:), alpha(:,:)
    real(kind=WP), allocatable :: W(:), beta(:)
    real(kind=WP) :: r1, ss
    integer :: la_rank, la_info, nObsThis, i, j, jj, M, N
    real(kind=DBL), allocatable :: fitted_B(:)

#ifdef _parallel_
    integer :: procRow, procCol, l_sp, m_sp, x_sp,y_sp, localRow, localCol
    integer, allocatable :: ipiv(:)
    integer :: rsrc_a, csrc_a, rsrc_b, csrc_b, M_, N_
    integer :: ia,ib, ja, jb, info, mb_a, mb_b, nb_a, nb_b, nrhs, gg, lwork
    integer, external :: numRoc, indxg2p
    character :: trans
    real(kind=DBL), allocatable :: work(:)
    integer :: iroffa, icoffa, iarow, iacol, MpA0, NqA0, ltau, lwf, lws
    integer :: iroffb, icoffb, ibrow, ibcol, Mpb0, Npb0, nrhsqb0
#else
    real(kind=WP), allocatable :: la_A(:,:), la_B(:), la_S(:)
    integer :: nObs
#endif


#ifdef _parallel_

    nObsThis = size ( dataIn )

    ! Upper left
    do i=1,nObs
        do j=1,nBFns

            procRow   = mod ( rowStartProc + (i-1)/rowBlockSize, npRow )
            procCol   = mod ( colStartProc + (j-1)/colBlockSize, npCol )

            if ( procRow==myRow .and. procCol==myCol ) then

                    l_sp = (i-1)/(npRow*rowBlockSize)
                    m_sp = (j-1)/(npCol*rowBlockSize)

                    x_sp = mod (i-1,rowBlockSize)+1
                    y_sp = mod (j-1,colBlockSize)+1

                    localRow = l_sp*rowBlockSize+x_sp
                    localCol = m_sp*colBlockSize+y_sp

                    la_A(localRow,localCol) = basis(i,j)%bTh
            endif

        enddo
    enddo

    ! Upper Right
    do i=1,nObs
        do j=nBFns+1,2*nBFns

            procRow   = mod ( rowStartProc + (i-1)/rowBlockSize, npRow )
            procCol   = mod ( colStartProc + (j-1)/colBlockSize, npCol )

            if ( procRow==myRow .and. procCol==myCol ) then

                    l_sp = (i-1)/(npRow*rowBlockSize)
                    m_sp = (j-1)/(npCol*rowBlockSize)

                    x_sp = mod (i-1,rowBlockSize)+1
                    y_sp = mod (j-1,colBlockSize)+1

                    localRow = l_sp*rowBlockSize+x_sp
                    localCol = m_sp*colBlockSize+y_sp

                    la_A(localRow,localCol) = basis(i,j-nBFns)%bPh
            endif

        enddo
    enddo

    ! Lower Left
    do i=nObs+1,2*nObs
        do j=1,nBFns

            procRow   = mod ( rowStartProc + (i-1)/rowBlockSize, npRow )
            procCol   = mod ( colStartProc + (j-1)/colBlockSize, npCol )

            if ( procRow==myRow .and. procCol==myCol ) then

                    l_sp = (i-1)/(npRow*rowBlockSize)
                    m_sp = (j-1)/(npCol*rowBlockSize)

                    x_sp = mod (i-1,rowBlockSize)+1
                    y_sp = mod (j-1,colBlockSize)+1

                    localRow = l_sp*rowBlockSize+x_sp
                    localCol = m_sp*colBlockSize+y_sp

                    la_A(localRow,localCol) = -basis(i-nObs,j)%bPh
            endif

        enddo
    enddo

    ! Lower Right
    do i=nObs+1,2*nObs
        do j=nBFns+1,2*nBFns

            procRow   = mod ( rowStartProc + (i-1)/rowBlockSize, npRow )
            procCol   = mod ( colStartProc + (j-1)/colBlockSize, npCol )

            if ( procRow==myRow .and. procCol==myCol ) then

                    l_sp = (i-1)/(npRow*rowBlockSize)
                    m_sp = (j-1)/(npCol*rowBlockSize)

                    x_sp = mod (i-1,rowBlockSize)+1
                    y_sp = mod (j-1,colBlockSize)+1

                    localRow = l_sp*rowBlockSize+x_sp
                    localCol = m_sp*colBlockSize+y_sp

                    la_A(localRow,localCol) = basis(i-nObs,j-nBFns)%bTh
            endif

        enddo
    enddo

    ! Top part of RHS
    do i=1,nObs

        procRow   = mod ( rowStartProc + (i-1)/rowBlockSize, npRow )

        if ( procRow==myRow .and. procCol==myCol ) then

                l_sp = (i-1)/(npRow*rowBlockSize)

                x_sp = mod (i-1,rowBlockSize)+1

                localRow = l_sp*rowBlockSize+x_sp

                la_B(localRow) = dataIn(i)%bT
        endif

    enddo

    ! Bottom part of RHS
    do i=nObs+1,2*nObs

        procRow   = mod ( rowStartProc + (i-1)/rowBlockSize, npRow )

        if ( procRow==myRow .and. procCol==myCol ) then

                l_sp = (i-1)/(npRow*rowBlockSize)

                x_sp = mod (i-1,rowBlockSize)+1

                localRow = l_sp*rowBlockSize+x_sp

                la_B(localRow) = dataIn(i-nObs)%bP
        endif

    enddo

    write(*,*) "DONE"

    m   = nRow
    n   = nCol
    nrhs    = 1

    ia  = 1!myRow * rowBlockSize + 1
    ja  = 1!myCol * colBlockSize + 1 
    ib  = 1!myRow * rowBlockSize + 1
    jb  = 1

    mb_a    = rowBlockSize
    nb_a    = colBlockSize

    mb_b    = rowBlockSize
    nb_b    = 1

    rsrc_a  = rowStartProc
    csrc_a  = colStartProc
    
    rsrc_b  = rowStartProc 
    csrc_b  = colStartProc

    iroffa  = mod( ia-1, mb_a )
    icoffa  = mod( ja-1, nb_a )
    iarow   = indxg2p( ia, mb_a, myRow, rsrc_a, npRow )
    iacol   = indxg2p( ja, nb_a, myCol, csrc_a, npCol )
    MpA0    = numroc( m+iroffa, mb_a, myrow, iarow, nprow )
    NqA0    = numroc( n+icoffa, nb_a, mycol, iacol, npcol )

    iroffb  = mod( ib-1, mb_b )
    icoffb  = mod( jb-1, nb_b )
    ibrow   = indxg2p( ib, mb_b, myRow, rsrc_b, npRow )
    ibcol   = indxg2p( jb, nb_b, myCol, csrc_b, npCol )
    Mpb0    = numroc( m+iroffb, mb_b, myRow, ibrow, npRow )
    Npb0    = numroc( n+iroffb, mb_b, myRow, ibrow, npRow )
    nrhsqb0 = numroc( nrhs+icoffb, nb_b, myCol, ibcol, npCol )

    ltau    = numroc( ja+min(m,n)-1, nb_a, myCol, csrc_a, npCol )
    lwf     = nb_a * ( Mpa0 + Nqa0 + nb_a )
    lws     = max( (nb_a*(nb_a-1))/2, (nrhsqb0 + Mpb0)*nb_a ) + &
                nb_a * nb_a

    lWork   = ltau + max ( lwf, lws )

    write(*,*) 'lwork: ', lwork
    
    allocate(work(lwork))

    trans = 'N'

    !allocate(ipiv(numroc ( M_, mb_a, myRow, rsrc_a, npRow ) + mb_a ))

    write(*,*) "Calling pzgesv ..."
    call pdgels ( trans, m, n, nrhs, la_A, ia, ja, descA,&
        la_B, ib, jb, descB, work, lwork, info ) 
    write(*,*) "DONE"

    allocate (coeffs(n))
    coeffs = 0

    if (myCol==0) then 
        do j=1,nColLocal

            !gg   =  myRow*rowBlockSize+1 &
            !        +mod(i-1,rowBlockSize) &
            !        +(i-1)/rowBlockSize * npRow*rowBlockSize
            gg   =  myRow*rowBlockSize+1 &
                    +mod(j-1,rowBlockSize) &
                    +(j-1)/rowBlockSize * npRow*rowBlockSize
 
            coeffs(gg) = la_B(j)

        enddo
    endif

    call cgSum2D ( iContext, 'All', ' ', nRow, 1, coeffs, nRow, -1, -1 )

    deallocate ( la_A, la_B )

#else

    nObs= size ( dataIn )

    write(*,*) '    Size of basis array: ', nBFns*2d0*nObs*2d0*16/(1024d0**2), 'MB'

    write(*,*) __FILE__, __LINE__
    write(*,*) "nObs: ", nObs, " nBFns: ", nBFns

    write(*,*) "Filling local matrix ..."

    allocate ( la_A(nObs*2, nBFns*2))

    ! M = number of data points
    ! N = number of basis functions

    N = size ( la_A, 2 ) 
    M = size ( la_A, 1 )

    write(*,*) '    N: ', N
    write(*,*) '    M: ', M

    la_A(1:nObs,1:nBFns) = basis%bTh
    la_A(nObs+1:2*nObs,nBFns+1:2*nBFns) = basis%bTh
    la_A(1:nObs,nBFns+1:2*nBFns) = basis%bPh
    la_A(nObs+1:2*nObs,1:nBFns) = -basis%bPh

    allocate ( la_B(M), la_S(N) )
    la_B(1:nObs) = dataIn%bT
    la_B(nObs+1:nObs*2) = dataIn%bP

    write(*,*) 'Solving Ax=B system ...'
  
     checkShape: &
    if (nObs*2 /= size(la_A, 1) ) then 
      err_stat = 11
      err_msg = 'ERROR: Array shape mismatch in ampFit_solve_svd'
      write(*,*) err_msg
      write(*,*) ' : ', shape(la_A), shape(dataIn)
      return 
    endif checkShape

    if(size(la_B)/=max(size(la_A,1),size(la_A,2))) then
        write(*,*) '*** ERROR: size B not right'
        stop 
    endif

    if(size(la_s)/=min(size(la_A,1),size(la_A,2))) then
        write(*,*) '*** ERROR: size S not right'
        stop
    endif


    write(*,*) '    Calling la_gelss ...'
    call la_gelss ( la_A, la_B, la_rank, la_s, 0.00001_WP, info = la_info ) 

    write(*,*) '    Rank: ', la_rank
    write(*,*) '    Info: ', la_info

    allocate (coeffs(N))
    coeffs  = la_B(1:N)
#endif
    write(*,*) 'DONE' 

end subroutine ampFit_solve_svd
!
! ------------------------------------------------------------------------------------
!
subroutine ampFit_sumBasis ( basis, dataIn, coeffs )

    type(ampBasis), intent(in) :: basis(:,:)
    type(ampData), intent(inout) :: dataIn(:)
    real(DBL), intent(in) :: coeffs(:)

    real(DBL), allocatable :: tmpA(:,:)
    integer :: nObs, M, N, j
    real(kind=DBL), allocatable :: fitted_B(:), jPar(:)

    nObs = size ( basis, 1 )

    if(nObs /= size(dataIn)) then
      err_stat = 12
      err_msg = 'Incorrect Array size in sumBasis'
      write(*,*) err_msg
      print*,'nObs=',nObs
      return
    endif

    if(size(coeffs) /= size(basis,2)*2) then
      err_stat = 13
      err_msg = 'Incorrect Array size in sumBasis'
      write(*,*) err_msg
      print*,'Coeffs array size=',size(coeffs)
      return
    endif

    allocate ( tmpA(nObs*2,nBFns*2) )

    tmpA(1:nObs,1:nBFns) = basis%bTh
    tmpA(nObs+1:2*nObs,nBFns+1:2*nBFns) = basis%bTh
    tmpA(1:nObs,nBFns+1:2*nBFns) = basis%bPh
    tmpA(nObs+1:2*nObs,1:nBFns) = -basis%bPh

    N = size ( tmpA, 2 ) 
    M = size ( tmpA, 1 )
    
    allocate(fitted_B(M))
    fitted_B = matMul ( tmpA, coeffs )

    allocate(jPar(size(dataIn)))

    jPar = matMul ( basis%Y, coeffs(nBFns+1:)*(-nArr*(nArr+1.0)) ) / (u0_*(rE+rSat)) * 1e-9 * 1e6! uAm^-2
    !jPar = matMul ( basis%Y, coeffs(1:nBFns)*(-nArr*(nArr+1.0)) ) / (u0_*(rE+rSat)) * 1e-9 * 1e6! uAm^-2

    dataIn%jPar = jPar
    dataIn%bT  = fitted_B(1:nObs)
    dataIn%bP  = fitted_B(nObs+1:nObs*2)

! Update the cartesian components from these new spherical ones
    call rtp_to_xyz_vec ( dataIn )

  end subroutine ampFit_sumBasis

end module ampFit_solve
!
! =====================================================================================

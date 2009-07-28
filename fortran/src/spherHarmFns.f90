module spherHarmFns
    use constants

    implicit none
    integer :: nTh

contains
    subroutine genLegendreFns ( nPts, theta_, m, nBFns, &
        lowEndBC_, highEndBC_, nVals_, &
        legendreFns, dLegendreFns, kVals_, minTheta, maxTheta )
        use evcrg_int
        use epirg_int
        use umach_int
        use wrcrn_int
        use sort_real_int
        use imsl_libraries
        use dislin
        use read_nameList

        implicit none

        !   In / out variables
        
        integer, intent ( in ) :: nPts
        integer, intent ( in ), optional :: lowEndBC_, highEndBC_
        real(kind=DBL), dimension ( nPts ), intent ( in ) :: theta_
        real(kind=DBL), intent ( in ) :: minTheta, maxTheta
        real(kind=DBL), allocatable, dimension ( : ), intent ( out ), optional :: nVals_, kVals_
        
        real(kind=DBL), allocatable, intent ( out ), dimension ( :, : ), optional :: &
            legendreFns, dLegendreFns
    
        integer, intent ( in ) :: m
        integer, intent ( out ) :: nBFns

        !   Internal variables

        real(kind=DBL), allocatable :: theta(:), coeffMatrix(:,:), &
          kVals(:), eigenValues(:), nValsTmp(:), nVals(:), &
          eigenVectors(:,:),eigenVectors_(:,:), dEigenVectors_(:,:)
       
        complex(kind=DBL), allocatable :: eigenValuesComplex(:), eigenVectorsComplex(:,:)

        integer, allocatable :: iiNVals(:)
       
        real(kind=DBL) :: dTh, term0, term1, term2, term3, term4
        

        integer :: lowEndBC, highEndBC, i, j

        !   Newtons method variables

        integer :: maxIterations, nIterations
        real(kind=DBL) :: maxError, nGuess, errVal, fN, dfNdN, nGuessOld

        !   Sorting variables

        integer, allocatable :: iiKeep(:)
        integer :: aStatus!, iiGoodK ( maxK )
        
!       ! laPack call variables
!       character :: jobVl = 'N', jobVr = 'V'
!       integer, parameter :: N = nTh
!       integer, parameter :: ldA = nTh
!       integer, parameter :: lWork = 4 * N
!       integer, parameter :: ldVr = N
!       integer, parameter :: ldVl = 1
!       integer :: info
!
!       real(kind=DBL) :: A ( ldA, N ), vL ( ldVl, N ), vR ( ldVr, N ), &
!           wI ( N ), work ( lWork ), wR ( N )
        
        !write (*,*) 
        !write (*,*) maxK, nPts
    !   read (*,*)

        nTh =  maxK * 10
        allocate ( theta(nTh), coeffMatrix(nTh,nTh), &
          kVals(nTh), eigenValues(nTh), nValsTmp(nTh), nVals(nTh), &
          eigenVectors(nTh,nTh),eigenVectors_(nTh,nTh), dEigenVectors_(nTh,nTh), &
          iiNVals(nTh), iiKeep(maxK), eigenValuesComplex(nTh), eigenVectorsComplex(nTh,nTh) )

        dTh = ( maxTheta - minTheta ) / ( nTh - 1.0 ) * degToRad
        theta   = (/ ( i, i = 0, nTh - 1 ) /) * dTh + minTheta * degToRad

        !write (*,*)
        !write (*,*) theta * radToDeg
        !write (*,*)
        !read (*,*)

        !   Initialize coeffMatrix

        coeffMatrix = 0.0

        if ( .not. present ( lowEndBC_ ) ) then
            
            lowEndBC    = 0
        
        else
            
            lowEndBC    = lowEndBC_

        end if

        if ( .not. present ( highEndBC_ ) ) then
            
            highEndBC   = 0
        
        else
            
            highEndBC   = highEndBC_

        end if


        fillMatrix: &
        do j = 1, nTh

        if ( j > 1 .and. j < nTh ) then

            term3   = 1.0 / ( dTh ** 2 ) &
                - cos ( theta(j) ) / sin ( theta(j) ) / ( 2.0 * dTh )!
            
            term1   = 1.0 / ( dTh ** 2 ) &
                + cos ( theta(j) ) / sin ( theta(j) ) / ( 2.0 * dTh )!

            term2   = - 2.0 / ( dTh ** 2 ) &
                - m ** 2 / ( sin ( theta(j) ) ** 2 )!

            coeffMatrix(j,j)    = term2
            coeffMatrix(j-1,j)  = term3
            coeffMatrix(j+1,j)  = term1

        end if

        if ( j == 1 ) then
    
            if ( highEndBC == 1 ) then 

                term3   = 1.0 / ( dTh ** 2 ) &
                    - cos ( theta(j) ) / sin ( theta(j) ) / ( 2.0 * dTh )!
                
                term1   = 1.0 / ( dTh ** 2 ) &
                    + cos ( theta(j) ) / sin ( theta(j) ) / ( 2.0 * dTh )!
    
                term2   = - 2.0 / ( dTh ** 2 ) &
                    - m ** 2 / ( sin ( theta(j) ) ** 2 )!
    
                coeffMatrix(j+1,j)  = term1 - term3 / 3.0
                coeffMatrix(j,j)    = term2 + term3 / 3.0 * 4.0
    
            end if

        !   if ( highEndBC == 2 .and. m == 0 ) then 

        !       term3   = 1.0 / ( dTh ** 2 ) &
        !           - cos ( theta(j) ) / sin ( theta(j) ) / ( 2.0 * dTh )!
        !       
        !       term1   = 1.0 / ( dTh ** 2 ) &
        !           + cos ( theta(j) ) / sin ( theta(j) ) / ( 2.0 * dTh )!
    
        !       term2   = - 2.0 / ( dTh ** 2 ) &
        !           - m ** 2 / ( sin ( theta(j) ) ** 2 )!
    
        !       coeffMatrix(j+1,j)  = term1 - term3 / 3.0
        !       coeffMatrix(j,j)    = term2 + term3 / 3.0 * 4.0
    
        !   end if


        end if
        
        if ( j == nTh ) then
    
            if ( lowEndBC == 1 ) then

                term3   = 1.0 / ( dTh ** 2 ) &
                    - cos ( theta(j) ) / sin ( theta(j) ) / ( 2.0 * dTh )!
                
                term1   = 1.0 / ( dTh ** 2 ) &
                    + cos ( theta(j) ) / sin ( theta(j) ) / ( 2.0 * dTh )!
    
                term2   = - 2.0 / ( dTh ** 2 ) &
                    - m ** 2 / ( sin ( theta(j) ) ** 2 )!
    
                coeffMatrix(j-1,j)  = term3 - term1 / 3.0
                coeffMatrix(j,j)    = term2 + term1 / 3.0 * 4.0

            end if

        end if


        !if ( j == 2 ) then 
        !   
        !   if ( highEndBC == 1 ) then

        !           term0   = -1.0 / dTh ** 2 &
        !               - cos ( theta(j) ) / ( 2.0 * sin ( theta(j) ) * dTh ) &
        !               - m ** 2 / sin ( theta(j) ) ** 2!
        !           term1   = 1.0 / dTh ** 2 &
        !               + cos ( theta(j) ) / ( 2.0 * sin ( theta(j) ) * dTh )!

        !           coeffMatrix(j,j)    = term0
        !           coeffMatrix(j+1,j)  = term1
        !           coeffMatrix(j-1,j)  = 0.0

        !   end if

        !   if ( highEndBC == 2 .and. m == 0 ) then 

        !           term0   = -1.0 / dTh ** 2 &
        !               - cos ( theta(j) ) / ( 2.0 * sin ( theta(j) ) * dTh ) &
        !               - m ** 2 / sin ( theta(j) ) ** 2!
        !           term1   = 1.0 / dTh ** 2 &
        !               + cos ( theta(j) ) / ( 2.0 * sin ( theta(j) ) * dTh )!

        !           coeffMatrix(j,j)    = term0
        !           coeffMatrix(j+1,j)  = term1
        !           coeffMatrix(j-1,j)  = 0.0

        !   end if

        !end if

        !if ( j == ( nTh - 1 ) ) then

        !   if ( lowEndBC == 1 ) then

        !       term0   = -1.0 / dTh ** 2 &
        !           + cos ( theta(j) ) / ( 2.0 * sin ( theta(j) ) * dTh ) &
        !           - m ** 2 / sin ( theta(j) ) ** 2
        !       term1   = 1.0 / dTh ** 2 &
        !           -   cos ( theta(j) ) / ( 2.0 * sin ( theta(j) ) * dTh )
    
        !       coeffMatrix(j,j)    = term0
        !       coeffMatrix(j-1,j)  = term1
        !       coeffMatrix(j+1,j)  = 0.0

        !   end if
    
        !end if 

        end do fillMatrix

!       A   = coeffMatrix
!
!       call sgeev ( jobVl, jobVr, N, A, lDa, wR, wI, vL, ldVl, vR, &
!           ldVr, Work, lWork, info )

        !eigenValues = 0.
        !eigenVectors    = 0.0

        call evcrg ( transpose ( coeffMatrix ), &
            eigenValuesComplex, eigenVectorsComplex, 2, nTh, nTh )

        eigenValues = real ( eigenValuesComplex )
        eigenVectors    = real ( eigenVectorsComplex )

        maxIterations = 200
        maxError = 0.001

!       calculateNVals: &
!       do i = 1, nTh
!
!           nGuess  = aimag ( sqrt ( cmplx ( eigenValues(i) ) ) )
!           nIterations = 0
!           errVal  = 1.0
!
!           solveNVal: &
!           do
!               fN  = nGuess ** 2 + nGuess + eigenValues(i)
!               dfNdn   = 2.0 * nGuess + 1.0
!               nGuessOld   = nGuess
!               nGuess  = nGuess - fN / dfNdN
!               errVal  = nGuess - nGuessOld
!               nIterations = nIterations + 1
!
!               if ( abs ( errVal ) < maxError .or. nIterations > maxIterations ) exit
!
!           end do solveNVal
!
!           if ( nGuess < 1 ) nGuess = 1d6 ! This is easier than removing these values later
!           nVals(i)    = nGuess
!
!       end do calculateNVals

        !write (*,*) eigenValues
        !write (*,*)
        !write (*,*) nVals
        !write (*,*)
        !write (*,*) abs ( -1.0 + sqrt ( 1.0 - 4.0 * eigenValues ) ) / 2.0  
        !write (*,*)
        !write (*,*) abs ( -1.0 - sqrt ( 1.0 - 4.0 * eigenValues ) ) / 2.0  

        nVals   = abs ( -1.0 + sqrt ( 1.0 - 4.0 * eigenValues ) ) / 2.0 


        !write (*,*) nVals
        !read (*,*)

        ! Before sort make sure the nVals that are less than or equal to 
        !   one do not get used. The extra 0.1 is to catch those nVals that
        ! are like 1.000000000769 etc.

        where ( nVals <= 1.1 ) nVals = 9999.9

        !   Sort to ascending order

        iiNVals = (/ ( i, i = 1, nTh ) /)
        call sort_real ( nVals, nValsTmp, iPerm = iiNVals )
        
        nVals   = nValsTmp
        eigenVectors    = eigenVectors(:,iiNVals)

    
        !   Remove the nVals less than or equal to 1
        
        !   Create k array

        if ( m == 0 ) then
    
            if ( lowEndBC == 1 ) then 
                kVals   = (/ ( i, i = 0, nTh - 1 ) /) * 2 + 2 
            else 
                kVals   = (/ ( i, i = 0, nTh - 1 ) /) * 2 + 1!
            end if
    
        else
    
            if ( lowEndBC == 1 ) then 
                kVals   = (/ ( i, i = 0, nTh - 1 ) /) * 2 + abs ( m ) 
            else 
                kVals   = (/ ( i, i = 0, nTh - 1 ) /) * 2 + abs ( m ) + 1!
            end if
    
        end if
    
        !write (*,*)
        !write (*,'(f8.2)') nVals
        !write (*,*)
        !write (*,'(f8.2)') eigenValues(iiNVals)
        !write (*,*)
        !write (*,'(f4.0)') kVals


        !   Select out maxK eigenvalues greater than m but less than maxK
       
        eigenValues = eigenValues (iiNVals)

        iiKeep  = 0
        nBFns   = 0
        !write (*,*)
        !write (*,*) m, maxK
        do i = 1, nTh
            if ( nVals(i) > abs ( m ) .and. nVals(i) > 1 .and. kVals(i) <= maxK ) then
                
                iiKeep(nBFns+1) = i
                nBFns   = nBFns + 1
                
                !write (*,*) nBFns, i, nVals(i), kVals(i), m, eigenValues(i)
            endif
        end do
    !   write (*,*) 
    !   write (*,*) kVals(1:10)
    !   write (*,*) nVals(1:10)
    !   write (*,*)
    
        !   Allocate memory for nBFns

        allocate ( eigenVectors_ ( nTh, nBFns ), dEigenVectors_ ( nTh, nBFns ), &
            stat = aStatus )

        if ( present ( legendreFns ) ) &    
            allocate ( legendreFns ( nPts, nBFns ), stat = aStatus )
        if ( present ( dLegendreFns ) ) &
            allocate ( dLegendreFns ( nPts, nBFns ), stat = aStatus )
        if ( present ( nVals_ ) ) &
            allocate ( nVals_ ( nBFns ), stat = aStatus )
        if ( present ( kVals_ ) ) &
            allocate ( kVals_ ( nBFns ), stat = aStatus )

        eigenVectors_   = eigenVectors(:,iiKeep(1:nBFns))

        !   Reconstruct the end points if m == 0 or lowEndBC  == 1

        if ( lowEndBC == 1 ) then
            do i = 1, nBFns

                eigenVectors_(nTh,i)    = 4.0 / 3.0 * eigenVectors_(nTh-1,i) &
                    - 1.0 / 3.0 * eigenVectors_(nTh-2,i)
    
            end do
        end if

        if ( highEndBC == 1 ) then
            do i = 1, nBFns

                eigenVectors_(1,i)  = 4.0 / 3.0 * eigenVectors_(2,i) &
                    - 1.0 / 3.0 * eigenVectors_(3,i)

            end do
        end if

        !if ( highEndBC == 2 .and. m == 0 ) then
        !   do i = 1, nBFns

        !       eigenVectors_(1,i)  = 4.0 / 3.0 * eigenVectors_(2,i) &
        !           - 1.0 / 3.0 * eigenVectors_(3,i)

        !   end do
        !end if

        !   Create the first derivative of the Legendre functions using
        !   an IMSL routine for quadratic polynomial derivative estimation.
        !   It does not seem to accept a vector of desired locations so it loops :(

        do i = 1, nBFns
            do j = 1, nTh

                !   qdder is an IMSL library function to estimate the derivative
                dEigenVectors_(j,i) = qdder ( 1, theta(j), theta, eigenVectors_(:,i) )

            end do
        end do

        !   Interpolate the derivative and eigenFunctions back to desired
        !   locations using IMSL routine

        do i = 1, nBFns
            do j = 1, nPts

                if ( present ( legendreFns ) ) &        
                    legendreFns(j,i)    = qdval ( theta_(j), theta, eigenVectors_(:,i) )

                if ( present ( dLegendreFns ) ) &
                    dLegendreFns(j,i)   = qdval ( theta_(j), theta, dEigenVectors_(:,i) )

            end do
        end do

!       write (*,*)
!       write (*,*) theta * radToDeg
!       write (*,*) 
!       write (*,*) theta_ * radToDeg
!       read (*,*)
        !   Return keyWord values

        if ( present ( nVals_ ) ) nVals_ = nVals(iiKeep(1:nBFns))
        if ( present ( kVals_ ) ) kVals_ = kVals(iiKeep(1:nBFns))   

        !write (*,*) nVals(iiKeep(1:nBFns))

        !   Plot legendre functions using dislin library

        !call metaFl ( 'XWIN' )
        !call page ( 3000, 3000 )
        !call winSiz ( 600, 600 )
        !call disIni () 
        !call errMod ( 'ALL', 'OFF' )
        !call axsPos ( 200, 1200 )
        !call axsLen ( 2600, 1000 )
        !call graf ( real ( minColat * 0.8 ), real ( maxColat * 1.2 ), real ( minColat ), 10.0, -0.5, 0.5, -0.5, 0.1 )
        !call incMrk ( 0 )
        !do i = 1, nBFns
        !   call color ( 'white' )
        !   if ( present ( legendreFns ) ) &
        !       call curve ( real ( theta_ * radToDeg ), real ( legendreFns(:,i) ), nPts )
        !   call color ( 'magenta' )
        !   call curve ( real ( theta * radToDeg ), real ( eigenVectors_(:,i) ), nTh )
        !end do
        !call endGrf ()

        !call axsPos ( 200, 2800 )
        !call axsLen ( 2600, 1000 )
        !call color ( 'white' )
        !call graf ( real ( minColat * 0.8 ), real ( maxColat * 1.2 ), real ( minColat ), 10.0, -1.5, 1.5, -1.5, 1.0 )
        !if ( present ( dLegendreFns ) ) then 
        !   do i = 1, nBFns
        !       call color ( 'white' )
        !       if ( present ( dLegendreFns ) ) &
        !           call curve ( real ( theta_ * radToDeg ), real ( dLegendreFns(:,i) ), nPts )
        !       call color ( 'magenta' )
        !       call curve ( real ( theta * radToDeg ), real ( dEigenVectors_(:,i) ), nTh )
        !   end do
        !end if
        !call disFin () 
    
        !   Deallocate the internal variables

        deallocate ( eigenVectors_, dEigenVectors_, stat = aStatus )
    
    end subroutine genLegendreFns


    integer function numberBFns ( )

        use read_nameList
        implicit none
        integer :: m, k

        numberBFns  = 0
        do m = -maxM, maxM
            do k = 1, maxK

                if ( mod ( k - abs ( m ), 2 ) == 0  .and. k >= abs ( m ) ) &
                    numberBFns = numberBFns + 1

                if ( mod ( k - abs ( m ), 2 ) == 1  .and. k >= abs ( m ) ) &
                    numberBFns = numberBFns + 1

            end do
        end do

        ! There are m more functions for the specific basis functions we use
        ! **[NOT]**
    
        !numberBFns = numberBFns + maxM 

        return

    end function numberBFns


    subroutine setupSHFns ( nPts, nBFns, r, coLat, lon, &
        brBFnArr, bThBFnArr, bPhBFnArr, YBFnArr, &
        minTheta, maxTheta, dLat, dLon )

        use read_nameList
        implicit none

        !   In / Out variable list

        integer, intent ( in ) :: nPts, nBFns
        real(kind=DBL), intent ( in ) :: r, coLat(nPts), lon(nPts)
        real(kind=DBL), optional, intent ( in ) :: dLat, dLon, minTheta, maxTheta
        real(kind=DBL), optional, intent ( out ) :: YBFnArr(nBFns,nPts), &
            brBFnArr(nBFns,nPts), bThBFnArr(nBFns,nPts), bPhBFnarr(nBFns,nPts)

        !   Internal variable list

        real(kind=DBL), allocatable, dimension ( :, : ) :: legendreFns_A, dLegendreFns_A
        real(kind=DBL), allocatable, dimension ( : ) :: nVals_A
        real(kind=DBL), allocatable, dimension ( :, : ) :: legendreFns_B, dLegendreFns_B
        real(kind=DBL), allocatable, dimension ( : ) :: nVals_B
        real(kind=DBL), allocatable, dimension ( :, : ) :: legendreFns_, dLegendreFns_
        real(kind=DBL), allocatable, dimension ( : ) :: nVals
        real(kind=DBL), dimension ( nBFns ) :: mArr, nArr !, kArr
        integer :: bFnCnt, m, i, nBFnsL_A = 0, nBFnsL_B = 0, nBFnsL
        real(kind=DBL) :: rE_, rDep, dYdr_rDep, orthoNormFactor

!       YBFnArr = 1d0
!       brBFnArr    = 1d0
!       bThBFnArr   = 0d0
!       bPhBFnArr   = 1d0

        nArr    = 0d0
        rE_ = 1.0
        bFnCnt  = 1
        do m = -maxM, maxM
        
            call genLegendreFns ( nPts, coLat, abs ( m ), nBFnsL_A, minTheta = minTheta, maxTheta = maxTheta, &
                lowEndBC_ = 1, highEndBC_ = 1, legendreFns = legendreFns_A, dLegendreFns = dLegendreFns_A, &
                nVals_  = nVals_A )
            call genLegendreFns ( nPts, coLat, abs ( m ), nBFnsL_B, minTheta = minTheta, maxTheta = maxTheta, &
                lowEndBC_ = 2, highEndBC_ = 1, legendreFns = legendreFns_B, dLegendreFns = dLegendreFns_B, &
                nVals_  = nVals_B )

            nBFnsL  = nBFnsL_A + nBFnsL_B

            !write (*,*) 
            !write (*,*) nBFns, nBFnsL, nBFnsL_A, nBFnsL_B
            !read (*,*)

            allocate ( nVals ( nBFnsL ), &
                legendreFns_ ( nPts, nBFnsL ), &
                dLegendreFns_ ( nPts, nBFnsL ) )

            legendreFns_(:,1:nBFnsL_A)  = legendreFns_A
            legendreFns_(:,nBFnsL_A+1:nBFnsL_A+nBFnsL_B)    = legendreFns_B

            dLegendreFns_(:,1:nBFnsL_A) = dLegendreFns_A
            dLegendreFns_(:,nBFnsL_A+1:nBFnsL_A+nBFnsL_B)   = dLegendreFns_B

            nVals(1:nBFnsL_A)   = nVals_A
            nVals(nBFnsL_A+1:nBFnsL_A+nBFnsL_B) = nVals_B

            !write (*,*)
            !write (*,*) nPts, maxM, maxK, m, nBFnsL
            !write (*,*)

            !   write (*,*) 
            !   write (*,*) nBFnsL
            !   write (*,*)

                do i = 1, nBFnsL

                    rDep    = r ** nVals(i) * &
                        ( 1.0 + nVals(i) / ( nVals(i) + 1.0 ) * ( rE_ / r ) ** ( 2.0 * nVals(i) + 1.0 ) )
                    dYdr_rDep   = nVals(i) * r ** ( nVals(i) - 1.0 ) * &
                        ( 1.0 - ( rE_ / r ) ** ( 2.0 * nVals(i) + 1.0 ) )

                    if ( m >= 0 ) then
                        
                        YBFnArr(bFnCnt,:)   = rDep * legendreFns_(:,i) * cos ( m * lon )!
                        brBFnArr(bFnCnt,:)  = dYdr_rDep * legendreFns_(:,i) * cos ( m * lon )!
                        bThBFnArr(bFnCnt,:) = 1.0 / r * rDep * dLegendreFns_(:,i) * cos ( m * lon )!
                        bPhBFnArr(bFnCnt,:) = -m * rDep / ( r * sin ( coLat ) ) * &
                            legendreFns_(:,i) * sin ( m * lon )!

                        orthoNormFactor = 1d0 / sqrt ( sum ( YBFnArr(bFnCnt,:) ** 2 * &
                            dLat * dLon * sin ( coLat ) ) )

                        YBFnArr(bFnCnt,:)   = orthoNormFactor * YBFnArr(bFnCnt,:)!
                        brBFnArr(bFnCnt,:)  = orthoNormFactor * brBFnArr(bFnCnt,:)!
                        bThBFnArr(bFnCnt,:) = orthoNormFactor * bThBFnArr(bFnCnt,:)!
                        bPhBFnArr(bFnCnt,:) = orthoNormFactor * bPhBFnArr(bFnCnt,:)!

                        mArr(bFnCnt)    = m
                        nArr(bFnCnt)    = nVals(i)
                        bFnCnt  = bFnCnt + 1

                        !write (*,*) bFnCnt, nVals(i), m, orthoNormFactor
                    else 

                        YBFnArr(bFnCnt,:)   = rDep * legendreFns_(:,i) * sin ( abs ( m ) * lon )!
                        brBFnArr(bFnCnt,:)  = dYdr_rDep * legendreFns_(:,i) * sin ( abs ( m ) * lon )!
                        bThBFnArr(bFnCnt,:) = 1.0 / r * rDep * dLegendreFns_(:,i) * sin ( abs ( m ) * lon )!
                        bPhBFnArr(bFnCnt,:) = abs(m) * rDep / ( r * sin ( coLat ) ) * &
                            legendreFns_(:,i) * cos ( abs(m) * lon )!
                    
                        orthoNormFactor = 1d0 / sqrt ( sum ( YBFnArr(bFnCnt,:) ** 2 * &
                            dLat * dLon * sin ( coLat ) ) ) 

                        YBFnArr(bFnCnt,:)   = orthoNormFactor * YBFnArr(bFnCnt,:)!
                        brBFnArr(bFnCnt,:)  = orthoNormFactor * brBFnArr(bFnCnt,:)!
                        bThBFnArr(bFnCnt,:) = orthoNormFactor * bThBFnArr(bFnCnt,:)!
                        bPhBFnArr(bFnCnt,:) = orthoNormFactor * bPhBFnArr(bFnCnt,:)!

                        mArr(bFnCnt)    = m
                        nArr(bFnCnt)    = nVals(i)
                        bFnCnt  = bFnCnt + 1
                        !write (*,*) bFnCnt, nVals(i), m, orthoNormFactor
                    end if

                end do
                !write (*,*) nArr   
                !write (*,*)
                deallocate ( legendreFns_, dLegendreFns_, nVals )
        end do

        !write (*,*) nArr
        !write (*,*)
        !write (*,*) mArr
        !read (*,*)

    end subroutine setupSHFns

end module spherHarmFns

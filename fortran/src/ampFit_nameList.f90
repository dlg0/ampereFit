module ampFit_nameList

        use constants

        implicit none

        integer :: maxK
        integer :: maxM
        real(kind=DBL) :: minCoLat
        real(kind=DBL) :: maxCoLat
        character(len=100) :: deltab_fileName
        real(kind=DBL) :: sHr = 0.0
        real(kind=DBL) :: eHr = 1.0
        character(len=4) :: basisType = 'FULL' ! FULL, CAP0, CAP1, CAPBOTH

      namelist / ampereFit / &
        maxK, &
        maxM, &
        minColat, &
        maxCoLat, &
        deltab_fileName, &
        sHr, &
        eHr, &
        basisType

contains
        subroutine init_nameList

        implicit none
        character(len=100) :: nml_fileName

        nml_fileName = 'ampereFit.nml'

        write(*,*) 'Reading namelist file: ', trim(nml_fileName), ' ...'

        open ( 7, file = nml_fileName, delim = 'APOSTROPHE', action = 'READ' )
        read ( unit = 7, nml = ampereFit )
        close ( 7 )

        write(*,*) 'DONE'

        end subroutine init_nameList

end module ampFit_nameList

module read_nameList
      use constants
      implicit none

      integer :: maxK
      integer :: maxM
      real(kind=DBL) :: minCoLat
      real(kind=DBL) :: maxCoLat
      character(len=100) :: deltab_fileName

      namelist / iridFit / &
        maxK, &
        maxM, &
        minColat, &
        maxCoLat, &
		deltab_fileName

contains
        subroutine init_nameList

        implicit none
        character(len=100) :: nml_fileName

        nml_fileName = 'iridFit.nml'

        open ( 7, file = nml_fileName, delim = 'APOSTROPHE' )
        read ( unit = 7, nml = iridFit )
        close ( 7 )

        end subroutine init_nameList

end module read_nameList

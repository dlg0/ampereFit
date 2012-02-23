! Fortran driver code for AMPERE data processing
! Calls all subroutines required to:
!  read AMPERE dB data,
!  calc spherical harmonic basis functions,
!  calc jPar
!  output to netCDF files
!
! C.L. Waters and D.L. Green
! Centre for Space Physics
! University of Newcastle
! Australia
!
! Ver: 201202
! 
subroutine ampFit ( sYear,sMonth,sDay,sHour,sMin,sSec, &
             avg_sec, south, &
             maxK, maxM, &
             latStep, lonStep, gridcoLat_deg, &
             base_name, data_dir, &
             deltab_filename, grd_name, ird_name)

  use constants
  use netcdf
  use cnvtime            ! YearSec calc - CLW Nov 2011
  use ampFit_data
  use ampFit_sort
  use basisFns
  use ampFit_solve
  use ampFit_write_output
  use ampFit_shift
  use ampFit_rotate
  use ISO_C_BINDING
  use f_aacgm
  use ampFit_aacgm

  implicit none

! Input variables
  integer, intent(in) :: sYear, sMonth, sDay, &
                         sHour, sMin, sSec, avg_sec, &
                         maxK, maxM, south
  real(kind=DBL), intent(in) :: latStep,lonStep,gridcoLat_deg
  character(len=100), intent(in) :: deltab_filename
  character(len=60), intent(in) :: data_dir
  character(len=40), intent(in) :: grd_name, ird_name
  character(len=20), intent(in) :: base_name

! GEI and GEO grid and data variables
  integer :: trk_order(6)

  type(ampData), allocatable :: &
        dataGEI_orig(:), &
        dataGEO_orig(:), &
        dataGEO_Shift(:), &
        dataGEO_Shift_half(:), &
        dataGEO_Shift_half_ghost(:), &
        dataGEO_Shift_half_ghost_fit(:), &
        dataHalfSphere_unshifted_fit(:)

  type(ampBasis), allocatable :: &
        dataBFn(:,:), gridBFn(:,:)

! AACGM Grid vars
  integer :: nLatGrid, nLonGrid, jGEOPACK
  type(ampData), allocatable :: &
            gridAACGM(:), gridGEO_Shift(:), &
            gridGEO(:)

  real(kind=DBL), allocatable :: dbThet_th(:), &
                                 dBPhi_th(:), &
                                 dBThet_ph(:), &
                                 dBPhi_ph(:), &
               aacgm_clat_rad(:), aacgm_lon_rad(:)

  integer :: i, s, np, extend, minut

! Vars for testing
  integer :: j, nbf, isat, ipln, isplice, ityp
  real :: utc, qual, px, py, pz, dbx, dby, dbz, &
              r_km, cl_deg,ln_deg, dbr, dbth, dbph

! Date/Time variables
  integer :: eYear, eMonth, eDay, eHour, eMin, eSec
  integer :: dayofyear
  integer(DBL) :: year_sec, s_epoc, e_epoc
  real(DBL) :: sHr, eHr, mlon, mlt, grid_res
  real(DBL) :: vGSEX, vGSEY, vGSEZ,ep

  real(DBL), allocatable :: coeffs(:)

! ***********************************************
! ---------------- Start of code ----------------
! ***********************************************

  grid_res = 180.0/maxK

  extend = 0
  s_epoc = cnv_YMDHMSToEpoch(sYear,sMonth,sDay,sHour,sMin,sSec)
  e_epoc = s_epoc + avg_sec
  call cnv_EpochToYMDHMS(e_epoc,eYear,eMonth,eDay,eHour,eMin,eSec)

! Decimal start and end hours of data time interval
  sHr = sHour + sMin/60.0 + sSec/3600.0
  eHr = (eDay-sDay)*24.0 + eHour + eMin/60.0 + eSec/3600.0
  if ( (eDay-sDay).ne.0) extend = 1
  
! Read AMPERE data file (netCDF), reads whole file
! Fill dataGEI_orig data structure based on selected time interval
  call ampFit_fill_structures (deltab_filename, &
              data_dir, base_name, &
              sHr, eHr, sYear, sMonth, sDay, extend, &
              dataGEI_orig)
  if (err_stat .ne. 0) return

  np=size(dataGEI_orig)
  print*,'Number of data points=',np

! Default solar wind vector (see GEOPACK documentation)
  vGSEX   =  -400.0
  vGSEY   =  0.0
  vGSEZ   =  0.0 

  minut = int ( mod( (sHr+eHr)/2.0, 1.0) * 60.0)
  year_sec = cnv_YMDHMSToYrsec ( sYear,sMonth,sDay,sHour,minut,sSec )
  dayOfYear = year_sec / (60*60*24) + 1

  call RECALC_08 (sYear,dayOfYear,sHour,minut,sSec,vGSEX,vGSEY,vGSEZ)

! Check input data against IDL
!    np=size(dataGEI_orig)
!    do i=1,np
!      write(1,*) dataGEI_orig(i)%x,dataGEI_orig(i)%y,dataGEI_orig(i)%z,&
!                 dataGEI_orig(i)%bx, dataGEI_orig(i)%by,dataGEI_orig(i)%bz, &
!                 dataGEI_orig(i)%ipln
!    enddo
!stop
! Checked Nov 2011 and Jan 2012 - Ok

! Convert from GEI to GEO coords
  print*,'Convert GEI->GEO...'
  call ampFit_conv_GEI_GEO ( dataGEI_orig, dataGEO_orig )

! cnvtime routine checked - Nov 2011 - Ok
! Fort and IDL outputs for GEI.r,th,ph and GEO.r,th,ph check
! Feb 2012 - Ok
!    np=size(dataGEI_orig)
!    do i=1,np
!      write(2,*) dataGEI_orig(i)%x, &
!                 dataGEI_orig(i)%y, &
!                 dataGEI_orig(i)%z, &
!                 dataGEI_orig(i)%r, &
!                 dataGEI_orig(i)%t*radtodeg, &
!                 dataGEI_orig(i)%p*radtodeg, &
!                 dataGEO_orig(i)%x, &
!                 dataGEO_orig(i)%y, &
!                 dataGEO_orig(i)%z, &
!                 dataGEO_orig(i)%r, &
!                 dataGEO_orig(i)%t*radtodeg, &
!                 dataGEO_orig(i)%p*radtodeg, &
!                 dataGEO_orig(i)%ipln
!    enddo
!stop

  print*,'Calc intersection point...'
  call calculate_intersection_point ( dataGEO_orig, south )
  if (err_stat .ne. 0) return

  print*,'Calc rotation matrix...'
  call calculate_shift_rotation_matrix (south)
  if (err_stat .ne. 0) return

  print*,'Shift GEO data...'
  call ShiftData ( dataGEO_orig, dataGEO_Shift )
  if (err_stat .ne. 0) return

! Check GEO and shifted GEO data
! Checked against IDL - Ok, Feb, 2012
!    np=size(dataGEO_orig)
!    do i=1,np
!      write(3,*) dataGEO_orig(i)%x, dataGEO_Shift(i)%x, &
!                 dataGEO_orig(i)%y, dataGEO_Shift(i)%y, &
!                 dataGEO_orig(i)%z, dataGEO_Shift(i)%z, &
!                 dataGEO_orig(i)%bx, dataGEO_Shift(i)%bx, &
!                 dataGEO_orig(i)%by, dataGEO_Shift(i)%by, &
!                 dataGEO_orig(i)%bz, dataGEO_Shift(i)%bz, &
!                 dataGEO_orig(i)%ipln
!    enddo
!stop

  print*,'Get 1/2 sphere data...'
  call create_dataHalfSphere ( dataGEO_Shift, dataGEO_Shift_half, south )
  if (err_stat .ne. 0) return

! Sort data along each Track
  print*,'Sort data along tracks...'
  call sort_struc ( dataGEO_Shift_half, trk_order, south )
  if (err_stat .ne. 0) return

! Identify any longitude stray satellites
  print*,'Tag longitude strays...'
  call tag_lon_strays ( dataGEO_Shift_half, south )

! Check half hemisphere, sorted data
!    np = size(dataGEO_Shift_half)
!    do i=1, np
!      write(4,*) dataGEO_Shift_half(i)%x,  &
!            dataGEO_Shift_half(i)%y, &
!            dataGEO_Shift_half(i)%z, &
!            dataGEO_Shift_half(i)%R, &
!            dataGEO_Shift_half(i)%T*radTodeg,&
!            dataGEO_Shift_half(i)%P*radTodeg, &
!            dataGEO_Shift_half(i)%bx, &
!            dataGEO_Shift_half(i)%by, &
!            dataGEO_Shift_half(i)%bz, &
!            dataGEO_Shift_half(i)%bR, &
!            dataGEO_Shift_half(i)%bT, &
!            dataGEO_Shift_half(i)%bP, &
!            dataGEO_Shift_half(i)%typ, &
!            dataGEO_Shift_half(i)%ipln
!    enddo
!stop
! Checked with IDL - Ok, Feb 2012

! Initialise AACGM before ghost routine
  s = f_AACGMInit(sYear)
  print*,'Ghost data in longitude...'

! read data from IDL code - TEST
! ####################################################
!  open(unit=1,file='/home/colinw/ampereFit/gitfork/ampereFit/fortran/aacgm_in.dat',status='old',action='read')
!  read(1,*) np
!  deallocate(dataGEO_Shift_half)
!  allocate ( dataGEO_Shift_half(np))
!  do i=1,np
!    read(1,*) utc, isat, ipln, qual, isplice, ityp, &
!              px, py, pz, dbx, dby, dbz, &
!              r_km, cl_deg,ln_deg, dbr, dbth, dbph
!    dataGEO_Shift_half(i)%utc=utc
!    dataGEO_Shift_half(i)%isat=isat
!    dataGEO_Shift_half(i)%ipln=ipln
!    dataGEO_Shift_half(i)%qual=qual
!    dataGEO_Shift_half(i)%splice=isplice
!    dataGEO_Shift_half(i)%typ=ityp
!    dataGEO_Shift_half(i)%x=px
!    dataGEO_Shift_half(i)%y=py
!    dataGEO_Shift_half(i)%z=pz
!    dataGEO_Shift_half(i)%bx=dbx
!    dataGEO_Shift_half(i)%by=dby
!    dataGEO_Shift_half(i)%bz=dbz
!    dataGEO_Shift_half(i)%R=r_km
!    dataGEO_Shift_half(i)%T=cl_deg*pi/180.0
!    dataGEO_Shift_half(i)%P=ln_deg*pi/180.0
!    dataGEO_Shift_half(i)%br=dbr
!    dataGEO_Shift_half(i)%bt=dbth
!    dataGEO_Shift_half(i)%bp=dbph
!  enddo
!  close(unit=1)
!####################################################

  call ampFit_ghosts_aacgm ( dataGEO_Shift_half, dataGEO_Shift_half_ghost, trk_order, south )
  if (err_stat .ne. 0) return

! Check ghosted data
! Checked Feb 2012 -> *********
!    np=size(dataGEO_Shift_Half_ghost)
!    do i=1,np
!      write(7,*) dataGEO_Shift_Half_ghost(i)%x, &
!                 dataGEO_Shift_Half_ghost(i)%y, &
!                 dataGEO_Shift_Half_ghost(i)%z, &
!                 dataGEO_Shift_Half_ghost(i)%R, &
!                 dataGEO_Shift_Half_ghost(i)%T*radToDeg, &
!                 dataGEO_Shift_Half_ghost(i)%P*radToDeg, &
!                 dataGEO_Shift_Half_ghost(i)%bR, &
!                 dataGEO_Shift_Half_ghost(i)%bT, &
!                 dataGEO_Shift_Half_ghost(i)%bP, &
!                 dataGEO_Shift_Half_ghost(i)%typ, &
!                 dataGEO_SHift_Half_ghost(i)%ipln
!    enddo
!stop

  print*,'Calc basis function values at data locations...'
  call create_bFns_at_data ( dataGEO_Shift_half_ghost, dataBFn, maxK, maxM )
  if (err_stat .ne. 0) return

! Check basis function values for each basis function at each data point
!  np = size(dataGEO_Shift_half_ghost)
!  nbf = size(mArr)
!  write(8,*) np, nbf
!  do i=1,np
!    do j=1,nbf
!      write(8,*) kArr(j),mArr(j), &
!        dataBFn(i,j)%br, dataBFn(i,j)%bTh, dataBFn(i,j)%bPh
!    enddo
!  enddo
! stop

! Apply weightings to basis functions
  do i=1,size(dataGEO_Shift_half_ghost)
    do j=1,size(kArr)
      dataBFn(i,j)%br = dataBFn(i,j)%br/dataGEO_Shift_half_ghost(i)%qual
      dataBFn(i,j)%bTh = dataBFn(i,j)%bTh/dataGEO_Shift_half_ghost(i)%qual
      dataBFn(i,j)%bPh = dataBFn(i,j)%bPh/dataGEO_Shift_half_ghost(i)%qual
    enddo
  enddo
  
  call ampFit_solve_svd ( dataBFn, dataGEO_Shift_half_ghost, coeffs )
  if (err_stat .ne. 0) return
! Check coeffs in fort.9
! do i=1,size(coeffs)
!  write(9,*) coeffs(i)
! enddo
! stop

  allocate(dataGEO_Shift_half_ghost_fit(size(dataGEO_Shift_half_ghost)))
  dataGEO_Shift_half_ghost_fit = dataGEO_Shift_half_ghost
! Use coeffs to calc dB at input data (shift+ghost) locations - dbR not done yet
  print*,'Calc dB using coeffs ...'
  call ampFit_sumBasis ( dataBFn, dataGEO_Shift_half_ghost_fit, coeffs )
  if (err_stat .ne. 0) return

! Check fitted data at input locations
!    np = size(dataGEO_Shift_half_ghost_fit)
!    do i=1, np
!      write(10,*) i, &
!        dataGEO_Shift_half_ghost(i)%T*radtodeg,  &
!        dataGEO_Shift_half_ghost_fit(i)%T*radtodeg,  &
!        dataGEO_Shift_half_ghost(i)%P*radtodeg,  &
!        dataGEO_Shift_half_ghost_fit(i)%P*radtodeg,  &
!        dataGEO_Shift_half_ghost(i)%bR,  &
!        dataGEO_Shift_half_ghost_fit(i)%bR,  &
!        dataGEO_Shift_half_ghost(i)%bT,  &
!        dataGEO_Shift_half_ghost_fit(i)%bT,  &
!        dataGEO_Shift_half_ghost(i)%bP,  &
!        dataGEO_Shift_half_ghost_fit(i)%bP
!    enddo
!stop

  write(*,*) 'Creating AACGM Grid ...'
  nLatGrid = int(gridcoLat_deg/latStep)
  nLonGrid = int(360.0/lonStep)
  allocate (gridGEO(nLatGrid*nLonGrid), gridAACGM(nLatGrid*nLonGrid) )
  call create_aacgm_grid (gridAACGM, gridGEO, &
             nLatGrid, nLonGrid, LatStep, LonStep, south)
  if (err_stat .ne. 0) return

!   call write_data ( gridGEO, fileName = 'output/ampData_gridGEO.nc' )
!   call write_data ( gridAACGM, fileName = 'output/ampData_gridAACGM.nc' )

  call ShiftData ( gridGEO, gridGEO_Shift )

! Check uniform grid
!    np = size(gridGEO)
!    do i=1, np
!      write(12,*) i, &
!        gridGEO(i)%T*radtodeg, gridGEO_Shift(i)%T*radtodeg, &
!        gridGEO(i)%P*radtodeg, gridGEO_Shift(i)%P*radtodeg, &
!        gridGEO(i)%R, gridGEO_Shift(i)%R
!    enddo
!stop
  print*,'Calc basis functions over grid...'
  call create_bFns_at_data ( gridGEO_Shift, gridBFn, maxK, maxM )
  if (err_stat .ne. 0) return

  print*,'Calc jPar on uniform GEO_Shift grid...'
  call ampFit_sumBasis ( gridBFn, gridGEO_Shift, coeffs )
  if (err_stat .ne. 0) return
  call UnShiftData ( gridGEO_Shift, gridGEO )

! Check jPar data - Ok, Feb 2012
! np = size(gridGEO)
! do i=1, np
!   write(11,*) i, &
!     gridGEO(i)%T*radtodeg, gridaacgm(i)%T*radtodeg, &
!     gridGEO(i)%P*radtodeg, gridaacgm(i)%P*radtodeg, &
!     gridGEO(i)%bT, gridGEO(i)%bP, &
!     gridGEO(i)%jpar
! enddo
!stop

  print*,'Calc dB and jPar on AACGM grid...'
  np = size(gridGEO)
  allocate(aacgm_clat_rad(np))
  allocate(aacgm_lon_rad(np))
  allocate (dbThet_th(np), dbPhi_th(np), &
            dbThet_ph(np), dbPhi_ph(np) )

  call geosph_to_aacgmvec (gridGEO, &
            aacgm_clat_rad, aacgm_lon_rad, &
                 dbThet_th, dBPhi_th, &
                 dBThet_ph, dBPhi_ph)
  if (err_stat .ne. 0) return

! Output uniform grid data
  call write_griddata ( grd_name, &
         sYear, sMonth, sDay, &
         sHour, sMin, sSec, &
         avg_sec, maxK, maxM, grid_res, nLatGrid, nLonGrid, &
         gridAACGM%T, gridAACGM%P, &
         dbThet_th, dBPhi_th, &
         dBThet_ph, dBPhi_ph, gridGEO%jPar )
  print*,'jPar min:max -> ',minval(gridGEO%jPar), &
                            maxval(gridGEO%jPar)
  deallocate(dbThet_th)
  deallocate(dBPhi_th)
  deallocate(dBThet_ph)
  deallocate(dBPhi_ph)
  deallocate(aacgm_lon_rad)
  deallocate(aacgm_clat_rad)

! Convert GEO_Orig data to AACGM
  np = size(dataGEO_orig)
  allocate (aacgm_clat_rad(np))
  allocate (aacgm_lon_rad(np))
  allocate (dbThet_th(np), dbPhi_th(np), &
            dbThet_ph(np), dbPhi_ph(np) )

  call geosph_to_aacgmvec ( dataGEO_orig, &
            aacgm_clat_rad, aacgm_lon_rad, &
                 dbThet_th, dBPhi_th, &
                 dBThet_ph, dBPhi_ph)
  if (err_stat .ne. 0) return
 
! Output IRD data -> as GEO_orig at present
  call write_irddata ( ird_name, &
         sYear, sMonth, sDay, &
         sHour, sMin, sSec, &
         avg_sec, &
         aacgm_clat_rad, aacgm_lon_rad, &
         dBThet_th, dBPhi_th, &
         dBThet_ph, dBPhi_ph, &
         dataGEO_orig%iPln, dataGEO_orig%isat, &
         dataGEO_orig%qual, dataGEO_orig%splice  )
    
  deallocate(dbThet_th)
  deallocate(dBPhi_th)
  deallocate(dBThet_ph)
  deallocate(dBPhi_ph)
  deallocate(aacgm_lon_rad)
  deallocate(aacgm_clat_rad)
    
end subroutine ampFit
!
! ==============================================================

program ampereFit_run
! Driver program for AmpFit.f90 and associated subroutines
! C.L. Waters, D.L. Green
! Space Physics Research Group
! University of NEwcastle
! New South Wales, Australia
! Jan 2012
!
  use cnvtime            ! YearSec calc - CLW Nov 2011
  use constants

  implicit none

! Date and time variables
  integer :: sYear, sMonth, sDay, sHour, sMin, sSec
  integer :: dayofyear, year_sec
  real(kind=DBL) :: st_run_time, en_run_time 
  integer, dimension(8) :: tme_vals
! Fit resolutions
  integer :: maxK, maxM, south, avg_sec
  real(DBL) :: latStep, lonStep, gridcoLat_deg, grid_res
! filename vars
  character(len=100) :: deltab_filename
  character(len=60) :: data_dir
  character(len=40) :: grd_name, ird_name
  character(len=20) :: base_name
  character(len=8) :: date_str
  character(len=4) :: yr_str, frmt
  character(len=3) :: hemis
  character(len=2) :: hr_str, mn_str
  integer :: date_int
!
! ---------------- Start of code ----------------
!
  maxK = 60
  maxM = 8
  south = 0
  gridcoLat_deg = 50.0
  latStep = 1.0
  lonStep = 15.0
  data_dir = '/data/ampere/'
  sYear = 2010
  sMonth = 8
  sDay = 24
  sHour = 5
  sMin = 20
  sSec = 0
  avg_sec = 600           ! int time in sec

  print*,'maxK, maxM : ',maxK,maxM
  print*,'south = ',south
  print*,'grid_coLat_deg = ',gridcoLat_deg
  print*,'dLat = ',latStep
  print*,'dLon = ',lonStep
  print*,'Data_dir = ',data_dir
  print*,'YY MM DD : ',sYear, sMonth, sDay
  print*,'HH:MM ',sHour,sMin
  print*,'Avg_Sec = ',avg_sec

  grid_res = 180.0/maxK

! Form input file name from the date information
  base_name = 'Amp_invert.ncdf'
  frmt = '(i8)'
  date_int = sYear*10000 + sMonth*100 + sDay
  write(date_str,frmt) date_int         ! format using internal file
  frmt = '(i4)'
  write(yr_str,frmt) sYear
  deltab_filename = trim(data_dir) // yr_str // '/' // trim(date_str) // trim(base_name)

  hr_str = char(ichar('0')+sHour/10) // &
           char(ichar('0')+(mod(sHour,10)/1))

  mn_str = char(ichar('0')+sMin/10) // &
           char(ichar('0')+(mod(sMin,10)/1))

  hemis = 'nth'
  if (south .eq. 1) hemis = 'sth'

  grd_name = 'output/' // trim(date_str) // '_' // hr_str &
             // mn_str // hemis // '_grd.ncdf'

  ird_name = 'output/' // trim(date_str) // '_' // hr_str &
             // mn_str // hemis // '_ird.ncdf'

  call date_and_time (VALUES = tme_vals)
  st_run_time = 86400.d0*tme_vals(3) + 3600.d0*tme_vals(5) &
       + 60.d0*tme_vals(6) + tme_vals(7) + 0.001d0*tme_vals(8)

  call ampFit (sYear,sMonth,sDay,sHour,sMin,sSec, &
         avg_sec, south, &
         maxK,maxM, &
         latStep,lonStep,gridcoLat_deg, &
         base_name, data_dir, &
         deltab_filename, grd_name, ird_name)

  call date_and_time (VALUES=tme_vals)
  en_run_time = 86400.d0*tme_vals(3) + 3600.d0*tme_vals(5) &
       + 60.d0*tme_vals(6) + tme_vals(7) + 0.001d0*tme_vals(8)

  print*,'Time for run : ',en_run_time-st_run_time,' sec'


end program ampereFit_run
!
! ========================================================================

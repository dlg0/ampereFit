module cnvtime
! This module contains:
! 1. function CNV_YMDHMST0Yrsec - converts year, month, day, hour, minut, sec to yrsec
! 2. subroutine CNV_YrsecT0YMDHMS - convert Yrsec to YMDHMS
! 3. function cnv_YMDHMSToEpoch (yr,mo,dy,hr,mn,sc)
! 4. subroutine cnv_EpochToYMDHMS (etme,yr,mo,dy,hr,mn,sc)
!
! CLW - Nov 2011
!

  use constants

  implicit none

  integer(kind=DBL), dimension(13) :: nday=(/0,31,59,90,120,151,181,212,243,273,304,334,365 /)
  integer(kind=DBL), dimension(13) :: lday=(/0,31,60,91,121,152,182,213,244,274,305,335,366 /)
  integer(kind=DBL), dimension(12) :: jday=(/0,31,59,90,120,151,181,212,243,273,304,334 /)
  integer(kind=DBL), dimension(12) :: mday=(/31,28,31,30,31,30,31,31,30,31,30,31 /)

  contains

  function cnv_YMDHMSToYrsec (yr,mo,dy,hr,mn,sc)

    implicit none
  
    integer, intent(in) :: yr,mo,dy,hr,mn,sc
    integer(kind=DBL) :: cnv_YMDHMSToYrsec, t

    t = jday(mo) + dy - 1
    if ((mo.gt.2).and.(mod(yr,4).eq.0).and.((mod(yr,100).ne.0).or.(mod(yr,400).eq.0))) t=t+1
    t = ( (t*24 + hr)*60 + mn)*60 + sc
    cnv_YMDHMSToYrsec = t

  end function cnv_YMDHMSToYrsec

!------------------------------------------------------------------
  
  subroutine cnv_YrsecToYMDHMS (yrsec,yr,mo,dy,hr,mn,sc)

    implicit none
  
    integer, intent(in) :: yr
    integer(kind=DBL), intent(in) :: yrsec
    integer, intent(out) :: mo,dy,hr,mn,sc

    integer(kind=DBL), dimension(13) :: jlday
    integer(kind=DBL) :: scday, dt
    integer :: yd, tmin, n

    if ((mod(yr,4).eq.0).and.((mod(yr,100).ne.0).or.(mod(yr,400).eq.0))) then
      do n=1,13
        jlday(n)=lday(n)
      enddo
    else
      do n=1,13
        jlday(n)=nday(n)
      enddo
    endif
 
    scday = 24*60*60
    yd = int(yrsec/scday)
    n=0
    do
      if ((n.ge.12).or.(yd.lt.jlday(n+1))) exit
      n=n+1
    end do
    mo=n
    if (n.gt.0) then
      dy=1+yd-jlday(n)
    else
      dy=yd+1
    endif
    
    dt = mod(yrsec,scday)
    hr = int(dt/(60*60))
    tmin = mod(dt,3600)
    mn = int(tmin/60)
    sc = int(mod(dt,60))

  end subroutine cnv_YrsecToYMDHMS

!------------------------------------------------------------------

  function cnv_YMDHMSToEpoch (yr,mo,dy,hr,mn,sc)

    implicit none
  
    integer, intent(in) :: yr,mo,dy,hr,mn,sc
    integer(kind=DBL) :: cnv_YMDHMSToEpoch, yrsec

    integer :: lpyear, ryear
    real(kind=DBL) :: year_sec, lyear_sec

    year_sec =365.0*24.0*3600.0
    lyear_sec=366.0*24.0*3600.0
    if (yr.lt.1970) then
      print *,'Error: Year LT 1970 in YMDHMS_to_Epoch'
      stop
    endif

    yrsec = cnv_YMDHMSToYrsec(yr,mo,dy,hr,mn,sc)
    lpyear=(yr-1969)/4
    ryear=(yr-1970)-lpyear
    cnv_YMDHMSToEpoch=yrsec+(lpyear*lyear_sec)+(ryear*year_sec)

  end function cnv_YMDHMSToEpoch

!------------------------------------------------------------------

  subroutine cnv_EpochToYMDHMS (etme,yr,mo,dy,hr,mn,sc)

    implicit none
  
    integer(kind=DBL), intent(in) :: etme
    integer, intent(out) :: yr,mo,dy,hr,mn,sc

    integer(kind=DBL) :: yrsec, tmptme
    integer :: i
    real(kind=DBL) :: year_sec, lyear_sec

    i=0
    yrsec=0
    year_sec =365.0*24.0*3600.0
    lyear_sec=366.0*24.0*3600.0
 
!    print*,'0. year_sec,lyear_sec=',year_sec,lyear_sec

    do
      if ((yrsec.gt.etme).or.(i.ge.10000)) exit
      if (mod(i,4).eq.2) then
        yrsec=yrsec+lyear_sec
      else
        yrsec=yrsec+year_sec
      endif
      i=i+1
    end do

!    print*,'1. yrsec=',yrsec

    if ( mod(i-1,4).eq.2) then
      tmptme=etme-(yrsec-lyear_sec)
    else
      tmptme=etme-(yrsec-year_sec)
    endif
    yr = i + 1969

!    print*,'2. tmptme=',tmptme

    call cnv_YrsecToYMDHMS (tmptme,yr,mo,dy,hr,mn,sc)

  end subroutine cnv_EpochToYMDHMS

end module cnvtime

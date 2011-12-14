program test_superDARN_time

use ISO_C_BINDING
use f_time
use constants, ONLY: SGL, DBL

implicit none

integer :: stat

real(C_DOUBLE) :: jd
real(DBL) :: epoch
real(DBL) :: jdIn, jdOut
real(C_DOUBLE), target :: sc
integer(C_INT), target :: yr,mo,dy,hr,mt


! f_TimeJulianToYMDHMS
! --------------------

write(*,*)
write(*,*) 'f_TimeJulianToYMDHMS'
write(*,*) '--------------------'

! take note of the various precision inputs ...

write(*,*) 'Should be: 2011-12-13  04:32:23.424006'

jdIn = 2455908.689160 ! this ends up being single precision so we lose much of the hour/minute/sec info
stat = f_TimeJulianToYMDHMS(jdIn,C_LOC(yr),C_LOC(mo),C_LOC(dy),C_LOC(hr),C_LOC(mt),C_LOC(sc))
write(*,'(12x,i4.4,a1,i2.2,a1,i2.2,a2,i2.2,a1,i2.2,a1,f9.6,a13)') yr,'-',mo,'-',dy,'  ',hr,':',mt,':',sc, '  [_SGL]'

jdIn = 2455908.689160_DBL ! this seems to work
stat = f_TimeJulianToYMDHMS(jdIn,C_LOC(yr),C_LOC(mo),C_LOC(dy),C_LOC(hr),C_LOC(mt),C_LOC(sc))
write(*,'(12x,i4.4,a1,i2.2,a1,i2.2,a2,i2.2,a1,i2.2,a1,f9.6,a13)') yr,'-',mo,'-',dy,'  ',hr,':',mt,':',sc, '  [_DBL]'

jd = 2455908.689160_C_DOUBLE ! this has to work
stat = f_TimeJulianToYMDHMS(jd,C_LOC(yr),C_LOC(mo),C_LOC(dy),C_LOC(hr),C_LOC(mt),C_LOC(sc))
write(*,'(12x,i4.4,a1,i2.2,a1,i2.2,a2,i2.2,a1,i2.2,a1,f9.6,a13)') yr,'-',mo,'-',dy,'  ',hr,':',mt,':',sc, '  [_C_DOUBLE]'

write(*,*)
write(*,*) 'f_TimeYMDHMSToJulian'
write(*,*) '--------------------'

write(*,*) 'Should be: 2455908.689160'

jdOut = f_TimeYMDHMSToJulian(yr,mo,dy,hr,mt,sc)
write(*,'(11x,f16.7)') jdOut

yr=2011_C_INT
mo=12_C_INT
dy=13_C_INT
hr=4_C_INT
mt=32_C_INT
sc=23.424006_C_DOUBLE
jdOut = f_TimeYMDHMSToJulian(yr,mo,dy,hr,mt,sc)
write(*,'(11x,f16.7)') jdOut


! f_TimeYMDHMSToEpoch
! -------------------

write(*,*)
write(*,*) 'f_TimeYMDHMSToEpoch'
write(*,*) '--------------------'

write(*,*) 'Should be: 1323750743.424006'

epoch = f_TimeYMDHMSToEPOCH(yr,mo,dy,hr,mt,sc)
write(*,'(10x,f19.6)') epoch

yr=2011_C_INT
mo=12_C_INT
dy=13_C_INT
hr=4_C_INT
mt=32_C_INT
sc=23.424006_C_DOUBLE

epoch = f_TimeYMDHMSToEPOCH(yr,mo,dy,hr,mt,sc)
write(*,'(10x,f19.6)') epoch


! f_TimeEpochToYMDHMST
! -------------------

write(*,*)
write(*,*) 'f_TimeEpochToYMDHMS'
write(*,*) '--------------------'

write(*,*) 'Should be: 2011-12-13  04:32:23.424006'

call f_TimeEpochToYMDHMS(epoch,C_LOC(yr),C_LOC(mo),C_LOC(dy),C_LOC(hr),C_LOC(mt),C_LOC(sc))
write(*,'(12x,i4.4,a1,i2.2,a1,i2.2,a2,i2.2,a1,i2.2,a1,f9.6,a13)') yr,'-',mo,'-',dy,'  ',hr,':',mt,':',sc, ''

yr=2011_C_INT
mo=12_C_INT
dy=13_C_INT
hr=4_C_INT
mt=32_C_INT
sc=23.748945_C_DOUBLE
call f_TimeEpochToYMDHMS(epoch,C_LOC(yr),C_LOC(mo),C_LOC(dy),C_LOC(hr),C_LOC(mt),C_LOC(sc))
write(*,'(12x,i4.4,a1,i2.2,a1,i2.2,a2,i2.2,a1,i2.2,a1,f9.6,a13)') yr,'-',mo,'-',dy,'  ',hr,':',mt,':',sc, ''


end program  test_superDARN_time

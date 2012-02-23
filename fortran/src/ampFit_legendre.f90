module ampFit_legendre
!
! C. L. Waters
! Space Physics
! University of Newcastle
! Australia
!
! Dec 2011
!
!  function gammln
!  function factorial
!  function plgndr
!  subroutine dLegendre
!
! Ver : 201202
!
  use constants
  
  implicit none

  contains

  FUNCTION gammln(xx)
! Calc ln[gamma(xx)] for xx > 0
! Used to calc factorial for large values
! algorithm is from Numerical Recipes (Press)

    real(DBL) :: gammln
    real(DBL), intent(in) :: xx
    real(DBL), dimension(6) :: cof
    real(DBL) :: ser, stp, tmp, x, y
    
    integer :: j

    cof(1) =  76.18009172947146d0
    cof(2) = -86.50532032941677d0
    cof(3) =  24.01409824083091d0
    cof(4) = -1.231739572450155d0
    cof(5) =  0.1208650973866179d-2
    cof(6) = -0.5395239384953d-5

    stp = 2.5066282746310005d0

    x = xx
    y = x
    tmp = x+5.5d0
    tmp = (x+0.5d0)*log(tmp)-tmp
    ser = 1.00000000190015d0
    do j=1,6
      y=y+1.0d0
      ser = ser+cof(j)/y
    enddo
    gammln = tmp+log(stp*ser/x)

    return

  end function gammln
!
! ------------------------------------------------------------------
!
  FUNCTION factorial(n)
! Calc the factorial of an integer
! If n > 32, then use the gamma function
!
! C.L. Waters
! Space Physics
! University of Newcastle
! Australia
! Dec 2011
!
    implicit none

    real(DBL) :: factorial

    integer(SGL), intent(in) :: n

    integer :: j
    real(DBL) :: fac

    if (n .lt. 0) then
      print*,'Error in factorial'
      print*,'n=',n
      stop
    endif
    fac = 1.0d0     ! for n=0 or n=1
    if ( (n .gt. 1) .and. (n .le. 32) ) then
      do j=2,n
        fac = j*fac
      enddo
    else
      fac = exp(gammln(n+1.0d0))
    endif

    factorial = fac

  end function factorial
!
! ---------------------------------------------------------------
! 
  FUNCTION plgndr ( L, M, x)
! Calculate the associated Legendre polynomial
!   0 <= m <= l
!  -1 <= x <= 1
! From Numerical Recipes, Press
!
    implicit none

    real(kind=DBL) :: plgndr
    real(kind=DBL), intent(in) :: x

    integer(kind=SGL), intent(in) :: L,M

    real(kind=DBL) :: fact, pll, pmm, pmmp1, somx2
    integer(kind=SGL) :: i, ll

    if ( (m .lt. 0) .or. (m .gt. L) .or. (abs(x) .gt. 1.0) ) then
      err_stat = 7
      err_msg = 'Bad parameter pass to PLGNDR'
      write(*,*) err_msg
      print*,'L,M,x : ',L,M,x
      return
    endif

    pmm = 1.0d0
    if (m .gt. 0) then
      somx2 = sqrt((1.0d0-x)*(1.0d0+x))
      fact = 1.0d0
      do i=1,m
        pmm = -pmm*fact*somx2
        fact = fact+2.0d0
      enddo
    endif

    if (L .eq. M) then
      plgndr = pmm
    else
      pmmp1 = x*(2.0d0*m+1.0d0)*pmm
      if (L .eq. m+1) then
        plgndr = pmmp1
      else
        do ll=m+2,L
          pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
          pmm=pmmp1
          pmmp1=pll
        enddo
        plgndr = pll
      endif
    endif
    return

  end function plgndr
!
! -----------------------------------------------------------------------
!
  subroutine dLegendre(L, M, theta, Plm, dPlm)
! Calc Schmidt semi-normalised P(l,m) and dP/d_theta
! Inputs are L,M, theta (coLat in radians)
!
! C.L. Waters
! Space Physics
! University of Newcastle
! Australia
! Dec 2011
!
    real(DBL), intent(in) :: theta
    integer(SGL) , intent(in) :: L, M
    real(DBL), intent(out) :: Plm, dPlm
  
    real(DBL) :: schmdt, x, lgp, mlt, trm, Ltrm, Rtrm, sr
    integer(SGL) :: kk, mm

    mm = abs(M)
    schmdt = 1.0d0
    x = cos(theta)

    if (mm .le. L) then
      if (m .ne. 0) schmdt=sqrt(2.0d0*factorial(L-M)/factorial(L+M))
      Plm = (-1.0)**mm*schmdt*plgndr(L,M,x)
    endif

! Do dPlm calcs
    if (mm .le. L) then
      if (m .eq. 0) then
        mlt = sqrt(0.5d0*L*(L+1))
        schmdt = sqrt(2.0d0*factorial(L-1)/factorial(L+1))    ! calc for M=1
        trm = schmdt*plgndr(L,1,x)
        dPlm = mlt*trm
      else
        schmdt = sqrt(2.0d0*factorial(L-mm+1)/factorial(L+mm-1))
        lgp = (-1.0)**(mm-1)*schmdt*plgndr(L,mm-1,x)
        sr = (L+M)*(L-mm+1)
        Ltrm = 0.5*sqrt(sr)*lgp
        if (L .eq. mm) then
          lgp = 0.0d0
        else
          schmdt = sqrt(2.0d0*factorial(L-mm-1)/factorial(L+mm+1))
          lgp = (-1.0)**(mm+1)*schmdt*plgndr(L,mm+1,x)
        endif
        sr = (L+mm+1)*(L-mm)
        Rtrm = -0.5d0*sqrt(sr)*lgp
        dPlm = Ltrm + Rtrm
      endif
    endif

  end subroutine dLegendre
!
! -----------------------------------------------------------------------
!
end module ampFit_legendre

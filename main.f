      program oscillator

      implicit real*8 (a-h,o-z)
      parameter (nt_max=10000)  ! maximum number of time slices
      parameter (nsteps_max=1000000) ! maximum number of MCMC steps
      parameter (pi = 3.141592653589793d0)
      integer nt                ! number of time slices
      integer nsteps
      integer istart  !! initial configuration flag: 0=cold, 1=hot
      integer j, idx, i
      integer idum, idum2, iv(32), iy
c      common /dasav/ idum, idum2, iv, iy
      real y(nt_max)            ! array of the discretized path
      real s                    ! euclidean action / h_bar
      real*8 eta                  ! adimensional_parameter : eta = a*omega
      real*8 gamma, alpha, mu, sigma  ! gaussian parameters for the heat bath

c     Namelist
      namelist /params/istart

c     Read input parameters
      read(5, params)
      
      nsteps = 100000

      eta = 0.5d0               ! set eta value
      alpha = (eta / 2.d0) + (1.d0 / eta)
      sigma = 1.d0 / sqrt(2.d0 * alpha)

c     Initialize ran2 RNG
      call ranstart()

c     Initialize path according to istart
      if (istart .eq. 0) then
          call cold_start(y, nt_max)
      else
          call hot_start(y, nt_max)
      end if

c     Main MCMC loop
      do j = 1, nsteps
        call total_update(y, nt, sigma, alpha, eta)
      end do

      end program oscillator


c     ============================
      subroutine cold_start(y, nt)
c     ============================
c     Initialize path to zero

      implicit none
      integer nt
      real*8 y(nt)
      integer i

      do i = 1, nt
          y(i) = 0.d0
      end do
      
      end subroutine cold_start


c     ===========================
      subroutine hot_start(y, nt)
c     ===========================
c     Initialize path with random values between -1 and 1

      implicit none
      integer nt
      real*8 y(nt)
      integer i

      do i = 1, nt
        y(i) = 2.d0 * ran2() - 1.d0
      end do
        
      end subroutine hot_start





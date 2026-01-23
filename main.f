      program oscillator

      implicit real*8 (a-h,o-z)
      parameter (bhw=5.d0)  ! beta*h_bar*omega
      parameter (n_nt=3)   ! number of different nt values
      parameter (nt_max=10000)  ! maximum number of time slices
      parameter (nsteps_max=1000000) ! maximum number of MCMC steps
      parameter (pi = 3.141592653589793d0)
      integer nt_vals(n_nt)                ! number of time slices
      integer nsteps, nt, therm_steps
      integer istart  !! initial configuration flag: 0=cold, 1=hot
      integer j, idx, i, k
      integer idum, idum2, iv(32), iy
      real*8 y(nt_max)            ! array of the discretized path
      real*8 ym(nsteps_max), y2m(nsteps_max), em(nsteps_max)  ! arrays for measurements
      real*8 s                    ! euclidean action / h_bar
      real*8 eta                  ! adimensional_parameter : eta = a*omega
      real*8 gamma, alpha, mu, sigma  ! gaussian parameters for the heat bath
      real*8 average_y(n_nt), average_y2(n_nt), average_e(n_nt)

c     Namelist
      namelist /params/istart

c     Read input parameters
      read(5, params)
      
      nsteps = 100000   ! set number of MCMC steps
      therm_steps = 10000  ! set number of thermalization steps
      nt_vals = (/5, 10, 20/) ! set different nt values
      
c     Initialize ran2 RNG
      call ranstart()

      do j = 1, n_nt
        nt = nt_vals(j)
        write(*,*) 'Running simulation with nt = ', nt
        eta = bhw / dble(nt)   ! set eta value
        alpha = (eta / 2.d0) + (1.d0 / eta)
        sigma = 1.d0 / sqrt(2.d0 * alpha)

c       Initialize path according to istart
        if (istart .eq. 0) then
          call cold_start(y, nt)
        else
          call hot_start(y, nt)
        end if

c       Thermalization
        do i = 1, therm_steps
          call total_update(y, nt, sigma, alpha, eta)
        end do

c       Main MCMC loop, with the steps consisting of 10 cicles of
c       1 Heat Bath sweep + 5 Microcanonical sweeps
        do k = 1, nsteps

          call total_update(y, nt, sigma, alpha, eta)

          call path_y(y, nt, ym(k))             !! Measurement
          call path_y2(y, nt, ym(k), y2m(k))
          call path_ene(y, nt, eta, em(k))
          
        end do

        ay = 0.d0
        ay2 = 0.d0
        ae = 0.d0
        do i = 1, nsteps
          ay = ay + ym(i)
          ay2 = ay2 + y2m(i)
          ae = ae + em(i)
        end do

        average_y(j) = ay / dble(nsteps)
        average_y2(j) = ay2 / dble(nsteps)
        average_e(j) = ae / dble(nsteps)

c     End of loop over nt values
      end do  

c     Save results to file
      open(unit=10, file='results.txt', status='unknown')
      write(10,*) '# nt    <y>    <y^2>    <E>'
      do j = 1, n_nt
        write(10,*) nt_vals(j), average_y(j), average_y2(j),
     &  average_e(j)
      end do
c     Close file
      close(10)

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
      real*8 y(nt), ran2
      integer i

      do i = 1, nt
        y(i) = 2.d0 * ran2() - 1.d0
      end do
        
      end subroutine hot_start





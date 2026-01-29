      program oscillator

      implicit real*8 (a-h,o-z)
      parameter (bhw=5.d0)  ! beta*h_bar*omega
      parameter (n_nt=1)   ! number of different nt values
      parameter (nt_max=1000)  ! maximum number of time slices
      parameter (nsteps_max=1000000) ! maximum number of MCMC steps
      parameter (pi = 3.141592653589793d0)
      integer ncorr
      common /corr_params/ ncorr
      integer nt_vals(n_nt)                ! number of time slices
      integer nsteps, nt, therm_steps
      integer istart, ground  !! initial configuration flags
      integer j, idx, i, k
      integer idum, idum2, iv(32), iy
      integer valid_n
      real*8 y(nt_max)            ! array of the discretized path
      real*8 ym(nsteps_max), y2m(nsteps_max), em(nsteps_max)  ! arrays for measurements
      real*8 ycm(nsteps_max), y3m(nsteps_max), Am(nsteps_max)
      real*8 y2cm(nsteps_max), y3cm(nsteps_max), Acm(nsteps_max)
      real*8 ay, ay2, ae, ayc, ay3, aA
      real*8 mean_jack, st_err
      real*8 s                    ! euclidean action / h_bar
      real*8 eta                  ! adimensional_parameter : eta = a*omega
      external y1, y2, y3, A
      external y1_corr, y2_corr, y3_corr, A_corr
      external connected_corr
      real*8 gamma, alpha, mu, sigma  ! gaussian parameters for the heat bath
      real*8 average_y(n_nt), average_y2(n_nt), average_e(n_nt)
      real*8 average_yc(n_nt), average_y3(n_nt), average_A(n_nt)

c     Namelist
      namelist /params/istart, ground

      istart = 0  ! cold start by default
      ground = 0  ! ground state simulation flag

c     Read input parameters
      read(5, params)
      
      nsteps = 100000   ! set number of MCMC steps
      therm_steps = 10000  ! set number of thermalization steps
      nt_vals = (/30/) ! set different nt values
      
c     Initialize ran2 RNG
      call ranstart()

      if (ground .eq. 1) then
        call ground_state()
      end if

c     Set correlator separation for the analysis
      ncorr = 1

      call set_corr_param(500)  ! Set ncorr to 500

      do j = 1, n_nt
        nt = nt_vals(j)
        write(*,*) 'Running simulation with nt = ', nt
        eta = bhw / dble(nt)   ! set eta value
        alpha = (eta / 2.d0) + (1.d0 / eta)
        sigma = 1.d0 / sqrt(2.d0 * alpha)
        
c       Set ncorr for correlation functions
        call set_corr_param(ncorr)

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
        do i = 1, nsteps

          call total_update(y, nt, sigma, alpha, eta)

          call path_observable(y, nt, y1, ym(i))    ! Blocking
          call path_observable(y, nt, y2, y2m(i))   ! Blocking
          call path_observable(y, nt, y3, y3m(i))   ! Blocking 
          call path_observable(y, nt, A, Am(i))     ! Blocking
          call path_ene(y, nt, eta, em(i))          ! Blocking
c         These correlators will only enter in the jackknife analysis
c         of the connected correlators
          call path_observable(y, nt, y1_corr, ycm(i))
          call path_observable(y, nt, y2_corr, y2cm(i))
          call path_observable(y, nt, y3_corr, y3cm(i))
          call path_observable(y, nt, A_corr, Acm(i))
          
        end do

        k = 1000  ! block size for jackknife analysis
        valid_n = (nsteps/k) * k

        call jackknife(ycm, ym, valid_n, k, connected_corr,
     &          mean_jack, st_err)
        write(*,*) 'Connected <y(0)y(nc)> - <y>^2 = ', mean_jack,
     &     ' +/- ', st_err

        call jackknife(y2cm, y2m, valid_n, k, connected_corr,
     &          mean_jack, st_err)
        write(*,*) 'Connected <y^2(0)y^2(nc)> - <y^2>^2 = ', mean_jack,
     &     ' +/- ', st_err

        call jackknife(y3cm, y3m, valid_n, k, connected_corr,
     &          mean_jack, st_err)
        write(*,*) 'Connected <y^3(0)y^3(nc)> - <y^3>^2 = ', mean_jack,
     &     ' +/- ', st_err

        call jackknife(Acm, Am, valid_n, k, connected_corr,
     &          mean_jack, st_err)
        write(*,*) 'Connected <A(0)A(nc)> - <A>^2 = ', mean_jack,
     &     ' +/- ', st_err

        ay = 0.d0
        ay2 = 0.d0
        ay3 = 0.d0
        ae = 0.d0
        ayc = 0.d0
        aA = 0.d0
        do i = 1, nsteps
          ay = ay + ym(i)
          ay2 = ay2 + y2m(i)
          ay3 = ay3 + y3m(i)
          ae = ae + em(i)
          ayc = ayc + ycm(i)
          aA = aA + Am(i)
        end do

        average_y(j) = ay / dble(nsteps)
        average_y2(j) = ay2 / dble(nsteps)
        average_e(j) = ae / dble(nsteps)
        average_yc(j) = ayc / dble(nsteps)
        average_y3(j) = ay3 / dble(nsteps)
        average_A(j) = aA / dble(nsteps)
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

      call ranfinish()

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





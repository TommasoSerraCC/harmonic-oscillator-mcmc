      program preliminary_run

      implicit real*8 (a-h,o-z)
      parameter (bhw=5.d0)  ! beta*h_bar*omega
      parameter (n_nt=3)   ! number of different nt values
      parameter (nt_max=200)  ! maximum number of time slices
      parameter (nsteps=100000) ! number of MCMC steps (10^5)
      parameter (ncorr_max=100) ! maximum nt/2 = 200/2 = 100
      integer ncorr
      common /corr_params/ ncorr
      integer nt_vals(n_nt)                ! number of time slices
      integer nt, ncorr_points
      integer istart  !! initial configuration flag
      integer j, idx, i, n
      real*8 y(nt_max)            ! array of the discretized path
      real*8 ym, y2m, y3m, Am, em  ! measurements at each step
      real*8 ycm(ncorr_max), y2cm(ncorr_max)
      real*8 y3cm(ncorr_max), Acm(ncorr_max)
      external y1, y2, y3, A
      external y1_corr, y2_corr, y3_corr, A_corr
      real*8 gamma, alpha, mu, sigma  ! gaussian parameters for the heat bath
      real*8 eta                  ! adimensional_parameter : eta = a*omega
      character*100 filename

c     Set nt values
      nt_vals(1) = 50
      nt_vals(2) = 100
      nt_vals(3) = 200

      istart = 0  ! cold start

c     Initialize ran2 RNG
      call ranstart()

c     Loop over different nt values
      do j = 1, n_nt
        nt = nt_vals(j)
        ncorr_points = nt / 2  ! number of correlator points

        write(*,*) 'Running simulation with nt = ', nt
        write(*,*) 'Number of correlator points: ', ncorr_points

        eta = bhw / dble(nt)   ! set eta value
        alpha = (eta / 2.d0) + (1.d0 / eta)
        sigma = 1.d0 / sqrt(2.d0 * alpha)

c       Initialize path (cold start)
        call cold_start(y, nt)

c       Open output file for this nt value
        write(filename, '(A,I0,A)') 
     &    'preliminary_data/raw_data_nt', nt, '.dat'
        open(unit=10, file=filename, status='unknown')

c       Main MCMC loop
        do i = 1, nsteps

          call total_update(y, nt, sigma, alpha, eta)

c         Compute observables
          call path_observable(y, nt, y1, ym)
          call path_observable(y, nt, y2, y2m)
          call path_observable(y, nt, y3, y3m)
          call path_observable(y, nt, A, Am)
          call path_ene(y, nt, eta, em)

c         Compute correlators for all n from 1 to nt/2
          do n = 1, ncorr_points
            call set_corr_param(n)
            call path_observable(y, nt, y1_corr, ycm(n))
            call path_observable(y, nt, y2_corr, y2cm(n))
            call path_observable(y, nt, y3_corr, y3cm(n))
            call path_observable(y, nt, A_corr, Acm(n))
          end do

c         Write all data to file: y, y2, y3, A, E, then correlators
          write(10, '(5(E20.12,1X))', advance='no') 
     &      ym, y2m, y3m, Am, em
          do n = 1, ncorr_points - 1
            write(10, '(4(E20.12,1X))', advance='no')
     &        ycm(n), y2cm(n), y3cm(n), Acm(n)
          end do
          write(10, '(4(E20.12,1X))')
     &      ycm(ncorr_points), y2cm(ncorr_points),
     &      y3cm(ncorr_points), Acm(ncorr_points)

c         Progress indicator every 10000 steps
          if (mod(i, 10000) .eq. 0) then
            write(*,*) '  Step: ', i, ' / ', nsteps
          end if

        end do

        close(10)
        write(*,*) 'Completed nt = ', nt
        write(*,*) 'Data saved to: ', trim(filename)
        write(*,*) ''

      end do

      call ranfinish()

      write(*,*) 'All simulations completed.'

      end program preliminary_run


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

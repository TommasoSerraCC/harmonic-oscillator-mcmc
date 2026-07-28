      program main_data_collection

      implicit real*8 (a-h,o-z)
      parameter (nt_max=200)  ! maximum number of time slices
      parameter (ncorr_max=100) ! maximum nt/2 = 200/2 = 100
      parameter (max_nt_list=20)
      real*8 bhw
      integer nsteps
      integer n_nt
      integer nt_vals(max_nt_list)
      integer energy_only
      namelist /params/ bhw, nsteps, n_nt, nt_vals, energy_only
      integer nt, ncorr_points
      integer j, i, n, nprog
      real*8 y(nt_max)
      real*8 ym, y2m, y3m, Am, em
      real*8 ycm(ncorr_max), y2cm(ncorr_max)
      real*8 y3cm(ncorr_max), Acm(ncorr_max)
      external y1, y2, y3, A
      external y1_corr, y2_corr, y3_corr, A_corr
      real*8 alpha, sigma
      real*8 eta
      character*200 filename
      character*200 datadir

c     Defaults
      energy_only = 0
      bhw = 10.d0
      nsteps = 1000000
      n_nt = 11
      nt_vals(1) = 4
      nt_vals(2) = 12
      nt_vals(3) = 24
      nt_vals(4) = 30
      nt_vals(5) = 36
      nt_vals(6) = 42
      nt_vals(7) = 50
      nt_vals(8) = 75
      nt_vals(9) = 100
      nt_vals(10) = 150
      nt_vals(11) = 200

c     Read namelist from stdin
      read(*, nml=params)

c     Build output directory
      write(datadir, '(A,I0,A,I0)')
     &  'data/bhw', nint(bhw), '_nstep', nsteps
      call system('mkdir "'//trim(datadir)//'" 2>nul')

      write(*,*) 'bhw =', bhw
      write(*,*) 'nsteps =', nsteps
      write(*,*) 'n_nt =', n_nt
      write(*,*) 'Output: ', trim(datadir)

      nprog = max(nsteps / 10, 1)

c     Initialize ran2 RNG
      call ranstart()

c     Loop over different nt values
      do j = 1, n_nt
        nt = nt_vals(j)
        ncorr_points = nt / 2

        write(*,*) 'nt =', nt, '  ncorr =', ncorr_points

        eta = bhw / dble(nt)
        alpha = (eta / 2.d0) + (1.d0 / eta)
        sigma = 1.d0 / sqrt(2.d0 * alpha)

c       Initialize path (cold start)
        call cold_start(y, nt)

c       Open output file
        if (energy_only .eq. 0) then
          write(filename, '(A,I0,A)')
     &      trim(datadir)//'/raw_data_nt', nt, '.dat'
        else
          write(filename, '(A,I0,A)')
     &      trim(datadir)//'/raw_energy_nt', nt, '.dat'
        end if
        open(unit=10, file=filename, status='unknown')

c       Main MCMC loop
        do i = 1, nsteps

          call total_update(y, nt, sigma, alpha, eta)

          if (energy_only .eq. 0) then
c           Compute observables
            call path_observable(y, nt, y1, ym)
            call path_observable(y, nt, y2, y2m)
            call path_observable(y, nt, y3, y3m)
            call path_observable(y, nt, A, Am)
            call path_ene(y, nt, eta, em)

c           Compute correlators for all n from 1 to nt/2
            do n = 1, ncorr_points
              call set_corr_param(n)
              call path_observable(y, nt, y1_corr, ycm(n))
              call path_observable(y, nt, y2_corr, y2cm(n))
              call path_observable(y, nt, y3_corr, y3cm(n))
              call path_observable(y, nt, A_corr, Acm(n))
            end do

c           Write all data to file: y, y2, y3, A, E, then correlators
            write(10, '(5(E20.12,1X))', advance='no') 
     &        ym, y2m, y3m, Am, em
            do n = 1, ncorr_points - 1
              write(10, '(4(E20.12,1X))', advance='no')
     &          ycm(n), y2cm(n), y3cm(n), Acm(n)
            end do
            write(10, '(4(E20.12,1X))')
     &        ycm(ncorr_points), y2cm(ncorr_points),
     &        y3cm(ncorr_points), Acm(ncorr_points)
          else
c           Energy-only mode: compute and write only energy
            call path_ene(y, nt, eta, em)
            write(10, '(E20.12)') em
          end if

          if (mod(i, nprog) .eq. 0) then
            write(*,*) '  Step: ', i, ' / ', nsteps
          end if

        end do

        close(10)
        write(*,*) 'Done nt =', nt

      end do

      call ranfinish()

      write(*,*) 'All simulations completed.'

      end program main_data_collection


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

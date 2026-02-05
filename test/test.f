c     MCMC simulation for a 1D Quantum Harmonic Oscillator

      program test
      implicit real*8 (a-h,o-z)
      parameter (nt_max=10000)  ! maximum number of time slices
      parameter (nsteps_max=1000000) ! maximum number of MCMC steps
      parameter (pi = 3.141592653589793d0)
      integer nt                ! number of time slices
      integer nsteps
      integer i
      integer j, idx
      integer idum, idum2, iv(32), iy
c      common /dasav/ idum, idum2, iv, iy
      real y(nt_max)            ! array of the discretized path
      real s                    ! euclidean action / h_bar
      real*8 eta                  ! adimensional_parameter : eta = a*omega
      real*8 gamma, alpha, mu, sigma  ! gaussian parameters for the heat bath
      
      nsteps = 100000

      eta = 0.5d0               ! set eta value
      alpha = (eta / 2.d0) + (1.d0 / eta)
      sigma = 1.d0 / sqrt(2.d0 * alpha)

c     Initialize ran2 RNG
      call ranstart()

c     Run tests
      call test_ran2()
      call test_box_muller()
      call test_get_indexes()
      call test_heat_bath(sigma, alpha, eta)
      call test_microcanonical(alpha)
      call test_total_update()

c     Finalize ran2 RNG
      call ranfinish()

      end program test

c     =========================================================
      subroutine test_get_indexes()
c     =========================================================
c     Test the get_indexes subroutine for periodic boundary conditions

      implicit real*8 (a-h,o-z)
      integer nt, idx, il, ir
      parameter (nt=4)
      integer e_il(nt), e_ir(nt)
      logical all_ok

      write(*,*) 'TEST - Get Indexes with Periodic Boundary Conditions'
c     Explicit test for nt=4 with hardcoded expected values

      e_il(1) = 4
      e_ir(1) = 2
      e_il(2) = 1
      e_ir(2) = 3
      e_il(3) = 2
      e_ir(3) = 4
      e_il(4) = 3
      e_ir(4) = 1

      all_ok = .true.
      do idx = 1, nt
        call get_indexes(idx, nt, il, ir)
c        write(*,*) 'Index: ', idx, ' Left: ', il, ' Right: ', ir
        if (il .ne. e_il(idx) .or. ir .ne. e_ir(idx)) then
          write(*,*) '  -> MISMATCH for idx=', idx,
     &      '-- expected (', e_il(idx), ',', e_ir(idx), ')'
          all_ok = .false.
        end if
      end do

        if (all_ok) then
              write(*,*) 'TEST COMPLETED'
        else
              write(*,*) 'TEST FAILED'
        end if

        write(*,*) ' '
      return
      end subroutine test_get_indexes

c     =========================================================
      subroutine test_ran2()
c     =========================================================
c     Test the ran2 random number generator

      implicit real*8 (a-h,o-z)
      real*4 ran2
      integer i
      parameter (npoints=50)
      real*8 x
      real*8 sum, sum2, mean, variance

      write(*,*) 'TEST - RAN2 Uniform RNG'
      sum = 0.d0
      sum2 = 0.d0
      do i = 1, npoints
          x = ran2()
          sum = sum + x
          sum2 = sum2 + x**2
      end do
      mean = sum / dble(npoints)
      variance = (sum2 / dble(npoints)) - mean**2
      write(*,*) 'RAN2 Test Results:'
      write(*,*) 'Mean: ', mean    
      write(*,*) 'Variance: ', variance

      write(*,*) 'TEST COMPLETED'
      write(*,*) ' '
      return
      end subroutine test_ran2



c     =========================================================
      subroutine test_box_muller()
c     Box-Muller Gaussian random number generator test

      implicit real*8 (a-h,o-z)
      real*4 ran2
      parameter (npoints=1000000)
      integer i
      real*8 x
      real*8 mean, sigma
      real*8 sum, sum2, sum4

      write(*,*) 'TEST - Box-Muller Gaussian RNG'

      mean = 0.d0
      sigma = 1.d0
      sum = 0.d0
      sum2 = 0.d0
      sum4 = 0.d0
      do i = 1, npoints
          call box_muller(x, mu, sigma)
          sum = sum + x
          sum2 = sum2 + x**2
          sum4 = sum4 + x**4
      end do
      mean = sum / dble(npoints)
      variance = (sum2 / dble(npoints)) - mean**2
      binder = (sum4 / dble(npoints)) / (variance**2)
      write(*,*) 'Box-Muller Test Results:'
      write(*,*) 'Mean: ', mean    
      write(*,*) 'Variance: ', variance
      write(*,*) 'Binder Cumulant: ', binder

      write(*,*) 'TEST COMPLETED'
      write(*,*) ' '
      return
      end subroutine test_box_muller


      subroutine test_heat_bath(sigma, alpha, eta)
c     Heat Bath MCMC Simulation test

      implicit real*8 (a-h,o-z)
      parameter (nt=100)
      real*8 y(nt)

      nsteps = 1000

      write(*,*) 'TEST - Heat Bath MCMC: cold start'

c     'COLD' Path initialization: y(i) = 0
      write(*,*) 'Initializing path to zero...'
      do i = 1, nt
          y(i) = 0.d0
      end do

      do j = 1, nsteps
        call heat_bath_sweep(y, nt, sigma, alpha, eta)
      end do

c      do i = 1, nt
c          write(*,*) y(i)
c      end do

      write(*,*) 'TEST COMPLETED'
      write(*,*) ' '
      return
      end subroutine test_heat_bath


c     =============================================
      subroutine test_microcanonical(alpha, eta)
c     =============================================

      implicit real*8 (a-h,o-z)
      parameter (nt=100)
      integer nsteps, j, i
      real*8 y(nt), alpha, eta
      nsteps = 1000
      write(*,*) 'TEST - Microcanonical MCMC: cold start'

c     'COLD' Path initialization: y(i) = 0
      write(*,*) 'Initializing path to zero...'
      do i = 1, nt
          y(i) = 0.d0
      end do

      y(1) = 1.d0  !! Introduce a small perturbation to avoid trivial path
      do j = 1, nsteps
        call microcanonical_sweep(y, nt, alpha, eta)
      end do

c      do i = 1, nt
c          write(*,*) y(i)
c      end do

      write(*,*) 'TEST COMPLETED'
      write(*,*) ' '
      return
      end subroutine test_microcanonical


      subroutine test_total_update()
c     Total Update MCMC Simulation test

      implicit real*8 (a-h,o-z)
      parameter (nt=100)
      integer nsteps, j, i
      real*8 y(nt), sigma, alpha, eta
      nsteps = 1000
      write(*,*) 'TEST - Total Update MCMC: cold start'
c     'COLD' Path initialization: y(i) = 0
      write(*,*) 'Initializing path to zero...'
      do i = 1, nt
          y(i) = 0.d0
      end do

      eta = 0.5d0               ! set eta value
      alpha = (eta / 2.d0) + (1.d0 / eta)
      sigma = 1.d0 / sqrt(2.d0 * alpha)

      do j = 1, nsteps
        call total_update(y, nt, sigma, alpha, eta)
      end do

      do i = 1, nt
          write(*,*) y(i)
      end do

      end subroutine test_total_update
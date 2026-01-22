c     MCMC simulation for a 1D Quantum Harmonic Oscillator

      program oscillator
      implicit real*8 (a-h,o-z)
      parameter (nt_max=10000)  ! maximum number of time slices
      parameter (nsteps_max=1000000) ! maximum number of MCMC steps
      parameter (pi = 3.141592653589793d0)
      integer nt                ! number of time slices
      integer nsteps
      integer i
      integer j, idx
      real y(nt_max)            ! array of the discretized path
      real s                    ! euclidean action / h_bar
      real eta                  ! adimensional_parameter : eta = a*omega
      real gamma, alpha, mu, sigma  ! gaussian parameters for the heat bath
      
      nsteps = 100000

      eta = 0.5d0               ! set eta value
      alpha = (eta / 2.d0) + (1.d0 / eta)
      sigma = 1.d0 / sqrt(2.d0 * alpha)

c     'COLD' Path initialization: y(i) = 0
      nt = 100
      do i = 1, nt
          y(i) = 0.0
      end do

      do j = 1, nsteps
        idx = mod(j-1, nt) + 1
        call heat_bath_step(y, nt, idx, sigma, gamma, alpha, eta)
      end do

      do i = 1, nt
          write(*,*) y(i)
      end do

      end program oscillator


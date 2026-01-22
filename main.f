c     MCMC simulation for a 1D Quantum Harmonic Oscillator

      program oscillator
c      implicit real*8 (a-h,o-z)
      parameter (nt_max=10000)  ! maximum number of time slices
      integer nt                ! number of time slices
      real y(nt_max)            ! array of the discretized path
      real s                    ! euclidean action / h_bar
      real eta                  ! adimensional_parameter : eta = a*omega
      
c     Path initialization: y(i) = 0
      nt = 1000
      do i = 1, nt
          y(i) = 0.0
      end do

      end program oscillator

      subroutine euclidean_action(s, y, nt, eta)
c     computes the euclidean action for a given path y

      eta_comb = (eta / 2.d0) + (1.d0 / eta)
      do i = 1, nt-1
c       skip the last term for periodic boundary conditions
        s = s + y(i)**2 * eta_comb - (1.d0/eta)*y(i) * y(i+1)
      end do
      s = s + y(nt)**2 * eta_comb - (1.d0/eta)*y(nt) * y(1)  ! periodic BC
      end subroutine euclidean_action


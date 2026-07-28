c     =========================================================
      subroutine microcanonical_sweep(y, nt, alpha, eta)
c     =========================================================
c     Performs a microcanonical sweep update on the entire path
c     Each site is updated deterministically to conserve the Euclidean action.
c     This is done by reflecting the current value about the mean of its neighbors,

      implicit real*8 (a-h,o-z)
      integer nt, idx, il, ir
      real*8 y(nt)
      real*8 gamma, mu, alpha

c     First element at the left boundary
      gamma = (y(nt) + y(2)) / eta
      mu = gamma / (2.d0 * alpha)
      y(1) = 2.d0 * mu - y(1) !! Reflection

c     Middle elements
      do idx = 2, nt-1
        gamma = (y(idx-1) + y(idx+1)) / eta
        mu = gamma / (2.d0 * alpha)
        y(idx) = 2.d0 * mu - y(idx) !! Reflection
      end do

c     Last element at the right boundary
      gamma = (y(nt-1) + y(1)) / eta
      mu = gamma / (2.d0 * alpha)
      y(nt) = 2.d0 * mu - y(nt) !! Reflection

      return
      end subroutine microcanonical_sweep

c     =========================================================
      subroutine heat_bath_sweep(y, nt, sigma, alpha, eta)
c     =========================================================
c     Performs a Heat Bath sweep update on the entire path.
c     Each site is updated by sampling from its full conditional probability,
c     fixing all other degrees of freedom. This distribution is derived from 
c     the local Euclidean action and is a Gaussian with mean and variance 
c     determined by neighboring sites, independent of the current site value.

      implicit real*8 (a-h,o-z)
      integer nt, idx, il, ir
      real*8 y(nt), sigma, gamma, alpha, mu, eta
      
c     First element at the left boundary
      gamma = (y(nt) + y(2)) / eta
      mu = gamma / (2.d0 * alpha)
      call box_muller(y(1), mu, sigma) !! Sample from Gaussian

c     Middle elements
      do idx = 2, nt-1
        gamma = (y(idx-1) + y(idx+1)) / eta
        mu = gamma / (2.d0 * alpha)
        call box_muller(y(idx), mu, sigma)
      end do

c     Last element at the right boundary
      gamma = (y(nt-1) + y(1)) / eta
      mu = gamma / (2.d0 * alpha)
      call box_muller(y(nt), mu, sigma)

      end subroutine heat_bath_sweep


c     ========================================================
      subroutine total_update(y, nt, sigma, alpha, eta)
c     ========================================================
c     Performs 10 times a combined update on the path:
c     first a Heat Bath sweep, followed by 5 Microcanonical sweeps.

      implicit real*8 (a-h,o-z)
      integer nt
      real*8 y(nt), sigma, alpha, eta
      integer i

      do j = 1, 10
        call heat_bath_sweep(y, nt, sigma, alpha, eta)
        do i = 1, 5
          call microcanonical_sweep(y, nt, alpha, eta)
        end do
      end do

      return
      end subroutine total_update
      

c     =========================================================
      subroutine get_indexes(idx, nt, il, ir)
c     =========================================================
c     Determines the left and right neighbor indices with periodic boundary conditions

      implicit none
      integer idx, nt, il, ir

      if (idx .eq. 1) then  !! left boundary
        il = nt
        ir = 2
      else if (idx .eq. nt) then  !! right boundary
        il = nt - 1
        ir = 1
      else
        il = idx - 1
        ir = idx + 1
      end if

      end subroutine get_indexes


c     =========================================================
      subroutine box_muller(x, mu, sigma)
c     =========================================================
c     Generates a Gaussian random number using the Box-Muller algorithm

      implicit real*8 (a-h,o-z)
      real*4 ran2
      parameter (pi = 3.141592653589793d0)
      real*8 x, mu, sigma
      real*8 u1, u2, z0

      u1 = dble(ran2())
      u2 = dble(ran2())

      z0 = sqrt(-2.d0 * log(u1)) * cos(2.d0 * pi * u2)
      x = mu + sigma * z0

      end subroutine box_muller


c     =========================================================
      subroutine euclidean_action(s, y, nt, eta, alpha)
c     =========================================================
c     Computes the discrete Euclidean action of a path y.
c     Used by the test suite to verify that the microcanonical
c     (over-relaxation) sweep leaves the action invariant.

      implicit none
      integer nt, i
      real*8 s, y(nt), eta, alpha

      s = 0.d0
      do i = 1, nt-1
        s = s + y(i)**2 * alpha - (1.d0/eta) * y(i) * y(i+1)
      end do
c     last term closes the periodic boundary condition
      s = s + y(nt)**2 * alpha - (1.d0/eta) * y(nt) * y(1)

      end subroutine euclidean_action


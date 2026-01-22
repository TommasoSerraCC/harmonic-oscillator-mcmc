c     =========================================================
      subroutine heat_bath_step(y, nt, idx, sigma, gamma, alpha,
     &                          eta)
c     =========================================================
c     Performs a heat bath update on the path at index idx
      implicit real*8 (a-h,o-z)
      integer nt, idx
      real*8 y(nt), sigma, gamma, alpha, mu, eta
      
c     Determine neighboring indices with periodic boundary conditions
      il = mod(idx - 2, nt) + 1
      ir = mod(idx, nt) + 1

c     Compute the mean (mu) for the Gaussian distribution
      gamma = (y(il) + y(ir)) / eta
      mu = gamma / (2.d0 * alpha)

      call box_muller(y(idx), mu, sigma)

      end subroutine heat_bath_step


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
c     Computes the euclidean action for a given path y
      implicit real*8 (a-h,o-z)

      s = 0.d0
      do i = 1, nt-1
c       skip the last term for periodic boundary conditions
        s = s + y(i)**2 * alpha - (1.d0/eta)*y(i) * y(i+1)
      end do
      s = s + y(nt)**2 * alpha - (1.d0/eta)*y(nt) * y(1)  ! periodic BC

      end subroutine euclidean_action



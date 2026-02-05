c ============================
c OBSERVABLES SUBROUTINES
c ============================


c     Local observable: y
      function y1(y, nt, i)
      implicit none
      integer nt, i
      real*8 y(nt), y1
      y1 = y(i)
      return
      end

c     Local observable: y^2
      function y2(y, nt, i)
      implicit none
      integer nt, i
      real*8 y(nt), y2
      y2 = y(i)**2
      return
      end

c     Local observable: y^3
      function y3(y, nt, i)
      implicit none
      integer nt, i
      real*8 y(nt), y3
      y3 = y(i)**3
      return
      end

c     Local observable: A = y^3 - 1.5*y
      function A(y, nt, i)
      implicit none
      integer nt, i
      real*8 y(nt), A
      A = y(i)**3 - 1.5d0 * y(i)
      return
      end


c =====================================================
c CORRELATION FUNCTIONS SUBROUTINES
c =====================================================

c     Common block for correlation parameter
      block data corr_params_init
      implicit none
      integer ncorr
      common /corr_params/ ncorr
      data ncorr /0/
      end


c     Correlator for a generic observable: C(i) = obs(i) * obs(i + ncorr)
      function correlator_gen(y, nt, i, obs_func)
      implicit none
      integer nt, i, ncorr_global, j
      real*8 y(nt), correlator_gen, obs_func
      external obs_func
      common /corr_params/ ncorr_global

      j = i + ncorr_global
      if (j > nt) j = j - nt

      correlator_gen = obs_func(y, nt, i) * obs_func(y, nt, j)
      return
      end

c     Correlator wrapper for y
      function y1_corr(y, nt, i)
      implicit none
      integer nt, i
      real*8 y(nt), y1_corr, correlator_gen, y1
      external correlator_gen, y1
      y1_corr = correlator_gen(y, nt, i, y1)
      return
      end

c     Correlator wrapper for y^2
      function y2_corr(y, nt, i)
      implicit none
      integer nt, i
      real*8 y(nt), y2_corr, correlator_gen, y2
      external correlator_gen, y2
      y2_corr = correlator_gen(y, nt, i, y2)
      return
      end

c     Correlator wrapper for y^3
      function y3_corr(y, nt, i)
      implicit none
      integer nt, i
      real*8 y(nt), y3_corr, correlator_gen, y3
      external correlator_gen, y3
      y3_corr = correlator_gen(y, nt, i, y3)
      return
      end

c     Correlator wrapper for A
      function A_corr(y, nt, i)
      implicit none
      integer nt, i
      real*8 y(nt), A_corr, correlator_gen, A
      external correlator_gen, A
      A_corr = correlator_gen(y, nt, i, A)
      return
      end


c     Compute path average of a generic observable
      subroutine path_observable(y, nt, obs_func, result)
      implicit none
      integer nt, i
      real*8 y(nt), result, obs_func
      external obs_func
      
      result = 0.d0
      do i = 1, nt
        result = result + obs_func(y, nt, i)
      end do
      result = result / dble(nt)

      end subroutine path_observable

c     Set correlation parameter
      subroutine set_corr_param(ncorr_val)
      implicit none
      integer ncorr_val, ncorr_global
      common /corr_params/ ncorr_global
      ncorr_global = ncorr_val
      end subroutine set_corr_param


c     Connected correlator: <O(0)O(nc)> - <O>^2
      subroutine connected_corr(corr_values, obs_values, n, result)
      implicit none
      integer n, i
      real*8 corr_values(n), obs_values(n), result
      real*8 mean_corr, mean_obs

      mean_corr = 0.d0
      mean_obs = 0.d0
      do i = 1, n
        mean_corr = mean_corr + corr_values(i)
        mean_obs = mean_obs + obs_values(i)
      end do
      mean_corr = mean_corr / dble(n)
      mean_obs = mean_obs / dble(n)

      result = mean_corr - mean_obs**2

      end subroutine connected_corr


c     Compute the energy over the path
      subroutine path_ene(y, nt, eta, energy)

      implicit none
      integer nt, i
      real*8 y(nt), eta, energy, kin, pot

      energy = 1.d0 / (2.d0 * eta) !! offset

c     Explicitly compute kinetic and potential energy
c     for the first point to handle periodic BC
      kin = (y(1) - y(nt))**2
      pot = y(1)**2

c     Sum over the rest of the path
      do i = 2, nt
        kin = kin + (y(i) - y(i-1))**2
        pot = pot + y(i)**2
      end do

      kin = kin / (2.d0 * nt * eta**2)
      pot = pot / (2.d0 * nt)

      energy = energy - kin + pot

      end subroutine path_ene
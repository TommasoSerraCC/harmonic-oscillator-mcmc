c     =========================================================
c     Unit tests for the MCMC quantum harmonic oscillator code
c     =========================================================
c     Build and run with:   make test
c     Exits with status 1 if any check fails.

      program test
      implicit real*8 (a-h,o-z)
      integer nfail
      real*8 eta, alpha, sigma

      eta   = 0.5d0
      alpha = (eta / 2.d0) + (1.d0 / eta)
      sigma = 1.d0 / sqrt(2.d0 * alpha)

      nfail = 0

      call ranstart()

      call test_ran2(nfail)
      call test_box_muller(nfail)
      call test_get_indexes(nfail)
      call test_action_invariance(nfail, sigma, alpha, eta)
      call test_heat_bath(nfail, sigma, alpha, eta)
      call test_total_update(nfail, sigma, alpha, eta)
      call test_path_ene(nfail, sigma, alpha, eta)
      call test_correlator(nfail, sigma, alpha, eta)

      call ranfinish()

      write(*,*) '==========================================='
      if (nfail .eq. 0) then
        write(*,*) 'ALL TESTS PASSED'
      else
        write(*,*) 'TESTS FAILED: ', nfail, ' check(s)'
      end if
      write(*,*) '==========================================='

      if (nfail .ne. 0) stop 1

      end program test


c     =========================================================
      subroutine check(nfail, value, expected, tol, label)
c     =========================================================
c     Compares value against expected within tol and reports.

      implicit none
      integer nfail
      real*8 value, expected, tol
      character*(*) label

      if (abs(value - expected) .le. tol) then
        write(*,*) '  [ OK ] ', label
      else
        write(*,*) '  [FAIL] ', label
        write(*,*) '         got ', value, ' expected ', expected,
     &             ' +/- ', tol
        nfail = nfail + 1
      end if

      end subroutine check


c     =========================================================
      subroutine test_ran2(nfail)
c     =========================================================
c     The ran2 generator must be uniform on [0,1):
c     mean = 1/2, variance = 1/12.

      implicit real*8 (a-h,o-z)
      integer nfail, i
      real*4 ran2
      parameter (npoints=100000)
      real*8 x, s1, s2, mean, variance

      write(*,*) 'TEST - ran2 uniform RNG'

      s1 = 0.d0
      s2 = 0.d0
      do i = 1, npoints
        x = dble(ran2())
        s1 = s1 + x
        s2 = s2 + x**2
      end do
      mean = s1 / dble(npoints)
      variance = (s2 / dble(npoints)) - mean**2

      write(*,*) '  mean     = ', mean
      write(*,*) '  variance = ', variance

      call check(nfail, mean, 0.5d0, 0.01d0, 'ran2 mean')
      call check(nfail, variance, 1.d0/12.d0, 0.005d0,
     &           'ran2 variance')
      write(*,*) ' '

      end subroutine test_ran2


c     =========================================================
      subroutine test_box_muller(nfail)
c     =========================================================
c     Box-Muller must produce a standard normal deviate:
c     mean = 0, variance = 1, <x^4>/sigma^4 = 3.

      implicit real*8 (a-h,o-z)
      integer nfail, i
      parameter (npoints=1000000)
      real*8 x, mu, sig, s1, s2, s4
      real*8 mean, variance, binder

      write(*,*) 'TEST - Box-Muller Gaussian RNG'

      mu  = 0.d0
      sig = 1.d0
      s1 = 0.d0
      s2 = 0.d0
      s4 = 0.d0
      do i = 1, npoints
        call box_muller(x, mu, sig)
        s1 = s1 + x
        s2 = s2 + x**2
        s4 = s4 + x**4
      end do
      mean = s1 / dble(npoints)
      variance = (s2 / dble(npoints)) - mean**2
      binder = (s4 / dble(npoints)) / variance**2

      write(*,*) '  mean     = ', mean
      write(*,*) '  variance = ', variance
      write(*,*) '  kurtosis = ', binder

      call check(nfail, mean, 0.d0, 0.01d0, 'box-muller mean')
      call check(nfail, variance, 1.d0, 0.01d0,
     &           'box-muller variance')
      call check(nfail, binder, 3.d0, 0.05d0,
     &           'box-muller <x^4>/sigma^4')
      write(*,*) ' '

      end subroutine test_box_muller


c     =========================================================
      subroutine test_get_indexes(nfail)
c     =========================================================
c     Neighbour indices must wrap around (periodic lattice).

      implicit none
      integer nfail, nt, idx, il, ir
      parameter (nt=4)
      integer e_il(nt), e_ir(nt)
      logical all_ok

      write(*,*) 'TEST - periodic neighbour indices'

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
        if (il .ne. e_il(idx) .or. ir .ne. e_ir(idx)) then
          write(*,*) '         mismatch at idx = ', idx
          write(*,*) '         got      (', il, ',', ir, ')'
          write(*,*) '         expected (', e_il(idx), ',',
     &               e_ir(idx), ')'
          all_ok = .false.
        end if
      end do

      if (all_ok) then
        write(*,*) '  [ OK ] get_indexes wrap-around'
      else
        write(*,*) '  [FAIL] get_indexes wrap-around'
        nfail = nfail + 1
      end if
      write(*,*) ' '

      end subroutine test_get_indexes


c     =========================================================
      subroutine test_action_invariance(nfail, sigma, alpha, eta)
c     =========================================================
c     The microcanonical (over-relaxation) sweep is a reflection
c     about the conditional mean, so it must leave the Euclidean
c     action exactly invariant.

      implicit real*8 (a-h,o-z)
      integer nfail, i, j
      parameter (nt=100)
      real*8 y(nt), sigma, alpha, eta
      real*8 s_before, s_after, drift

      write(*,*) 'TEST - microcanonical sweep conserves the action'

c     Generate a non-trivial configuration with the heat bath
      do i = 1, nt
        y(i) = 0.d0
      end do
      do j = 1, 200
        call heat_bath_sweep(y, nt, sigma, alpha, eta)
      end do

      call euclidean_action(s_before, y, nt, eta, alpha)
      do j = 1, 20
        call microcanonical_sweep(y, nt, alpha, eta)
      end do
      call euclidean_action(s_after, y, nt, eta, alpha)

      drift = abs(s_after - s_before) / abs(s_before)
      write(*,*) '  S before      = ', s_before
      write(*,*) '  S after       = ', s_after
      write(*,*) '  relative drift= ', drift

      call check(nfail, drift, 0.d0, 1.d-8,
     &           'action invariance under over-relaxation')
      write(*,*) ' '

      end subroutine test_action_invariance


c     =========================================================
      subroutine test_heat_bath(nfail, sigma, alpha, eta)
c     =========================================================
c     Heat-bath sampling alone must already reproduce the
c     symmetric equilibrium distribution: <y> = 0.

      implicit real*8 (a-h,o-z)
      integer nfail, i, j
      parameter (nt=100)
      parameter (ntherm=1000)
      parameter (nmeas=5000)
      real*8 y(nt), sigma, alpha, eta
      real*8 s1, s2, mean, mean2

      write(*,*) 'TEST - heat-bath sweep equilibrium'

      do i = 1, nt
        y(i) = 0.d0
      end do
      do j = 1, ntherm
        call heat_bath_sweep(y, nt, sigma, alpha, eta)
      end do

      s1 = 0.d0
      s2 = 0.d0
      do j = 1, nmeas
        call heat_bath_sweep(y, nt, sigma, alpha, eta)
        do i = 1, nt
          s1 = s1 + y(i)
          s2 = s2 + y(i)**2
        end do
      end do
      mean  = s1 / dble(nmeas * nt)
      mean2 = s2 / dble(nmeas * nt)

      write(*,*) '  <y>   = ', mean
      write(*,*) '  <y^2> = ', mean2

      call check(nfail, mean, 0.d0, 0.05d0, 'heat-bath <y> = 0')
      call check(nfail, mean2, 0.5d0, 0.06d0,
     &           'heat-bath <y^2> near ground state')
      write(*,*) ' '

      end subroutine test_heat_bath


c     =========================================================
      subroutine test_total_update(nfail, sigma, alpha, eta)
c     =========================================================
c     The production update (10 x [1 heat bath + 5 over-relax])
c     must sample the same distribution, with <y^2> close to the
c     continuum ground-state value 1/2 at this low temperature.

      implicit real*8 (a-h,o-z)
      integer nfail, i, j
      parameter (nt=100)
      parameter (ntherm=200)
      parameter (nmeas=2000)
      real*8 y(nt), sigma, alpha, eta
      real*8 s1, s2, mean, mean2

      write(*,*) 'TEST - total_update equilibrium'

      do i = 1, nt
        y(i) = 0.d0
      end do
      do j = 1, ntherm
        call total_update(y, nt, sigma, alpha, eta)
      end do

      s1 = 0.d0
      s2 = 0.d0
      do j = 1, nmeas
        call total_update(y, nt, sigma, alpha, eta)
        do i = 1, nt
          s1 = s1 + y(i)
          s2 = s2 + y(i)**2
        end do
      end do
      mean  = s1 / dble(nmeas * nt)
      mean2 = s2 / dble(nmeas * nt)

      write(*,*) '  <y>   = ', mean
      write(*,*) '  <y^2> = ', mean2

      call check(nfail, mean, 0.d0, 0.05d0, 'total_update <y> = 0')
      call check(nfail, mean2, 0.5d0, 0.06d0,
     &           'total_update <y^2> near ground state')
      write(*,*) ' '

      end subroutine test_total_update


c     =========================================================
      subroutine test_path_ene(nfail, sigma, alpha, eta)
c     =========================================================
c     At beta*hbar*omega = nt*eta = 50 the system is frozen in
c     the ground state, so the energy estimator must return
c     approximately E0 = 1/2 (up to O(eta^2) lattice effects).

      implicit real*8 (a-h,o-z)
      integer nfail, i, j
      parameter (nt=100)
      parameter (ntherm=200)
      parameter (nmeas=2000)
      real*8 y(nt), sigma, alpha, eta
      real*8 em, s1, mean_e

      write(*,*) 'TEST - energy estimator (path_ene)'

      do i = 1, nt
        y(i) = 0.d0
      end do
      do j = 1, ntherm
        call total_update(y, nt, sigma, alpha, eta)
      end do

      s1 = 0.d0
      do j = 1, nmeas
        call total_update(y, nt, sigma, alpha, eta)
        call path_ene(y, nt, eta, em)
        s1 = s1 + em
      end do
      mean_e = s1 / dble(nmeas)

      write(*,*) '  <E> = ', mean_e

      call check(nfail, mean_e, 0.5d0, 0.06d0,
     &           'ground state energy ~ 1/2')
      write(*,*) ' '

      end subroutine test_path_ene


c     =========================================================
      subroutine test_correlator(nfail, sigma, alpha, eta)
c     =========================================================
c     Consistency of the correlator machinery: at separation
c     n = 0 the correlator of y reduces exactly to y^2.

      implicit real*8 (a-h,o-z)
      integer nfail, i, j
      parameter (nt=100)
      real*8 y(nt), sigma, alpha, eta
      real*8 c0, y2m
      external y1_corr, y2

      write(*,*) 'TEST - correlator at zero separation'

      do i = 1, nt
        y(i) = 0.d0
      end do
      do j = 1, 200
        call total_update(y, nt, sigma, alpha, eta)
      end do

      call set_corr_param(0)
      call path_observable(y, nt, y1_corr, c0)
      call path_observable(y, nt, y2, y2m)

      write(*,*) '  C_y(0) = ', c0
      write(*,*) '  <y^2>  = ', y2m

      call check(nfail, c0, y2m, 1.d-12, 'C_y(0) equals <y^2>')
      write(*,*) ' '

      end subroutine test_correlator

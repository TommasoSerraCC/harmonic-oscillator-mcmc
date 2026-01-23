c     =============================================
      subroutine pos_mean(y, nt, mean)
c     =============================================
c     Compute the mean of the path positions

      implicit none
      integer nt
      real*8 y(nt), mean

      mean = 0.d0

      do i = 1, nt
        mean = mean + y(i)
      end do

      mean = mean / dble(nt)

      end subroutine pos_mean


c     =============================================
      subroutine pos_variance(y, nt, mean, variance)
c     =============================================

      implicit none
      integer nt, i
      real*8 y(nt), mean, variance

      variance = 0.d0

      do i = 1, nt
        variance = variance + (y(i) - mean)**2
      end do

      variance = variance / dble(nt - 1.d0)
    
      end subroutine pos_variance


c     =============================================
      subroutine ene(y, nt, eta, energy)
c     =============================================

      implicit none
      integer nt, i
      real*8 y(nt), eta, energy, kin, pot

      gr_ene = 1.d0 / (2.d0 * eta) !! Ground state energy

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

      energy = gr_ene + kin + pot

      end subroutine ene


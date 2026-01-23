c     ESTIMATES FOR THE COMMON OBSERVABLES OF THE HARMONIC OSCILLATOR  

c     =============================================
      subroutine path_y(y, nt, mean)
c     =============================================
c     Compute the position mean over the path

      implicit none
      integer nt, i
      real*8 y(nt), mean

      mean = 0.d0
      do i = 1, nt
        mean = mean + y(i)
      end do
      mean = mean / dble(nt)

      end subroutine path_y


c     =============================================
      subroutine path_y2(y, nt, mean, y2mean)
c     =============================================
c     Compute the mean of y^2 over the path

      implicit none
      integer nt, i
      real*8 y(nt), mean, y2mean

      y2mean = 0.d0
      do i = 1, nt
        y2mean = y2mean + y(i)**2
      end do
      y2mean = y2mean / dble(nt)
    
      end subroutine path_y2


c     =============================================
      subroutine path_ene(y, nt, eta, energy)
c     =============================================
c     Compute the energy over the path

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




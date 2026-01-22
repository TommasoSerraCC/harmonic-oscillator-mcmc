c     =========================================================
      subroutine heat_bath_step(y, nt, idx, sigma, gamma, alpha,
     &                          eta)
c     =========================================================
c     Performs a heat bath update on the path at index idx
      integer nt, idx
      real y(nt), sigma, gamma, alpha, mu
      
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

      parameter (pi = 3.141592653589793d0)
      real x, mu, sigma
      real u1, u2, z0

      u1 = ran2()
      u2 = ran2()

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


c============================================================================
c  RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
c============================================================================
      function ran2()
      implicit real*4 (a-h,o-z)
      implicit integer*4 (i-n)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real ran2,am,eps,rnmx
      parameter(im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &          ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &          ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,
     &          rnmx=1.-eps)
      integer idum2,j,k,iv,iy
      common /dasav/ idum,idum2,iv(ntab),iy
c      save iv,iy,idum2
c      data idum2/123456789/, iv/NTAB*0/, iy/0/

      if(idum.le.0) then
         idum=max0(-idum,1)
         idum2=idum
         do j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if(idum.lt.0) idum=idum+im1
            if(j.le.ntab) iv(j)=idum
         enddo
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)

      return
      end
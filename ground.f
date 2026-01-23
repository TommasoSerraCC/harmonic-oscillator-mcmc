      subroutine ground_state()

      implicit real*8 (a-h,o-z)
      parameter (bhw=10.d0)  ! beta*h_bar*omega
      parameter (nt=50)      ! number of time slices  
      parameter (nsteps=300000) ! number of MCMC steps for ground state
      parameter (n_conf=30000)   ! number of path configurations to save
      
      integer therm_steps
      integer istart  ! initial configuration flag: 0=cold, 1=hot
      integer i, k, conf_idx, step_interval
      integer idum, idum2, iv(32), iy
      real*8 y(nt)                    ! array of the discretized path
      real*8 path_matrix(n_conf, nt) ! matrix to store path configurations
      real*8 eta                      ! adimensional parameter: eta = a*omega
      real*8 alpha, sigma             ! gaussian parameters for the heat bath
      
      therm_steps = 100000              ! thermalization steps
      step_interval = nsteps / n_conf  ! interval between saved configurations
      
      write(*,*) 'Ground State Histogram Generation'
      write(*,*) 'bhw =', bhw
      write(*,*) 'nt =', nt
      write(*,*) 'nsteps =', nsteps
      write(*,*) 'n_conf =', n_conf
      write(*,*) 'step_interval =', step_interval

      ! Set parameters
      eta = bhw / dble(nt)
      alpha = (eta / 2.d0) + (1.d0 / eta)
      sigma = 1.d0 / sqrt(2.d0 * alpha)

c     Initialize path to zero
      call cold_start(y, nt)

c     Thermalization
      do i = 1, therm_steps
        call total_update(y, nt, sigma, alpha, eta)
      end do

c     Main MCMC loop - collect path configurations
      conf_idx = 1
      do k = 1, nsteps
        call total_update(y, nt, sigma, alpha, eta)
        
c       Save configuration every step_interval steps
        if (mod(k, step_interval) .eq. 0 
     &     .and. conf_idx .le. n_conf) then
          do i = 1, nt
            path_matrix(conf_idx, i) = y(i)
          end do
          conf_idx = conf_idx + 1
        end if
      end do


c     Save path configurations to file
      open(unit=11, file='ground_state_paths.txt', status='unknown')
      write(11,*) '# Ground state path configurations'
      write(11,*) '# n_conf =', n_conf, ' nt =', nt  
      write(11,*) '# bhw =', bhw
      do i = 1, n_conf
        do k = 1, nt
          write(11,*) path_matrix(i, k)
        end do
      end do
      close(11)
      
      end subroutine ground_state


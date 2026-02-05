
      subroutine blocking(x, nt, k, st_dev)

      implicit none
      integer nt, i, j, k, nblocks, valid_nt, is, ie
      real*8 x(nt), block_sum(nt), st_dev, mean_block
      nblocks = nt / k
      valid_nt = nblocks * k  ! Calculate the largest multiple of k less than or equal to nt

      do i = 1, nblocks
        is = (i - 1) * k + 1
        ie = is + k -1
        block_sum(i) = 0.d0
        do j = is, ie
          block_sum(i) = block_sum(i) + x(j)
        end do
      end do

      st_dev = 0.d0
      mean_block = 0.d0
      do i = 1, nblocks
        mean_block = mean_block + block_sum(i)
      end do

      mean_block = mean_block / dble(nblocks)
      do i = 1, nblocks
        st_dev = st_dev + (block_sum(i) - mean_block)**2
      end do

      st_dev = sqrt(st_dev / dble(nblocks - 1))

      end subroutine blocking



      subroutine mean_skip_block(x, nt, k, nblocks,
     &                          idx_start, mean_val)
      implicit none
      integer nt, i, k, idx_start, nblocks, valid_nt
      real*8 x(nt), mean_val

      valid_nt = nblocks * k

      do i = 1, idx_start - 1
        mean_val = mean_val + x(i)
      end do

      do i = idx_start + k, valid_nt
        mean_val = mean_val + x(i)
      end do

      mean_val = mean_val / dble(valid_nt - k)

      end subroutine mean_skip_block



      subroutine compute_mean(x, n, mean_val)
c     ======================
c     Compute mean of array
c     ======================

      implicit none
      integer n, i
      real*8 x(n), mean_val

      mean_val = 0.d0
      do i = 1, n
        mean_val = mean_val + x(i)
      end do
      mean_val = mean_val / dble(n)

      end subroutine compute_mean


c     Generic jackknife that skips one block and calls observable_func
      subroutine jackknife(x1, x2, valid_n, k,
     &                  observable_func, mean_jack, st_dev)

      implicit none
      integer valid_n, i, j, nblocks, idx_start, idx_end, k
      real*8 x1(valid_n), x2(valid_n), st_dev, mean_jack
      real*8 x1_without_block(valid_n - k)
      real*8 x2_without_block(valid_n - k)
      real*8 obs_values(valid_n)
      external observable_func

      nblocks = valid_n / k
      mean_jack = 0.d0

c     For each block, compute observable excluding that block
      do i = 1, nblocks
        idx_start = (i - 1) * k + 1
        idx_end = idx_start + k - 1

c       Copy data excluding current block
        j = 1
        do j = 1, idx_start - 1
          x1_without_block(j) = x1(j)
          x2_without_block(j) = x2(j)
        end do
        do j = idx_end + 1, valid_n
          x1_without_block(j - k) = x1(j)
          x2_without_block(j - k) = x2(j)
        end do

c       Call observable function on data without block
        call observable_func(x1_without_block, x2_without_block,
     &                       valid_n - k, obs_values(i))
        mean_jack = mean_jack + obs_values(i)
      end do

      mean_jack = mean_jack / dble(nblocks)

c     Compute standard deviation
      st_dev = 0.d0
      do i = 1, nblocks
        st_dev = st_dev + (obs_values(i) - mean_jack)**2
      end do

      st_dev = sqrt(st_dev * dble(nblocks - 1) / dble(nblocks))

      end subroutine jackknife
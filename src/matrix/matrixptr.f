      subroutine matrixptr(a, n)
      implicit none
      integer, intent(in) :: n
      integer, intent(inout) :: a(n)
      integer i
      do i = 1, n
        write(*,*) 'a[i] = ', i, a(i)
      end do

      do i = 1, n
        a(i) = a(i) + 2
      end do

      return
      end subroutine matrixptr
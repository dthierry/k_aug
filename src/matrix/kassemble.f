      subroutine kassemble(aij, acol, arow, nza, nvar)
      implicit none
      integer, intent(in) :: acol(nza), arow(nza), nza, nvar
      real*8, intent(in) :: aij(nza)
      integer i, nnumber
      
      write(*, *) "Jacobian matrix with ", nza, "Elements"
      write(*, *) "Check zeroes in the main d "
      

      nnumber = 0
      do i = 1, nza
        if (acol(i) == arow(i)) then
            nnumber = nnumber + 1
        end if
      end do
      write(*,*) "Found zeores", nnumber
      write(*,*) "n_var ", nvar
      if (nnumber == nvar) then
        write(*,*) "There are zeroes in the whole main diag"
      else
        write(*,*) "Not enough zeroes"
      end if

      return
      end subroutine kassemble

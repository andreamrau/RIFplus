SUBROUTINE test(x, n)

      integer i
      integer n
      double precision x(n)

      do i = 1, n
          x(i) = x(i) ** 2
    enddo

end
program testf90
    use umatutils, only: m33eigvalsh
    implicit none
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: f(3, 3), c(3, 3)
    integer :: i
    f = reshape([2.53680432, 0.73945601, -0.4530953 ,  0.86524763,  1.45320696,&
        0.75649414, -0.21189672,  0.48545667,  2.69808608], [3, 3])      
    c = matmul(transpose(f), f)
    write (*, *) sqrt(m33eigvalsh(c))

contains
end program testf90
program testf90
    use umatutils, only: m33inv
    implicit none
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: m(3, 3)
    m = reshape([1, 2, 3, 3, 2, 1, 2, 5, 7], [3, 3])
    write (*, *) m33inv(m)

contains
end program testf90
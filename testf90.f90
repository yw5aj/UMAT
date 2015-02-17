program testf90
    implicit none
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: a(2)
    a = [1._dp, 2._dp]
    write (*, *) a(1::2)
contains
        
end program testf90
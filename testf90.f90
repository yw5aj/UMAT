program testf90
    implicit none
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: a(6)=[0, 1, 2, 3, 4, 5]
    call alterslice(a)
    write(*, *) ('asdf'=='asdf').eqv..true.

contains
    subroutine alterslice(a)
        real(dp), intent(inout) :: a(:)
        a(1) = a(1) + 1
    end subroutine alterslice
        
end program testf90
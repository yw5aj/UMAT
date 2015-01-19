program testf90
    implicit none
    integer, parameter :: dp=kind(0.d0)
    real(dp) :: a
    integer :: i, j
    do i = 1, 3
        do j = 1, i
            write(*, *) i, j
            if(i /= j) then
                write(*, *) j, i
            end if
        end do
    end do

contains
        
end program testf90
module umatutils
    !!! Utility module to be used for UMATs
    implicit none
    private
    public dp, delta, m31tensorprod, m33det, m33tensorprod, mapnotation
    integer, parameter :: dp=kind(0.d0)
    real(dp), parameter :: delta(3, 3)=reshape([1, 0, 0, 0, 1, 0, 0, 0, 1],&
        [3, 3])
        
contains
    subroutine mapnotation(sigma, ccj, ntens, stress, ddsdde)
        !! Map the matrix notation to vector notation
        real(dp), intent(in) :: sigma(3, 3), ccj(3, 3, 3, 3)
        integer, intent(in) :: ntens
        real(dp), intent(out) :: stress(ntens), ddsdde(ntens, ntens)
        integer, parameter :: notationmap(6, 2)=reshape([1, 2, 3, 1, 1, 2, 1,&
            2, 3, 2, 3, 3], [6, 2])
        integer :: k1, k2
        do k1 = 1, 6
            stress(k1) = sigma(notationmap(k1, 1), notationmap(k1, 2))    
            do k2 = 1, 6
                ddsdde(k1, k2) = ccj(notationmap(k1, 1), notationmap(k1, 2),&
                    notationmap(k2, 1), notationmap(k2, 2))
            end do
        end do
    end subroutine mapnotation
    
    function m31tensorprod(a, b)
        !! Return the tensor product from two 3x1 vectors
        real(dp), dimension(3), intent(in) :: a, b
        real(dp) :: m31tensorprod(3, 3)
        integer :: k1, k2
        forall(k1=1:3, k2=1:3)
            m31tensorprod(k1, k2) = a(k1) * b(k2)
        end forall
        return
    end function m31tensorprod
    
    function m33tensorprod(a, b)
        !! Return the tensor product from two 3x3 matrices
        real(dp), dimension(3, 3), intent(in) :: a, b
        real(dp) :: m33tensorprod(3, 3, 3, 3)
        integer :: k1, k2, k3, k4
        forall(k1=1:3, k2=1:3, k3=1:3, k4=1:3)
            m33tensorprod(k1, k2, k3, k4) = a(k1, k2) * b(k3, k4)
        end forall
        return
    end function m33tensorprod
    
    real(dp) function m33det(m33)
        !! Return determinant of a 3 by 3 matrix
        real(dp), intent(in) :: m33(3, 3)
        m33det = m33(1, 1) * m33(2, 2) * m33(3, 3)&
            - m33(1, 1) * m33(2, 3) * m33(3, 2)&
            - m33(1, 2) * m33(2, 1) * m33(3, 3)&
            + m33(1, 2) * m33(2, 3) * m33(3, 1)&
            + m33(1, 3) * m33(2, 1) * m33(3, 2)&
            - m33(1, 3) * m33(2, 2) * m33(3, 1)
        return
    end function m33det
end module umatutils

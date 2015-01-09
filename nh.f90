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


module numerichyper
    !!! Numeric method for hyperelastic materials
    !!! Module to be embedded in the UMAT
    !!! Paramters
    !!! ---------
    !!! det : Jacobian
    !!! dfgrd : deformation gradient
    !!! fbar : modified deformation gradient
    !!! rcg : right Cauchy-Green tensor
    !!! rcgbar : modified right Cauchy-Green tensor
    !!! ibar1 : first invariant of rcgbar
    !!! sigma : Cauchy stress
    !!! ccj : tangent modulus for Cauchy stress in Jaumann rate    
    use umatutils, only: dp, delta, m33det, m31tensorprod, mapnotation
    implicit none
    private
    public update_stress_ddsdde

contains
    function getpsi(rcg, props) result(psi)
        !! Return strain energy density given C and material properties
        real(dp), intent(in) :: rcg(3, 3), props(:)
        real(dp) :: rcgbar(3, 3), det, ibar1, psi
        real(dp) :: c10, d1 ! Specific to Neo-Hookean model
        integer :: k1
        c10 = props(1)
        d1 = props(2)
        det = sqrt(m33det(rcg))
        rcgbar = det**(-2._dp/3) * rcg
        ibar1 = sum([(rcgbar(k1, k1), k1=1, 3)])
        psi = c10 * (ibar1 - 3) + (det - 1)**2 / d1
        return
    end function getpsi

    subroutine update_stress_ddsdde(props, dfgrd, ntens, stress, ddsdde)
        !! Update the stress and ddsdde for Neo-Hookean
        real(dp), intent(in) :: props(:), dfgrd(3, 3)
        integer, intent(in) :: ntens
        real(dp), intent(out) :: stress(ntens), ddsdde(ntens, ntens)
        real(dp), parameter :: eps_s=1d-4, eps_c=1d-6
        real(dp) :: sigma(3, 3), ccj(3, 3, 3, 3)
        call get_sigma_ccj(dfgrd, props, eps_c, eps_s, sigma, ccj)
        call mapnotation(sigma, ccj, ntens, stress, ddsdde)
    end subroutine update_stress_ddsdde
        
    function gettau(dfgrd, props, eps_s) result(tau)
        !! Return the Cauchy stress given deformation gradient
        real(dp), intent(in) :: dfgrd(3, 3), props(:), eps_s
        real(dp) :: tau(3, 3), rcg(3, 3), psi, rcgptb(3, 3), psiptb, pk2(3, 3)
        integer :: k1, k2
        rcg = matmul(transpose(dfgrd), dfgrd)
        psi = getpsi(rcg, props)
        ! Outer two layers are indices for the 2nd PK stress
        do k1 = 1, 3
        do k2 = 1, 3
            rcgptb = rcg + eps_s * (&
                m31tensorprod(delta(:, k1), delta(:, k2)) +&
                m31tensorprod(delta(:, k2), delta(:, k1)))
            psiptb = getpsi(rcgptb, props)
            pk2(k1, k2) = (psiptb - psi) / eps_s
        end do
        end do
        ! Calculate Kirchoff stress from 2nd PK stress
        tau = matmul(matmul(dfgrd, pk2), transpose(dfgrd))
        return
        end function gettau

    subroutine get_sigma_ccj(dfgrd, props, eps_c, eps_s, sigma, ccj)
        real(dp), intent(in) :: dfgrd(3, 3), props(:), eps_c, eps_s
        real(dp), intent(out) :: sigma(3, 3), ccj(3, 3, 3, 3)
        real(dp) :: tau(3, 3), tauptb(3, 3), fptb(3, 3), det
        integer :: k3, k4
        det = m33det(dfgrd)
        tau = gettau(dfgrd, props, eps_s)
        sigma = tau / det ! Pass out sigma
        ! Use k3 & k4 rather than k1 & k2 to denote that it's the k, l 
        ! component of the elasticity tensor being calculated
        do k3 = 1, 3
        do k4 = 1, 3
            fptb = dfgrd + eps_c / 2 * (&
                matmul(m31tensorprod(delta(:, k3), delta(:, k4)), dfgrd) +&
                matmul(m31tensorprod(delta(:, k4), delta(:, k3)), dfgrd))
            tauptb = gettau(fptb, props, eps_s)
            ccj(:, :, k3, k4) = (tauptb - tau) / eps_c / det
        end do
        end do
        end subroutine get_sigma_ccj
end module numerichyper


program nh
    !!! Testing numerical method for NH model
    use numerichyper, only: update_stress_ddsdde
    use umatutils, only: dp
    implicit none
    real(dp) :: d1, c10, props(2)
    integer, parameter :: ntens=6
    real(dp) :: stress(ntens), ddsdde(ntens, ntens), dfgrd(3, 3)
    ! Assign initial parameters
    dfgrd(1, 1) = 1.5488135
    dfgrd(1, 2) = 0.71518937
    dfgrd(1, 3) = 0.60276338
    dfgrd(2, 1) = 0.54488318
    dfgrd(2, 2) = 1.4236548
    dfgrd(2, 3) = 0.64589411
    dfgrd(3, 1) = 0.43758721
    dfgrd(3, 2) = 0.891773
    dfgrd(3, 3) = 1.96366276    
    ! dfgrd = reshape([1._dp, 0._dp, 0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 0._dp, 1._dp], [3, 3])
    ! dfgrd = dfgrd * 1.6d0
    c10 = 8d4
    d1 = 2d-1
    props = [c10, d1]
    ! Calculate
    call update_stress_ddsdde(props, dfgrd, ntens, stress, ddsdde)
    write(*, *) 'Stress: ', stress
    
contains
    !!! Does not contain any subroutine / functions
end program nh






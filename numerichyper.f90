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
    use modpsi, only: getpsi
    implicit none
    private
    public update_umat

contains
    subroutine update_umat(props, dfgrd, ntens, stress, ddsdde, &
        statev)
        !! Update the stress and ddsdde for Neo-Hookean
        real(dp), intent(in) :: props(:), dfgrd(3, 3)
        integer, intent(in) :: ntens
        real(dp), intent(out) :: stress(ntens), ddsdde(ntens, ntens)
        real(dp), intent(inout) :: statev(:)
        real(dp), parameter :: eps_s=1d-6, eps_c=1d-4
        real(dp) :: sigma(3, 3), ccj(3, 3, 3, 3)
        call get_sigma_ccj(dfgrd, props, eps_c, eps_s, sigma, ccj, statev)
        call mapnotation(sigma, ccj, ntens, stress, ddsdde)
    end subroutine update_umat
        
    function gettau(dfgrd, props, statev, eps_s) result(tau)
        !! Return the Cauchy stress given deformation gradient
        real(dp), intent(in) :: dfgrd(3, 3), props(:), statev(:), eps_s
        real(dp) :: tau(3, 3), rcg(3, 3), psi, rcgptb(3, 3), psiptb, pk2(3, 3)
        integer :: k1, k2
        rcg = matmul(transpose(dfgrd), dfgrd)
        psi = getpsi(rcg, props, statev)
        ! Outer two layers are indices for the 2nd PK stress
        do k1 = 1, 3
        do k2 = 1, 3
            rcgptb = rcg + eps_s * (&
                m31tensorprod(delta(:, k1), delta(:, k2)) +&
                m31tensorprod(delta(:, k2), delta(:, k1)))
            psiptb = getpsi(rcgptb, props, statev)
            pk2(k1, k2) = (psiptb - psi) / eps_s
        end do
        end do
        ! Calculate Kirchoff stress from 2nd PK stress
        tau = matmul(matmul(dfgrd, pk2), transpose(dfgrd))
        return
        end function gettau

    subroutine get_sigma_ccj(dfgrd, props, eps_c, eps_s, sigma, ccj, statev)
        real(dp), intent(in) :: dfgrd(3, 3), props(:), eps_c, eps_s
        real(dp), intent(out) :: sigma(3, 3), ccj(3, 3, 3, 3)
        real(dp), intent(inout) :: statev(:)
        real(dp) :: tau(3, 3), tauptb(3, 3), fptb(3, 3), det
        integer :: k3, k4
        det = m33det(dfgrd)
        tau = gettau(dfgrd, props, statev, eps_s)
        sigma = tau / det ! Pass out sigma
        ! Use k3 & k4 rather than k1 & k2 to denote that it's the i, j 
        ! component of the elasticity tensor being calculated
        do k3 = 1, 3
        do k4 = 1, 3
            fptb = dfgrd + eps_c / 2 * (&
                matmul(m31tensorprod(delta(:, k3), delta(:, k4)), dfgrd) +&
                matmul(m31tensorprod(delta(:, k4), delta(:, k3)), dfgrd))
            tauptb = gettau(fptb, props, statev, eps_s)
            ccj(:, :, k3, k4) = (tauptb - tau) / eps_c / det
        end do
        end do
        end subroutine get_sigma_ccj
end module numerichyper

module modpsi
    !!! Module for the function getpsi for Holzapfel-Gasser-Ogden model
    use umatutils, only: dp, m33det
    implicit none
    private
    public getpsi

contains
    function getpsi(rcg, props, statev) result(psi)
        !! Return strain energy density given C and material properties
        real(dp), intent(in) :: rcg(3, 3), props(:), statev(:)
        real(dp) :: rcgbar(3, 3), det, ibar1, psi, psi_iso, psi_ani, psi_vol,&
            a1(3), a2(3), ebar1, ebar2, ibar411, ibar422
        real(dp) :: c10, d1, k1, k2, kappa !Material properties
        integer :: i
        ! Assign variable to more legible names
        c10 = props(1)
        d1 = props(2)
        k1 = props(3)
        k2 = props(4)
        kappa = props(5)
        a1 = statev(1:3)
        a2 = statev(4:6)
        ! Get isotropic isochoric response
        det = sqrt(m33det(rcg))
        rcgbar = det**(-2._dp/3) * rcg
        ibar1 = sum([(rcgbar(i, i), i=1, 3)])
        psi_iso = c10 * (ibar1 - 3)
        ! Get isotropic volumetric response
        psi_vol = ((det**2 - 1) / 2 - log(det)) / d1
        ! Get anisotropic response
        ibar411 = dot_product(a1, matmul(rcgbar, a1))
        ibar422 = dot_product(a2, matmul(rcgbar, a2))
        ebar1 = kappa * (ibar1 - 3) + (1 - 3*kappa) * (ibar411 - 1)
        ebar2 = kappa * (ibar1 - 3) + (1 - 3*kappa) * (ibar422 - 1)
        psi_ani = k1 / 2 / k2 * (&
            (exp(k2 * (ebar1 + abs(ebar1))**2 / 4) - 1) + &
            (exp(k2 * (ebar2 + abs(ebar2))**2 / 4) - 1))
        ! Summate to get psi
        psi = psi_iso + psi_ani + psi_vol
        return
    end function getpsi
end module modpsi
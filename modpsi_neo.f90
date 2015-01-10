module modpsi
    !!! Module for the function getpsi for Neo-Hookean material
    use umatutils, only: dp, m33det
    implicit none
    private
    public getpsi

contains
    function getpsi(rcg, props, statev) result(psi)
        !! Return strain energy density given C and material properties
        real(dp), intent(in) :: rcg(3, 3), props(:), statev(:)
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
end module modpsi
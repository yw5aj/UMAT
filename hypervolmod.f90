module hypervolmod
    !!! Module for the volumetric part hyperelastic material
    !!! Psi_vol = Sum(i=1, N){1/Di*(J-1)**(2i)}
    use umatutils, only: dp, delta, m33det, ii, ccc2ccj
    implicit none
    private
    public hypervol

contains
    subroutine hypervol(f, d, nterms, sigma, ccj)
        integer, intent(in) :: nterms
        real(dp), intent(in) :: f(3, 3), d(nterms)
        real(dp), intent(out) :: sigma(3, 3), ccj(3, 3, 3, 3)
        integer :: i, j, k, l, n
        real(dp) :: p = 0, ptilde = 0, det, ccc(3, 3, 3, 3)
        det = m33det(f)
        do n = 1, nterms
            p = p + 2 * n * (det - 1) ** (2 * n - 1) / d(n)
            ptilde = ptilde + 2 * n * (det - 1) ** (2 * n - 2) * &
                (2 * n * det - 1) / d(n)
        end do
        do i = 1, 3
            do j = i, 3
                sigma(i, j) = p * delta(i, j)
                ! Fill symmetric part
                if (i /= j) then
                    sigma(j, i) = sigma(i, j)
                end if
                do k = 1, 3
                    do l = k, 3
                        ccc(i, j, k, l) = ptilde * delta(i, j) * delta(k, l) &
                            - 2 * p * ii(i, j, k, l)
                        ! Fill symmetric part k, l interchange
                        if (k /= l) then
                            ccc(i, j, l, k) = ccc(i, j, k, l)
                        end if
                    end do
                ! Fill symmetric part i, j interchange
                if (i /= j) then
                    ccc(j, i, :, :) = ccc(i, j, :, :)
                end if
                end do
            end do
        end do
        ! Switch to Jaumann rate
        ccj = ccc2ccj(ccc, sigma)
        end subroutine hypervol
end module hypervolmod

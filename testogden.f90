program testogden
    use umatutils, only: dp, delta, m33det, m33eigval, ii, eps
    implicit none
    real(dp) :: dfgrd1(3, 3), props(2), sigma(3, 3), ccj(3, 3, 3, 3)
    dfgrd1 = reshape([1., 0., 0., 0., 1., 0., .45, 0., 1.], [3, 3])
    props = [160e3_dp, 2._dp]
    call hyperiso(dfgrd1, props, size(props)/2, sigma, ccj)

contains
    subroutine hyperiso(f, props, nterms, sigma, ccj)
        integer, intent(in) :: nterms
        real(dp), intent(in) :: f(3, 3), props(:)
        real(dp), intent(out) :: sigma(3, 3), ccj(3, 3, 3, 3)
        integer :: i, j, k, l, n
        real(dp) :: mu(nterms), alpha(nterms), b(3, 3), lam2(3), lambar(3),&
            lbpow(3, nterms), det, lam(3), beta(3), gamma(3, 3), m(3, 3, 3),&
            d(3), dprime(3), i1, i3, ib(3, 3, 3, 3), dmterm1, dmterm2, dmterm3,&
            dm(3, 3, 3, 3, 3)
        mu = props(1::2)
        alpha = props(2::2)
        det = m33det(f)
        b = matmul(f, transpose(f))
        lam2 = m33eigval(b)
        lam = sqrt(lam2)
        lambar = lam * det**(-1._dp/3)
        ! I1, I3, and Ib
        i1 = b(1, 1) + b(2, 2) + b(3, 3)
        i3 = det ** 2
        do i = 1, 3
            do j = i, 3
                do k = 1, 3
                    do l = k, 3
                        ib(i, j, k, l) = (b(i, k) * b(j, l) +&
                            b(i, l) * b(j, k)) / 2
                        if (k /= l) then
                            ib(i, j, l, k) = ib(i, j, k, l)
                        end if
                    end do
                end do
                if (i /= j) then
                    ib(j, i, :, :) = ib(i, j, :, :)
                end if
            end do
        end do
        ! Di, D'i, mi
        d = 2 * lam2**2 - i1 * lam2 + i3 / lam2
        dprime = 8 * lam**3 - 2 * i1 * lam - 2 * i3 / lam**3
        do i = 1, 3
            m(:, :, i) = (matmul(b, b) - (i1 - lam2(i)) * b +&
                i3 / lam2(i) * delta) / d(i)
        end do
        ! dm_i
        do n = 1, 3
            do i = 1, 3
                do j = i, 3
                    do k = 1, 3
                        do l = k, 3
                            dmterm1 = ib(i, j, k, l) - b(i, j) * b(k, l)&
                                + i3 / lam2(n) * (delta(i, j) * delta(k, l)&
                                - ii(i, j, k, l))
                            dmterm2 = lam2(n) * (b(i, j) * m(k, l, n)&
                                + m(i, j, n) * b(k, l)) - dprime(n)/2 * lam(n)&
                                * m(i, j, n) * m(k, l, n)
                            dmterm3 = i3 / lam2(n) * (delta(i, j) * m(k, l, n)&
                                + m(i, j, n) * delta(k, l))
                            dm(i, j, k, l, n) = (dmterm1+dmterm2-dmterm3)/d(n)
                            if (k /= l) then
                                dm(i, j, l, k, n) = dm(i, j, k, l, n)
                            end if
                        end do
                    end do
                    if (i /= j) then
                        dm(j, i, :, :, n) = dm(i, j, :, :, n)
                    end if
                end do
            end do        
        end do
        ! lambdabar_alpha, 3 x nterms
        do n = 1, nterms
            lbpow(:, n) = lambar**(alpha(n))
        end do
        ! beta_i
        do i = 1, 3
            beta(i) = 0._dp
            do n = 1, nterms
                beta(i) = beta(i) + mu(n) * (lbpow(i, n) - sum(lbpow(:, n))/3)
            end do
        end do
        ! gamma_ij
        do i = 1, 3
            do j = i, 3
                gamma(i, j) = 0._dp
                do n = 1, nterms
                    if (i == j) then
                        gamma(i, j) = gamma(i, j) + mu(n) * alpha(n) * (&
                            lbpow(i, n)/3 + sum(lbpow(:, n))/9)
                    elseif (i /= j) then
                        gamma(i, j) = gamma(i, j) + mu(n) * alpha(n) * (&
                            -lbpow(i, n)/3 - lbpow(j, n)/3 +&
                            sum(lbpow(:, n))/9)
                    end if
                end do
                if (i /= j) then
                    gamma(j, i) = gamma(i, j)
                end if
            end do
        end do
        
        ! Output for testing
        write (*, *) dm(:, :, :, :, 1)
    end subroutine hyperiso
end program testogden
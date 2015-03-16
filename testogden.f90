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


module hyperisomod
    !!! Module for the isochoric part of Ogden hyperelastic material
    !!! Used definition for Holzapfel's book
    use umatutils, only: dp, delta, m33det, m33eigval, ii, eps, ccc2ccj, m33inv
    implicit none
    private
    public hyperiso

contains
    subroutine hyperiso(f, isoprops, nterms, sigma, ccj)
        !! isoprops : mu1, alpha1, mu2, alpha2, etc.
        integer, intent(in) :: nterms
        real(dp), intent(in) :: f(3, 3), isoprops(:)
        real(dp), intent(out) :: sigma(3, 3), ccj(3, 3, 3, 3)
        integer :: i, j, k, l, n, k1, k2, a, c
        real(dp) :: mu(nterms), alpha(nterms), b(3, 3), lam2(3), lambar(3),&
            lbpow(3, nterms), det, lam(3), beta(3), gamma(3, 3), m(3, 3, 3),&
            d(3), dprime(3), i1, i3, ib(3, 3, 3, 3), dmterm1, dmterm2, dmterm3,&
            dm(3, 3, 3, 3, 3), ccc(3, 3, 3, 3)
        mu = isoprops(1::2)
        alpha = isoprops(2::2)
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
                        ! Fill symmetric part, k, l interchange
                        if (k /= l) then
                            ib(i, j, l, k) = ib(i, j, k, l)
                        end if
                    end do
                end do
                ! Fill symmetric part, i, j interchange
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
                            ! Fill symmetric part k, l interchange
                            if (k /= l) then
                                dm(i, j, l, k, n) = dm(i, j, k, l, n)
                            end if
                        end do
                    end do
                    ! Fill symmetric part i, j interchange
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
        beta = 0
        do i = 1, 3
            do n = 1, nterms
                beta(i) = beta(i) + mu(n) * (lbpow(i, n) - sum(lbpow(:, n))/3)
            end do
        end do
        ! gamma_ij
        gamma = 0
        do i = 1, 3
            do j = i, 3
                do n = 1, nterms
                    if (i == j) then
                        gamma(i, j) = gamma(i, j) + mu(n) * alpha(n) * (&
                            lbpow(i, n)/3 + sum(lbpow(:, n))/9)
                    else if (i /= j) then
                        gamma(i, j) = gamma(i, j) + mu(n) * alpha(n) * (&
                            -lbpow(i, n)/3 - lbpow(j, n)/3 +&
                            sum(lbpow(:, n))/9)
                    end if
                end do
                ! Fill symmetric part
                if (i /= j) then
                    gamma(j, i) = gamma(i, j)
                end if
            end do
        end do
        ! Put everything together, for three distinct cases
        sigma = 0
        ccc = 0
        if ((abs(lam(1) - lam(2)) >= eps) .and. (abs(lam(1) - lam(3)) >= eps)&
            .and. (abs(lam(2) - lam(3)) >= eps)) then
            ! If all three principal stretches are different
            ! Get stress first
            do i = 1, 3
                do j = i, 3
                    do k = 1, 3
                        sigma(i, j) = sigma(i, j) + beta(k) * m(i, j, k) / det
                    end do
                    ! Fill symmetric part
                    if (i /= j) then
                        sigma(j, i) = sigma(i, j)
                    end if
                end do
            end do
            ! Get elasticity tensor for Cauchy stress, convective rate
            do i = 1, 3
                do j = i, 3
                    do k = 1, 3
                        do l = k, 3
                            do k1 = 1, 3
                                ccc(i, j, k, l) = ccc(i, j, k, l) +&
                                    2 * beta(k1) * dm(i, j, k, l, k1) / det
                                do k2 = 1, 3
                                    ccc(i, j, k, l) = ccc(i, j, k, l) +&
                                        gamma(k1, k2) * m(i, j, k1) *&
                                        m(k, l, k2) / det
                                end do
                            end do
                            ! Fill symmetric part, k, l interchange
                            if (k /= l) then
                                ccc(i, j, l, k) = ccc(i, j, k, l)
                            end if
                        end do    
                    end do
                    ! Fill symmetric part, i, j interchange
                    if (i /= j) then
                        ccc(j, i, :, :) = ccc(i, j, :, :)
                    end if
                end do
            end do
        else if (abs(lam(1) - lam(2)) < eps .and. abs(lam(1) - lam(3)) < eps)&
            then
            ! All three are the same
            do i = 1, 3
                do j = i, 3
                    do k = 1, 3
                        do l = k, 3
                            ccc(i, j, k, l) = gamma(1, 1) * 1.5_dp / det * (&
                                ii(i, j, k, l) - delta(i, j) *&
                                delta(k, l) / 3)
                            if (k /= l) then
                                ccc(i, j, l, k) = ccc(i, j, k, l)
                            end if
                        end do
                    end do
                    if (i /= j) then
                        ccc(j, i, :, :) = ccc(i, j, :, :)
                    end if
                end do
            end do
        else
            ! Two are the same
            ! First simplify into one case
            if (abs(lam(1) - lam(2)) < eps) then
                c = 3
            else if (abs(lam(1) - lam(3)) < eps) then
                c = 2
            else if (abs(lam(2) - lam(3)) < eps) then
                c = 1
            end if
            a = mod(c + 1, 3)
            write (*, *) lam
            ! Plug in
            do i = 1, 3
                do j = i, 3
                    sigma(i, j) = (beta(a) * delta(i, j) + (beta(c) - beta(a))&
                        * m(i, j, c)) / det
                    if (i /= j) then
                        sigma(j, i) = sigma(i, j)
                    end if
                    do k = 1, 3
                        do l = k, 3
                            ccc(i, j, k, l) = (&
                                -beta(a) * 2 * ii(i, j, k, l)&
                                + (beta(c) - beta(a)) * 2 * dm(i, j, k, l, c)&
                                + gamma(a, a) * (&
                                    delta(i, j) - m(i, j, c)) * (&
                                    delta(k, l) - m(k, l, c))&
                                + gamma(c, c) * m(i, j, c) * m(k, l, c)&
                                + gamma(a, c) * (&
                                    m(i, j, c) * (delta(k, l) - m(k, l, c))&
                                    + (delta(i, j) - m(i, j, c)) * m(k, l, c)&
                                    / det))&
                                / det
                            if (k /= l) then
                                ccc(i, j, l, k) = ccc(i, j, k, l)
                            end if
                        end do
                    end do
                    if (i /= j) then
                        ccc(j, i, :, :) = ccc(i, j, :, :)
                    end if
                end do
            end do
        end if
        ! Switch to Jaumann rate
        ccj = ccc2ccj(ccc, sigma)
    end subroutine hyperiso
end module hyperisomod


module ogdenmod
    use umatutils, only: dp, delta, m33det, m33inv
    use hyperisomod, only: hyperiso
    use hypervolmod, only: hypervol
    implicit none
    private
    public ogden, ogdenpk2

contains
    subroutine ogden (f, props, nterms, sigma, ccj)
        real(dp), intent(in) :: f(3, 3), props(:)
        integer, intent(in) :: nterms
        real(dp), intent(out) :: sigma(3, 3), ccj(3, 3, 3, 3)
        real(dp) :: sigmaiso(3, 3), sigmavol(3, 3), propsiso(2*nterms),&
            propsvol(nterms), ccjiso(3, 3, 3, 3), ccjvol(3, 3, 3, 3)
        ! Assign material properties
        propsiso = props(:2*nterms)
        propsvol = props(2*nterms+1:)
        ! Calculate both parts separately
        call hyperiso(f, propsiso, nterms, sigmaiso, ccjiso)
        call hypervol(f, propsvol, nterms, sigmavol, ccjvol)
        ! Add them up
        sigma = sigmaiso + sigmavol
        ccj = ccjiso + ccjvol
    end subroutine ogden
    
    subroutine ogdenpk2 (f, props, nterms, sigma, ccj, siso, svol)
        !! subroutine to include output of pk2 stress
        real(dp), intent(in) :: f(3, 3), props(:)
        integer, intent(in) :: nterms
        real(dp), intent(out) :: sigma(3, 3), ccj(3, 3, 3, 3), siso(3, 3),&
            svol(3, 3)
        real(dp) :: sigmaiso(3, 3), sigmavol(3, 3), det, finv(3, 3)
        integer :: i
        ! Get the spatial representation
        call ogden (f, props, nterms, sigma, ccj)
        sigmavol = sum([(sigma(i, i), i = 1, 3)]) / 3 * delta
        sigmaiso = sigma - sigmavol
        write (*, *) sigmaiso
        ! Convert to 2nd PK stress
        det = m33det(f)
        finv = m33inv(f)
        siso = det * matmul(matmul(finv, sigmaiso), transpose(finv))
        svol = det * matmul(matmul(finv, sigmavol), transpose(finv))
    end subroutine ogdenpk2
end module ogdenmod

program testogden
    use umatutils, only: dp
    use hyperisomod, only: hyperiso
    use hypervolmod, only: hypervol
    use ogdenmod, only: ogden, ogdenpk2
    implicit none
    real(dp) :: dfgrd1(3, 3), props(3), sigma(3, 3), ccj(3, 3, 3, 3),&
        siso(3, 3), svol(3, 3), f2same(3, 3), f3same(3, 3)
    ! dfgrd1 = reshape([1., 0., 0., 0., 1., 0., .45, 0., 1.], [3, 3])
    ! dfgrd1 = reshape([2., 0., 0., 0., 2., 0., 0., 0., 2.], [3, 3])
    dfgrd1 = reshape([1.5488135, 0.54488318, 0.43758721, 0.71518937,&
        1.4236548, 0.891773, 0.60276338, 0.64589411, 1.96366276], [3, 3])
    f2same = reshape([2.53680432, 0.73945601, -0.4530953 ,  0.86524763,  1.45320696,&
        0.75649414, -0.21189672,  0.48545667,  2.69808608], [3, 3])        
    props = [160e3_dp, 2._dp, .2_dp]
    ! call hyperiso(dfgrd1, props(:size(props)*2/3), size(props)/3, sigma, ccj)
    ! call hypervol(dfgrd1, reshape([.2_dp], [1]), 1, sigma, ccj)
    call ogdenpk2(f2same, props, size(props)/3, sigma, ccj, siso, svol)
end program testogden
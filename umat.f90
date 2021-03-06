module umatutils
    !!! Utility module to be used for UMATs
    implicit none
    private
    public dp, delta, m31tensorprod, m33det, m33tensorprod, mapnotation, pi,&
        eps, m33eigvalsh, m33eigvect, ii, ccc2ccj, m33inv
    integer, parameter :: dp=kind(0.d0)
    real(dp), parameter :: delta(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1],&
        [3, 3]), eps = 1e-6_dp, pi = 4*atan(1._dp), ii(3, 3, 3, 3) = reshape([&
        1._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp,&
        0._dp, 0.5_dp, 0._dp, 0.5_dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp,&
        0._dp, 0._dp, 0.5_dp, 0._dp, 0._dp, 0._dp, 0.5_dp, 0._dp, 0._dp,&
        0._dp, 0.5_dp, 0._dp, 0.5_dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp,&
        0._dp, 0._dp, 0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 0._dp, 0._dp,&
        0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0.5_dp, 0._dp, 0.5_dp, 0._dp,&
        0._dp, 0._dp, 0.5_dp, 0._dp, 0._dp, 0._dp, 0.5_dp, 0._dp, 0._dp,&
        0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0.5_dp, 0._dp, 0.5_dp, 0._dp,&
        0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 1._dp],&
        [3, 3, 3, 3])
        
contains
    function m33inv(a) result(ainv)
        !! Inverse the real non-symmetric 3x3 matrix
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: ainv(3, 3), cofactor(3, 3), det
        det = m33det(a)
        cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
        cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
        cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
        cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
        cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
        cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
        cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
        cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
        cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))
        ainv = transpose(cofactor) / det
    end function m33inv
    
    function ccc2ccj(ccc, sigma) result(ccj)
        !! Convert ccc to ccj given sigma
        real(dp), intent(in) :: ccc(3, 3, 3, 3), sigma(3, 3)
        real(dp) :: ccj(3, 3, 3, 3)
        integer :: i, j, k, l
        do i = 1, 3
            do j = i, 3
                do k = 1, 3
                    do l = k, 3
                        ccj(i, j, k, l) = ccc(i, j, k, l) + (&
                            delta(i, k) * sigma(j, l) + &
                            delta(i, l) * sigma(j, k) + &
                            delta(j, k) * sigma(i, l) + &
                            delta(j, l) * sigma(i, k)) / 2
                        ! Fill l/k symmetry
                        if (l /= k) then
                            ccj(i, j, l, k) = ccj(i, j, k, l)
                        end if
                    end do
                end do
                ! Fill i/j symmetry
                if (i /= j) then
                    ccj(j, i, :, :) = ccj(i, j, :, :)
                end if
            end do
        end do
        return
    end function ccc2ccj

    function m33eigvalsh(a) result(w)
        !! Returns the eigenvalues of a 3x3 real symmetric matrix
        !! Unless lambda1 = lambda2 = lambda3, lambda1 > lambda2
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: w(3), m, p, q, phi, a_mi(3, 3), c, s, sqp, sq3, p3_q2
        m = (a(1, 1) + a(2, 2) + a(3, 3)) / 3
        a_mi = a - m * delta
        q = m33det(a_mi) / 2
        p = sum(a_mi ** 2) / 6
        p3_q2 = p ** 3 - q ** 2
        if (abs(p - q) < eps) then
            w = m
        else if (abs(p3_q2) < eps) then
            sqp = sqrt(p)
            w(1) = m - 2 * sqp
            w(2:3) = m + sqp
        else
            phi = atan(sqrt(p3_q2) / q) / 3
            if (phi <= 0) then
                phi = pi + phi
            end if
            c = cos(phi)
            s = sin(phi)
            sqp = sqrt(p)
            sq3 = sqrt(3._dp)
            w(1) = m + 2 * sqp * c
            w(2) = m - sqp * (c + sq3 * s)
            w(3) = m - sqp * (c - sq3 * s)
        end if
        return
    end function m33eigvalsh
    
    function m31cross(a, b) result(c)
        !! Computes the cross product between two 3x1 vectors
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
        return
    end function m31cross
    
    function m31colinear(a, b) result(mu)
        !! If a and b are colinear, mu = |a|/|b|; otherwise 0
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: det, mu
        det = a(1) * b(2) - a(2) * b(1)
        if(abs(det) < eps) then
            mu = norm2(a) / norm2(b)
        else
            mu = 0
        end if
        return
    end function m31colinear
    
    subroutine m33eigvect(a, eigval, eigvect)
        !! Computes the eigenvalues and eigenvectors for real symmetric 3x3
        !! matrix a. Eigenvectors are normalized to unit length.
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: eigval(3), eigvect(3, 3)
        real(dp) :: mu, cross1(3), cross2(3)
        integer :: i
        eigval = m33eigvalsh(a)
        ! Get the 1st and 2nd eigenvector, assume x1 /= x2
        ! Deal with x1=x2(=x3) and the end of first loop
        do i = 1, 2
            ! Construct the 1st and 2nd column to do cross product
            ! Deal with x1=x2(=x3), for i == 2 only
            if((abs(eigval(1) - eigval(2)) < eps) .and. (i == 2)) then
                cross1 = eigvect(:, 1)
            else
                cross1 = a(:, 1) - eigval(i) * delta(:, 1)
            end if
            cross2 = a(:, 2) - eigval(i) * delta(:, 2)
            ! If cross1 is zero, then [1, 0, 0] is eigvect, similar for cross2
            ! Doesn't matter if both cross1 and cross2 are zero, that means
            ! x(2) is a degenerate eigenvalue, but x(1) is not ... ?
            if(norm2(cross1) < eps) then
                eigvect(:, i) = [1, 0, 0]
            elseif(norm2(cross2) < eps) then
                eigvect(:, i) = [0, 1, 0]
            else
                mu = m31colinear(cross1, cross2)
                if(mu == 0._dp) then
                    eigvect(:, i) = m31cross(cross1, cross2)
                else
                    ! Will be normalized at the end
                    eigvect(:, i) = [1._dp, -mu, 0._dp]
                end if
            end if
        end do
        ! Get the 3rd eigenvector
        eigvect(:, 3) = m31cross(eigvect(:, 1), eigvect(:, 2))
        ! Normalize all
        do i = 1, 3
            eigvect(:, i) = eigvect(:, i) / norm2(eigvect(:, i))
        end do
    end subroutine m33eigvect
    
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


module psimod
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
end module psimod

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
    use psimod, only: getpsi
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
        ! Only calculate the symmetric part
        do k1 = 1, 3
            do k2 = 1, k1
                rcgptb = rcg + eps_s * (&
                    m31tensorprod(delta(:, k1), delta(:, k2)) +&
                    m31tensorprod(delta(:, k2), delta(:, k1)))
                psiptb = getpsi(rcgptb, props, statev)
                pk2(k1, k2) = (psiptb - psi) / eps_s
                ! Fill the other part from symmetry
                if(k1 /= k2) then
                    pk2(k2, k1) = pk2(k1, k2)
                end if
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
        ! Only calculate the symmetric part
        do k3 = 1, 3
        do k4 = 1, k3
            fptb = dfgrd + eps_c / 2 * (&
                matmul(m31tensorprod(delta(:, k3), delta(:, k4)), dfgrd) +&
                matmul(m31tensorprod(delta(:, k4), delta(:, k3)), dfgrd))
            tauptb = gettau(fptb, props, statev, eps_s)
            ccj(:, :, k3, k4) = (tauptb - tau) / eps_c / det
            ! Fill the other part from symmetry
            if(k3 /= k4) then
                ccj(:, :, k4, k3) = ccj(:, :, k3, k4)
            end if
        end do
        end do
        end subroutine get_sigma_ccj
end module numerichyper


subroutine umat(stress,statev,ddsdde,sse,spd,scd,&
    rpl,ddsddt,drplde,drpldt,&
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
    ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
    use numerichyper, only: update_umat
    use umatutils, only: dp
    implicit none
    ! This is a hack. The content of the 'aba_param.inc' is simply the
    ! Following two lines. I commented the first one, and modified the
    ! second one.
    ! include 'aba_param.inc'
    ! implicit real*8(a-h,o-z)
    ! parameter (nprecd=2)
    integer, parameter :: nprecd=2
    character*80, intent(in) :: cmname
    real(dp), intent(in) :: stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
        props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3), dtime, temp, &
        dtemp, celent
    integer, intent(in) :: ndi, nshr, ntens, nstatv, nprops, noel, npt, layer,&
        kspt, kstep, kinc
    real(dp), intent(inout) :: stress(ntens), statev(nstatv), sse, spd, scd,&
        rpl, ddsdde(ntens, ntens), ddsddt(ntens), drplde(ntens), drpldt, pnewdt
    ! Rotate SDVs from cylindrical coordinate to cartesian
    if((kstep==1).and.(kinc==1)) then
        statev(1:3) = cyl2cart(statev(1:3), coords)
        statev(4:6) = cyl2cart(statev(4:6), coords)
    end if
    ! Update the stress and ddsdde
    call update_umat(props, dfgrd1, ntens, stress, ddsdde, statev)
contains
    function cyl2cart(cyl, coords) result(cart)
        real(dp), intent(in) :: cyl(3), coords(3)
        real(dp) :: cart(3), r, c, a, x, y, z, theta
        r = cyl(1)
        c = cyl(2)
        a = cyl(3)
        theta = atan(coords(2)/coords(1))
        x = r * cos(theta) - c * sin(theta)
        y = r * sin(theta) + c * cos(theta)
        z = a
        cart = [x, y, z]
        return
    end function cyl2cart
end subroutine umat


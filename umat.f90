module umatutils
    !!! Utility module to be used for UMATs
    implicit none
    private
    public dp, delta, m31tensorprod, m33det, m33tensorprod, mapnotation,&
        rotm31
    integer, parameter :: dp=kind(0.d0)
    real(dp), parameter :: delta(3, 3)=reshape([1, 0, 0, 0, 1, 0, 0, 0, 1],&
        [3, 3])
        
contains
    subroutine rotm31(rot, m31)
        !! Rotate a vector by matrix rot
        real(dp), intent(in) :: rot(3, 3)
        real(dp), intent(inout) :: m31(3)
        m31 = matmul(rot, m31)
    end subroutine rotm31
    
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
    public update_stress_ddsdde

contains
    subroutine update_stress_ddsdde(props, dfgrd, ntens, stress, ddsdde, &
        statev)
        !! Update the stress and ddsdde for Neo-Hookean
        real(dp), intent(in) :: props(:), dfgrd(3, 3)
        integer, intent(in) :: ntens
        real(dp), intent(out) :: stress(ntens), ddsdde(ntens, ntens)
        real(dp), intent(inout) :: statev(:)
        real(dp), parameter :: eps_s=1d-4, eps_c=1d-6
        real(dp) :: sigma(3, 3), ccj(3, 3, 3, 3)
        call get_sigma_ccj(dfgrd, props, eps_c, eps_s, sigma, ccj, statev)
        call mapnotation(sigma, ccj, ntens, stress, ddsdde)
    end subroutine update_stress_ddsdde
        
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
        ! Use k3 & k4 rather than k1 & k2 to denote that it's the k, l 
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


subroutine umat(stress,statev,ddsdde,sse,spd,scd,&
    rpl,ddsddt,drplde,drpldt,&
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
    ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
    use numerichyper, only: update_stress_ddsdde
    use umatutils, only: dp, rotm31
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
    real(dp) :: c10, d1
    ! Update the stress and ddsdde
    call update_stress_ddsdde(props, dfgrd1, ntens, stress, ddsdde, statev)
    ! Update fiber orientation if it is two families of fibers
    if(nstatv==6) then
        call rotm31(drot, statev(1:3))
        call rotm31(drot, statev(4:6))
    end if
end subroutine umat


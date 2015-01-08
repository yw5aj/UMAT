module numeric_nh
    !!! Module to be embedded in the UMAT
    !!! Paramters
    !!! ---------
    !!! det : Jacobian
    !!! dfgrd : deformation gradient
    !!! fbar : modified deformation gradient
    !!! rcg : right Cauchy-Green tensor
    !!! rcgbar : modified right Cauchy-Green tensor
    !!! ibar1 : first invariant of rcgbar
    !!! tau : Kirchoff stress
    !!! ccj : tangent modulus for Cauchy stress in Jaumann rate
    implicit none
    integer, parameter :: dp=kind(1.d0)
    integer, parameter :: notationmap(6, 2)=reshape((/1, 2, 3, 1, 1, 2, 1, 2,&
        3, 2, 3, 3/), (/6, 2/))
    real(dp), parameter :: delta(3, 3)=reshape((/1, 0, 0, 0, 1, 0, 0, 0, 1/),&
        (/3, 3/))

contains
    !!! Utility functions
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
    
    !!! Neo-Hookean specific functions
    real(dp) function getpsi(rcg, c10, d1)
        !! Return strain energy density given C and material properties
        real(dp), intent(in) :: rcg(3, 3), c10, d1
        real(dp) :: rcgbar(3, 3), det, ibar1
        integer :: k1, k2
        det = (m33det(rcg))**(1.d0/2)
        rcgbar = det**(-2.d0/3) * rcg
        ibar1 = sum((/(rcgbar(k1, k1), k1=1, 3)/))
        getpsi = c10 * (ibar1 - 3.d0) + 1.d0 / d1 * (det - 1.d0)**2.d0
        return
    end function getpsi
        
    function gettau(dfgrd, c10, d1, epss)
        !! Return the Kirchoff stress given deformation gradient
        real(dp), intent(in) :: dfgrd(3, 3), c10, d1, epss
        real(dp) :: gettau(3, 3)
        real(dp) :: det, rcg(3, 3), psi, rcgptb(3, 3), psiptb, pk2(3, 3)
        integer :: k1, k2, k3, k4
        det = m33det(dfgrd)
        rcg = matmul(transpose(dfgrd), dfgrd)
        psi = getpsi(rcg, c10, d1)
        ! Outer two layers are indices for the 2nd PK stress
        do k1 = 1, 3
            do k2 = 1, 3
                rcgptb = rcg + epss * (&
                    m31tensorprod(delta(:, k1), delta(:, k2)) +&
                    m31tensorprod(delta(:, k2), delta(:, k1)))
                psiptb = getpsi(rcgptb, c10, d1)
                pk2(k1, k2) = (psiptb - psi) / epss
            end do
        end do
        ! Calculate Kirchoff stress from 2nd PK stress
        gettau = matmul(matmul(dfgrd, pk2), transpose(dfgrd))
        return
        end function gettau

    function getccj(dfgrd, c10, d1, epsc, epss)
        real(dp), intent(in) :: dfgrd(3, 3), c10, d1, epsc, epss
        real(dp) :: getccj(3, 3, 3, 3)
        real(dp) :: tau(3, 3), tauptb(3, 3), fptb(3, 3), det
        integer :: k3, k4
        det = m33det(dfgrd)
        tau = gettau(dfgrd, c10, d1, epss)
        ! Use k3 & k4 rather than k1 & k2 to denote that it's the k, l 
        ! component of the elasticity tensor being calculated
        do k3 = 1, 3
        do k4 = 1, 3
            fptb = dfgrd + epsc / 2.d0 * (&
                matmul(m31tensorprod(delta(:, k3), delta(:, k4)), dfgrd) +&
                matmul(m31tensorprod(delta(:, k4), delta(:, k3)), dfgrd))
            tauptb = gettau(fptb, c10, d1, epss)
            getccj(:, :, k3, k4) = 1.d0 / det / epsc * (tauptb - tau)
        end do
        end do
        return
        end function getccj
end module numeric_nh


subroutine umat(stress,statev,ddsdde,sse,spd,scd,&
    rpl,ddsddt,drplde,drpldt,&
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
    ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
    use numeric_nh, only: gettau, getccj, m33det, notationmap
    implicit none
    ! This is a hack. The content of the 'aba_param.inc' is simply the
    ! Following two lines. I commented the first one, and modified the
    ! second one.
    ! include 'aba_param.inc'
    ! implicit real*8(a-h,o-z)
    ! parameter (nprecd=2)
    integer, parameter :: nprecd=2
    character*80 :: cmname
    integer, parameter :: dp=kind(1.d0)
    real(dp), intent(in) :: stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
        props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3), dtime, temp, &
        dtemp, celent
    integer, intent(in) :: ndi, nshr, ntens, nstatv, nprops, noel, npt, layer,&
        kspt, kstep, kinc
    real(dp), intent(inout) :: stress(ntens), statev(nstatv), sse, spd, scd,&
        rpl, ddsdde(ntens, ntens), ddsddt(ntens), drplde(ntens), drpldt, pnewdt
    real(dp) :: sigma(3, 3), ccj(3, 3, 3, 3), det, epss, epsc, c10, d1
    integer :: k1, k2
    ! Assign the epsilons
    epss = 1d-4
    epsc = 1d-6
    ! Get material properties
    c10 = props(1)
    d1 = props(2)
    ! Get the standard matrix form of Cauchy stress and tangent modulus
    det = m33det(dfgrd1)
    sigma = gettau(dfgrd1, c10, d1, epss) / det
    ccj = getccj(dfgrd1, c10, d1, epsc, epss)
    ! Fill in stress and ddsdde
    do k1 = 1, ntens
        stress(k1) = sigma(notationmap(k1, 1), notationmap(k1, 2))    
        do k2 = 1, ntens
            ddsdde(k1, k2) = ccj(notationmap(k1, 1), notationmap(k1, 2),&
                notationmap(k2, 1), notationmap(k2, 2))
        end do
    end do
    return
    end subroutine umat
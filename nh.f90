module numeric_nh
    !!! Module to be embedded in the UMAT
    implicit none
    integer, parameter :: dp=kind(1.d0)
    real(dp) :: delta(3, 3)
    data delta /1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/

contains
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
                ! Inner two layers are for perturbing the rcg
                do k3 = 1, 3
                    do k4 = 1, 3
                        rcgptb(k3, k4) = rcg(k3, k4)&
                            + epss * delta(k1, k3) * delta(k2, k4)&
                            + epss * delta(k2, k3) * delta(k1, k4)
                    end do
                end do
                psiptb = getpsi(rcgptb, c10, d1)
                pk2(k1, k2) = (psiptb - psi) / epss
            end do
        end do
        ! Calculate Kirchoff stress from 2nd PK stress
        gettau = matmul(matmul(dfgrd, pk2), transpose(dfgrd))
        return
        end

    function getccj(dfgrd, c10, d1, epsc, epss)
        real(dp), intent(in) :: dfgrd(3, 3), c10, d1, epsc, epss
        real(dp) :: getccj(3, 3, 3, 3)
        real(dp) :: tau(3, 3), tauptb(3, 3), fptb(3, 3), det
        integer :: k1, k2, k3, k4
        det = m33det(dfgrd)
        tau = gettau(dfgrd, c10, d1, epss)
        do k3 = 1, 3
        do k4 = 1, 3
            do k1 = 1, 3
            do k2 = 1, 3
                fptb(k1, k2) = dfgrd(k1, k2) + epsc / 2.d0 * (&
                    (delta(k1,k3)*delta(1,k4)+delta(1,k3)*delta(k1,k4))*dfgrd(1,k2)&
                    +(delta(k1,k3)*delta(2,k4)+delta(2,k3)*delta(k1,k4))*dfgrd(2,k2)&
                    +(delta(k1,k3)*delta(3,k4)+delta(3,k3)*delta(k1,k4))*dfgrd(3,k2))
            end do
            end do
            tauptb = gettau(fptb, c10, d1, epss)
            do k1 = 1, 3
            do k2 = 1, 3
                getccj(k1, k2, k3, k4) = 1.d0 / det / epsc * &
                    (tauptb(k1, k2) - tau(k1, k2))
            end do
            end do
        end do
        end do
        return
        end function getccj
end module numeric_nh


program nh
    !!! Testing numerical method for NH model
    use numeric_nh, only: gettau, getccj
    implicit none
    ! det : Jacobian
    ! dfgrd : deformation gradient
    ! fbar : modified deformation gradient
    ! rcg : right Cauchy-Green tensor
    ! rcgbar : modified right Cauchy-Green tensor
    ! ibar1 : first invariant of rcgbar
    ! tau : Kirchoff stress
    ! ccj : tangent modulus for Cauchy stress in Jaumann rate
    integer, parameter :: dp = kind(0.d0)
    real(dp) :: ibar1, epss, epsc, getdet, d1, c10
    real(dp), dimension(3, 3) :: dfgrd, fbar, rcg, rcgbar, tau
    real(dp), dimension(3, 3, 3, 3) :: ccj 
    integer :: k1, k2, k3, k4, fid, delta
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
    ! dfgrd = reshape((/2.d0, 0.d0, 0.d0, 0.d0, 2.d0, 0.d0, 0.d0, 0.d0, 2.d0/), (/3, 3/))
    c10 = 8d4
    d1 = 2d-1
    epss = 1d-4
    epsc = 1d-6
    ! Calculate
    tau = gettau(dfgrd, c10, d1, epss)
    ccj = getccj(dfgrd, c10, d1, epsc, epss)
    ! Write to a file
    ! fid = 1
    ! open (fid, file='ccj.txt', status='new')
    ! write (fid, *) ccj
    ! fid = 2
    ! open (fid, file='tau.txt', status='new')
    ! write (fid, *) tau
    write (*, *) ccj
    write (*, *) tau

contains
    !!! Does not contain any subroutine / functions
end program nh






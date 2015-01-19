program nhtest
    !!! Testing numerical method (Neo-Hookean model used herein)
    use numerichyper, only: update_umat
    use umatutils, only: dp
    implicit none
    real(dp) :: d1, c10, props(2)
    integer, parameter :: ntens=6
    real(dp) :: stress(ntens), ddsdde(ntens, ntens), dfgrd(3, 3), statev(6)
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
    ! dfgrd = reshape([1._dp, 0._dp, 0._dp, 0._dp, 1._dp, 0._dp, 0._dp, 0._dp, 1._dp], [3, 3])
    ! dfgrd = dfgrd * 1.6d0
    c10 = 8e4_dp
    d1 = 2e-1_dp
    props = [c10, d1]
    ! Calculate
    call update_umat(props, dfgrd, ntens, stress, ddsdde, statev)
    write(*, *) 'Stress: ', stress
    write(*, *) 'DDSDDE: ', ddsdde
end program nhtest
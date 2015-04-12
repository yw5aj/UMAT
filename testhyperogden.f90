program testogden
    use umatutils, only: dp
    use hypermod, only: hyper, hyperpk2
    implicit none
    real(dp) :: dfgrd1(3, 3), props(3), sigma(3, 3), ccj(3, 3, 3, 3),&
        siso(3, 3), svol(3, 3)
    dfgrd1 = reshape([1., 0., 0., 0., 1., 0., .45, 0., 1.], [3, 3])
    props = [160e3_dp, 2._dp, .2_dp]
    call hyper(dfgrd1, props, size(props)/3, sigma, ccj)
    call hyperpk2(dfgrd1, props, size(props)/3, sigma, ccj, siso, svol)
    write (*, *) siso, svol
end program testogden
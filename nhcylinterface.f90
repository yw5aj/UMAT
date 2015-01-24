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
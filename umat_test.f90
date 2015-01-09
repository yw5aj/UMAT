program umat_test
    use umatutils, only: dp
    implicit none
    character*80 :: cmname='test'
    integer, parameter :: nprecd=2
    real(dp) :: stran(6),dstran(6),time(2),predef(1),dpred(1),&
        props(2),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3), dtime, temp, &
        dtemp, celent
    real(dp) :: stress(6), statev(1), sse, spd, scd,&
        rpl, ddsdde(6, 6), ddsddt(6), drplde(6), drpldt, pnewdt
    dfgrd1 = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
    props = [8000., 0.2]
    call umat(stress,statev,ddsdde,sse,spd,scd,&
        rpl,ddsddt,drplde,drpldt,&
        stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
        3,3,6,1,props,2,coords,drot,pnewdt,&
        celent,dfgrd0,dfgrd1,1,1,1,1,1,1)
    write(*, *) stress
end program umat_test

      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c
      include 'aba_param.inc'
c
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
c    local arrays
c ----------------------------------------------------------------
c det : Jacobian
c dfgrd1 : deformation gradient
c fbar : modified deformation gradient
c rcg : right Cauchy-Green tensor
c rcgbar : modified right Cauchy-Green tensor
c ibar1 : first invariant of rcgbar
c tau : Kirchoff stress
c sigma : Cauchy stress
c ccj : tangent modulus for Cauchy stress in Jaumann rate
c ----------------------------------------------------------------
c
      double precision det, fbar(3, 3), rcg(3, 3), 
     1    rcgbar(3, 3), ibar1, tau(3, 3), ccj(3, 3, 3, 3),
     2    epss, epsc, getdet, d1, c10, sigma(3, 3)
      integer k1, k2, fid, delta
c
c ----------------------------------------------------------------
c umat for compressible neo-hookean hyperelasticity
c cannot be used for plane stress
c ----------------------------------------------------------------
c props(1) : c10
c props(2) - d1
c ----------------------------------------------------------------
      c10=props(1)
      d1=props(2)
c First calculate stress      
      call gettau(dfgrd1, c10, d1, epss, tau)
      det = getdet(dfgrd1)
      do k1 = 1, 3
      do k2 = 1, 3
          sigma(k1, k2) = tau(k1, k2) / det
      end do
      end do
c Then calculate ccj
      call getccj(dfgrd1, c10, d1, epsc, epss, ccj)      
c Convert to abaqus format
c First for stress
      stress(1) = sigma(1, 1)
      stress(2) = sigma(2, 2)
      stress(3) = sigma(3, 3)
      stress(4) = sigma(1, 2)
      if (nshr.eq.3) then
          stress(5) = sigma(1, 3)
          stress(6) = sigma(2, 3)
      end if
c Then for tangent modulus
      ddsdde(1, 1) = ccj(1, 1, 1, 1)
      ddsdde(2, 2) = ccj(2, 2, 2, 2)
      ddsdde(3, 3) = ccj(3, 3, 3, 3)
      ddsdde(1, 2) = ccj(1, 1, 2, 2)
      ddsdde(1, 3) = ccj(1, 1, 3, 3)
      ddsdde(2, 3) = ccj(2, 2, 3, 3)
      ddsdde(1, 4) = ccj(1, 1, 1, 2)
      ddsdde(2, 4) = ccj(2, 2, 1, 2)
      ddsdde(3, 4) = ccj(3, 3, 1, 2)
      ddsdde(4, 4) = ccj(1, 2, 1, 2)
      if (nshr.eq.3) then
          ddsdde(1, 5) = ccj(1, 1, 1, 3)
          ddsdde(2, 5) = ccj(2, 2, 1, 3)
          ddsdde(3, 5) = ccj(3, 3, 1, 3)
          ddsdde(1, 6) = ccj(1, 1, 2, 3)
          ddsdde(2, 6) = ccj(2, 2, 2, 3)
          ddsdde(3, 6) = ccj(3, 3, 2, 3)
          ddsdde(5, 5) = ccj(1, 3, 1, 3)
          ddsdde(6, 6) = ccj(2, 3, 2, 3)
          ddsdde(4, 5) = ccj(1, 2, 1, 3)
          ddsdde(4, 6) = ccj(1, 2, 2, 3)
          ddsdde(5, 6) = ccj(1, 3, 2, 3)
      end if
      do k1 = 1, ntens
          do k2 = 1, k1 - 1
              ddsdde(k1, k2) = ddsdde(k2, k1)
          end do
      end do
      return
      end      
c
c
c The subroutine to calculate strain energy given C
      subroutine getpsi (rcg, c10, d1, psi)
      double precision rcg(3, 3), c10, d1, rcgbar(3, 3), psi, ibar1,
     c    getdet, det
      integer k1, k2
      det = (getdet(rcg))**(1./2.)
      do k1 = 1, 3
      do k2 = 1, 3
          rcgbar(k1, k2) = det**(-2./3.) * rcg(k1, k2)
      end do
      end do
      ibar1 = rcgbar(1, 1) + rcgbar(2, 2) + rcgbar(3, 3)
      psi = c10 * (ibar1 - 3.) + 1. / d1 * (det - 1.)**2
      return
      end
c
c The subroutine to calculate the Kirchoff stress given F
      subroutine gettau (dfgrd1, c10, d1, epss, tau)
      double precision dfgrd1(3, 3), c10, d1, epss, rcg(3, 3), 
     1    pk2(3, 3), rcgptb(3, 3), tau(3, 3), psi, psiptb, det,
     2    getdet
      integer k1, k2, k3, k4, delta
      det = getdet(dfgrd1)
      do k1 = 1, 3
      do k2 = 1, 3
          rcg(k1, k2) = dfgrd1(1, k1) * dfgrd1(1, k2)
     1        + dfgrd1(2, k1) * dfgrd1(2, k2)
     2        + dfgrd1(3, k1) * dfgrd1(3, k2)
      end do
      end do
      call getpsi(rcg, c10, d1, psi)
      do k1 = 1, 3
      do k2 = 1, 3
          do k3 = 1, 3
          do k4 = 1, 3
              rcgptb(k3, k4) = rcg(k3, k4) 
     1            + epss * delta(k1, k3) * delta(k2, k4)
     2            + epss * delta(k2, k3) * delta(k1, k4)
          end do
          end do
          call getpsi(rcgptb, c10, d1, psiptb)
          pk2(k1, k2) = (psiptb - psi) / epss
      end do
      end do
      do k1 = 1, 3
      do k2 = 1, 3
          tau(k1, k2) = dfgrd1(k1, 1) * pk2(1, 1) * dfgrd1(k2, 1)
     1        + dfgrd1(k1, 2) * pk2(2, 1) * dfgrd1(k2, 1)
     2        + dfgrd1(k1, 3) * pk2(3, 1) * dfgrd1(k2, 1)
     3        + dfgrd1(k1, 1) * pk2(1, 2) * dfgrd1(k2, 2)
     4        + dfgrd1(k1, 2) * pk2(2, 2) * dfgrd1(k2, 2)
     5        + dfgrd1(k1, 3) * pk2(3, 2) * dfgrd1(k2, 2)
     6        + dfgrd1(k1, 1) * pk2(1, 3) * dfgrd1(k2, 3)
     7        + dfgrd1(k1, 2) * pk2(2, 3) * dfgrd1(k2, 3)
     8        + dfgrd1(k1, 3) * pk2(3, 3) * dfgrd1(k2, 3)
      end do
      end do
      return
      end
c
c The subroutine to calculate ccj
      subroutine getccj (dfgrd1, c10, d1, epsc, epss, ccj)
      double precision dfgrd1(3, 3), c10, d1, epsc, epss, 
     1    getdet, det, tau(3, 3), tauptb(3, 3), fptb(3, 3),
     2    ccj(3, 3, 3, 3)
      integer k1, k2, k3, k4, k5, k6, delta
      det = getdet(dfgrd1)
      call gettau(dfgrd1, c10, d1, epss, tau)
      do k3 = 1, 3
      do k4 = 1, 3
          do k1 = 1, 3
          do k2 = 1, 3
              fptb(k1, k2) = dfgrd1(k1, k2) + epsc / 2. * (
     1  (delta(k1,k3)*delta(1,k4)+delta(1,k3)*delta(k1,k4))*dfgrd1(1,k2)
     2 +(delta(k1,k3)*delta(2,k4)+delta(2,k3)*delta(k1,k4))*dfgrd1(2,k2)
     3 +(delta(k1,k3)*delta(3,k4)+delta(3,k3)*delta(k1,k4))*dfgrd1(3,k2)
     4 )
          end do
          end do
          call gettau(fptb, c10, d1, epss, tauptb)
          do k1 = 1, 3
          do k2 = 1, 3
              ccj(k1, k2, k3, k4) = 1. / det / epsc * 
     1        (tauptb(k1, k2) - tau(k1, k2))
          end do
          end do
      end do
      end do
      return
      end
c   
c
c A function to calculate determinant
      double precision function getdet(dfgrd1)
      double precision dfgrd1(3, 3)
      getdet = dfgrd1(1, 1) * dfgrd1(2, 2) * dfgrd1(3, 3)
     1    + dfgrd1(1, 2) * dfgrd1(2, 3) * dfgrd1(3, 1)
     2    + dfgrd1(1, 3) * dfgrd1(2, 1) * dfgrd1(3, 2)
     3    - dfgrd1(1, 1) * dfgrd1(2, 3) * dfgrd1(3, 2)
     4    - dfgrd1(1, 2) * dfgrd1(2, 1) * dfgrd1(3, 3)
     5    - dfgrd1(1, 3) * dfgrd1(2, 2) * dfgrd1(3, 1)      
      return
      end
c
c A function to calculate delta
      integer function delta(k1, k2)
      integer k1, k2
      if (k1.eq.k2) then
          delta = 1
      else
          delta = 0
      end if
      return
      end
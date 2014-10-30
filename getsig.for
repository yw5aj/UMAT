      subroutine getsig (c10, d1, f, sigma)
      real c10, d1
      real f(3, 3), sigma(3, 3)
c 
c local variables
c 
      real j
      integer k1, k2
      real fbar(3, 3), bbar(3, 3), sigma(3, 3)      
c 
c jacobian
c
      j = f(1, 1) * f(2, 2) * f(3, 3)
     1    + f(1, 2) * f(2, 3) * f(3, 1)
     2    + f(1, 3) * f(3, 2) * f(2, 1)
     3    - f(1, 2) * f(2, 1) * f(3, 3)
     4    - f(1, 3) * f(3, 1) * f(2, 2)
     5    - f(2, 3) * f(3, 2) * f(1, 1)
      end if
c 
c modified f and b
c 
      do k1 = 1, 3
          do k2 = 1, 3
              fbar(k1, k2) = j**(-1./3.) * f(k1, k2)
          end do
      end do
      bbar(1, 1) = fbar(1, 1)**2 + fbar(1, 2)**2 + fbar(1, 3)**2
      bbar(2, 2) = fbar(2, 1)**2 + fbar(2, 2)**2 + fbar(2, 3)**2
      bbar(3, 3) = fbar(3, 1)**2 + fbar(3, 2)**2 + fbar(3, 3)**2
      bbar(1, 2) = fbar(1, 1) * fbar(2, 1) + fbar(1, 2) * fbar(2, 2)
     1    + fbar(1, 3) * fbar(2, 3)
      bbar(1, 3) = fbar(1, 1) * fbar(3, 1) + fbar(1, 2) * fbar(3, 2)
     1    + fbar(1, 3) * fbar(3, 3)
      bbar(2, 3) = fbar(2, 1) * fbar(3, 1) + fbar(2, 2) * fbar(3, 2)
     1    + fbar(2, 3) * fbar(3, 3)     
      bbar(2, 1) = bbar(1, 2)
      bbar(3, 1) = bbar(1, 3)
      bbar(3, 2) = bbar(2, 3)
c 
c calculate the cauchy stress tensor
c 
      sigma(1, 1) = 2. / j * c10 * (bbar(1, 1) - 1./3. * (bbar(1, 1) 
     1    + bbar(2, 2) + bbar(3, 3))) + 2. / d1 * (j - 1.)
      sigma(2, 2) = 2. / j * c10 * (bbar(2, 2) - 1./3. * (bbar(1, 1) 
     1    + bbar(2, 2) + bbar(3, 3))) + 2. / d1 * (j - 1.)
      sigma(3, 3) = 2. / j * c10 * (bbar(3, 3) - 1./3. * (bbar(1, 1) 
     1    + bbar(2, 2) + bbar(3, 3))) + 2. / d1 * (j - 1.)
      sigma(1, 2) = 2. / j * c10 * bbar(1, 2)
      sigma(1, 3) = 2. / j * c10 * bbar(1, 3)
      sigma(2, 3) = 2. / j * c10 * bbar(2, 3)
      sigma(2, 1) = sigma(1, 2)
      sigma(3, 1) = sigma(1, 3)
      sigma(3, 2) = sigma(2, 3)
      return
      end      
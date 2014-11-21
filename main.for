      program main
c The main program used for testing      
      real det delta pertf
      real f(3, 3), j, fbar(3, 3), bbar(3, 3), sigma(3, 3), c10, d1
      integer k1, k2, k3, k4
c      
      do 10 k1 = 1, 3
      do 20 k2 = 1, 3
          f(k1, k2) = 0.5 * delta(k1, k2)
20    continue
10    continue
      call leftcg(f, j, fbar, bbar)
      c10 = 1e4
      d1 = .1
      call nh(c10, d1, j, bbar, sigma)
c Output result      
      write (*, *) sigma
      stop
      end
      
      
      real function det(m)
      real m(3, 3)
c Determinant of 3-by-3 matrix m      
      det = m(1, 1) * m(2, 2) * m(3, 3)
     1    + m(1, 2) * m(2, 3) * m(3, 1)
     2    + m(1, 3) * m(3, 2) * m(2, 1)
     3    - m(1, 2) * m(2, 1) * m(3, 3)
     4    - m(1, 3) * m(3, 1) * m(2, 2)
     5    - m(2, 3) * m(3, 2) * m(1, 1)
      return
      end
      
      
      real function delta(i, j)
      integer i, j
c Return delta(i, j) as the Kronecker delta      
      if(i .eq. j) then 
          delta = 1
      else
          delta = 0
      end if
      return
      end


      real function pertf(f, i, j)
      real f(3, 3)
      integer i, j
c Return the perturbation of F on its (i, j) component
      do 
      return
      end


      subroutine leftcg (f, j, fbar, bbar)
      real f(3, 3), j, fbar(3, 3), bbar(3, 3)
c Jacobian, the modified deformation gradient and left Cauchy Green tensor
      integer k1, k2
      real det
c
      j = det(f)
      do 30 k1 = 1, 3
      do 40 k2 = 1, 3
          fbar(k1, k2) = j**(-1./3.) * f(k1, k2)
40    continue
30    continue
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
      return
      end
      
      
      subroutine nh (c10, d1, j, bbar, sigma)
      real c10, d1, j, bbar(3, 3), sigma(3, 3)
c Cauchy stress
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
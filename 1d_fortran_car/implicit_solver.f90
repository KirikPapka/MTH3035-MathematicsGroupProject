      SUBROUTINE implicit_solver(chi,k_diff,z,zn,dt,b_1,c_1,r_1,r_kkp,chi_inter)

      ! Subroutine to calculated turbulent diffusion 
      ! Using an implicit method

      IMPLICIT none
      INCLUDE 'params.inc'

      ! In The time level n value of the variable
      REAL, INTENT(IN), DIMENSION(kkp) :: chi
      ! Out the intermediate n+1 value after turbulent diffusion
      REAL, INTENT(OUT), DIMENSION(kkp) :: chi_inter

      ! The turbulent diffusion
      REAL, INTENT(IN), DIMENSION(kkp) :: k_diff 
      ! The vertical coordinates
      REAL, INTENT(IN), DIMENSION(kkp) :: z, zn
      
      ! Coefficients which change for u,v or theta
      REAL, INTENT(IN) :: b_1, c_1, r_1, r_kkp
     
      ! Locals 
      ! Coefficients for tri-diagonal solver
      REAL, DIMENSION(kkp) :: a, b, c, r
      INTEGER :: k

      REAL, DIMENSION(kkp) :: dz, dzn

      ! Define grid increments
      DO k=2, kkp
        dz(k) =z(k) - z(k-1)
        dzn(k)=zn(k) - zn(k-1)
      ENDDO

      ! Initialise a, b, c, r at k=1 and kkp
      a(1)=0.0
      b(1)=b_1
      c(1)=c_1
      r(1)=r_1

      a(kkp)=0.0
      b(kkp)=1
      c(kkp)=0.0
      r(kkp)=r_kkp       

      DO k=2, kkp-1
        a(k)=-dt*gamma*k_diff(k-1)/(dz(k)*dzn(k))      
        b(k)=1.0 + dt*gamma*k_diff(k)/(dz(k)*dzn(k+1))  &
             + dt*gamma*k_diff(k-1)/(dz(k)*dzn(k))
        c(k)=-dt*gamma*k_diff(k)/(dz(k)*dzn(k+1))
        r(k)= dt*(1.0-gamma)*k_diff(k)*chi(k+1)/(dz(k)*dzn(k+1))         &
              + (1 -  dt*(1.0-gamma)*k_diff(k)/(dz(k)*dzn(k+1))          &
              -  dt*(1.0-gamma)*k_diff(k-1)/(dz(k)*dzn(k)) )*chi(k)      &
              + dt*(1.0-gamma)*k_diff(k-1)*chi(k-1)/(dz(k)*dzn(k))
      ENDDO  

      DO k=1,kkp
        chi_inter(k)=chi(k)
      ENDDO

      ! Call tridiagonal solver
      CALL tridag(a,b,c,r,chi_inter,kkp)

      
      END

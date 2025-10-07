      SUBROUTINE implicit_solver_flux(chi,u,v,k_diff,z,zn,dt,r_kkp  &
                                      ,chi_inter,h_0)

      ! Subroutine to calculated turbulent diffusion 
      ! using an implicit method
      ! for applied surface heat flux.

      IMPLICIT none
      INCLUDE 'params.inc'

      ! In The time level n value of the variable
      REAL, INTENT(IN), DIMENSION(kkp) :: chi

      ! u and v
      REAL, INTENT(IN), DIMENSION(kkp) :: u, v

       

      ! Array for implicit solution
      ! chi_inter(2:kkp) = chi_temp(1:kkp-1)
      ! chi_inter(1) is determined by extrapolation
      ! using a neutral exchange coefficient
     
      REAL, DIMENSION(kkp-1) :: chi_temp
      ! Out the intermediate n+1 value after turbulent diffusion
      REAL, INTENT(OUT), DIMENSION(kkp) :: chi_inter

      ! The turbulent diffusion
      REAL, INTENT(IN), DIMENSION(kkp) :: k_diff 
      ! The vertical coordinates
      REAL, INTENT(IN), DIMENSION(kkp) :: z, zn

      ! Surface heat flux
      REAL, INTENT(IN) :: h_0 
      ! Surface heat flux at times n and n+1
      REAL :: H_n, H_n_pl_one

      ! Surface neutral exhange coefficient
      REAL :: kh_1_neut

      ! wind speed at level 2
      REAL :: wind_2
      
      ! Coefficients which change for u,v or theta
      REAL :: b_1, c_1, r_1, r_kkp
     
      ! Locals 
      ! Coefficients for tri-diagonal solver
      REAL, DIMENSION(kkp-1) :: a, b, c, r
      INTEGER :: k

      REAL, DIMENSION(kkp) :: dz, dzn
      
      ! h_0 initially fixed in time 
      H_n=h_0
      H_n_pl_one=h_0

      ! Define grid increments
      DO k=2, kkp
        dz(k) =z(k) - z(k-1)
        dzn(k)=zn(k) - zn(k-1)
      ENDDO

      ! Initialise a, b, c, r at k=1 and kkp

      a(kkp-1)=0.0
      b(kkp-1)=1
      c(kkp-1)=0.0
      r(kkp-1)=r_kkp       

      
      a(1)=0
      b(1)=1.0 + dt*gamma*k_diff(2)/(dz(2)*dzn(3))
      c(1)=-dt*gamma*k_diff(2)/(dz(2)*dzn(3))

! Implicit surface fluxes
!      r(1)= dt*(1.0-gamma)*k_diff(2)*chi(3)/(dz(2)*dzn(3))               &
!              + (1 -  dt*(1.0-gamma)*k_diff(2)/(dz(2)*dzn(3)) )*chi(2)   &
!              +  dt*(gamma*H_n_pl_one+(1-gamma)*H_n)/(dz(2)*rho_0*c_p)

! Explicit surface fluxes gamma=0 in h_n specification
      r(1)= dt*(1.0-gamma)*k_diff(2)*chi(3)/(dz(2)*dzn(3))               &
              + (1 -  dt*(1.0-gamma)*k_diff(2)/(dz(2)*dzn(3)) )*chi(2)   &
              +  dt*H_n/(dz(2)*rho_0*c_p)

      DO k=3, kkp-1
        a(k-1)=-dt*gamma*k_diff(k-1)/(dz(k)*dzn(k))      
        b(k-1)=1.0 + dt*gamma*k_diff(k)/(dz(k)*dzn(k+1))  &
             + dt*gamma*k_diff(k-1)/(dz(k)*dzn(k))
        c(k-1)=-dt*gamma*k_diff(k)/(dz(k)*dzn(k+1))
        r(k-1)= dt*(1.0-gamma)*k_diff(k)*chi(k+1)/(dz(k)*dzn(k+1))         &
              + (1 -  dt*(1.0-gamma)*k_diff(k)/(dz(k)*dzn(k+1))          &
              -  dt*(1.0-gamma)*k_diff(k-1)/(dz(k)*dzn(k)) )*chi(k)      &
              + dt*(1.0-gamma)*k_diff(k-1)*chi(k-1)/(dz(k)*dzn(k))
      ENDDO  

      DO k=1,kkp-1
        chi_temp(k)=chi(k+1)
      ENDDO

      ! Call tridiagonal solver
      CALL tridag(a,b,c,r,chi_temp,kkp-1)
      
      DO k=2,kkp
        chi_inter(k)=chi_temp(k-1)
      ENDDO

      ! Calculate surface neutral exchange coefficient
      wind_2=sqrt(u(2)**2 + v(2)**2)
      kh_1_neut= wind_2 * (kar**2)*z_0/alog( (zn(2)+z_0)/z_0 ) 
      
      ! Surface extrapolation

      chi_inter(1) = (H_n_pl_one/(rho_0*c_p)) * dzn(2)/kh_1_neut   & 
                     + chi_inter(2)     

      END

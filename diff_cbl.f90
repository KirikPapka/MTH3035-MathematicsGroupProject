      SUBROUTINE diff_cbl(u_star,w_star,h_cbl,z,zn,L_obuk,k_m,k_h)
      ! Calculate diffusion for the CBL
      IMPLICIT none
      INCLUDE 'params.inc'
      ! Inputs
      REAL, INTENT(IN) :: u_star, w_star, h_cbl, L_obuk
      REAL, INTENT(IN), DIMENSION(kkp) :: z, zn
      INTEGER :: k
      ! Outputs
      REAL, INTENT(INOUT), DIMENSION(kkp) :: k_m, k_h
      ! Locals      
      ! Velocity scales for momentum and heat
      REAL :: w_m, w_t
     
      ! Demensionless velocity and heat gradients
      REAL :: phi_m, phi_h, phi_h_on_phi_m

      ! Prandtl no.
      REAL :: pr
      
      ! Calculate Prandtl no.
      ! Independent of height 
      ! Calculated at 0.1h_cbl
      phi_h=(1 - 15.0*0.1*h_cbl/L_obuk)**(-0.5)
      phi_m=(1 - 15.0*0.1*h_cbl/L_obuk)**(-1.0/3.0)
      phi_h_on_phi_m=phi_h/phi_m
      w_m= (u_star**3 + c_1_hb*w_star**3)**(1.0/3.0)
      pr=   phi_h_on_phi_m              &
            +  a_hb*kar* 0.1 * w_star/w_m
      DO k=2, kkp-1
        
        phi_h=(1 - 15.0*zn(k+1)/L_obuk)**(-0.5)
        phi_m=(1 - 15.0*z(k)/L_obuk)**(-1.0/3.0)

        ! Calculate velocity scales

        IF (z(k)/h_cbl <= 0.1) THEN
          ! In surface layer
          w_m = u_star/phi_m
        ELSE
          ! Above surface layer
          w_m= (u_star**3 + c_1_hb*w_star**3)**(1.0/3.0)
        ENDIF

        w_t=w_m/pr 
        IF (z(k) <= h_cbl) THEN
           ! Diffusion coefficients
          k_m(k)=w_m*kar*z(k)*(1-z(k)/h_cbl)**2
        ENDIF 
        IF (zn(k+1) <= h_cbl) THEN
           ! Diffusion coefficients
          k_h(k)=w_t*kar*zn(k+1)*(1-zn(k+1)/h_cbl)**2 
        ENDIF 
      ENDDO

      END
  

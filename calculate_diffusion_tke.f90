      SUBROUTINE calculate_diffusion_tke(u,v,th,z,zn,k_m,k_h,ri)
      ! A subroutine to calculate the diffusion using
      ! TKE closure
      IMPLICIT none
      INCLUDE 'params.inc'
      ! Inputs
      REAL, INTENT(IN), DIMENSION(kkp) :: u, v, th, z, zn
      INTEGER :: k
      ! Outputs
      REAL, INTENT(INOUT), DIMENSION(kkp) :: k_m, k_h,ri
      ! Locals      
      REAL, DIMENSION(kkp) :: shear,f_m,f_h      
      REAL, PARAMETER :: ri_c=0.25
      REAL :: lambda 
      ! Calculate gradient richardson number
      ! shear 
      ! and stability functions


      DO k=1, kkp-1
        ri(k)= (g/th_ref)*                                     &         
               (th(k+1)-th(k))*(zn(k+1)-zn(k))/                &
               ( (u(k+1)-u(k))**2 + (v(k+1)-v(k))**2 +1e-9)

        shear(k)=( (u(k+1)-u(k))**2 + (v(k+1)-v(k))**2 )/      &
                 ((zn(k+1)-zn(k))**2)

        IF ((ri(k) <= ri_c) .and. (ri(k) > 0.0)) THEN  
          f_m(k)=1.0-ri(k)/ri_c
          f_h(k)=1.0-ri(k)/ri_c
           
           
        ELSEIF (ri(k) <= 0.0) THEN
          f_m(k)=1.0
          f_h(k)=1.0
        ELSE
          f_m(k)=0.0
          f_h(k)=0.0
        END IF
         
!        ! Long tails 
!        IF (ri(k) > 0.0) THEN
!          f_m(k)=1/(1+10.0*ri(k))
!          f_h(k)=1/(1+10.0*ri(k))
!        ELSE
!          f_m(k)=1.0
!          f_h(k)=1.0
!        ENDIF
        
      ENDDO

      ! top boundary condition
      ri(kkp)=ri(kkp-1)
      shear(kkp)=shear(kkp-1)
      f_m(kkp)=f_m(kkp-1)
      f_h(kkp)=f_h(kkp-1)
     
         
      DO k=1, kkp
        lambda=1.0/(1.0/(kar*(z(k)+0.1)) + 1.0/lambda_const)
!        lambda=lambda_const
        k_m(k)=(lambda**2)*shear(k)
!*f_m(k)
        k_h(k)=(lambda**2)*shear(k)
!*f_h(k)
      ENDDO

      ! Artifically over-ride level 1 k's
!      k_m(1)=0.01 
!      k_h(1)=k_m(1)

      END
  

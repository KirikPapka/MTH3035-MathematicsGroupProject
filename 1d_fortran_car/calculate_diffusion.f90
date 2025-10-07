      SUBROUTINE calculate_diffusion(u,v,th,z,zn,k_m,k_h,ri_m,ri_h)
      ! Calculate local diffusion     
      IMPLICIT none
      INCLUDE 'params.inc'
      ! Inputs
      REAL, INTENT(IN), DIMENSION(kkp) :: u, v, th, z, zn
      INTEGER :: k
      ! Outputs
      REAL, INTENT(INOUT), DIMENSION(kkp) :: k_m, k_h, ri_m, ri_h
      ! Locals      
      REAL, DIMENSION(kkp) :: shear_m, shear_h, f_m,f_h  
      REAL, DIMENSION(kkp) :: delta_s, delta_th   
      REAL, PARAMETER :: ri_c=0.25
      REAL :: lambda_m, lambda_h 
      REAL :: a_sub, b_sub, c_sub
      REAL :: g_sub, pr_n
      REAL, DIMENSION(kkp) :: dz, dzn

      ! Define grid increments
      DO k=2, kkp
        dz(k) =z(k) - z(k-1)
        dzn(k)=zn(k) - zn(k-1)
      ENDDO
      dz(1)=dz(2)
      dzn(1)=dzn(2)
      ! Calculate gradient richardson number
      ! shear 
      ! and stability functions
      pr_n=0.7
      a_sub=1/pr_n
      b_sub=40.0
      c_sub=16.0
      g_sub=1.2

      DO k=2, kkp
        delta_s(k)= (u(k)-u(k-1))**2 + (v(k)-v(k-1)+1e-9)**2  
        delta_th(k)= th(k)-th(k-1)
      ENDDO
      delta_th(1)=delta_th(2)
      delta_s(1)=delta_s(2)
      ! delta_s(1) and delta_th(1)
!      PRINT*,'k,ri_m(k)'
      DO k=1, kkp-1
        ri_m(k)= (g/th_ref)*                                           &       
               (dzn(k+1)/delta_s(k+1))*                                &
               (delta_th(k+1)*dz(k)+ delta_th(k)*dz(k+1))/             &
               (dz(k)+dz(k+1))
        shear_m(k)=sqrt(( (u(k+1)-u(k))**2 + (v(k+1)-v(k))**2 )/        & 
                 ((zn(k+1)-zn(k))**2))             
!        PRINT*,k,ri_m(k)
      ENDDO

      DO k=1, kkp-2
        ri_h(k)= (g/th_ref)*                                     &         
                 delta_th(k+1)*dz(k+1)*(dzn(k+1)+dzn(k+2))       &
                 /(delta_s(k+1)*dzn(k+2)+delta_s(k+2)*dzn(k+1))
        shear_h(k)=(1/dz(k+1)) * sqrt(                              &
                   (delta_s(k+1)*dzn(k+2)+delta_s(k+2)*dzn(k+1))  &
                    /(dzn(k+1)+dzn(k+2)))
      ENDDO         

      ! top boundary condition
      ri_m(kkp)=ri_m(kkp-1)
      ri_h(kkp-1)=ri_h(kkp-2)
      ri_h(kkp)=ri_h(kkp-1) 
      shear_m(kkp)=shear_m(kkp-1)

      DO k=1, kkp 

        ! Sharp tails: make a decision re Prandtl no.
        ! at the moment its just the neutral value
        ! LEM CBL functions
        IF ( (ri_m(k) < 0.1) .and. (ri_m(k) >= 0.0) ) THEN  
          f_m(k)=(1.0-5.0*ri_m(k))**2
        ELSEIF (ri_m(k) < 0.0) THEN
          f_m(k)=sqrt(1.0-c_sub*ri_m(k))
        ELSE
          f_m(k)=( 1/(20*ri_m(k)) )**2
        END IF
         
      ENDDO
      DO k=1, kkp 

        ! Sharp tails: make a decision re Prandtl no.
        ! at the moment its just the neutral value
        ! LEM CBL functions
        IF ( (ri_h(k) < 0.1) .and. (ri_h(k) >= 0.0) ) THEN  
          f_h(k)=a_sub*(1.0-5.0*ri_h(k))**2
        ELSEIF (ri_h(k) < 0.0) THEN
          f_h(k)=a_sub*sqrt(1-b_sub*ri_h(k))
        ELSE
          f_h(k)=a_sub*(1/(20*ri_h(k)))**2
        END IF
         
!        ! Long tails 
!        IF (ri(k) > 0.0) THEN
!          f_m(k)=1/(1+10.0*ri(k))
!          f_h(k)=1/(1+10.0*ri(k))
!        ELSE
!          f_m(k)=1.0
!          f_h(k)=1.0
!        ENDIF
        ! Below required to remove noise from boundary layer diffusion.
        f_h(k)=f_m(k)
      ENDDO

      f_m(kkp)=f_m(kkp-1)
      f_h(kkp)=f_h(kkp-1)
     
         
      DO k=1, kkp
        lambda_m=1.0/(1.0/(kar*(z(k)+z_0)) + 1.0/lambda_const)
         
!        lambda=lambda_const
        k_m(k)=(lambda_m**2)*shear_m(k)*f_m(k)
      ENDDO

      DO k=1, kkp-1
        lambda_h=1.0/(1.0/(kar*(zn(k+1)+z_0)) + 1.0/lambda_const)
         
!        lambda=lambda_const
        k_h(k)=(lambda_h**2)*shear_h(k)*f_h(k)
      ENDDO
      k_h(kkp)=k_h(kkp-1)
      ! Artifically over-ride level 1 k's
!      k_m(1)=0.01 
!      k_h(1)=k_m(1)

      END
  

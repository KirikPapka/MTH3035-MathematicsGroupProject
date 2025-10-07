      SUBROUTINE non_gradient( u_star, w_star, wth_0, wqv_0, h_cbl, z, k_m, k_h,&
       du_dt_ng, dv_dt_ng, dth_dt_ng, dqv_dt_ng, wth_ng, wqv_ng, uw_ng, vw_ng)

      ! Calculates the non-gradient contribution
      ! a la Holstlag and Boville 1993 for theta
      ! and following Brown and Grant 1997 for momentum
      IMPLICIT none
      INCLUDE 'params.inc'

      ! Declare Inputs
      REAL :: w_star, wth_0, wqv_0, u_star
      REAL :: h_cbl
      REAL, DIMENSION(kkp) :: k_m, k_h, z
      ! Non-gradient flux
      REAL, DIMENSION(kkp) :: wth_ng, wqv_ng, uw_ng, vw_ng 
      ! Declare Outputs
      ! Explicit tendencies due to non-gradient terms
      REAL, DIMENSION(kkp) :: du_dt_ng, dv_dt_ng, dth_dt_ng, dqv_dt_ng       

      ! Declare locals
      REAL :: gamma_th, gamma_q, w_m
      REAL :: gamma_m
      INTEGER :: k

      w_m= (u_star**3 + c_1_hb*w_star**3)**(1.0/3.0)
      DO k=1,kkp
        wth_ng(k)=0.0
        uw_ng(k)=0.0
        vw_ng(k)=0.0
      ENDDO

      ! Thermal non-gradient term
      DO k=2, kkp 
        IF (z(k) < h_cbl) THEN
          gamma_th= a_hb * w_star * wth_0 /(h_cbl * w_m**2) 
          ! Non-gradient heat flux
          wth_ng(k)= k_h(k) * gamma_th
          dth_dt_ng(k)= -(wth_ng(k)-wth_ng(k-1))/(z(k)-z(k-1))
        ENDIF
      ENDDO
      
      !Water vapour mixing ration non-local term
      ! Thermal non-gradient term
      DO k=2, kkp 
        IF (z(k) < h_cbl) THEN
          gamma_q= a_hb * w_star * wqv_0 /(h_cbl * w_m**2) 
      !    ! Non-gradient heat flux
          wqv_ng(k)= k_h(k) * gamma_q
          dqv_dt_ng(k)= -(wqv_ng(k)-wqv_ng(k-1))/(z(k)-z(k-1))
        ENDIF
      ENDDO

      ! Momentum non-gradient term
      DO k=2, kkp 
        IF (z(k) < h_cbl) THEN
          gamma_m= a_brgr * u_star**2 *  w_star**3 /(h_cbl * w_m**4) 
          ! Non-gradient heat flux
          uw_ng(k)= -k_m(k) * gamma_m
          du_dt_ng(k)= -(uw_ng(k)-uw_ng(k-1))/(z(k)-z(k-1))

          vw_ng(k)=uw_ng(k)
          dv_dt_ng(k)=du_dt_ng(k)
          
        ENDIF
      ENDDO


      END





      

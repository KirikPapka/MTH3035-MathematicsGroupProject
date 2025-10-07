      SUBROUTINE calc_coriolis(u, v, u_g, v_g, du_dt_cor, dv_dt_cor)
     
      IMPLICIT none
      INCLUDE 'params.inc'
      ! Inputs
      REAL, INTENT(IN), DIMENSION(kkp) :: u, v

      INTEGER :: k

      ! Output Coriolis tendency
      REAL, INTENT(OUT), DIMENSION(kkp) :: du_dt_cor, dv_dt_cor


      DO k=1, kkp 
        du_dt_cor(k)=   f_0*(v(k)- V_g)
        dv_dt_cor(k)= - f_0*( u(k) - U_g )
      ENDDO

      END
  

      SUBROUTINE calculate_diffusion(u,v,th,z,zn,k_m,k_h)
     
      IMPLICIT none
      INCLUDE 'params.inc'
      ! Inputs
      REAL, INTENT(IN), DIMENSION(kkp) :: u, v, th
      INTEGER :: k
      ! Outputs
      REAL, INTENT(INOUT), DIMENSION(kkp) :: k_m, k_h

      ! Test values
      DO k=1, kkp
        k_m(k)=1.0
        k_h(k)=1.0
      ENDDO

      END

      SUBROUTINE calc_fluxes(u,v,th,z,zn,k_m,k_h,uw,vw,wth,n_in)
     
      IMPLICIT none
      INCLUDE 'params.inc'
      ! Inputs
      REAL, INTENT(IN), DIMENSION(kkp,2) :: u, v, th
      ! Vertical coordinates
      REAL, INTENT(IN), DIMENSION(kkp) :: z, zn, k_m, k_h

      INTEGER, INTENT(IN) :: n_in
      INTEGER :: k
      ! Outputs fluxes
      REAL, INTENT(INOUT), DIMENSION(kkp) :: uw, vw, wth

      REAL :: dzn, dzz

      DO k=1, kkp-1
        dzn=zn(k+1)-zn(k)
        uw(k)=   -k_m(k)*(u(k+1,n_in)-u(k,n_in))/dzn
        vw(k)=   -k_m(k)*(v(k+1,n_in)-v(k,n_in))/dzn
        wth(k)=  -k_h(k)*(th(k+1,n_in)-th(k,n_in))/dzn  
      ENDDO

      ! Upper boundary condition ensures no turbulent
      ! increments

      uw(kkp)=uw(kkp-1)
      vw(kkp)=vw(kkp-1)
      wth(kkp)=wth(kkp-1)

      END
  

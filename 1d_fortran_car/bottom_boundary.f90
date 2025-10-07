      SUBROUTINE bottom_boundary(u,v,th,wth,du_dt, dv_dt, dth_dt,grid_factor)

      ! Apply bottom boundary condition to u,v, theta

      IMPLICIT none

      INCLUDE 'params.inc'      
      ! Declare Inputs
      ! winds, potential temperature
      REAL, DIMENSION(kkp,2) :: u, v, th
      ! Vertical pot. temp. fluxes
      REAL, DIMENSION(kkp) :: wth
      ! Declare Outputs
      REAL, DIMENSION(kkp) :: du_dt, dv_dt, dth_dt
      ! surface grid factor
      REAL :: grid_factor

      ! Declare locals
      REAL :: dth_s_dt
      INTEGER :: k
   
      ! u=0, v=0 surface boundary condition at surface

      du_dt(1)=-du_dt(2) * grid_factor/(1-grid_factor)
      dv_dt(1)=-dv_dt(2) * grid_factor/(1-grid_factor)

! The lines below should be redundant
!      u(1,1)=-u(2,1)
!      v(1,1)=-v(2,1)

      ! tendency for surface potential temperature
      ! Impose a simple cooling
!      dth_s_dt=-0.25/3600.0
      ! Simple surface warming
       dth_s_dt=1.5/3600.0
 
!      dth_s_dt= (-IR/(rho_0*c_p) - wth(1))/SED



!      PRINT*,'dth_s_dt*3600.0 ',dth_s_dt*3600.0
!      PRINT*,'grid_factor ', grid_factor
      dth_dt(1)=dth_s_dt/(1-grid_factor)   & 
                - dth_dt(2)*grid_factor/(1-grid_factor)
      
!      dth_dt(1)=dth_s_dt
!      PRINT*,'dth_dt(1)*3600.0 ',dth_dt(1)*3600.0
      END






       

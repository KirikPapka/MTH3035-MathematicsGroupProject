      ! ***************************************
      ! 1DBL 
      ! A program to simulate
      ! 1D boundary layers
      ! Author: Bob Beare 
      ! Started: Sept 2007
      ! ***************************************
      ! Converted to Charney-Phillips Mar 2008
      ! ****************************************
      ! Further Developments: Georgios Efstathiou Nov 2020
      ! 1DBLv1.1
      ! * Set linear piecewise initial profiles of potential temperature 
      ! and vapour mixing ratio.
      ! * Specified varying surface sensible/latent heat fluxes
      ! * Vapour mixing ratio diffusion with a non-local term
      ! * Moist convection diagnostics

      PROGRAM main
      IMPLICIT none
      ! Parameters
      INCLUDE 'params.inc'

      INTEGER :: k, n, ntim, ntemp, n_out, n_freq
      
      ! Vertical coordinates
      REAL, DIMENSION(kkp) :: z, zn
      ! Extra coordinate for C-P grid.
      REAL, DIMENSION(kkp-1) :: z_wt
      REAL :: dz_2, dzn_2
      ! Time
      REAL:: t
      ! Lid value of theta
      REAL:: th_lid, qv_lid   

      ! Wind and potential temperature 
      ! on different time levels
      REAL, DIMENSION(kkp,2) :: u, v, th, qv, ql
      
      ! Pressure as a diagnostic
      REAL, DIMENSION(kkp) :: p, rh
      
      ! Initial theta and qv profiles
      real, dimension(kkp) :: th_init, qv_init

      ! Explicit tendencies due to Coriolis
      REAL, DIMENSION(kkp) :: du_dt_cor, dv_dt_cor

      ! Explicit tendencies due to non-gradient terms
      REAL, DIMENSION(kkp) :: du_dt_ng, dv_dt_ng, dth_dt_ng, dqv_dt_ng
      
      ! Explicit tendencies due to moist convection
      REAL, DIMENSION(kkp) :: dth_dt_cu, dqv_dt_cu, dql_dt_cu

      ! Intermediate values after turbulent diffusion
      ! has been applied
      REAL, DIMENSION(kkp) :: u_inter, v_inter, th_inter, qv_inter
      
      ! Vertical diffusion coefficients
      REAL, DIMENSION(kkp) :: k_m, k_h, K_q

      ! Viscous Courant no.
      REAL, DIMENSION(kkp) :: c_vis
       
      ! Outputs fluxes
      REAL, DIMENSION(kkp) :: uw, vw, wth, wqv
      ! Non-gradient flux
      REAL, DIMENSION(kkp) :: wth_ng, wqv_ng, uw_ng, vw_ng       
      ! Gradient Richardson number on flux levels
!      REAL, DIMENSION(kkp) :: ri 
      REAL, DIMENSION(kkp) :: ri_h, ri_m
      ! surface grid factor
      REAL :: grid_factor
      
      ! surface potential temperature
      REAL :: th_surf

      ! Friction velocity
      REAL :: u_star

      ! Obukhov length
      REAL :: L_obuk

      ! boundary layer depth
      REAL :: h_bl

      ! Convective velocity scale
      REAL :: w_star
     
      ! Suface layer non-dimensional u gradient
      REAL :: phi_m
      REAL :: w_m
 
      ! Parameters for implicit solver
       
      REAL :: b_1, c_1, r_1

      ! Surface sensible heat flux and latent heat flux
      REAL :: h_0, q_0

      REAL :: pi
      
      ! Moist dynamics variables
      real :: z_cb, z_lcl
      real :: z_cloud_base, z_cloud_top
      
      ! String for multiple output
      character(len=6) :: str

      pi=4.0*atan(1.0)

      ! Read in namelist
      OPEN(5,file='nml_1d_bl')
      READ(5,nml=run)
      READ(5,nml=geowind)
      CLOSE(5) 

      PRINT*,'run_length (h): ',run_length/3600.0
      PRINT*,'freq_out (mins): ',freq_out/60.0
      PRINT*,'time_step (s):  ',dt
      PRINT*,'Geostrophic u wind ',U_g
      PRINT*,'Geostrophic v wind', V_g
!      h_0=h_sbl
!      h_0=0
      
      ! Initialise 
      ! Call routine to define grid  
      CALL vert_levs(z,zn,z_wt)

      ! Initialise the wind

      CALL init_wind(u(1,1), v(1,1), zn, u_g, v_g)

      ! Initialise potential temperature
      call init_theta(z, th_init)
        th(:,1)=th_init
        th_lid=th(kkp,1)
        
      call init_qv(z, qv_init)
        qv(:,1)=qv_init
        qv_lid=qv(kkp,1)

      !Initialise bl depth
      h_bl=10.

      ! Surface potential temperature
      th_surf=th_ref !- th_i


      grid_factor= (z(1)-zn(1)) / (zn(2)-zn(1))


      ! u=0 at z=0
      u(1,1)=-u(2,1) * grid_factor/(1-grid_factor)
      ! v=0 at z=0
      v(1,1)=-v(2,1) * grid_factor/(1-grid_factor)

      ! th=th0_init at z=0
    !  th(1,1)=th0_init
        
      ! Calculate local diffusion
      CALL calculate_diffusion(u(1,1),v(1,1),th(1,1),z,zn,k_m  &
      ,k_h,ri_m,ri_h)

      ! Calculate viscous courant no.
      DO k=1,kkp-1
        c_vis(k)=abs( k_m(k)*dt/ ( (zn(k+1)-zn(k))**2 ) )
      ENDDO

      CALL calc_fluxes(u,v,th,z,zn,k_m,k_h,uw,vw,wth,1)

      u_star= (uw(1)**2 + vw(1)**2)**0.25

      IF (surf_opt==1) THEN
        call set_fluxes_l(0.0 ,h_0, q_0)
        wth(1)=h_0/(rho_0*c_p)
      ENDIF
      L_obuk= - th_ref*u_star**3/( kar * g * wth(1) ) 

      IF ((wth(1) > 0) .and. (surf_opt==1)) THEN
        ! phi_m from Holtslag and Boville '93
        phi_m=(1- 15*zn(2)/L_obuk)**(-1.0/3.0)
        w_m=u_star/phi_m
        th_surf=th(k_sl,1) + b_hb*wth(1)/w_m 
 ! above not required for high res grid
      ELSEIF ((wth(1) < 0) .and. (surf_opt==1)) THEN
        th_surf=th(k_sl,1)
      ENDIF
             
      ! no. of timesteps
      ntim=INT(run_length/dt)
      
      !no. of output fragments
      n_out=INT(run_length/freq_out)+1

      ! No. of timesteps for each output
      n_freq=INT(freq_out/dt)

      ! Write important integers to file
      ! Open initialisation output file
      OPEN(10,file='output/arm_shcu_init.dat')  
     ! WRITE(10,'(2i5)') n_out, kkp
     ! WRITE(10,'(F10.4)') freq_out

      ! Write single level vars
      WRITE(10,'(5F18.7)') u_star, w_star, l_obuk, h_bl, th_surf
      ! Write initial profile to output
      DO k=1, kkp
        WRITE(10,'(13F18.7)') zn(k),z(k),u(k,1),v(k,1),th(k,1),1000*qv(k,1) &
       ,uw(k), vw(k), wth(k),k_m(k),k_h(k),c_vis(k)
      ENDDO
      WRITE(10,*) '*********Next time ****************'
      close (10)


      PRINT*,' Executing ',ntim,' timesteps' 
      t=0.0
        
      ! Time integration
      DO n=1,ntim

        ! Calculate diffusion coefficients for timestep n
        CALL calculate_diffusion(u(1,1),v(1,1),th(1,1),z,zn,k_m, &
        k_h,ri_m,ri_h)
 
        ! Calculate viscous courant no.
        DO k=1,kkp-1
          c_vis(k)=abs( k_m(k)*dt/ ( (zn(k+1)-zn(k))**2 ) )
        ENDDO        
 
        ! Calculate explicit fluxes for diagnosis 
        CALL calc_fluxes(u,v,th,z,zn,k_m,k_h,uw,vw,wth,1)

        u_star= (uw(1)**2 + vw(1)**2)**0.25

        IF (surf_opt==1) THEN
          wth(1)=h_0/(rho_0*c_p)
        ENDIF

        L_obuk= - th_ref*u_star**3/( kar * g * wth(1) )       

        IF ((wth(1) > 0) .and. (surf_opt==1)) THEN
          ! phi_m from Holtslag and Boville '93
          phi_m=(1- 15*zn(2)/L_obuk)**(-1.0/3.0)
          w_m=u_star/phi_m
          th_surf=th(k_sl,1) + b_hb*wth(1)/w_m

        ELSEIF ((wth(1) < 0) .and. (surf_opt==1)) THEN
          th_surf=th(k_sl,1)
        ENDIF

!        IF (surf_opt==1) THEN
!          th_surf=th(2,1)
!        ENDIF

        ! diagnose bl depth
        CALL bl_depth(u(1,1), v(1,1), th(1,1), th_surf, zn, z, &
                      u_star, L_obuk, wth(1), ri_m, h_bl)

        ! Diagnose w_star
        IF (wth(1) > 0 ) THEN
          w_star=(h_bl*g*wth(1)/th_ref)**(1.0/3.0)
          ! zero coefficients
!          DO k=2,kkp
!            k_m(k)=0.0
!            k_h(k)=0.0
!          ENDDO
          CALL diff_cbl(u_star,w_star,h_bl,z,zn,L_obuk,k_m,k_h)
          ! Recalculate fluxes
          CALL calc_fluxes(u,v,th,z,zn,k_m,k_h,uw,vw,wth,1)

          u_star= (uw(1)**2 + vw(1)**2)**0.25

          IF (surf_opt==1) THEN
            call set_fluxes_l(t, h_0, q_0)
            wth(1)=h_0/(rho_0*c_p)
          ENDIF
          CALL non_gradient( u_star, w_star, wth(1), q_0, h_bl, z, k_m, k_h,&
             du_dt_ng, dv_dt_ng, dth_dt_ng, dqv_dt_ng, wth_ng, wqv_ng, uw_ng, vw_ng)

          DO k=1,kkp
            wth(k)=wth(k)+wth_ng(k)
            wqv(k)=wqv(k)+wqv_ng(k)
            uw(k)=uw(k)+uw_ng(k)
            vw(k)=vw(k)+vw_ng(k)
          ENDDO
        ELSE
          w_star=0.0
          DO k=1, kkp
            dth_dt_ng(k)=0.0
            dqv_dt_ng(k)=0.0
            du_dt_ng(k)=0.0
            dv_dt_ng(k)=0.0
          ENDDO
        ENDIF
 
        ! Time integration 
        CALL implicit_solver(u(1,1),k_m,z,zn,dt,0.5,0.5,0.0,u_g,u_inter)
        CALL implicit_solver(v(1,1),k_m,z,zn,dt,0.5,0.5,0.0,v_g,v_inter)
       
        IF (surf_opt==0) THEN      
          ! Prescribed surface temperature

          !Increment th_surf to n+1 value
          th_surf=th_surf+ heating*dt/3600.0

          ! Calculate intermediate values after turbulent diffusion  

          CALL implicit_solver_th(th(1,1),k_h,z,zn,dt,1.0   &
               ,0.0,th_surf &
               ,th_lid,th_inter) 
        ELSEIF (surf_opt==1) THEN
          call set_fluxes_l(t, h_0, q_0)
          !for comparison with LES
!          h_0=100.0*t/(6.0*3600.0)
          CALL implicit_solver_flux_th(th(1,1),u(1,1),v(1,1)  &
               ,k_h,z,zn,dt,th_lid,th_inter,h_0)
          CALL implicit_solver_flux_qv(qv(1,1),u(1,1),v(1,1)  &
               ,k_h,z,zn,dt,qv_lid,qv_inter,q_0)

        ENDIF                           
       
      ! Calculate Coriolis
      CALL calc_coriolis(u(1,1), v(1,1), u_g, v_g, du_dt_cor, dv_dt_cor)
      
     ! call pressure_from_theta(z, th(:,1), p)
     ! call shallow_conv_adj(z, p, th(:,1), qv(:,1), ql(:,1),  &
     !               dth_dt_cu, dqv_dt_cu, dql_dt_cu, z_cloud_base, z_cloud_top)
  
      DO k=1, kkp      
        u(k,2) =u_inter(k)  + dt*du_dt_cor(k) + dt*du_dt_ng(k)
        v(k,2) =v_inter(k)  + dt*dv_dt_cor(k) + dt*dv_dt_ng(k)
        th(k,2)=th_inter(k)                   + dt*dth_dt_ng(k) !+ dt*dth_dt_cu(k)
        qv(k,2)=qv_inter(k)                   + dt*dqv_dt_ng(k) !+ dt*dqv_dt_cu(k)
        ql(k,2)=  0.0                                           !   + dt*dql_dt_cu(k)
        ql(k,2)=max(1e-7, ql(k,2))
      ENDDO 

      ! swap variables
      DO k=1, kkp
        u(k,1) =u(k,2) 
        v(k,1) =v(k,2)  
        th(k,1)=th(k,2)  
        qv(k,1)=qv(k,2) 
        ql(k,1)=ql(k,2)       
      ENDDO
      
      ! Moist dynamics diagnostics
        if (moist_dyn) then
        call pressure_from_theta(z, th(:,1), p)
        call lcl_from_Tq(z, p(1), th(1,1), qv(1,1), z_lcl, z_cb)
        call rh_z(p, th(:,1),qv(:,1),rh)
        endif
 
        ! Update time
        t=t+dt
        
        IF ( MOD(n,n_freq) == 0) THEN

          ! Write single level vars
          write(str , '(i5)') int(t)
          open ( 10,file='output/arm_shcu_t_'//trim(adjustl(str))//'.dat')       
          if (moist_dyn) then
          WRITE(10,'(10F18.7)') u_star, w_star, l_obuk, h_bl, th_surf, wth(1), q_0, z_cb, &
            z_cloud_base, z_cloud_top
           ! Write vertical profiles
           DO k=1, kkp
            WRITE(10,'(15F18.7)') zn(k),z(k), p(k), u(k,1),v(k,1),th(k,1),1000*qv(k,1), &
            1000*ql(k,1),rh(k),uw(k), vw(k), wth(k),k_m(k),k_h(k),c_vis(k)
          ENDDO
        
          else  
          WRITE(10,'(7F18.7)') u_star, w_star, l_obuk, h_bl, th_surf, wth(1), q_0
          ! Write vertical profiles
          DO k=1, kkp
            WRITE(10,'(13F18.7)') zn(k),z(k),u(k,1),v(k,1),th(k,1),1000*qv(k,1) &
            ,uw(k), vw(k), wth(k),k_m(k),k_h(k),c_vis(k)
          ENDDO
          endif
      
          WRITE(10,*) '*********Next time ****************'
          close(10)
       !   write(*,*) 'cloud base = ', z_cloud_base, '  cloud top = ', z_cloud_top
        !  call flush(6)
          
        ENDIF

      ENDDO ! loop over n

!      CLOSE(10)  ! Close output file 

      END  !End of program
      INCLUDE 'calculate_diffusion.f90'
      INCLUDE 'calc_fluxes.f90'
      INCLUDE 'calc_coriolis.f90'
      INCLUDE 'vert_levs.f90'
      INCLUDE 'implicit_solver.f90'
      INCLUDE 'implicit_solver_th.f90'
      INCLUDE 'implicit_solver_flux_th.f90'
      INCLUDE 'implicit_solver_flux_qv.f90'
      INCLUDE 'tridag.f90'
      INCLUDE 'init_wind.f90'
      INCLUDE 'bl_depth.f90'
      INCLUDE 'diff_cbl.f90'
      INCLUDE 'non_gradient.f90'
      include 'set_fluxes.f90'
      include 'set_fluxes_l.f90'
      include 'init_theta.f90'
      include 'init_qv.f90'
      include 'pressure_from_theta.f90'
      include 'lcl_from_Tq.f90'
      include 'shallow_conv_adj.f90'
      include 'rh_z.f90'
     
      

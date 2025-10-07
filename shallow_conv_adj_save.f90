!==============================================================================
! Shallow convective adjustment on a w-grid (interfaces)
! Inputs on w-levels: theta_w, qv_w, qc_w, p_w, rho_w, z_w
! Non-precipitating: all condensation goes to qc (no autoconversion)
! Outputs tendencies at w-levels and cloud base/top heights (m).
!------------------------------------------------------------------------------
subroutine shallow_conv_adj(z_w, p_w, theta_w, qv_w, qc_w,  &
                           dtheta_dt, dqv_dt, dqc_dt, z_cloud_base, z_cloud_top)
  implicit none

  INCLUDE 'params.inc'
  !-------------------- arguments --------------------
  REAL, DIMENSION(kkp) , INTENT(in) :: z_w, theta_w, p_w, qv_w, qc_w   
  REAL, DIMENSION(kkp) , INTENT(out) :: dtheta_dt, dqv_dt, dqc_dt
  real, intent(out) :: z_cloud_base, z_cloud_top
  !-------------------- tunables ---------------------
  real, parameter :: tau_adj  = 1800.d0     ! relaxation timescale [s]
  real, parameter :: tol_stab = 1.0d-8      ! stability tolerance
  real, parameter :: dz_min   = 0.5d0       ! min layer thickness [m]
  real, parameter :: qeps     = 1.0d-9      ! moisture floor
  real, parameter :: qc_cloud = 1.0d-7      ! threshold for cloud detection [kg/kg]
  !-------------------- locals -----------------------
  REAL, DIMENSION(kkp) :: rho_w
  integer :: k, k1, k2, kb
  real :: kappa, alpha
  real, allocatable :: Pi(:), T(:), qv(:), qc(:), tl(:), qt(:), thv(:)
  real, allocatable :: theta_new(:), qv_new(:), qc_new(:)
  real, allocatable :: dz(:)
  integer :: icb, ict

  rho_w = 1.0
  dtheta_dt = 0.d0; dqv_dt = 0.d0; dqc_dt = 0.d0
  z_cloud_base = -1.d0; z_cloud_top = -1.d0
  if (kkp <= 1) return

  kappa = Rd/c_p
  alpha = 1.d0 - exp( -max(0.d0, dt) / max(1.d0, tau_adj) )

  allocate(Pi(kkp), T(kkp), qv(kkp), qc(kkp), tl(kkp), qt(kkp), thv(kkp))
  allocate(theta_new(kkp), qv_new(kkp), qc_new(kkp))
  allocate(dz(kkp-1))

  do k = 1, kkp
     Pi(k) = (p_w(k)/p0)**kappa
  end do
  do k = 1, kkp-1
     dz(k) = max(dz_min, z_w(k+1) - z_w(k))
  end do

  ! Copy inputs
  qv = max(qeps, qv_w); qc = max(0.d0, qc_w)
  T  = theta_w * Pi      ! T = θ Π at fixed level

  !---------------------------------------------------------------
  ! 1) Local saturation adjustment (no rain): trade qv <-> qc, adjust T
  !---------------------------------------------------------------
  do k = 1, kkp
     call sat_adjust_no_rain_T(T(k), qv(k), qc(k), p_w(k), c_p, Lv, Rd, Rv)
  end do

  ! 2) Build near-conserved variables and virtual potential temperature
  do k = 1, kkp
     qt(k) = qv(k) + qc(k)
     tl(k) = theta_l_from_Tqc(T(k), qc(k), p_w(k), p0, c_p, Rd, Lv)
     thv(k)= theta_from_T(T(k), p_w(k), p0, Rd, c_p) * (1.d0 + 0.61d0*qv(k) - qc(k))
  end do

  !---------------------------------------------------------------
  ! 3) Find unstable blocks in thv (decreasing upward) and mix in (theta_l, qt)
  !    Blocks are defined on w-level indices k1..k2, k2>k1.
  !---------------------------------------------------------------
  kb = 0
  do k = 1, kkp-1
     if (thv(k+1) + tol_stab < thv(k)) then
        if (kb == 0) kb = k
     else
        if (kb /= 0) then
           k1 = kb; k2 = k
           call mix_block_w(k1, k2, z_w, rho_w, tl, qt)
           kb = 0
        end if
     end if
     if (k == kkp-1 .and. thv(k+1) + tol_stab < thv(k)) then
        k1 = kb; k2 = k+1
        call mix_block_w(k1, k2, z_w, rho_w, tl, qt)
        kb = 0
     end if
  end do

  !---------------------------------------------------------------
  ! 4) Reconstruct (T, qv, qc) from (θ_l, q_t, p) at each w-level
  !---------------------------------------------------------------
  do k = 1, kkp
     call recover_T_qv_qc_from_tl_qt(tl(k), qt(k), p_w(k), p0, c_p, Rd, Rv, Lv, &
                                     T(k), qv(k), qc(k))
  end do
  do k = 1, kkp
     theta_new(k) = theta_from_T(T(k), p_w(k), p0, Rd, c_p)
     qv_new(k)    = max(qeps, qv(k))
     qc_new(k)    = max(0.d0, qc(k))
  end do

  ! Optional partial relaxation toward adjusted state
  do k = 1, kkp
     theta_new(k) = theta_w(k) + alpha * (theta_new(k) - theta_w(k))
     qv_new(k)    = qv_w(k)    + alpha * (qv_new(k)    - qv_w(k))
     qc_new(k)    = qc_w(k)    + alpha * (qc_new(k)    - qc_w(k))
  end do

  !---------------------------------------------------------------
  ! 5) Tendencies and cloud base/top diagnosis
  !---------------------------------------------------------------
  
    ! ---- LIMITERS to prevent negative qv/qc after caller's update ----
  do k = 1, kkp
     ! raw tendencies from adjusted -> original
     dtheta_dt(k) = (theta_new(k) - theta_w(k)) / max(1.d-9, dt)

     ! Limit moisture tendencies so next-step state stays >= floors
     ! qv floor = qeps; qc floor = 0
     dqv_dt(k) = (qv_new(k) - qv_w(k)) / max(1.d-9, dt)
     dqc_dt(k) = (qc_new(k) - qc_w(k)) / max(1.d-9, dt)

     dqv_dt(k) = max( (qeps - qv_w(k)) / max(1.d-9,dt), dqv_dt(k) )
     dqc_dt(k) = max( (0.d0 - qc_w(k)) / max(1.d-9,dt), dqc_dt(k) )
  end do

  do k = 1, kkp
     dtheta_dt(k) = (theta_new(k) - theta_w(k)) / max(1.d-9, dt)
     dqv_dt(k)    = (qv_new(k)    - qv_w(k))    / max(1.d-9, dt)
     dqc_dt(k)    = (qc_new(k)    - qc_w(k))    / max(1.d-9, dt)
  end do

  ! Cloud base/top from qc_new on w-levels
  icb = 0; ict = 0
  do k = 1, kkp
     if (qc_new(k) > qc_cloud) then
        icb = k; exit
     end if
  end do
  if (icb > 0) then
     do k = kkp, 1, -1
        if (qc_new(k) > qc_cloud) then
           ict = k; exit
        end if
     end do
  end if

  if (icb > 0 .and. ict >= icb) then
     z_cloud_base = z_w(icb)
     z_cloud_top  = z_w(ict)
  else
     z_cloud_base = -1.d0
     z_cloud_top  = -1.d0
  end if

contains
  ! ---------- helpers (local scope) ----------------------------------

  pure function theta_from_T(T, p, p0, Rd, c_p) result(theta)
    implicit none
    real, intent(in) :: T, p, p0, Rd, c_p
    real :: theta, kappa
    kappa = Rd/c_p
    theta = T * (p0/max(1.d0,p))**kappa
  end function theta_from_T

  pure function qsat_ew(T, p, Rd, Rv) result(qs)
    implicit none
    real, intent(in) :: T, p, Rd, Rv
    real :: es, eps, qs
    eps = 0.622d0
    es  = 611.2d0 * exp(17.67d0*(T-273.15d0)/(T-29.65d0))
    qs  = eps*es / max(1.d-6, (p - es))
  end function qsat_ew

  ! Local saturation adjustment exchanging qv <-> qc, updating T. No rain.
  subroutine sat_adjust_no_rain_T(T, qv, qc, p, c_p, Lv, Rd, Rv)
    implicit none
    real, intent(inout) :: T, qv, qc
    real, intent(in)    :: p, c_p, Lv, Rd, Rv
    real :: qs, dq
    qs = qsat_ew(T, p, Rd, Rv)
    if (qv > qs) then
       dq = qv - qs
       qv = qs
       qc = max(0.d0, qc + dq)
       T  = T + (Lv/c_p) * dq
    elseif (qv < qs .and. qc > 0.d0) then
       dq = min(qs - qv, qc)
       qv = qv + dq
       qc = qc - dq
       T  = T - (Lv/c_p) * dq
    end if
  end subroutine sat_adjust_no_rain_T

  ! Liquid-water potential temperature from (T, qc, p)
  pure function theta_l_from_Tqc(T, qc, p, p0, c_p, Rd, Lv) result(tl)
    implicit none
    real, intent(in) :: T, qc, p, p0, c_p, Rd, Lv
    real :: tl, kappa, theta
    kappa = Rd/c_p
    theta = T * (p0/max(1.d0,p))**kappa
    tl = theta * exp( -Lv*max(0.d0,qc) / max(1.d-6, c_p*T) )
  end function theta_l_from_Tqc

  ! Mix a block [k1..k2] on w-levels by mass-weighting θ_l and q_t.
  ! Layer masses use interface densities: m_k ≈ 0.5*(ρ_k+ρ_{k+1})*Δz_k for k=k1..k2-1
  subroutine mix_block_w(k1, k2, z_w, rho_w, tl, qt)
    implicit none
    integer, intent(in) :: k1, k2
    real, intent(in) :: z_w(:), rho_w(:)
    real, intent(inout) :: tl(:), qt(:)
    integer :: k
    real :: sum_m, sum_tl, sum_qt, dzloc, mlay, tl_layer, qt_layer

    sum_m = 0.d0; sum_tl = 0.d0; sum_qt = 0.d0
    do k = k1, k2-1
       dzloc    = max(0.5d0, z_w(k+1) - z_w(k))
       mlay     = 0.5d0 * (rho_w(k) + rho_w(k+1)) * dzloc
       tl_layer = 0.5d0 * (tl(k) + tl(k+1))
       qt_layer = 0.5d0 * (qt(k) + qt(k+1))
       sum_m  = sum_m  + mlay
       sum_tl = sum_tl + mlay * tl_layer
       sum_qt = sum_qt + mlay * qt_layer
    end do
    if (sum_m <= 0.d0) return
    sum_tl = sum_tl / sum_m
    sum_qt = max(0.d0, sum_qt / sum_m)
    do k = k1, k2
       tl(k) = sum_tl
       qt(k) = sum_qt
    end do
  end subroutine mix_block_w

  ! Recover (T, qv, qc) from given (θ_l, q_t, p) via a small Newton solve in T.
  subroutine recover_T_qv_qc_from_tl_qt(tl, qt, p, p0, c_p, Rd, Rv, Lv, T, qv, qc)
    implicit none
    real, intent(in)  :: tl, qt, p, p0, c_p, Rd, Rv, Lv
    real, intent(out) :: T, qv, qc
    real :: kappa, theta, qs, f, df, expfac, dqs
    integer :: it

    kappa = Rd/c_p
    ! Start from unsaturated guess: θ = θ_l
    theta = tl
    T     = theta * (p/max(1.d0,p0))**kappa

    do it = 1, 8
       qs = qsat_ew(T, p, Rd, Rv)
       if (qt <= qs) then
          qv = qt; qc = 0.d0
          return
       else
          ! θ_l = θ * exp(-Lv*(qt-qs)/(c_p*T))
          expfac = exp( -Lv*(qt-qs) / max(1.d-6, c_p*T) )
          f  = tl - ( T * (p/max(1.d0,p0))**kappa ) * expfac

          ! dθ/dT = θ/T; d(expfac)/dT = expfac * [ Lv*(qt-qs)/(c_p*T^2) - (Lv/c_p)*( -dqs/dT )/T ]
          dqs = dqs_dT(T, p, Rd, Rv)
          df = - ( (p/max(1.d0,p0))**kappa ) * ( expfac + (T*expfac) * &
               ( + Lv*(qt-qs)/(c_p*T*T) - (Lv/c_p)*( -dqs ) / max(1.d-12, T) ) )

          if (abs(df) < 1.d-10) then
             T = max(150.d0, min(330.d0, T - 0.2d0*f))
          else
             T = max(150.d0, min(330.d0, T - f/df))
          end if
       end if
    end do

    ! finalize saturation state
    qs = qsat_ew(T, p, Rd, Rv)
    if (qt <= qs) then
       qv = qt; qc = 0.d0
    else
      qv = qs; qc = max(0.d0, qt - qs)
    end if
  end subroutine recover_T_qv_qc_from_tl_qt

  pure function dqs_dT(T, p, Rd, Rv) result(dq)
    implicit none
    real, intent(in) :: T, p, Rd, Rv
    real :: es, desdT, eps
    real :: dq
    eps   = 0.622d0
    es    = 611.2d0 * exp(17.67d0*(T-273.15d0)/(T-29.65d0))
    desdT = es * 17.67d0 * 29.65d0 / ( (T-29.65d0)**2 )
    dq    = eps * ( desdT*(p - es) + es*desdT ) / max(1.d-12, (p - es)**2)
  end function dqs_dT

end subroutine shallow_conv_adj
!==============================================================================

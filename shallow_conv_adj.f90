subroutine shallow_conv_adj(z_w, p_w, theta_w, qv_w, qc_w,  &
                            dtheta_dt, dqv_dt, dqc_dt, z_cloud_base, z_cloud_top)
  implicit none
  INCLUDE 'params.inc'   ! expects: kkp, p0, Rd, Rv, c_p, g, Lv, dt

  ! arguments
  real, dimension(kkp), intent(in)  :: z_w, p_w, theta_w, qv_w, qc_w
  real, dimension(kkp), intent(out) :: dtheta_dt, dqv_dt, dqc_dt
  real,               intent(out)   :: z_cloud_base, z_cloud_top

  ! tunables / guards
  real, parameter :: tau_adj   = 1800.0
  real, parameter :: buoy_tol  = 0.0        ! K of θv parcel minus env to count as buoyant
  real, parameter :: qc_thr    = 1.0e-6
  real, parameter :: qeps      = 1.0e-8
  real, parameter :: pmin      = 10000.0    ! Pa; safer top-floor to avoid fake saturation aloft
  real, parameter :: pmax      = 2000000.0
  real, parameter :: Tmin      = 150.0
  real, parameter :: Tmax      = 340.0
  real, parameter :: bl_tol_thv = 0.3       ! K; θv uniformity tolerance for BL mixing
  integer, parameter :: nmin_cloud_levels = 2

  ! locals
  integer :: k, k_lcl, k1, k2, k_bl_top
  real :: kappa, alpha, denom_dt
  real :: z_h, z_LCL, p_LCL, T_src, Td_src, T_LCL
  real :: tl_src, qt_src

  ! environment state
  real, allocatable :: Pi(:), T(:), qv(:), qc(:), qt(:), thv_env(:)
  real, allocatable :: dpm(:)

  ! parcel reference from LCL upward
  real, allocatable :: T_ref(:), qv_ref(:), qc_ref(:), theta_ref(:), thv_ref(:)

  ! outputs/targets
  real, allocatable :: theta_new(:), qv_new(:), qc_new(:)

  logical :: in_cloud
  integer :: runtt

  ! init outputs
  dtheta_dt = 0.0; dqv_dt = 0.0; dqc_dt = 0.0
  z_cloud_base = -1.0; z_cloud_top = -1.0
  if (kkp <= 1) return

  ! coefficients
  kappa = Rd / max(1.0e-12, c_p)
  if (tau_adj <= 0.0) then
     alpha    = 1.0
     denom_dt = max(1.0, dt)
  else
     if (dt > 0.0) then
        alpha    = 1.0 - exp(-dt/tau_adj)
        denom_dt = dt
     else
        alpha    = 0.0
        denom_dt = 1.0
     end if
  end if

  ! allocate
  allocate(Pi(kkp), T(kkp), qv(kkp), qc(kkp), qt(kkp), thv_env(kkp))
  allocate(T_ref(kkp), qv_ref(kkp), qc_ref(kkp), theta_ref(kkp), thv_ref(kkp))
  allocate(theta_new(kkp), qv_new(kkp), qc_new(kkp))
  allocate(dpm(kkp-1))

  ! Exner and layer masses (with guards)
  do k = 1, kkp
     Pi(k) = ( max(pmin, min(pmax, p_w(k))) / max(pmin, p0) )**kappa
  end do
  do k = 1, kkp-1
     dpm(k) = abs( max(pmin,min(pmax,p_w(k))) - max(pmin,min(pmax,p_w(k+1))) ) / max(1.0e-12, g)
  end do

  ! environment (guarded)
  qv = max(qeps, qv_w)
  qc = max(0.0,  qc_w)
  T  = max(Tmin, min(Tmax, theta_w * Pi))
  qt = qv + qc
  do k = 1, kkp
     thv_env(k) = theta_w(k) * (1.0 + 0.61*qv(k) - qc(k))
  end do

  ! ---- diagnose BL height z_h from θv uniformity ----
  call diagnose_bl_height_mixed(z_w, thv_env, bl_tol_thv, k_bl_top, z_h)

  ! ---- BL-mean parcel and its LCL (Bolton-style) ----
  call bl_means(theta_w, T, qv, dpm, k_bl_top, T_src, qt_src, tl_src)
  ! dewpoint and TLCL from BL-mean parcel at near-sfc pressure
  call dewpoint_from_qv_p(qt_src, max(pmin,min(pmax,p_w(1))), Td_src)
  call tlcl_bolton(T_src, Td_src, T_LCL)
  call lcl_pressure_from_T_qt(T_LCL, qt_src, p_LCL)       ! uses eps and es(T_LCL)

  ! map p_LCL -> z_LCL, and find index k_lcl
  call interp_z_at_p(z_w, p_w, max(pmin,min(pmax,p_LCL)), z_LCL, k_lcl)

  ! ---- TRIGGER: only proceed if BL reaches LCL ----
  if (z_h + 1.0e-6 < z_LCL) then
     ! BL top below LCL => no convection
     theta_new = theta_w
     qv_new    = qv
     qc_new    = qc
     ! tendencies are already zero; exit cleanly
     deallocate(Pi,T,qv,qc,qt,thv_env,T_ref,qv_ref,qc_ref,theta_ref,thv_ref,theta_new,qv_new,qc_new,dpm)
     return
  end if

  ! ---- build parcel reference from LCL upward (reversible moist adiabat) ----
  ! invariants: tl_src and qt_src (from BL means). Below LCL, ref = env.
  do k = 1, kkp
     if (k < k_lcl) then
        T_ref(k)     = T(k)
        qv_ref(k)    = qv(k)
        qc_ref(k)    = 0.0
        theta_ref(k) = theta_w(k)
        thv_ref(k)   = thv_env(k)
     else
        call recover_T_qv_qc_from_tl_qt( tl_src, qt_src,                     &
             max(pmin,min(pmax,p_w(k))), p0, c_p, Rd, Rv, Lv, Tmin, Tmax,    &
             T_ref(k), qv_ref(k), qc_ref(k) )
        theta_ref(k) = theta_from_T( T_ref(k), max(pmin,min(pmax,p_w(k))), p0, Rd, c_p )
        thv_ref(k)   = theta_ref(k) * (1.0 + 0.61*qv_ref(k) - qc_ref(k))
     end if
  end do

  ! ---- determine saturated buoyant cloud layer using Δθv > buoy_tol ----
  k1 = -1; k2 = -1; in_cloud = .false.; runtt = 0
  do k = 1, kkp
     if ( (qc_ref(k) > qc_thr) .and. (thv_ref(k) - thv_env(k) > buoy_tol) .and. (k >= k_lcl) ) then
        if (.not. in_cloud) then
           k1 = k; in_cloud = .true.; runtt = 1
        else
           runtt = runtt + 1
        end if
     else
        if (in_cloud) then
           if (runtt >= nmin_cloud_levels) then
              k2 = k-1
              exit
           else
              in_cloud = .false.; runtt = 0; k1 = -1
           end if
        end if
     end if
  end do
  if (in_cloud .and. runtt >= nmin_cloud_levels .and. k2 < 0) k2 = kkp

  ! defaults
  theta_new = theta_w
  qv_new    = qv
  qc_new    = qc

  if (k1 > 0 .and. k2 >= k1) then
     do k = k1, k2
        theta_new(k) = theta_ref(k)
        qv_new(k)    = qv_ref(k)
        qc_new(k)    = qc_ref(k)
     end do
     z_cloud_base = z_w(k1)
     z_cloud_top  = z_w(k2)
  else
     z_cloud_base = -1.0
     z_cloud_top  = -1.0
  end if

  ! ---- tendencies ----
  if (alpha > 0.0) then
     do k = 1, kkp
        dtheta_dt(k) = alpha * (theta_new(k) - theta_w(k)) / denom_dt
        dqv_dt(k)    = alpha * (qv_new(k)    - qv_w(k))     / denom_dt
        dqc_dt(k)    = alpha * (qc_new(k)    - qc_w(k))     / denom_dt
     end do
  else
     dtheta_dt = 0.0; dqv_dt = 0.0; dqc_dt = 0.0
  end if

  ! tidy
  deallocate(Pi,T,qv,qc,qt,thv_env,T_ref,qv_ref,qc_ref,theta_ref,thv_ref,theta_new,qv_new,qc_new,dpm)

contains
  !=================== thermodynamic helpers ===================
  pure function theta_from_T(T, p, p0, Rd, c_p) result(theta)
    real, intent(in) :: T, p, p0, Rd, c_p
    real :: theta
    theta = T * ( p0/max(1.0e-12,p) )**( Rd/max(1.0e-12,c_p) )
  end function theta_from_T

  pure function qsat_ew(T, p, Rd, Rv) result(qs)
    ! Tetens (liquid) with guards
    real, intent(in) :: T, p, Rd, Rv
    real :: es, eps, qs
    eps = Rd/max(1.0e-12,Rv)
    es  = 611.2 * exp( max(-60.0, min(60.0, 17.67*(T-273.15)/(T-29.65))) )
    es  = min(0.99*max(50.0,p), es)
    qs  = eps * es / max(1.0e-6, p - es)
  end function qsat_ew

  subroutine recover_T_qv_qc_from_tl_qt(tl, qt, p, p0, c_p, Rd, Rv, Lv, Tmin, Tmax, T, qv, qc)
    ! Robust saturation adjustment (bisection with guards)
    real, intent(in)  :: tl, qt, p, p0, c_p, Rd, Rv, Lv, Tmin, Tmax
    real, intent(out) :: T, qv, qc
    real :: kappa, Pi, Tdry, qs_dry, f_low, f_high, Tlow, Thigh, Tmid, fmid
    integer :: it, itmax

    kappa = Rd/max(1.0e-12,c_p)
    Pi    = (p/max(1.0e-12,p0))**kappa
    Tdry  = max(Tmin, min(Tmax, tl*Pi))
    qs_dry = qsat_ew(Tdry, p, Rd, Rv)

    if (qt <= qs_dry) then
       T  = Tdry
       qc = 0.0
       qv = max(0.0, qt)
       return
    end if

    itmax = 50
    Tlow  = max(Tmin, Tdry - min(80.0, Lv*max(0.0,qt)/max(1.0e-6,c_p)))
    Thigh = min(Tmax, Tdry + 5.0)

    f_low  = theta_l_from_Tqc(Tlow,  max(0.0, qt - qsat_ew(Tlow,  p, Rd, Rv)), p, p0, c_p, Rd, Lv) - tl
    f_high = theta_l_from_Tqc(Thigh, max(0.0, qt - qsat_ew(Thigh, p, Rd, Rv)), p, p0, c_p, Rd, Lv) - tl

    if (f_low*f_high > 0.0) then
       Tlow  = max(Tmin, Tlow  - 15.0)
       Thigh = min(Tmax, Thigh + 15.0)
       f_low  = theta_l_from_Tqc(Tlow,  max(0.0, qt - qsat_ew(Tlow,  p, Rd, Rv)), p, p0, c_p, Rd, Lv) - tl
       f_high = theta_l_from_Tqc(Thigh, max(0.0, qt - qsat_ew(Thigh, p, Rd, Rv)), p, p0, c_p, Rd, Lv) - tl
    end if

    if (f_low*f_high > 0.0) then
       T = Tdry
       do it = 1, itmax
          T = 0.5*T + 0.5*tl*Pi * &
              exp( min(40.0, Lv*max(0.0,qt-qsat_ew(T,p,Rd,Rv))/max(1.0e-6,c_p*T)) )
          T = max(Tmin, min(Tmax, T))
       end do
    else
       do it = 1, itmax
          Tmid = 0.5*(Tlow + Thigh)
          fmid = theta_l_from_Tqc(Tmid, max(0.0, qt - qsat_ew(Tmid, p, Rd, Rv)), p, p0, c_p, Rd, Lv) - tl
          if (abs(fmid) < 1.0e-6 .or. (Thigh-Tlow) < 1.0e-4) exit
          if (f_low*fmid <= 0.0) then
             Thigh = Tmid; f_high = fmid
          else
             Tlow = Tmid;  f_low  = fmid
          end if
       end do
       T = 0.5*(Tlow + Thigh)
    end if

    qc = max(0.0, qt - qsat_ew(T, p, Rd, Rv))
    qv = max(0.0, qt - qc)
  end subroutine recover_T_qv_qc_from_tl_qt

  pure function theta_l_from_Tqc(T, qc, p, p0, c_p, Rd, Lv) result(tl)
    real, intent(in) :: T, qc, p, p0, c_p, Rd, Lv
    real :: tl, expo
    expo = -min(60.0, Lv*max(0.0,qc)/max(1.0e-6,c_p*T))
    tl   = theta_from_T(T, p, p0, Rd, c_p) * exp(expo)
  end function theta_l_from_Tqc

  !=================== BL + LCL helpers ===================
  subroutine diagnose_bl_height_mixed(z, thv, tolK, k_top, z_h)
    real, intent(in)  :: z(:), thv(:), tolK
    integer, intent(out) :: k_top
    real, intent(out) :: z_h
    integer :: j
    real :: thv0
    thv0 = thv(1)
    k_top = 1
    do j = 2, size(thv)
       if ( abs(thv(j) - thv0) <= tolK ) then
          k_top = j
       else
          exit
       end if
    end do
    ! simple interface placement
    if (k_top < size(thv)) then
       z_h = z(k_top)
    else
       z_h = z(size(thv))
    end if
  end subroutine diagnose_bl_height_mixed

  subroutine bl_means(theta, T, qv, dpm, k_top, T_bar, qt_bar, tl_bar)
    real, intent(in) :: theta(:), T(:), qv(:), dpm(:)
    integer, intent(in) :: k_top
    real, intent(out) :: T_bar, qt_bar, tl_bar
    integer :: i
    real :: M, sum_theta, sum_T, sum_q
    ! mass-weight BB average using layers 1..k_top-1 (interfaces variables)
    M = 0.0; sum_theta = 0.0; sum_T = 0.0; sum_q = 0.0
    do i = 1, max(1,k_top-1)
       M         = M + dpm(i)
       sum_theta = sum_theta + 0.5*(theta(i)+theta(i+1)) * dpm(i)
       sum_T     = sum_T     + 0.5*(T(i)    +T(i+1))     * dpm(i)
       sum_q     = sum_q     + 0.5*(qv(i)   +qv(i+1))    * dpm(i)
    end do
    if (M > 0.0) then
       T_bar  = sum_T / M
       qt_bar = max(qeps, sum_q / M)
       ! unsaturated BL ⇒ θ_l ≈ θ; approximate with θ from mean T at sfc p
       tl_bar = sum_theta / M
    else
       T_bar  = T(1)
       qt_bar = qv(1)
       tl_bar = theta(1)
    end if
  end subroutine bl_means

  subroutine dewpoint_from_qv_p(qv_in, p, TdK)
    ! convert mixing ratio + p to dewpoint (Magnus inversion)
    real, intent(in)  :: qv_in, p
    real, intent(out) :: TdK
    real :: eps, e, gamma
    eps = Rd/max(1.0e-12,Rv)
    e   = p * qv_in / max(1.0e-12, eps + qv_in)           ! Pa
    gamma = log( max(1.0e-12, e/611.2) )
    TdK = (243.5*gamma)/(17.67 - gamma) + 273.15
  end subroutine dewpoint_from_qv_p

  subroutine tlcl_bolton(TK, TdK, TLCL)
    real, intent(in)  :: TK, TdK
    real, intent(out) :: TLCL
    TLCL = 1.0 / ( 1.0/(TdK - 56.0) + log(TK/TdK)/800.0 ) + 56.0
    TLCL = max(150.0, min(340.0, TLCL))
  end subroutine tlcl_bolton

  subroutine lcl_pressure_from_T_qt(TLCL, qt, pLCL)
    ! at LCL: qsat(TLCL, pLCL) = qt ⇒ pLCL = es(TLCL)*(eps+qt)/qt
    real, intent(in)  :: TLCL, qt
    real, intent(out) :: pLCL
    real :: eps, es
    eps = Rd/max(1.0e-12,Rv)
    es  = 611.2 * exp( max(-60.0, min(60.0, 17.67*(TLCL-273.15)/(TLCL-29.65))) )
    pLCL = es * (eps + max(1.0e-12, qt)) / max(1.0e-12, qt)
    pLCL = max(pmin, min(pmax, pLCL))
  end subroutine lcl_pressure_from_T_qt

  subroutine interp_z_at_p(z, p, p_target, z_out, k_index)
    real, intent(in)  :: z(:), p(:), p_target
    real, intent(out) :: z_out
    integer, intent(out) :: k_index
    integer :: j, n
    n = size(p)
    ! assume p decreases with z
    if (p_target >= p(1)) then
       z_out = z(1); k_index = 1; return
    else if (p_target <= p(n)) then
       z_out = z(n); k_index = n; return
    end if
    do j = 1, n-1
       if ( (p(j) >= p_target) .and. (p_target >= p(j+1)) ) then
          z_out = z(j) + (z(j+1)-z(j)) * ( (p(j)-p_target)/max(1.0e-12, p(j)-p(j+1)) )
          k_index = j+1
          return
       end if
    end do
    z_out = z(n); k_index = n
  end subroutine interp_z_at_p

end subroutine shallow_conv_adj

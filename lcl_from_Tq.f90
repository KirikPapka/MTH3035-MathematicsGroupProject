subroutine lcl_from_Tq(z, p, theta_1, q_1, z_lcl, z_cb)
  implicit none
  INCLUDE 'params.inc'              ! expects: Rd, c_p, g, p0  (SI units)
  real, intent(in) :: theta_1       ! K (potential temperature)
  real, intent(in) :: q_1           ! mixing ratio r [kg/kg dry air]
  real, dimension(kkp), intent(in) :: z, p   ! z [m], p [Pa]
  real, intent(out) :: z_lcl, z_cb

  real :: RH, es_sfc, es_lcl, eps, Tlcl, p_lcl, T_sfc, kappa
  real :: r, q, r_s_lcl, Tv_sfc, Tv_lcl, Tv_bar, theta_v
  integer :: k, kcb

  kappa = Rd/c_p
  eps   = 0.622d0

  ! interpret input humidity as MIXING RATIO r
  r = q_1
  q = r / (1.d0 + r)                         ! specific humidity

  ! theta -> T at surface pressure p(1)
  T_sfc = theta_1 * (p(1)/p0)**kappa

  ! surface saturation vapor pressure (Pa)
  es_sfc = es_pa(T_sfc)

  ! RH from mixing ratio: e = p r / (eps + r)
  RH = (p(1) * r / (eps + r)) / max(1.d-6, es_sfc)
  RH = max(1.d-6, min(0.999999d0, RH))

  ! Bolton-like LCL temperature (K)
  Tlcl = 1.d0 / ( (1.d0/(T_sfc - 55.d0)) - log(RH)/2840.d0 ) + 55.d0

  ! LCL pressure via Poisson (cp/Rd exponent)
  p_lcl = p(1) * (Tlcl/T_sfc)**(c_p/Rd)

  ! Virtual temperatures for hypsometric
  Tv_sfc = T_sfc * (1.d0 + r/eps) / (1.d0 + r)

  es_lcl  = es_pa(Tlcl)
  r_s_lcl = eps * es_lcl / max(1.d0, (p_lcl - es_lcl))   ! saturated r at LCL
  Tv_lcl  = Tlcl * (1.d0 + r_s_lcl/eps) / (1.d0 + r_s_lcl)

  Tv_bar = 0.5d0 * (Tv_sfc + Tv_lcl)

  ! LCL height (m)
  z_lcl = max(0.d0, (Rd*Tv_bar)/g * log( p(1) / max(1.d0, p_lcl) ))

  ! optional: virtual potential temperature at surface
  theta_v = theta_1 * (1.d0 + r/eps) / (1.d0 + r)
  ! (return via argument if you need it)

  ! first model level with p(k) <= p_lcl
  kcb = 0
  do k = 1, kkp
     if (p(k) <= p_lcl) then
        kcb = k
        z_cb = z(kcb)
        exit
     end if
  end do
  if (kcb == 0) z_cb = 0.0   ! consider a missing value instead

contains
  real function es_pa(T) result(es)
    real, intent(in) :: T
    ! Buck over water, in Pa
    es = 611.2d0 * exp( 17.67d0*(T - 273.15d0) / (T - 29.65d0) )
  end function es_pa
end subroutine lcl_from_Tq

subroutine rh_z(p, theta, q_v, rh)
  implicit none
  INCLUDE 'params.inc'              ! expects: Rd, c_p, g, p0  (SI units)
  real, dimension(kkp), intent(in) :: theta      ! K (potential temperature)
  real, dimension(kkp), intent(in) :: q_v        ! mixing ratio r [kg/kg dry air]
  real, dimension(kkp), intent(in) :: p          ! p [Pa]
  real, dimension(kkp), intent(out) :: rh

  real :: e_sat, eps, T_env, kappa
  real :: theta_v
  integer :: k

  kappa = Rd/c_p
  eps   = 0.622d0

  ! interpret input humidity as MIXING RATIO r
  
   do k = 1, kkp

    ! Temperature
    T_env = theta(k) * (p(k)/p0)**kappa

   ! vapor pressure (Pa)
    e_sat = es_pa(T_env)

   ! RH from mixing ratio: e = p r / (eps + r)
   rh(k) = (p(k) * q_v(k) / (eps + q_v(k))) / max(1.d-6, e_sat)
   rh(k) = max(1.d-6, min(1.00000d0, rh(k)))

  ! optional: virtual potential temperature at surface
  theta_v = theta(k) * (1.d0 + 0.61*q_v(k))
  
  end do

contains
  real function es_pa(T) result(es)
    real, intent(in) :: T
    ! Buck over water, in Pa
    es = 611.2d0 * exp( 17.67d0*(T - 273.15d0) / (T - 29.65d0) )
  end function es_pa
  
end subroutine rh_z

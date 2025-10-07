 subroutine pressure_from_theta(z, theta, p)
    implicit none
    INCLUDE 'params.inc'
    
    REAL, DIMENSION(kkp) , INTENT(in) :: z, theta
    REAL, DIMENSION(kkp) , INTENT(out) :: p

    integer :: k
    real :: kappa, Pi, dPi, dz, inv_theta_k, inv_theta_kp1

    if (kkp < 1) return

    kappa = Rd / c_p

    ! Exner at surface
    Pi = max(1.0d-12, (psfc / p0)**kappa)
    p(1) = p0 * Pi**(1.0d0/kappa)

    do k = 1, kkp-1
       dz = max(1.0d0, z(k+1) - z(k))
       inv_theta_k   = 1.0d0 / max(1.0d-6, theta(k))
       inv_theta_kp1 = 1.0d0 / max(1.0d-6, theta(k+1))

       ! Trapezoidal step for Exner: dPi/dz = - g / (cp * theta)
       dPi = - (g / c_p) * 0.5d0 * (inv_theta_k + inv_theta_kp1) * dz
       Pi  = max(1.0d-12, Pi + dPi)

       p(k+1) = p0 * Pi**(1.0d0/kappa)
    end do
  end subroutine pressure_from_theta
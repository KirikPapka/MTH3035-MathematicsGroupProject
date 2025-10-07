subroutine set_fluxes_l(t, h, q)
  implicit none
  real, intent(in)  :: t            ! seconds since start
  real, intent(out) :: h, q         ! SHF (W/m^2), moisture flux (kg/m^2/s)

  integer :: i, n
  real :: tt, Lv, period

  ! 
  real :: shf(14), lhf(14), times(14)

  shf = (/-10.0,  90.0, 140.0, 140.0, 100.0,  1.0, -10.0,  &
           -10.0,  90.0, 140.0, 140.0, 100.0,  1.0, -10.0 /)

  lhf = (/  5.0, 250.0, 450.0, 500.0, 420.0, 180.0,   0.0,  &
             5.0, 250.0, 450.0, 500.0, 420.0, 180.0,   0.0 /)

  times = (/    0.0, 14400.0, 23400.0, 27000.0, 36000.0, 45000.0, 52200.0,  &
            86400.0,100800.0,109800.0,113400.0,122400.0,131400.0,138600.0 /)

  ! --- CONSTANTS ---
  Lv     = 2.5e6         ! J/kg
  period = 86400.0       ! seconds in a day
  n      = size(times)

 
  if (t < times(1)) then
     tt = times(1)
  else
     ! Modulo relative to times(1) avoids tiny numerical gaps at 0.0
     tt = mod(t - times(1), period) + times(1)
     if (tt > times(n)) tt = times(n)
  end if

  ! ---- Piecewise-linear interpolation on [times(i), times(i+1)) ----
  if (tt >= times(n)) then
     h = shf(n)
     q = lhf(n)/Lv
     return
  end if

  do i = 1, n-1
     if ( times(i) <= tt .and. tt < times(i+1) ) then
        h = shf(i) + (shf(i+1)-shf(i)) * (tt-times(i)) / (times(i+1)-times(i))
        q = ( lhf(i) + (lhf(i+1)-lhf(i)) * (tt-times(i)) / (times(i+1)-times(i)) ) / Lv
        return
     end if
  end do

  ! Fallback (shouldn't hit if times() is strictly increasing)
  h = shf(n)
  q = lhf(n)/Lv
end subroutine set_fluxes_l

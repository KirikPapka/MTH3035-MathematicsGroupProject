      subroutine init_qv(zn, qv)
     
      IMPLICIT none
      INCLUDE 'params.inc'
      ! Inputs
      REAL,  INTENT(IN), DIMENSION(kkp) :: zn
      INTEGER :: k, i
      
      real :: z_lev(27)
      real :: qv_n(27)
     
      ! Outputs a linear piece-wise profile of qv
      REAL, INTENT(OUT), DIMENSION(kkp) :: qv
      
     z_lev = (/0., 15., 48., 145., 181., 258., 550., 654., 712., 717., 780., &
                824., 1036., &
                1207., 1364., 1401., 1469., 2401., 2599., 3150.,&
                3500., 3599., 3690., 3736., 3784., 3840., 4001./)
     ! qv_n  =  (/0.011, 0.010, 0.0095, 0.009, 0.0088, 0.006, 0.003, 0.001/)
      qv_n = (/3.2e-3, 4.11e-3, 4.33e-3, 4.37e-3, 4.43e-3, 4.6e-3, 1.9e-3, 0.63e-3,&
      0.64e-3, 0.26e-3,0.51e-3, 0.52e-3, 0.73e-3, 0.78e-3, 0.94e-3, 1.61e-3, 1.73e-3, &
      1.93e-3, 1.49e-3, 1.45e-3, 1.27e-3, 1.12e-3, 1.49e-3, 1.49e-3, 0.73e-3, 0.63e-3, 0.63e-3 /)
      
      do k=1, kkp
       do i=1,size(z_lev)-1
         if (z_lev(i) <= zn(k) .and. zn(k) < z_lev(i+1)) then
           qv(k) = qv_n(i)+(qv_n(i+1)-qv_n(i)) * (zn(k)-z_lev(i))/(z_lev(i+1)-z_lev(i))
         end if
       end do
            
       if (zn(k) >=z_lev(size(z_lev))) then
       qv(k) = qv_n(size(z_lev))
       end if
       
      end do
      
      do k=1, kkp
      qv(k) = qv(k)/(1-qv(k))
      end do

      end
  

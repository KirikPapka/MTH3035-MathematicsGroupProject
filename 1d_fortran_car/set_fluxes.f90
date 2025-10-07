      subroutine set_fluxes(t,h,q)
     
      IMPLICIT none
      INCLUDE 'params.inc'
      ! Inputs
      REAL, INTENT(IN) :: t
      INTEGER :: k, i
      
      real :: times(26)
      real :: shf(26)
      real :: lhf(26)
     
      ! Outputs fluxes
      REAL, INTENT(OUT) :: h ,q
      
       times = (/0., 1800., 3600., 5400., 7200., 9000., 10800., 12600., 14400., 16200.,& 
       18000., 19800., 21600., 23400., 25200., 27000., 28800., 30600., 32400., &
       34200., 36000., 37800., 39600., 41400., 43200., 45000./)

       shf = (/-1.68, -5.70, 0.309,0.508, 13.0566, 34.7814, 49.9949, 77.3823, 118.247, &
               126.543, 126.05,&
                 144.196, 144.631, 138.641, 159.381, 159.913, &
                173.733, 144.033, 130.996, 103.84, 93.1718, &
                 52.5289, 51.4246, 23.937, 11.7942, -6.83469 /)
       lhf = (/-0.37, -1.48, 0.15, -0.06, 9.53, 20.8192, 41.2206, &
                          45.9119, 87.0906, 92.8874, 94.4597, &
                        97.5644, 124.068, 126.202, 125.416, 137.699, &
                         166.794, 163.311, 149.353, 148.503, 126.988, &
                     67.6469, 96.78, 68.4137, 90.3206, 40.8392 /)

      ! Surface fluxes (W / m^2)
       do i=1,size(times)-1
         if (times(i) <= t .and. t < times(i+1)) then
           h = shf(i)+(shf(i+1)-shf(i)) * (t-times(i))/(times(i+1)-times(i))
           q = ( lhf(i)+(lhf(i+1)-lhf(i)) * (t-times(i))/(times(i+1)-times(i)) )/2.5e6
         end if
       end do
       
       if (t >=times(size(times))) then
       h = shf(size(times))
       q = lhf(size(times))/2.6e6
       end if

      end
  

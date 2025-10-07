      SUBROUTINE init_wind(u, v, zn, u_g, v_g)
      IMPLICIT none

      ! Parameters
      INCLUDE 'params.inc'

      ! Initial wind
      REAL, DIMENSION(kkp) :: u, v
      REAL, DIMENSION(kkp) :: zn

      INTEGER :: k

      ! Parameter in jet profile
      REAL :: h_1
      REAL :: pi
      pi=4.0*atan(1.0)      

      ! Option for wind initialisation wind_opt
      ! 0 = geostrophic initialisation
      ! 1 = exponential in u, 0 in v
      ! 2 = Beare '06 jet
      ! 3 = Ekman spiral

      ! Read in namelist
      OPEN(5,file='nml_1d_bl')
      READ(5,nml=wind)
      CLOSE(5)      

      SELECT case(wind_opt)  
      CASE(0)
        PRINT*,'Geostrophic wind initialisation'
        ! geostrophic initialisation .
        DO k=1, kkp
          u(k)=u_g
          v(k)=0.0
        ENDDO
      CASE(1)
        ! Exponential u
        PRINT*,'Exponential wind initialisation'
        DO k=1, kkp
          u(k)=u_g*(1-exp(-zn(k)/h_m))
          v(k)=v_g
        ENDDO 
      CASE(2)
        ! Beare '06 jet
        h_1=h_beare*(1+sqrt(cc))
        PRINT*,'Beare 06 jet initialisation'
        DO k=1, kkp           
          IF (zn(k)<= h_1) THEN
            u(k)=u_g + u_g*(cc-(1-zn(k)/h_beare)**2)
          ELSE
            u(k)=u_g
          ENDIF
  
          IF (zn(k)<= h_2) THEN 
            v(k)=u_g*dd*exp(-zn(k)/h_beare)*sin(2*pi*zn(k)/h_2)
          ELSE
            v(k)=0.0
          ENDIF     
        ENDDO
      CASE(3)
        ! Ekman spiral
        DO k=1, kkp  
          u(k)=u_g*( 1.0 - exp(-pi*zn(k)/hh)*cos(pi*zn(k)/hh) )
          v(k)=u_g*exp(-pi*zn(k)/hh)*sin(pi*zn(k)/hh)
        ENDDO
      END SELECT ! wind_opt



      END

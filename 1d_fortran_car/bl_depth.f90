      SUBROUTINE bl_depth(u, v, th, th_surf, zn, z, u_star, &
                          L_obuk, wth_0, ri, h_bl)

      ! Routine to calculate the boundary layer depth 
      IMPLICIT none
      INCLUDE 'params.inc'

      ! Declare Inputs
      REAL, DIMENSION(kkp) :: u, v, th, zn, z, ri
      REAL :: th_surf
      REAL :: u_star, L_obuk, wth_0
      ! Declare Outputs
      REAL :: h_bl
 
      ! Declare locals
      INTEGER :: k
      ! Variables for calculating bl depth
      REAL :: min_diff, diff, h_temp
      REAL :: delta_th, rib
      REAL :: h_cbl, h_stablebl
      h_bl=0.0

! h_bl is on theta levels
!      PRINT*,'th_surf: ',th_surf

      ! Bulk Richardson no. method.
!      DO k=k_sl,kkp
!        h_temp = ri_cr * th_ref * ( u(k)**2 + v(k)**2 )      &
!                 / ( g * (th(k) - th_surf) + 1e-9 )
!        diff = abs(h_temp - zn(k))
!        IF (diff < min_diff) THEN
!          h_bl=zn(k)
!          min_diff=diff
!        ENDIF
!      ENDDO

      ! Bulk Richardson no. method.
      ! when rib > ri_cr
      !   ------------ z(k)
      !   ------------ zn(k)
      !   ------------ z(k-1)
      ! Calculate h on z levels as this will be used 
      ! for calculating k_m, k_h and later explicit
      ! entrainment

!      h_bl=z(2)
      ! Bulk Ri diagnosis
!!$      DO k=2,kkp
!!$
!!$        rib= (g/th_ref) * zn(k) * (th(k) - th_surf)     &
!!$              / ( u(k)**2 + v(k)**2 + 1e-9)
!!$       h_bl=z(k-1)
!!$         IF (rib > ri_cr) EXIT
!!$ 
!!$      ENDDO

!!$
      
      IF (wth_0 >= 0) THEN
        min_diff=1e8
        ! Parcel ascent for CBL
        DO k=k_sl,kkp-1
          ! minimise delta th
          delta_th=abs(th(k)-th_surf)
          IF (delta_th < min_diff) THEN
            h_bl=zn(k+1)
            min_diff=delta_th
          ENDIF
        ENDDO
      ELSE 
!        ! Zilitinkevich 72 for SBL
!!$        h_bl=c_zil*sqrt(u_star*L_obuk/f_0)
!!$
!!$        ! Parcel ascent for CBL
!!$        min_diff=0.1
!!$        h_cbl=z(2)
!!$        DO k=k_sl,kkp
!!$          ! minimise delta th
!!$          
!!$          delta_th=abs(th(k)-th_surf)
!!$          IF (delta_th < min_diff) THEN
!!$            h_cbl=z(k-1)
!!$            min_diff=delta_th
!!$          ENDIF
!!$        ENDDO
!!$         ! Gradient Richardson
        DO k=2, kkp-1
          h_bl=zn(k+1)
          IF (ri(k) > ri_cr) EXIT

        ENDDO  

!!$        h_bl=MAX(h_stablebl,h_cbl)       
      ENDIF

!      PRINT*,'h_bl: ',h_bl
      END





      

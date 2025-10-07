      SUBROUTINE vert_levs(z,zn,z_wt)
 
      ! Subroutine to calculate vertical levels
      ! according to the Weng and Taylor '03 distribution
      ! big_z= ln((z+z_0)/z_0) + z/b_0
      ! where big_z is uniformly distributed
      ! but z has many more points in the log layer.
      ! Distributing levels like this 
      ! means that the  you don't need to
      ! apply bulk surface exchange.

      IMPLICIT none        
      INCLUDE 'params.inc' 

      REAL, INTENT(INOUT), DIMENSION(kkp) :: z, zn

      REAL, INTENT(INOUT), DIMENSION(kkp-1) :: z_wt
    

      ! locals are below

      ! Upper case Z variable in eq. 33 in Weng and Taylor
      REAL, DIMENSION(kkp) :: big_z  

      ! Maximum value of big_z
      REAL :: max_big_z

      ! Parameter in Weng and Taylor
!      REAL, PARAMETER :: b_0=67.5

      REAL, DIMENSION(n_fine) :: z_fine
      
      ! Variables used in solution
      REAL:: diff, min_diff, big_z_temp, dz_grid      

      ! Loop integers
      INTEGER :: n, k
      
      if (stretch_grid)then 

        max_big_z = alog((zztop+b_grid)/b_grid) + zztop/a_grid

       ! Define the fine grid
       DO n=1, n_fine
         z_fine(n)=zztop*REAL(n-1)/REAL(n_fine-1)
       ENDDO
        
       ! Define upper case z
       ! evenly spaced up to max_big_z
       DO k=1, kkp 
         big_z(k)=max_big_z*REAL(k-1)/REAL(kkp-1)
       ENDDO
     
       DO k=1, kkp 
         min_diff=1.0e6
         DO n=1, n_fine

           big_z_temp=alog((z_fine(n)+b_grid)/b_grid) + z_fine(n)/a_grid
           diff=ABS(big_z(k)-big_z_temp)

           IF (diff < min_diff) THEN
             min_diff=diff
             z(k)=z_fine(n)
           ENDIF
 
         ENDDO
       ENDDO
       
       ELSE
       
       dz_grid = zztop  / (kkp - 1)
       
       ! Define z points (w - theta levels)
      DO k = 1, kkp
         z(k) = (k-1) * dz_grid
      END DO
       
       ENDIF

      ! Calculate zn from z
      DO k=2,kkp
        zn(k)=(z(k)+z(k-1))/2.0
      ENDDO
      zn(1)=-zn(2) 

      ! Calculate z_wt
      DO k=1, kkp-1
        z_wt(k)=zn(k+1)
      ENDDO



!!$      PRINT*,'k  z(k)  zn(k)'
!!$      DO k=1,kkp
!!$        PRINT*,k,z(k),zn(k)
!!$      ENDDO
!!$      PRINT*,k,'zn(k)-zn(k-1)'
!!$      DO k=2,kkp
!!$        PRINT*,k,zn(k)-zn(k-1)
!!$      ENDDO
      END

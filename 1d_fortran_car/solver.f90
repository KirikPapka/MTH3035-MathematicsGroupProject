     PROGRAM main
     INTEGER n,i
     REAL, DIMENSION(20) :: a,b,c,r,u
     n=10

     DO i=1,n
       a(i)=-1.0
       b(i)=2.0
       c(i)=-1.0
       r(i)=0.0
     ENDDO
     r(1)=1
     CALL tridag(a,b,c,r,u,n)
     DO i=1,n

       PRINT*,i,u(i)
     ENDDO

     END
     INCLUDE 'tridag.f90'

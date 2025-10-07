      SUBROUTINE tridag(a,b,c,r,u,n)
      ! Tri-diagonal solver from numerical recipes.
      INTEGER:: n
      REAL,DIMENSION(n) :: a,b,c,r,u,gam
      INTEGER:: j
      REAL:: bet
      if(b(1).eq.0.)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
      enddo
      do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
      enddo
      return
      END

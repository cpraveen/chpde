c
c
c =========================================================
      subroutine ode1dif(maxmx,meqn,mbc,mx,q,t,dt)
c =========================================================
      implicit real*8(a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
c
      parameter (maxmode = 502)
      dimension x(-1:maxmode)
      common /comeps/  eps
      common /grid/  x
      common /combc/ g0,g1  !# boundary values at end of step
c
      parameter (msize = 500)
      dimension c(msize),d(msize),e(msize),b(msize)
c
      dx = x(2) - x(1)
      r = 0.5d0 * eps * dt / dx**2
      gamma1 = 1.d0 + 2.d0*r
      gamma2 = 1.d0 - 2.d0*r
      do 30 i=1,mx
	 c(i) = -r
	 d(i) = gamma1
	 e(i) = -r
	 b(i) = r*q(i-1,1) + gamma2*q(i,1) + r*q(i+1,1)
   30    continue
c
c     # modify RHS using boundary conditions:
      g0 = q(0,1)
      g1 = q(mx,1)
      b(1) = b(1) + r*g0
      b(mx) = b(mx) + r*g1
c
  40  continue
      call dgtsl(mx,c,d,e,b,info)
      if (info.ne.0) write(6,*) 'info =',info
c
      do 50 i=1,mx
	 q(i,1) = b(i)
   50    continue
c
      return 
      end

c
c    ================================
      subroutine dgtsl(n,c,d,e,b,info)
c    ================================
      integer n,info
      double precision c(1),d(1),e(1),b(1)
c
c     from LINPACK
c     dgtsl given a general tridiagonal matrix and a right hand
c     side will find the solution.
c
c     on entry
c
c        n       integer
c                is the order of the tridiagonal matrix.
c
c        c       double precision(n)
c                is the subdiagonal of the tridiagonal matrix.
c                c(2) through c(n) should contain the subdiagonal.
c                on output c is destroyed.
c
c        d       double precision(n)
c                is the diagonal of the tridiagonal matrix.
c                on output d is destroyed.
c
c        e       double precision(n)
c                is the superdiagonal of the tridiagonal matrix.
c                e(1) through e(n-1) should contain the superdiagonal.
c                on output e is destroyed.
c
c        b       double precision(n)
c                is the right hand side vector.
c
c     on return
c
c        b       is the solution vector.
c
c        info    integer
c                = 0 normal value.
c                = k if the k-th element of the diagonal becomes
c                    exactly zero.  the subroutine returns when
c                    this is detected.
c
c     linpack. this version dated 08/14/78 .
c     jack dongarra, argonne national laboratory.
c
c     no externals
c     fortran dabs
c
c     internal variables
c
      integer k,kb,kp1,nm1,nm2
      double precision t
c     begin block permitting ...exits to 100
c
         info = 0
         c(1) = d(1)
         nm1 = n - 1
         if (nm1 .lt. 1) go to 40
            d(1) = e(1)
            e(1) = 0.0d0
            e(n) = 0.0d0
c
            do 30 k = 1, nm1
               kp1 = k + 1
c
c              find the largest of the two rows
c
               if (dabs(c(kp1)) .lt. dabs(c(k))) go to 10
c
c                 interchange row
c
                  t = c(kp1)
                  c(kp1) = c(k)
                  c(k) = t
                  t = d(kp1)
                  d(kp1) = d(k)
                  d(k) = t
                  t = e(kp1)
                  e(kp1) = e(k)
                  e(k) = t
                  t = b(kp1)
                  b(kp1) = b(k)
                  b(k) = t
   10          continue
c
c              zero elements
c
               if (c(k) .ne. 0.0d0) go to 20
                  info = k
c     ............exit
                  go to 100
   20          continue
               t = -c(kp1)/c(k)
               c(kp1) = d(kp1) + t*d(k)
               d(kp1) = e(kp1) + t*e(k)
               e(kp1) = 0.0d0
               b(kp1) = b(kp1) + t*b(k)
   30       continue
   40    continue
         if (c(n) .ne. 0.0d0) go to 50
            info = n
         go to 90
   50    continue
c
c           back solve
c
            nm2 = n - 2
            b(n) = b(n)/c(n)
            if (n .eq. 1) go to 80
               b(nm1) = (b(nm1) - d(nm1)*b(n))/c(nm1)
               if (nm2 .lt. 1) go to 70
               do 60 kb = 1, nm2
                  k = nm2 - kb + 1
                  b(k) = (b(k) - d(k)*b(k+1) - e(k)*b(k+2))/c(k)
   60          continue
   70          continue
   80       continue
   90    continue
  100 continue
c
      return
      end

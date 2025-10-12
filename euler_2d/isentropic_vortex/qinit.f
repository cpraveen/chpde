
c
c
c
c     =====================================================
       subroutine qinit(meqn,mbc,mx,my,xlower,ylower,
     &                   dx,dy,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
c      # Data is piecewise constant with 4 values in 4 quadrants
c      # 2D Riemann problem from Figure 4 of
c        @article{csr-col-glaz,
c          author="C. W. Schulz-Rinne and J. P. Collins and H. M. Glaz",
c          title="Numerical Solution of the {R}iemann Problem for 
c                 Two-Dimensional Gas Dynamics",
c          journal="SIAM J. Sci. Comput.",
c          volume="14",
c          year="1993",
c          pages="1394-1414" }

c
       implicit double precision (a-h,o-z)
       dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
       dimension aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
       common /cparam/  gamma
       common /problem/ dmach, alpha, beta 
c
       gamma1 = gamma - 1.d0
       pi  = 4.0d0 * datan(1.0d0)
       vx0 = dmach * dcos(alpha * pi / 180.0d0)
       vy0 = dmach * dsin(alpha * pi / 180.0d0)
       a   = gamma1 * beta**2 / (8.0d0 * gamma * pi**2)
       b   = beta/(2.0d0*pi)
c
       do 15 i=1,mx
          xcell = xlower + (i-0.5d0)*dx
          do 15 j=1,my
             ycell    = ylower + (j-0.5d0)*dy
             r2       = xcell**2 + ycell**2
             temp     = 1.0d0 - a * dexp(1.0d0 - r2)
             rho      = temp**(1.0d0/gamma1)
             pre      = rho**gamma
             vx       = vx0 - b * ycell * dexp(0.5d0*(1.0d0-r2))
             vy       = vy0 + b * xcell * dexp(0.5d0*(1.0d0-r2))
             q(1,i,j) = rho
             q(2,i,j) = rho * vx
             q(3,i,j) = rho * vy
             q(4,i,j) = pre/gamma1 + 0.5d0*rho*(vx**2 + vy**2)
   15        continue
       return
       end

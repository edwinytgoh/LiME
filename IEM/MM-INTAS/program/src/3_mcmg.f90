!#1 subroutine mcmg

!--- tools
!#2 subroutine geterf1
!#3 subroutine inter1d
!#4 subroutine sortip


!#1
!=======================================================================
!       subroutine mcmg
!=======================================================================
!  function :
!    routine to advance scalars for one time step according to the
!    mapping closure models for which the evolution of the scalar pdf
!    is calculated by a time dependent mapping of a standard Gaussian
!    Pope (1991), Valino, Ros & Dopazo (1991).
!
!       - Pope (1991). Mapping closures for turbulent mixing and
!         reaction.  Theoret. Comput. Fluid Dynamics, 2, 255-270.
!
!       - Valino, L., Ros, J. and Dopazo, C. (1991). Monte Carlo
!         implementation and analytic solution of an inert-scalar
!         turbulent-mixing test problem using a mapping closure.
!         Physics of Fluids A, 3, 2191-2198.
!
!  input:
!    dt            : timestep
!    cphi          : constant in IEM model
!    ncell         : number of cells in the domain
!    omean_cell(l) : mean turbulent frequency in cell l
!    np            : number of particles
!    npd           : first dimension size of array f
!    nsv           : number of independent scalars
!    nf(l)         : first particle in cell l
!    nl(l)         : last particle in cell l
!    f_kwt(np)     : particle weights
!
!  input/output:
!    f_kf(npd,nsv) : particle independent scalars
!
!  includes:
!    param
!    logunits.cmn
!
!  calls:
!    sortip
!    inter1d
!    geterf1
!-----------------------------------------------------------------------
      SUBROUTINE mcmg(loutf, &
                      dt, cphi, &
                      ncell, &
                      omean_cell, &
                      np, npd, nsv, nf, nl, &
                      f_kwt, f_kf)

!-----  used modules
      USE cmn_param, ONLY : zero, half, one
      USE blas_lapack_module, ONLY : dpbsv

      IMPLICIT NONE

!=====  VARIABLES

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: ncell
      INTEGER (KIND=4), INTENT(IN) :: np, npd, nsv
      INTEGER (KIND=4), INTENT(IN) :: nf(ncell)
      INTEGER (KIND=4), INTENT(IN) :: nl(ncell)
      REAL    (KIND=8), INTENT(IN) :: dt, cphi
      REAL    (KIND=8), INTENT(IN) :: omean_cell(ncell)
      REAL    (KIND=8), INTENT(IN) :: f_kwt(np)

!-----  input/output variables
      REAL    (KIND=8), INTENT(INOUT) :: f_kf(npd,nsv)

!-----  local variables
      INTEGER (KIND=4) :: i, ig, igadd, k, l, nfl, nll, npt
      REAL    (KIND=8) :: hcphomdt, inc_wt, omega, sumwt, &
                          fold, fnew, ffold, ffnew, &
                          etai, etaph, etap, gph, fac

      INTEGER (KIND=4), ALLOCATABLE :: isp(:)                 ! (npt)
      REAL    (KIND=8), ALLOCATABLE :: cdf(:)                 ! (npt+1)
      REAL    (KIND=8), ALLOCATABLE :: eta(:)                 ! (npt+1)
      REAL    (KIND=8), ALLOCATABLE :: phi(:)                 ! (npt)
      REAL    (KIND=8), ALLOCATABLE :: bph(:)                 ! (npt-1)

!-----  variables needed for call of LAPACK routine DPBSV
      INTEGER (KIND=4), PARAMETER :: kd = 1, nrhs = 1, ldab = kd+1

      INTEGER (KIND=4) :: info
      REAL    (KIND=8), ALLOCATABLE :: ab(:,:)                ! (ldab,npt)
      CHARACTER(LEN=1) :: uplo

      REAL    (KIND=8), PARAMETER :: rtpi = 0.39894228d0
      REAL    (KIND=8), PARAMETER :: phimin = 0.d0, phimax = 1.d0
      REAL    (KIND=8), PARAMETER :: tiny = 1.E-30, small = 1.E-10


!=====  MAIN ACTION

      !--- set variables needed for DSBEV
      uplo = "U"


      !--- loop over cells
      CELL_LOOP: &
      DO l = 1, ncell

         !--- number of particles in this cell
         nfl = nf(l)
         nll = nl(l)
         npt = nll - nfl + 1
         IF (npt <= 1) CYCLE CELL_LOOP

         !--- mean scalar dissipation frequency in this cell
         omega = omean_cell(l)
         hcphomdt = half * cphi * omega * dt

         !--- total weight in this cell
         sumwt = zero
         DO i = nfl, nll
            sumwt = sumwt + f_kwt(i)
         END DO
         IF (sumwt <= tiny) CYCLE CELL_LOOP

         ALLOCATE(isp(npt))
         ALLOCATE(cdf(npt+1))
         ALLOCATE(eta(npt+1))
         ALLOCATE(phi(npt))
         ALLOCATE(bph(npt-1))

         ALLOCATE(ab(ldab,npt))

         !--- loop over independent scalars
         SCALAR_LOOP: &
         DO k = 1, nsv

            !--- calculate f mean
            fold  = zero
            ffold = zero
            DO i = nfl,nll
               fold  = fold  + f_kwt(i) * f_kf(i,k)
               ffold = ffold + f_kwt(i) * f_kf(i,k)**2
            END DO
            fold  = fold / sumwt
            ffold = (ffold / sumwt) - fold**2

            !--- skip mixing if variance is very small
            IF (ffold < small) THEN
               CYCLE SCALAR_LOOP
            END IF

            !--- order particles on scalar value
            igadd = nfl - 1
            DO i = 1, npt
               ig     = i + igadd
               isp(i) = ig
               phi(i) = f_kf(ig,k)
            END DO
            CALL sortip(npt, phi, isp)

            !--- probabilities
            inc_wt = zero
            DO i = 1, npt-1
               ig     = isp(i)
               inc_wt = inc_wt + f_kwt(ig)
               cdf(i) = inc_wt / sumwt
            END DO

            !--- to construct eta(1)
            cdf(npt) = half * cdf(1)

            !--- to construct eta(npt)
            cdf(npt+1) = half * (cdf(npt-1) + one)

            !--- mapped particles eta(1.5), eta(2.5),.... eta(1), eta(npt)
            CALL geterf1(npt+1, cdf , eta)


            !--- construct matrix coefficients

            !--- B_1+1/2
            etai   = eta(npt)
            etaph  = eta(1)
            etap   = half * (eta(1) + eta(2))
            gph    = rtpi * EXP(-half * etaph**2)
            bph(1) = DBLE(npt) * gph / (etap - etai)

            !--- B_i+1/2
            DO i = 2, npt-2
               etai   = half * (eta(i-1) + eta(i  ))
               etaph  = eta(i)
               etap   = half * (eta(i  ) + eta(i+1))
               gph    = rtpi * EXP(-half * etaph**2)
               bph(i) = DBLE(npt) * gph / (etap - etai)
            END DO

            !--- B_npt-1/2
            etai       = half * (eta(npt-2) + eta(npt-1))
            etaph      = eta(npt-1)
            etap       = eta(npt+1)
            gph        = rtpi * EXP(-half * etaph**2)
            bph(npt-1) = DBLE(npt) * gph / (etap - etai)

            !--- construct matrix (I-dt*A)
            ab(2,1) = one + bph(1)*hcphomdt
            DO i = 2, npt-1
               ab(2,i) = one + (bph(i-1) + bph(i))*hcphomdt
            END DO
            ab(2,npt) = one + bph(npt-1)*hcphomdt
 
            DO i = 1, npt-1
               ab(1,i+1) = -bph(i)*hcphomdt
            END DO

            !--- solve system
            CALL dpbsv(uplo, npt, kd, nrhs, ab, ldab, phi, npt, info)
            IF (info /= 0) THEN
               WRITE (loutf,*) "Error in subroutine MCMG"
               WRITE (loutf,*) "On return from DPBSV: INFO = ", info
               STOP 'MCMG: see output file.'
            END IF

            !--- new scalar values
            DO i = 1, npt
               f_kf(isp(i),k) = phi(i)
            END DO

            !--- calculate f mean
            fnew  = zero
            ffnew = zero
            DO i = nfl, nll
               fnew  = fnew  + f_kwt(i) * f_kf(i,k)
               ffnew = ffnew + f_kwt(i) * f_kf(i,k)**2
            END DO
            fnew  = fnew / sumwt
            ffnew = (ffnew / sumwt) - fnew**2

            !--- perform correction
            fac = SQRT(ffold/(ffnew + tiny)) * EXP(-hcphomdt)
            DO i = nfl, nll
               f_kf(i,k) = fold + (f_kf(i,k) - fnew) * fac
               f_kf(i,k) = MAX(MIN(f_kf(i,k),phimax),phimin)
            END DO
 
         END DO SCALAR_LOOP

         DEALLOCATE(isp)
         DEALLOCATE(cdf)
         DEALLOCATE(eta)
         DEALLOCATE(phi)
         DEALLOCATE(bph)

         DEALLOCATE(ab)

      END DO CELL_LOOP

      END SUBROUTINE mcmg


!#2
!=======================================================================
!  subroutine geterf1( n, x, y )
!=======================================================================
!  method:
!    approximation to the inverse normal cumulative density function
!    according to the approximation by Abramowitz & Stegun formula
!    26.2.23 approximation error < 4.5E-04
!
!  limits :
!    x           [zero,one]
!    y           [-big,big] , y = -big for x < small
!                             y =  big for x > one - small
!
!    small       [2.8D-17]    chosen smallest number possible
!    big         [1000   ]    chosen arbitrary large number one order
!                             of magnitude larger than
!                             ERF1[ONE-SMALL+E]
!
!  input :
!    n           : number of function values
!    x           : probability
!
!  output :
!    f           : inverse error function of y
!----------------------------------------------------------------------
      SUBROUTINE geterf1(n, x, y)

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: n
      REAL    (KIND=8), INTENT(IN) :: x(n)                    !(n)

!----- output variables
      REAL    (KIND=8), INTENT(OUT) :: y(n)                   !(n)

!----- local variables
      INTEGER (KIND=4) :: i
      REAL    (KIND=8) :: xi, x1, x2, x3, yi, y1, t, t1, t2, f, f1, f2, s

      REAL    (KIND=8), PARAMETER :: zero = 0.d0, &
                                     half = 0.5d0, &
                                     one = 1.d0, &
                                     two = 2.d0, &
                                     big = 1.D+03, &
                                     small = 2.8D-17

      REAL    (KIND=8), PARAMETER :: c0 = 2.515517d0, &
                                     c1 = 0.802853d0, &
                                     c2 = 0.010328d0, &
                                     d1 = 1.432788d0, &
                                     d2 = 0.189269d0, &
                                     d3 = 0.001308d0

!-----  calculate inverse error function
      DO i = 1,n

!--------  endpoints clamped to -BIG, BIG
         IF     (x(i).LT.(zero + small)) THEN
            y(i) = -big
         ELSE IF (x(i).GT.(one  - small)) THEN
            y(i) =  big

!--------  normal calulation
         ELSE
            x1   = x(i)
            x2   = half - x1
            x3   = ABS(x2)
            xi   = half - x3
            s    = -SIGN(one, x2)

            t1   = LOG(xi)
            t2   = -two * t1
            t    = SQRT(t2)

            f1   = c0  + c1 * t + c2 * t2
            f2   = one + d1 * t + d2 * t2 + d3 * t*t2
            f    = f1 / f2

            y1   = t - f
            yi   = s * y1
            y(i) = yi
         END IF
      END DO

      END SUBROUTINE geterf1


!#3
!=======================================================================
!  subroutine inter1d
!=======================================================================
!  function:
!    interpolates in an array of function values and arguments
!    it is assumed that the arguments are sorted
!    extra endpoints can be passed using XIN, XMAX, FMIN, FMAX
!    values outsside [XMIN,XMAX] are clamped to FMIN or FMAX
!
!  input:
!    nin           : number of function values
!    xin           : list of arguments
!    xmin          : minimum argument value
!    xmax          : maximum argument value
!    fin           : list of function values
!    fmin          : minimum function value
!    fmax          : maximum function value
!    nout          : number of output values
!    xout          : list of arguments for interpolation
!
!  output:
!    fout          : new interpolated values
!-----------------------------------------------------------------------
      SUBROUTINE inter1d (nin, xin, xmin, xmax, fin, fmin, fmax, &
                          nout, xout, fout)

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: nin, nout
      REAL    (KIND=8), INTENT(IN) :: xmin, xmax, fmin, fmax
      REAL    (KIND=8), INTENT(IN) :: xin(nin)                ! (nin)
      REAL    (KIND=8), INTENT(IN) :: fin(nin)                ! (nin)
      REAL    (KIND=8), INTENT(IN) :: xout(nout)              ! (nout)
                                     

!-----  output variables
      REAL    (KIND=8), INTENT(INOUT) :: fout(nout)           ! (nout)

!-----  local variables
      INTEGER (KIND=4) :: i, j
      REAL    (KIND=8) :: f1, f2, fac1, fac2, x1, x2, dx, xmm, xxx, xouti

!=====  MAIN ACTION

!-----  limits
      xmm = xin(1  )
      xxx = xin(nin)

!-----  loop over points
      j = 2
      DO i = 1,nout
         xouti = xout(i)
         IF (xouti.LT.xmin) THEN
            fout(i) = fmin
         ELSE IF (xouti.GT.xmax) THEN
            fout(i) = fmax
         ELSE
            IF (xouti.LT.xmm) THEN
               x1 = xmin
               x2 = xmm
               f1 = fmin
               f2 = fin(1)
            ELSE IF (xouti.GT.xxx) THEN
               x1 = xxx
               x2 = xmax
               f1 = fin(nin)
               f2 = fmax
            ELSE
               DO
                  IF (xin(j).GT.xouti) THEN
                     x1 = xin(j-1)
                     x2 = xin(j  )
                     f1 = fin(j-1)
                     f2 = fin(j  )
                     EXIT
                  ELSE
                     j = j + 1
                  END IF
               END DO
               dx      =  x2    - x1
               fac1    = (x2    - xouti) / dx
               fac2    = (xouti - x1   ) / dx
               fout(i) = fac1 * f1 + fac2 * f2
            END IF
         END IF
      END DO
!-----  end of loop over points

      END SUBROUTINE inter1d


!#4
!=======================================================================
!  subroutine sortip
!=======================================================================
!  function:
!    routine to sort array RA into ascending order, and to rearrange
!    the array JB accordingly.  Modified version of numerical recipes
!    routine sort2.
!
!  input:
!    n      : number of array elements
!
!  input/output:
!    ra     : array to be sorted
!    jb     : integer array to be sorted accordingly
!-----------------------------------------------------------------------
      SUBROUTINE sortip(n, ra, jb)

      IMPLICIT NONE

!=====  VARIABLES

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: n

!-----  input/output variables
      REAL    (KIND=8), INTENT(INOUT) :: ra(n)                ! (n)
      INTEGER (KIND=4), INTENT(INOUT) :: jb(n)                ! (n)

!-----  local variables
      INTEGER (KIND=4) :: i, j, ir, l, jjb
      REAL    (KIND=8) :: rra

!=====  MAIN ACTION

      l  = n/2 + 1
      ir = n

      DO

         IF (l > 1) THEN
            l   = l - 1
            rra = ra(l)
            jjb = jb(l)
         ELSE
            rra    = ra(ir)
            jjb    = jb(ir)
            ra(ir) = ra(1 )
            jb(ir) = jb(1 )
            ir     = ir - 1

            IF (ir == 1) THEN
               ra(1) = rra
               jb(1) = jjb
               EXIT
            END IF
         END IF

         i = l
         j = l + l

         DO
            IF (j <= ir) THEN
               IF (j < ir) THEN
                  IF (ra(j) < ra(j+1)) j = j + 1
               END IF

               IF (rra < ra(j)) THEN
                  ra(i) = ra(j)
                  jb(i) = jb(j)
                  i     = j
                  j     = j  + j
               ELSE
                  j     = ir + 1
               END IF

            ELSE
               EXIT
            END IF
         END DO

         ra(i) = rra
         jb(i) = jjb

      END DO

      END SUBROUTINE sortip

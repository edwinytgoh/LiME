!#1 subroutine codi

!--- tools
!#2 subroutine rnp82b
!#3 subroutine rnp82c
!#4 subroutine rnp82e

!#1
!=======================================================================
!  subroutine codi
!=======================================================================
!  function : 
!    routine to advance scalars for one time step
!    according to coalescence dispersion model
!    Reference: S.B. Pope, An improved turbulent mixing model,
!               Comb. Sci. Tech. 28: 131-145 (1982)
!
!  input :
!    mixext        : determines the extend of mixing
!                    1-full mixing
!                    2-mixing uniformly distributed between 0 and 1
!                    3-distribution (b) from Pope (1982, Comb.Sci.Tech.)
!                    4-distribution (c) ""
!                    5-distribution (d) ""
!                    6-distribution (e) ""
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
!    f_kf(npd,nsv) : particle independent scalars
!
!  output :
!    f_kf(npd,nsv) : particle independent scalars
!
!  includes :
!    param
!-----------------------------------------------------------------------
      SUBROUTINE codi (loutf, &
                       mixext, dt, cphi, &
                       ncell, &
                       omean_cell, &
                       np, npd, nsv, nf, nl, &
                       f_kwt, f_kf)

!-----  used modules
      USE cmn_param, ONLY : one, two, three, four, ten, PI, &
                            half, d3d2, d4d3
      USE rndgen_module, ONLY : rnu

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: mixext
      INTEGER (KIND=4), INTENT(IN) :: ncell
      INTEGER (KIND=4), INTENT(IN) :: np, npd, nsv
      INTEGER (KIND=4), INTENT(IN) :: nf(ncell)
      INTEGER (KIND=4), INTENT(IN) :: nl(ncell)
      REAL    (KIND=8), INTENT(IN) :: dt, cphi
      REAL    (KIND=8), INTENT(IN) :: omean_cell(ncell)
      REAL    (KIND=8), INTENT(IN) :: f_kwt(np)

!-----  input/output variable
      REAL    (KIND=8), INTENT(INOUT) :: f_kf(npd,nsv)

!----- local variables
      INTEGER (KIND=4) :: i, k, lc, nmixed, totnum
      REAL    (KIND=8) :: sumwt, omean, tauwt, prwt1, prmix, newscv

      INTEGER (KIND=4) :: mixnum, ip, iq, j
      REAL    (KIND=8) :: wghtmax, &
                          corrmax, corrwt, corrpr, corrext, &
                          realpq(2), rannum(1), amext(1)

      INTEGER (KIND=4), SAVE :: icall_codi = 0
      REAL    (KIND=8), SAVE :: corrext2, corrext3, corrext4, &
                                corrext5, corrext6

!------ MAIN ACTION

      IF (icall_codi == 0) THEN
         icall_codi = 1
         corrext2 = one / (one - one/three)
         corrext3 = one / (d3d2 - (ten/6.d0 - 30.d0/28.d0))
         corrext4 = one / (one + four/PI**2 - two/PI**2 - one/three)
         corrext5 = one / (d4d3 - half)
         corrext6 = one / (1.25d0 - 9.d0/20.d0)
      END IF

!-----  calculate the average inverse factor in which the variance
!-----  decay is increased by the inclusion of a variable extend of mixing

      SELECT CASE (mixext)
         CASE (2); corrext = corrext2
         CASE (3); corrext = corrext3
         CASE (4); corrext = corrext4
         CASE (5); corrext = corrext5
         CASE (6); corrext = corrext6
         CASE DEFAULT; corrext = 1.d0
      END SELECT


      !--- loop over cells 
      DO lc = 1, ncell 

         !--- number of particles in cell
         totnum = nl(lc) - nf(lc) + 1

         !--- check number of particles
         !    no mixing if totnum = 0 (empty cell) or 1 (only one particle)
         IF ((totnum == 0) .OR. (totnum == 1)) CYCLE

         !--- calculate total weight in cell and the maximum particle weight
         sumwt   = f_kwt(nf(lc))
         wghtmax = f_kwt(nf(lc))
         DO i = nf(lc) + 1, nl(lc)
            sumwt = sumwt + f_kwt(i)
            IF (f_kwt(i) > wghtmax) wghtmax = f_kwt(i)
         END DO

         !--- calculate the maximum correction factor possible and the
         !    part of the correction factor that is the same for all pairs
         !    in this cell
         corrwt  = half * (totnum / sumwt)
         corrmax = corrext * corrwt * (wghtmax + wghtmax)
       
         !--- calculate mean omega in cell
         omean = omean_cell(lc)
      
         !--- calculate time scale for one pair
         tauwt = one / (omean * cphi * totnum)

         !--- calculate would-be mixing probability for one pair
         !    with weight 1 for both particles
         !    this probability can be > 1 
         prwt1  = dt / tauwt 
!        prwt1  = REAL(totnum) * (1._rwp - EXP(-omean * cphi * dt))

         !--- calculate the number of times one pair has to mix with
         !    the same probability (try mixing at least two times)
         mixnum = 2.d0 * INT(one + prwt1 * corrmax)

         prwt1  = corrext * corrwt * prwt1 / DBLE(mixnum)


         !--- mixing loop
         nmixed = 0
         DO j = 1, mixnum

            !--- create a mixing pair and calculate the mixing probability for
            !    this pair

            DO
               CALL rnu(realpq, 2)
               ip = nf(lc) + INT(DBLE(totnum) * realpq(1))
               iq = nf(lc) + INT(DBLE(totnum) * realpq(2))
               IF (ip /= iq) EXIT
            END DO

            corrpr = f_kwt(ip) + f_kwt(iq)

            prmix = prwt1 * corrpr

            IF (prmix > one) THEN
               WRITE (loutf,*) 'prmix too large, prmix = ', prmix
               WRITE (loutf,*) '                 C     = ', corrwt * corrpr
            END IF

            !--- determine whether pair should mix 
            CALL rnu(rannum(1), 1)
            IF (rannum(1) < prmix) THEN
               nmixed = nmixed + 1
               SELECT CASE (mixext)
                  CASE (2); CALL rnu   (amext(1), 1)
                  CASE (3); CALL rnp82b(amext(1), 1)
                  CASE (4); CALL rnp82c(amext(1), 1)
                  CASE (5); CALL rnu   (amext(1), 1); amext(1) = SQRT(amext(1))
                  CASE (6); CALL rnp82e(amext(1), 1)
                  CASE DEFAULT; amext(1) = one
               END SELECT

               DO k = 1, nsv
                  newscv  = (f_kwt(ip)*f_kf(ip,k) + f_kwt(iq)*f_kf(iq,k)) &
                            / (f_kwt(ip) + f_kwt(iq))
                  f_kf(ip,k) = (one - amext(1))*f_kf(ip,k) + amext(1)*newscv
                  f_kf(iq,k) = (one - amext(1))*f_kf(iq,k) + amext(1)*newscv
               END DO

            END IF

         END DO

      END DO

      END SUBROUTINE codi


!#2
!=======================================================================
!  subroutine rnp82b( x, n )
!=======================================================================
!  routine to generate pseudo-random numbers which have a pdf
!  A(a) = 10a^3 * (1-3a/4)
!  this pdf is used in Pope (1982, Comb.Sci.Tech. Vol.28 pp.131-135)
!
!  method: first guess with A(a) = a^3, followed by two step Newton-
!          Raphson
!  110294: new method for generating random numbers implemented in all
!          subroutines. method comparable to hit-or-miss monte-carlo
!          integration technique. Old subroutines moved to
!          'rnfuncs_old.f'
!
!  input :
!    n           : length of the array of random numbers
!
!  output :
!    x           : array of random numbers
!----------------------------------------------------------------------
      SUBROUTINE rnp82b(x,n)

      USE rndgen_module, ONLY : rnu

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: n

!----- output variables
      REAL    (KIND=8), INTENT(OUT) :: x(n)

!----- local variables
      INTEGER (KIND=4) :: i
      REAL    (KIND=8) :: y(1), c(2), p

!----- generate array of n uniformly distriduted random numbers
      CALL rnu(x, n)

!----- make n random numbers
      DO i = 1, n
         CALL rnu(y, 1)
         DO
            p = x(i)**3 * (4.d0 - 3.d0 * x(i))
            IF (y(1) > p) THEN
               CALL rnu(c, 2)
               x(i) = c(1)
               y(1) = c(2)
               CYCLE
            ELSE
               EXIT
            END IF
         END DO
      END DO

      END SUBROUTINE rnp82b


!#3
!=======================================================================
!  subroutine rnp82c( x, n )
!=======================================================================
!  routine to generate pseudo-random numbers which have a pdf
!  A(a) = 1 - cos(pi * a)
!  this pdf is used in Pope (1982, Comb.Sci.Tech. Vol.28 pp.131-135)
!
!  method: first guess with A(a) = a^2.4, followed by two steps Newton-
!          Raphson
!  110294: new method for generating random numbers implemented in all
!          subroutines. method comparable to hit-or-miss monte-carlo
!          integration technique. Old subroutines moved to
!          'rnfuncs_old.f'
!
!  input :
!    n           : length of the array of random numbers
!
!  output :
!    x           : array of random numbers
!----------------------------------------------------------------------
      SUBROUTINE rnp82c(x,n)

      USE cmn_param, ONLY : PI
      USE rndgen_module, ONLY : rnu

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: n

!----- output variables
      REAL    (KIND=8), INTENT(OUT) :: x(n)

!----- local variables
      INTEGER (KIND=4) :: i
      REAL    (KIND=8) ::  y(1), c(2), p

!----- generate array of n uniformly distriduted random numbers
      CALL rnu(x, n)

!----- make n random numbers
      DO i = 1, n
         CALL rnu(y, 1)
         DO
            p = 1.d0 - DCOS(PI * x(i))
            IF (2.d0 * y(1) > p) THEN
               CALL rnu(c ,2)
               x(i) = c(1)
               y(1) = c(2)
               CYCLE
            ELSE
               EXIT
            END IF
         END DO
      END DO

      END SUBROUTINE rnp82c


!#4
!=======================================================================
!  subroutine rnp82e( x, n )
!=======================================================================
!  routine to generate pseudo-random numbers which have a pdf
!  A(a) = 3a(1 - 0.5a)
!  this pdf is used in Pope (1982, Comb.Sci.Tech. Vol.28 pp.131-135)
!
!  method: first guess with A(a) = a^1.7, followed by two steps Newton-
!          Raphson
!  110294: new method for generating random numbers implemented in all
!          subroutines. method comparable to hit-or-miss monte-carlo
!          integration technique. Old subroutines moved to
!          'rnfuncs_old.f'
!
!  input :
!    n           : length of the array of random numbers
!
!  output :
!    x           : array of random numbers
!----------------------------------------------------------------------
      SUBROUTINE rnp82e(x,n)

      USE rndgen_module, ONLY : rnu

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: n

!----- output variables
      REAL    (KIND=8), INTENT(OUT) :: x(n)

!----- local variables
      INTEGER (KIND=4) :: i
      REAL    (KIND=8) :: y(1), c(2), p

!----- generate array of n uniformly distriduted random numbers
      CALL rnu(x, n)

!----- make n random numbers
      DO i = 1, n
         CALL rnu(y, 1)
         DO
            p = 2.d0 * x(i) * (1.d0 - 0.5d0 * x(i))
            IF (y(1) > p) THEN
               CALL rnu(c, 2)
               x(i) = c(1)
               y(1) = c(2)
               CYCLE
            ELSE
               EXIT
            END IF
         END DO
      END DO

      END SUBROUTINE rnp82e

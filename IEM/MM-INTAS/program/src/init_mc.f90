!#1 subroutine init_mc
!#2 subroutine initpp
!#3 function genbet


!#1
!======================================================================= 
!  subroutine init_mc
!======================================================================= 
!  function :
!    Routine to initialise the Monte Carlo procedure. Allocates and
!    initialises the property array f
!-----------------------------------------------------------------------
      SUBROUTINE init_mc (loutf, llogf, icasestudied, nindv)

!-----  used modules
      USE cmn_gaspdf, ONLY : np, nkd, nsv, ndv, f, kf, kdf, kvol, &
                             relvol, relvol_init
      USE cmn_lookup, ONLY : case_methane

      IMPLICIT NONE

!-----  VARIABLES

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf, llogf
      INTEGER (KIND=4), INTENT(IN) :: icasestudied
      INTEGER (KIND=4), INTENT(IN) :: nindv

!-----  local variables
      INTEGER (KIND=4) :: i, npd


!-----  MAIN ACTION

      !---  allocate particle property array
      ALLOCATE (f(np,nkd))

      !--- fresh start: initialise particle properties
      WRITE(llogf,'(A)') 'INIT_MC: calling initpp'
      CALL initpp (np, nkd, nindv, f)

      !--- at this point, only the independent scalars are initialised
      !    initialise dependent scalars if necessary
      SELECT CASE (icasestudied)

      CASE (0)
         relvol      = 1.d0
         relvol_init = relvol

      CASE (case_methane)
         npd = np
         WRITE (llogf,'(A)') 'INIT_MC: calling lookup'
         CALL lookup (loutf, llogf, 0, &
                      np, npd, nsv, ndv, &
                      f(1,kf), f(1,kdf), f(1,kvol))

         !--- compute total relative volume
         relvol = 0.d0
         DO i = 1, np
            relvol = relvol + f(i,kvol)
         END DO
         relvol_init = relvol

      END SELECT

      END SUBROUTINE init_mc


!#2
!======================================================================= 
!      subroutine initpp
!======================================================================= 
!  function : 
!    routine to initialize particle properties
!
!  input : 
!    np         : number of particles
!    nkd        : number of particle properties
!    nsv        : number of independent scalars
!
!  input / output :
!    f(npd,nkd) : array containing particle properties
!-----------------------------------------------------------------------
      SUBROUTINE initpp (np, nkd, nsv, f)

!-----  used modules
      USE cmn_param, ONLY :  zero, one
      USE cmn_gaspdf, ONLY : initscal, &
                             mean_mixf,var_mixf,min_mixf,max_mixf, &
                             kf

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: np, nkd, nsv

!-----  input/output variables
      REAL    (KIND=8), INTENT(INOUT) :: f(np,nkd)


!-----  local variables
      INTEGER (KIND=4) :: i, half_np

      !--- for initialisation of scalars according to beta-PDF
      REAL    (KIND=8) :: var_max, aa, bb
      REAL    (KIND=8) :: sum_mixf, mixf_init

      REAL    (KIND=8), PARAMETER :: var_mixf_thres = 1.d-10
      REAL    (KIND=8), PARAMETER :: small = 1.d-20

      REAL    (KIND=8), EXTERNAL :: genbet

!-----  MAIN ACTION


!-----  INITIALISE INDEPENDENT SCALAR FIELDS
!       three possibilities: a. one peak
!                            b. two equal peaks
!                            c. beta function PDF for mixture fraction
!
!       This subroutine is written only for the case when 1 independent
!       scalar (mixture fraction) is considered

      IF (nsv > 1) THEN
         WRITE(*,*) 'INITPP: you have to modify this subroutine.'
         WRITE(*,*) 'initialisation is made only for one independent scalar.'
         STOP
      END IF


      SELECT CASE(initscal) 

      !--- a. one peak
      CASE(0)

         !--- loop over particles
         DO i = 1, np
            f(i,kf) = mean_mixf
         END DO


      !--- b. two equal peaks
      CASE(1)

         half_np = np / 2
         IF (2*half_np /= np) THEN
            WRITE(*,*) 'INITPP: To initialise two equal peaks, we want an even'
            WRITE(*,*) 'number of particles.'
            STOP
         END IF

         !--- loop over particles
         DO i = 1, half_np
            f(i,kf) = min_mixf
         END DO

         DO i = half_np+1, np
            f(i,kf) = max_mixf
         END DO


      !--- c. initialisation of the indepent scalars according to an
      !       assumed-shape PDF for mixture fraction: beta-function
      CASE(2)

         !--- maximum mixture fraction variance allowed
         var_max = mean_mixf * (one - mean_mixf)


         !--- var_mixf = 0 : degenerate delta function
         !                   (no turbulent fluctuations)
         IF (var_mixf < var_mixf_thres) THEN
            WRITE(*,*) 'INITPP: mixture fraction variance is too small.'
            WRITE(*,*) '        --> initialize with one peak'
            STOP

         !--- var_mixf large : degenerate double Dirac delta function
         !                     (intermittency)
         ELSE IF (var_mixf+small >= var_max) THEN
            WRITE(*,*) 'INITPP: mixture fraction variance is too large.'
            WRITE(*,*) '        --> decrease its value'
            STOP

         !--- general case
         ELSE

            DO i = 1, np

               aa = mean_mixf * mean_mixf * (one - mean_mixf) &
                     / var_mixf - mean_mixf
               bb = aa * (one / mean_mixf - one)

               !---  mixture fraction according to beta function
               f(i,kf) = genbet(aa,bb)

            END DO

         END IF


         !--- Check the beta-initialisation
         sum_mixf = zero
         DO i = 1, np
            sum_mixf = sum_mixf + f(i,kf)
         END DO
         mixf_init = sum_mixf / DBLE(np)

         WRITE(*,'(a)')
         WRITE(*,*) 'INITPP: beta-initialisation.'
         WRITE(*,*) '  mean_mixf (input) = ', mean_mixf
         WRITE(*,*) '  mean_mixf (init)  = ', mixf_init
         WRITE(*,'(a)')

         sum_mixf = zero
         DO i = 1, np
            sum_mixf = sum_mixf + (f(i,kf)-mixf_init)*(f(i,kf)-mixf_init)
         END DO
         mixf_init = sum_mixf / DBLE(np)

         WRITE(*,*) '   var_mixf (input) = ', var_mixf
         WRITE(*,*) '   var_mixf (init)  = ', mixf_init
         WRITE(*,'(a)')

      END SELECT

      END SUBROUTINE initpp


!#3
!======================================================================
!     REAL*8 FUNCTION GENBET( A, B )
!======================================================================
!               GeNerate BETa random deviate
!
!
!                              Function
!
!
!     Returns a single random deviate from the beta distribution with
!     parameters A and B.  The density of the beta is
!               x^(a-1) * (1-x)^(b-1) / B(a,b) for 0 < x < 1
!
!
!                              Arguments
!
!
!     A --> First parameter of the beta distribution
!                         REAL*8 A
!     JJV                 (A > 1.0E-37)
!
!     B --> Second parameter of the beta distribution
!                         REAL*8 B
!     JJV                 (B > 1.0E-37)
!
!
!                              Method
!
!
!     R. C. H. Cheng
!     Generating Beta Variatew with Nonintegral Shape Parameters
!     Communications of the ACM, 21:317-322  (1978)
!     (Algorithms BB and BC)
!
!----------------------------------------------------------------------
      FUNCTION genbet(aa,bb)

      USE rndgen_module, ONLY : rnu_scalar

      IMPLICIT NONE

      REAL (KIND=8) genbet

!----- Scalar Arguments
      REAL    (KIND=8), INTENT(IN) :: aa, bb

!----- Parameters

!     Close to the largest number that can be exponentiated
!     JJV changed this - 89 was too high, and LOG(1.0E38) = 87.49823
!     REAL    (KIND=4), PARAMETER :: expmax = 87.49823d0
!     double precision: LOG(1.0d308) = 709.196208642
      REAL    (KIND=8), PARAMETER :: expmax = 700.0d0

!     Close to the largest representable single precision number
!     REAL    (KIND=4), PARAMETER :: infnty = 1.E38
      REAL    (KIND=8), PARAMETER :: infnty = 1.d308

!     JJV added the parameter minlog
!     Close to the smallest number of which a LOG can be taken.
!     REAL    (KIND=4), PARAMETER :: minlog = 1.0E-37
      REAL    (KIND=8), PARAMETER :: minlog = 1.0d-308

!----- Local Scalars
!     REAL    (KIND=8) :: a, b, delta, r, s, t, u1, u2, v, w, y, z
      REAL    (KIND=8) :: delta, r, s, t, u1, u2, v, w, y, z
      LOGICAL :: qsame

!----- Intrinsic Functions
!     INTRINSIC EXP,LOG,MAX,MIN,SQRT

!     JJV added a,b
      REAL    (KIND=8), SAVE :: a, b

!     JJV changed these to ridiculous values
!     REAL    (KIND=4), SAVE :: olda = -1.0E37
!     REAL    (KIND=4), SAVE :: oldb = -1.0E37
      REAL    (KIND=8), SAVE :: olda = -1.0d307
      REAL    (KIND=8), SAVE :: oldb = -1.0d307
      REAL    (KIND=8), SAVE :: alpha, beta, gamma, k1, k2

!     ..
!     .. Executable Statements ..
!     ON REAL*8 UNDERFLOW IGNORE


!=====  MAIN ACTION

      qsame = (olda == aa) .AND. (oldb == bb)

      IF (.NOT. qsame) THEN
!     JJV added small minimum for small log problem in calc of W
         IF (aa < minlog .OR. bb < minlog) THEN
            WRITE (*,*) ' AA or BB <', minlog, ' in GENBET - Abort!'
            WRITE (*,*) ' AA: ', aa, ' BB ', bb
            STOP ' AA or BB too small in GENBET - Abort!'
         END IF

         olda = aa
         oldb = bb
      END IF


      IF (MIN(aa,bb) > 1.0d0) THEN

!     Algorithm BB

!
!     Initialize
!
         IF (.NOT. qsame) THEN
            a = MIN(aa,bb)
            b = MAX(aa,bb)
            alpha = a + b
            beta = SQRT((alpha-2.0d0) / (2.0d0*a*b-alpha))
            gamma = a + 1.0d0/beta
         END IF

         LOOP_ALGORITHM_BB : &
         DO 
!
!     Step 1
!
            CALL rnu_scalar(u1)
            CALL rnu_scalar(u2)

            v = beta*LOG(u1/MAX(minlog,(1.0d0-u1)))
            IF (v >= expmax) THEN
               w = infnty
            ELSE
!     JJV added checker to see if a*exp(v) will overflow
!     JJV 50 _was_ w = a*exp(v); also note here a > 1.0
               w = EXP(v)
               IF (w > infnty/a) THEN
                  w = infnty
               ELSE
                  w = a*w
               END IF
            END IF
            z = u1*u1 *u2
            r = gamma*v - 1.3862944d0
            s = a + r - w
!
!     Step 2
!
            IF ((s+2.609438d0) >= (5.0d0*z)) EXIT LOOP_ALGORITHM_BB
!
!     Step 3
!
            t = LOG(z)
            IF (s > t) EXIT LOOP_ALGORITHM_BB
!
!     Step 4
!
!     JJV added checker to see if log(alpha/(b+w)) will
!     JJV overflow.  If so, we count the log as -INF, and
!     JJV consequently evaluate conditional as true, i.e.
!     JJV the algorithm rejects the trial and starts over
!     JJV May not need this here since ALPHA > 2.0
            IF (alpha/(b+w) < minlog) CYCLE LOOP_ALGORITHM_BB

            IF ((r+alpha*LOG(alpha/ (b+w))) < t) THEN
               CYCLE LOOP_ALGORITHM_BB
            ELSE
               EXIT LOOP_ALGORITHM_BB
            END IF

         END DO LOOP_ALGORITHM_BB
!
!     Step 5
!
         IF (aa == a) THEN
            genbet = w/ (b+w)
         ELSE
            genbet = b/ (b+w)
         END IF


      ELSE

!     Algorithm BC

!
!     Initialize
!
         IF (.NOT. qsame) THEN
            a = MAX(aa,bb)
            b = MIN(aa,bb)
            alpha = a + b
            beta  = 1.0d0 / b
            delta = 1.0d0 + a - b
            k1 = delta* (0.0138889d0+0.0416667d0*b)/ (a*beta-0.777778d0)
            k2 = 0.25d0 + (0.5d0+0.25d0/delta)*b
         END IF

         LOOP_ALGORITHM_BC : &
         DO 
!
!     Step 1
!
            CALL rnu_scalar(u1)
            CALL rnu_scalar(u2)
!
!     Step 2
!
            IF (u1 < 0.5d0) THEN
               y = u1*u2
               z = u1*y
               IF ((0.25d0*u2+z-y) >= k1) THEN
                  CYCLE LOOP_ALGORITHM_BC
               ELSE
                  GO TO 170
               END IF
            END IF
!
!     Step 3
!
            z = u1*u1 *u2
            IF (z <= 0.25d0) THEN
               v = beta*LOG(u1/ MAX(minlog,(1.0d0-u1)))
!               IF (v > expmax) THEN
!                  w = infnty
!               ELSE
!                  w = a*EXP(v)
!               END IF
!
!     JJV instead of checking v > expmax at top, I will check
!     JJV if a < 1, then check the appropriate values
               IF (a <= 1.0d0) THEN
!     JJV A < 1 so it can help out if EXP(V) would overflow
                  IF (v <= expmax) THEN
                     IF (v <= -expmax) THEN
                        w = 0.d0
                     ELSE
                        w = a*EXP(v)
                     END IF
                  ELSE
                     w = v + LOG(a)
                     IF (w <= expmax) THEN
                        IF (w <= -expmax) THEN
                           w = 0.d0
                        ELSE
                           w = EXP(w)
                        END IF
                     ELSE
                        w = infnty
                     END IF
                  END IF

!     JJV in this case A > 1
               ELSE IF (v <= expmax) THEN
                  IF (v <= -expmax) THEN
                     w = 0.d0
                  ELSE
                     w = EXP(v)
                  END IF
                  IF (w <= infnty/a) THEN
                     w = a*w
                  ELSE
                     w = infnty
                  END IF

               ELSE
                  w = infnty
               END IF


               EXIT LOOP_ALGORITHM_BC
            END IF

            IF (z >= k2) CYCLE LOOP_ALGORITHM_BC
!
!     Step 4
!
!
!     Step 5
!

  170       v = beta*LOG(u1/ (1.0d0-u1))
!     JJV same kind of checking as above
            IF (a <= 1.0d0) THEN
!     JJV A < 1 so it can help out if EXP(V) would overflow
               IF (v <= expmax) THEN
                  IF (v <= -expmax) THEN
                     w = 0.d0
                  ELSE
                     w = a*EXP(v)
                  END IF
               ELSE
                  w = v + LOG(a)
                  IF (w <= expmax) THEN
                     IF (w <= -expmax) THEN
                        w = 0.d0
                     ELSE
                        w = EXP(w)
                     END IF
                  ELSE
                     w = infnty
                  END IF
               END IF

!     JJV in this case A > 1
            ELSE IF (v <= expmax) THEN
               IF (v <= -expmax) THEN
                  w = 0.d0
               ELSE
                  w = EXP(v)
               END IF
               IF (w <= infnty/a) THEN
                  w = a*w
               ELSE
                  w = infnty
               END IF

            ELSE
               w = infnty
            END IF

!     JJV here we also check to see if log overlows; if so, we treat it
!     JJV as -INF, which means condition is true, i.e. restart
            IF (alpha/(b+w) < minlog) CYCLE LOOP_ALGORITHM_BC

            IF ((alpha* (LOG(alpha/ (b+w))+v)-1.3862944d0) < LOG(z)) THEN
               CYCLE LOOP_ALGORITHM_BC
            ELSE
               EXIT LOOP_ALGORITHM_BC
            END IF

         END DO LOOP_ALGORITHM_BC
!
!     Step 6
!
         IF (a == aa) THEN
            genbet = w/ (b+w)
         ELSE
            genbet = b/ (b+w)
         END IF

      END IF

      RETURN
      END FUNCTION genbet

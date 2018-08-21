!#1  subroutine rnu
!#2  subroutine rnu_scalar
!#3  subroutine rnuget
!#4  subroutine rnuput
!#5  subroutine rndg
!#6  subroutine rndgget
!#7  subroutine rndgput


      MODULE rndgen_module

!----- set parameter PI
      REAL    (KIND=8), PARAMETER, PRIVATE :: PI = 3.14159265358979d0

!----- seeds for rnu routine
      INTEGER (KIND=4), SAVE, PRIVATE :: i1 = 12345
      INTEGER (KIND=4), SAVE, PRIVATE :: i2 = 67890

!-----  variables for rndg routine
      INTEGER (KIND=4), PARAMETER, PRIVATE :: n1 = 8
      INTEGER (KIND=4), PARAMETER, PRIVATE :: n2 = 9
      INTEGER (KIND=4), PARAMETER, PRIVATE :: n3 = 10
      INTEGER (KIND=4), PARAMETER, PRIVATE :: n4 = 12
      INTEGER (KIND=4), SAVE, PRIVATE :: ifst = 0
      INTEGER (KIND=4), SAVE, PRIVATE :: nd
      REAL    (KIND=8), SAVE, PRIVATE :: g(n1*n2*n3*n4)


      CONTAINS


!#1
!***********************************************************************
!     subroutine rnu
!
!  function:
!    Routine to generate pseudo-random numbers uniformly in the
!    exclusive interval [0,1]. Method, see f. james, computer physics
!    communications, vol.60, p.329, 1990.
!    The entries  rnuget  and  rnuput  can be used to get and set
!    the seeds i1 and i2.
!
!  input:
!    n        : number of values requested
!
!  output:
!    x        : array of n random numbers
!
!  For generic purposes, two version of rnu are defined below:
!  rnu_array : (original rnu)
!              x is declared as a real array.
!  rnu_scalar: (new)
!              x is declared as a scalar.
!              in principle, n should be equal to 1
!              if not, this is overruled by rnu2, and n is replaced by 1.
!***********************************************************************
      SUBROUTINE rnu(x, n)

      IMPLICIT NONE

!----- array declaration
      INTEGER (KIND=4), INTENT(IN)  :: n
      REAL    (KIND=8), INTENT(OUT) :: x(n)

      INTEGER (KIND=4) :: i, ik, ix

      DO i = 1, n

         ik = i1 / 53668
         i1 = 40014 * ( i1 - ik * 53668 ) - ik * 12211
         IF (i1 < 0) i1 = i1 + 2147483563

         ik = i2 / 52774
         i2 = 40692 * ( i2 - ik * 52774 ) - ik * 3791
         IF (i2 < 0) i2 = i2 + 2147483399

         ix = i1 - i2
         IF (ix < 1) ix = ix + 2147483562

         x(i) = DBLE(ix) * 4.656612d-10
      END DO

      END SUBROUTINE rnu


!#2
!***********************************************************************
      SUBROUTINE rnu_scalar(x)
!***********************************************************************
!  Scalar version of rnu1.
!  Uses the same random seeds i1 and i2 as rnu1.
!-----------------------------------------------------------------------

      IMPLICIT NONE

!----- array declaration
      REAL    (KIND=8), INTENT(OUT) :: x

      INTEGER (KIND=4) :: ik, ix

      ik = i1 / 53668
      i1 = 40014 * ( i1 - ik * 53668 ) - ik * 12211
      IF (i1 < 0) i1 = i1 + 2147483563

      ik = i2 / 52774
      i2 = 40692 * ( i2 - ik * 52774 ) - ik * 3791
      IF (i2 < 0) i2 = i2 + 2147483399

      ix = i1 - i2
      IF (ix < 1) ix = ix + 2147483562

      x = DBLE(ix) * 4.656612d-10

      END SUBROUTINE rnu_scalar


!#3
!***********************************************************************
      SUBROUTINE rnuget(is1, is2)
!***********************************************************************
!  Get parameters i1 and i2 used by subroutine rnu.
!-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER (KIND=4), INTENT(OUT) :: is1, is2

      is1 = i1
      is2 = i2

      END SUBROUTINE rnuget


!#4
!***********************************************************************
      SUBROUTINE rnuput(is1, is2)
!***********************************************************************
!  Set parameters i1 and i2 used by subroutine rnu.
!-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER (KIND=4), INTENT(IN) :: is1, is2

      i1 = is1
      i2 = is2

      END SUBROUTINE rnuput


!#5
!***********************************************************************
!  subroutine rndg
!
!  function:
!    routine to generate discrete random numbers that approximate
!    a gaussian distribution
!
!  method:
!    the random number returned, x, is a discrete random
!    variable, with equally-probable values g(i), i=1, nd.
!    once the discrete values g(i) have been formed, sampling
!    is trivial (see the 100 loop).
!    the values of g(i) are formed as all possible sums of the
!    values of 4 other discrete variables, whose (positive) values,
!    given in data statements, are created by the program dggen.
!    the odd moments are zero by symmetry.  to machine accuracy
!    the variance is 1, and the flatness is 3.  with the
!    parameters set to n1=8, n2=9, n3=10, n4=12, the 6-th, 8-th
!    and 10-th moments are: 14.73, 98.30, and 811.6  (cf 15, 105,
!    945 for a gaussian).  for all x, the difference between the
!    cdf and the gaussian cdf (i.e. the error in the cdf) is less
!    than  0.005.
!
!  input:
!    n          : number of values requested
!
!  output:
!    x          : array of random numbers
!
!  includes:
!
!  calls:
!    rnu
!
!***********************************************************************
      SUBROUTINE rndg( x, n )

      IMPLICIT NONE

      INTEGER (KIND=4), INTENT(IN) :: n
      REAL    (KIND=8), INTENT(OUT) :: x(n)

      INTEGER (KIND=4) :: nc(4), amom(10)
      REAL    (KIND=8) :: gc(n4,4)

      INTEGER (KIND=4) :: i, ig, nh, j1, j2, j3, j4, icheck, im, j
      REAL    (KIND=8) :: sum, g1, g2, g3

      DATA ( gc(i,1), i = 1, 4 ) / 0.0791315883d0,  0.3056624234d0, &
     &                             0.6851569414d0,  1.8522603512d0  /

      DATA ( gc(i,2), i = 1, 4 ) / 0.1839743704d0,  0.4354168475d0, &
     &                             0.8178389072d0,  1.8993961811d0  /

      DATA ( gc(i,3), i = 1, 5 ) / 0.0856366232d0,  0.2853691578d0, &
     &                             0.5440011621d0,  0.9214153290d0, &
     &                                            1.9406925440d0  /

      DATA ( gc(i,4), i = 1, 6 ) / 0.0801580697d0, 0.2553057969d0, &
     &                             0.4600485563d0, 0.7161571383d0, &
     &                             1.0778864622d0, 2.0104796886d0  /

!----- on first call, determine discrete values
      IF (ifst == 0) THEN
         ifst = 1
         nd = n1*n2*n3*n4
         nc(1) = n1
         nc(2) = n2
         nc(3) = n3
         nc(4) = n4

!-------- complete component discrete gaussians
         DO ig = 1, 4
            nh = nc(ig) / 2
            gc(nh+1,ig) = 0.d0

            DO i = 1, nh
               gc(nc(ig)+1-i, ig) = - gc(i, ig)
            END DO
         END DO

!-------- form compound discrete gaussian
         i = 0
         DO j1 = 1, n1
            g1 = gc(j1, 1)

            DO j2 = 1, n2
               g2 = g1 + gc(j2, 2)

               DO j3 = 1, n3
                  g3 = g2 + gc(j3, 3)

                  DO j4 = 1, n4
                     i = i + 1
                     g(i) = 0.5d0*(g3 + gc(j4, 4))
                  END DO
               END DO
            END DO
         END DO

!-------- check g
         icheck = 0
         IF (icheck /= 0) THEN

            DO im = 1, 10
               sum = 0.d0

               DO i = 1, nd
                  sum = sum + g(i) ** im
               END DO

               amom(im) = sum / nd
            END DO

         END IF

      END IF

!----- end of set-up
!----- randomly sample from discrete values

      CALL rnu(x, n)

      DO j = 1, n
         i = 1 + x(j) * nd
         x(j) = g(i)
      END DO

      END SUBROUTINE rndg


!#6
!***********************************************************************
      SUBROUTINE rndgget(ifsts, nds, gs)
!***********************************************************************
!  Get parameters ifst, nd, g used by subroutine rndg.
!-----------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER (KIND=4), INTENT(OUT) :: ifsts, nds
      REAL    (KIND=8), INTENT(OUT) :: gs(*)
      INTEGER (KIND=4) :: i

      ifsts = ifst
      nds   = nd
      DO i = 1, n1*n2*n3*n4
         gs(i) = g(i)
      END DO

      END SUBROUTINE rndgget


!#7
!***********************************************************************
      SUBROUTINE rndgput(ifsts, nds, gs)
!***********************************************************************
!  Set parameters ifst, nd, g used by subroutine rndg.
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER (KIND=4), INTENT(IN) :: ifsts, nds
      REAL    (KIND=8), INTENT(IN) :: gs(*)

      INTEGER (KIND=4) :: i

      ifst = ifsts
      nd   = nds
      DO i = 1, n1*n2*n3*n4
         g(i) = gs(i)
      END DO

      END SUBROUTINE rndgput


      END MODULE rndgen_module

!#1 subroutine check_input


!#1
!=======================================================================
!  SUBROUTINE check_input
!=======================================================================
!  function :
!    routine to check input variables 
!-----------------------------------------------------------------------
      SUBROUTINE check_input (loutf, icasestudied, nstep, dt)

!-----  used modules
      USE cmn_gaspdf, ONLY : np
      USE cmn_lookup, ONLY : ilut, case_methane
      USE cmn_pdfoutp, ONLY : ipdfop
      USE cmn_mixing, ONLY : mixmod, mixext, cphi

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: icasestudied, nstep
      REAL    (KIND=8), INTENT(IN) :: dt

!-----  local variables
      INTEGER (KIND=4) :: nerr, nwar
      REAL    (KIND=8), PARAMETER :: small = 1.d-10


!-----  MAIN ACTION

      nerr = 0
      nwar = 0

      WRITE(loutf,'(/A,/A,/A)') &
         '****************************', &
         'SUBROUTINE CHECK_INPUT', &
         '****************************'


!-----  number of particles
      IF (np <= 0) THEN
         WRITE(loutf,*)  'Error : np <= 0'
         nerr = nerr + 1
      END IF


!-----  look-up
      IF (icasestudied < 0 .OR. icasestudied > case_methane) THEN
         WRITE(loutf,*) 'Error : icasestudied <> 0,1'
         nerr = nerr + 1
      END IF

      IF (ilut < 0) THEN
         WRITE(loutf,*)  'Error : ilut < 0'
         nerr = nerr + 1
      END IF


!-----  iterations
      IF (nstep < 0) THEN
         WRITE(loutf,*)  'Error : nstep < 0'
         nerr = nerr + 1
      END IF
      IF (dt <= 0.d0) THEN
         WRITE(loutf,*)  'Error : dt <= 0'
         nerr = nerr + 1
      END IF


!-----  output
      IF (ipdfop < 0) THEN
         WRITE(loutf,*)  'Error : ipdfop < 0'
         nerr = nerr + 1
      END IF


!-----  mixing
      IF (mixmod < 0 .OR. mixmod > 3) THEN
         WRITE(loutf,*)  'Error : mixmod <> 0 - 3'
         nerr = nerr + 1
      END IF

      IF (mixext < 1) THEN
         WRITE(loutf,*)  'Error : mixext < 1'
         nerr = nerr + 1
      END IF
      IF (mixext > 6) THEN
         WRITE(loutf,*)  'Error : mixext > 6'
         nerr = nerr + 1
      END IF

      !--- model constants
      IF (cphi < 0.) THEN
         WRITE(loutf,*) 'Error : cphi < 0.'
         nerr = nerr + 1
      ELSE IF (ABS(cphi - 2.0) >= 1E-6) THEN
         WRITE(loutf,*)  'Warning : ABS(cphi - 2.0) >= 1E-6'
         nwar = nwar + 1
      END IF


!-----  totals

      WRITE(loutf,'(/A,I4)') 'Number of errors   : ', nerr
      WRITE(loutf,'( A,I4)') 'Number of warnings : ', nwar
      IF (nerr > 0) THEN
         WRITE(loutf,'(/A/)')   'STOPPED IN CHECK_INPUT'
         STOP 'CHECK_INPUT: see output file'
      END IF
     
      END SUBROUTINE check_input

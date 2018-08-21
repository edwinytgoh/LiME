!#1 subroutine output
!#2  subroutine catnum


!#1
!======================================================================= 
!       subroutine output
!
!  function :
!    generates outputs
!
!  input :
!    istep   : current iteration for the current calculation
!              ( 1 < . < nstep)
!    nstep   : number of iterations for the current calculation
!    time    : time corresponding to istep
!-----------------------------------------------------------------------
      SUBROUTINE output (loutf, icasestudied, istep, nstep, time)


!-----  used modules
      USE cmn_param, ONLY : fileln
      USE cmn_gaspdf, ONLY : f, iprlab, np, nsv, ndv, kf, kl, kdf, kdl
      USE cmn_pdfoutp, ONLY : ipdfop
      USE files_module, ONLY : lpdfout, pdfout

!      USE interface_module, this_subroutine => output

      IMPLICIT NONE

!-----  VARIABLES

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: icasestudied
      INTEGER (KIND=4), INTENT(IN) :: istep, nstep
      REAL    (KIND=8), INTENT(IN) :: time

!----- local variables
      INTEGER (KIND=4) :: i, k
      INTEGER (KIND=4) :: lformt
      CHARACTER(LEN=32) :: formt
      CHARACTER(LEN=fileln+5) :: pdfout_it


!-----  MAIN ACTION

      IF (ipdfop == 0) RETURN
      IF (istep < nstep) THEN
         IF (MOD(istep, ipdfop) /= 0) RETURN
      END IF
 

      !---  construct datafile name
      pdfout_it = TRIM(pdfout)
      IF (istep < 10000) THEN
         CALL catnum(loutf, pdfout_it,0,pdfout_it)
         IF (istep < 1000) THEN
            CALL catnum(loutf, pdfout_it,0,pdfout_it)
            IF (istep < 100) THEN
               CALL catnum(loutf, pdfout_it,0,pdfout_it)
               IF (istep < 10) THEN
                  CALL catnum(loutf, pdfout_it,0,pdfout_it)
               END IF
            END IF
         END IF
      END IF
      CALL catnum(loutf, pdfout_it, istep, pdfout_it)
 
      !--- open file
      OPEN(UNIT=lpdfout, &
           FILE=pdfout_it, &
           FORM='FORMATTED', &
           ACCESS='SEQUENTIAL')

      !--- write time
      WRITE(lpdfout,'((A7),(F12.5),(A9))') '# After', time, ' seconds.'

      !--- write property names
      WRITE(lpdfout,'(A22)') '# Independent scalars:'
      DO k = kf, kl
         WRITE(lpdfout,'((A3),(A))') '#  ', TRIM(iprlab(k))
      END DO
      IF (icasestudied > 0) THEN
         WRITE(lpdfout,'(A20)') '# Dependent scalars:'
         DO k = kdf, kdl
            WRITE(lpdfout,'((A3),(A))') '#  ', TRIM(iprlab(k))
         END DO
      ELSE
         WRITE(lpdfout,'(A33)') '# There are no dependent scalars.'
      END IF

      !--- write data
      IF (icasestudied > 0) THEN
         formt(1:32) = '((I8),   (E15.4E3),   (E15.4E3))'
         lformt = 32
         WRITE(formt(7:9),'(I3)') nsv
         WRITE(formt(20:22),'(I3)') ndv
         DO i = 1, np
            WRITE(lpdfout,FMT=formt(1:lformt)) &
              i, (f(i,k), k=kf,kl), (f(i,k), k=kdf,kdl)
         END DO
      ELSE
         formt(1:19) = '((I8),   (E15.4E3))'
         lformt = 19
         WRITE(formt(7:9),'(I3)') nsv
         DO i = 1, np
            WRITE(lpdfout,FMT=formt(1:lformt)) i, (f(i,k), k=kf,kl)
         END DO
      END IF

      !--- close file
      CLOSE(lpdfout)

      END SUBROUTINE output


!#2
!======================================================================= 
!     subroutine catnum
!======================================================================= 
!  function :
!    construct output string from input string and integer:
!    
!  <output-string> = <input_string><integer>
!
!  input :
!    strin   : input string
!    intin   : input integer
!    
!  output :
!    strout  : output string
!----------------------------------------------------------------------- 
      SUBROUTINE catnum(loutf, strin, intin, strout)

      IMPLICIT NONE

!-----  input variables
      CHARACTER(LEN=*), INTENT(IN) :: strin
      INTEGER (KIND=4), INTENT(IN) :: loutf, intin

!-----  output variables
      CHARACTER(LEN=*), INTENT(OUT) :: strout

!-----  local variables
      INTEGER (KIND=4) :: nrstrlen, strlen
      CHARACTER(LEN=6) :: nrstr

!-----  MAIN ACTION

!-----  find length of input string
      strlen = LEN_TRIM(strin)

!-----  construct number string
      WRITE(nrstr(1:6), '(I6)') intin

!-----  strip leading spaces from number string
      nrstr = ADJUSTL(nrstr)
      nrstrlen = LEN_TRIM(nrstr)

!-----  check combined length
      IF ((strlen + nrstrlen) > LEN(strout)) THEN
         WRITE(loutf,*) 'Error in subroutine catnum:'
         WRITE(loutf,*) 'String too large'
         STOP
      END IF

!-----  construct combined string
      strout(1:strlen) = strin(1:strlen)
      strout(strlen + 1:strlen + nrstrlen) = nrstr(1:nrstrlen)
      strout(strlen + nrstrlen + 1:) = ' '

      END SUBROUTINE catnum

!#1 subroutine mainloop
!#2 subroutine logrew

!#1
!=======================================================================        
!       subroutine mainloop
!=======================================================================        
!  input : 
!    nstep            : number of iterations to be made
!    dt               : basic time step
!
!  output :
!    global_iteration : global iteration number
!    time             : global time
!-----------------------------------------------------------------------        
      SUBROUTINE mainloop (loutf, llogf, icasestudied, nstep, dt)

!-----  used modules
      USE cmn_gaspdf, ONLY : np, nsv, ndv, kvol, kf, kdf, &
                             f, relvol, relvol_init
      USE cmn_mixing, ONLY : omega

!      USE interface_module, this_subroutine => mainloop

      IMPLICIT NONE

!-----  VARIABLES

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf, llogf
      INTEGER (KIND=4), INTENT(IN) :: icasestudied
      INTEGER (KIND=4), INTENT(IN) :: nstep
      REAL    (KIND=8), INTENT(IN) :: dt

!-----  local variables
      INTEGER (KIND=4) :: i
      INTEGER (KIND=4) :: istep
      INTEGER (KIND=4) :: npd
      REAL    (KIND=8) :: time


!=====  MAIN ACTION

      !---  write outputs before mixing occured
      istep = 0
      time = 0.d0
      WRITE(llogf,'(A)') 'MAINLOOP: calling outputs (after 0 iteration)'
      CALL output (loutf, icasestudied, istep, nstep, time)


!-----  MAIN LOOP

      npd = np ! maximum number of particles
               ! = size of the first dimension of array f
               ! (useful in a full Monte Carlo simulation
               !  where the total number of particles varies)

      MAIN_LOOP: &
      DO istep = 1, nstep

         WRITE(llogf,'(A,I6,A,I6)') &
        'MAINLOOP: starting iteration ', istep, ' out of ', nstep


         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         !   MICRO-MIXING   MICRO-MIXING   MICRO-MIXING   MICRO-MIXING
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         WRITE(llogf,'(A)') 'MAINLOOP: calling stepf'
         CALL stepf (loutf, llogf, &
                     dt, omega, np, npd, nsv, &
                     f(1,kf))
 
         !---  increment time
         time = time + dt

         IF (icasestudied > 0) THEN

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !  NEW PARTICLE COMPOSITION AND SPECIFIC VOLUME (= 1/DENSITY)
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            WRITE (llogf,'(A)') 'MAINLOOP: calling lookup'
            CALL lookup (loutf, llogf, istep, &
                         np, npd, nsv, ndv, &
                         f(1,kf), f(1,kdf), f(1,kvol))
 

            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !  NEW TOTAL RELATIVE VOLUME
            !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            relvol = 0.d0
            DO i = 1, np
               relvol = relvol + f(i,kvol)
            END DO
            WRITE(*,'((I5),(A6),(F18.15),(A18),(F18.15))') &
              istep, '   t =', time, &
                     '   vol/vol_init = ', relvol/relvol_init

         ELSE

            WRITE(*,'((I5),(A6),(F18.15))') istep, '   t =', time

         END IF


         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         !  OUTPUTS
         !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         !---  write outputs if necessary
         WRITE(llogf,'(A)') 'MAINLOOP: calling outputs'
         CALL output (loutf, icasestudied, istep, nstep, time)

         !---  rewind logfile if necessary
         IF (istep /= nstep) CALL logrew

      END DO MAIN_LOOP

      END SUBROUTINE mainloop


!#2
!=======================================================================
!  subroutine logrew
!=======================================================================
!  function : rewind logfile
!-----------------------------------------------------------------------
      SUBROUTINE logrew

!-----  used modules
      USE files_module, ONLY : llog, loutput, logfile
      
      IMPLICIT NONE

!-----  local variables
      CHARACTER(LEN=31) :: words

!-----  MAIN ACTION

!-----  check if action is to be taken
      IF ((llog /= 6) .AND. (llog /= loutput)  ) THEN
      
!--------  close logfile to clear buffer
         CLOSE(llog)
      
!--------  open logfile
         OPEN (UNIT=llog, &
               FILE=logfile, &
               FORM='FORMATTED', &
               ACCESS='SEQUENTIAL')

!--------  find string 
         DO
            READ(llog, '(A31)', END = 200, ERR = 200) words
            IF (words == 'MAINLOOP: starting iteration') EXIT
         END DO

!--------  backspace over line just read
         BACKSPACE(llog)

 200     CONTINUE

      END IF

      END SUBROUTINE logrew

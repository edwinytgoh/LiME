!#1 subroutine init_files
!#1 subroutine close_infofiles

!#1
!======================================================================= 
!       subroutine init_files
!======================================================================= 
!  function :
!    read the file names from the standard input
!    call the subroutines that initialise the files used in each module
!-----------------------------------------------------------------------
      SUBROUTINE init_files (loutf, llogf, linpt, inpt)

!-----  used modules
      USE cmn_param, ONLY : fileln
      USE files_module, ONLY : loutput, llog, linpt_in, inpt_in, &
                               lfldat, fldat, &
                               read_filenames, open_infofiles
      USE fllookup_module, ONLY : init_fllookup_files

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4) , INTENT(OUT) :: loutf, llogf

!----- output variables
      INTEGER (KIND=4) , INTENT(OUT) :: linpt
      CHARACTER(LEN=fileln) , INTENT(OUT) :: inpt


!-----  MAIN ACTION

!------  read names of the files: inputs, output and loginf
      CALL read_filenames (6, 0)

      CALL init_fllookup_files (loutput, llog, lfldat, fldat)

!------  open files: output and loginf
      CALL open_infofiles

      loutf = loutput
      llogf = llog
      linpt = linpt_in
       inpt =  inpt_in

      END SUBROUTINE init_files

!#2
!=======================================================================
!     subroutine close_infofiles
!=======================================================================
!  function :
!    close PDFD output files depending on logical unit numbers
!-----------------------------------------------------------------------
      SUBROUTINE close_infofiles (loutf, llogf)

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf, llogf

!-----  local variables
      INTEGER (KIND=4) :: i

!-----  MAIN ACTION

!-----  write tails and close files
      WRITE(llogf, '(72A1)') ('*', i=1, 72)
      WRITE(llogf, '(A)') 'E N D   O F   L O G   O U T P U T'
      WRITE(llogf, '(72A1)') ('*', i=1, 72)

      WRITE(loutf, '(72A1)') ('*', i=1, 72)
      WRITE(loutf, '(A)') 'E N D   O F   O U T P U T   C A L C U L A T I O N'
      WRITE(loutf, '(72A1)') ('*', i=1, 72)

      IF ((llogf /= 6) .AND. (llogf /= loutf)) THEN
         CLOSE(llogf)
      END IF

      CLOSE(loutf)

      END SUBROUTINE close_infofiles

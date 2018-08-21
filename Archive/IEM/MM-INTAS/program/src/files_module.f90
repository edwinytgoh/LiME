!#1   subroutine read_filenames
!#2   subroutine open_infofiles



      MODULE files_module

      USE cmn_param, ONLY : fileln

!-----------------------------------------------------------------------
!
!     LOGICAL UNIT NUMBERS
!
!-----------------------------------------------------------------------

      INTEGER(KIND=4), PARAMETER :: linfile   =  5 ! argument file
                                                   ! (standard input)

      !--- input
      INTEGER(KIND=4), PARAMETER :: linpt_in  = 11 ! input file

      !--- output
      INTEGER(KIND=4), PARAMETER :: loutput   = 21 ! output file
      INTEGER(KIND=4), PARAMETER :: llog      = 22 ! logfile

      INTEGER(KIND=4), PARAMETER :: lpdfout   = 31 ! PDF output file

      !--- lookup
      INTEGER(KIND=4), PARAMETER :: lfldat    = 41 ! look-up datafile


!-----------------------------------------------------------------------
!
!     FILE NAMES
!
!-----------------------------------------------------------------------

      !--- input
      CHARACTER(LEN=fileln), SAVE :: inpt_in  = ' ' ! PDFD input file

      !--- output
      CHARACTER(LEN=fileln), SAVE :: output   = ' ' ! output file
      CHARACTER(LEN=fileln), SAVE :: logfile  = ' ' ! logfile

      CHARACTER(LEN=fileln), SAVE :: pdfout   = ' ' ! PDF output file

      !--- lookup
      CHARACTER(LEN=fileln), SAVE :: fldat    = ' ' ! look-up datafile


      CONTAINS


!#1
!======================================================================= 
!     subroutine read_filenames
!
!  function :
!    read names of files from standard input
!
!  variables in the argument list
!
!  i  lout     : unit number of output file
!  i  imon     : monitor flag
!----------------------------------------------------------------------
      SUBROUTINE read_filenames (lout, imon)

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: lout, imon

!-----  local variables
!-----  searched file strings
      INTEGER (KIND=4), PARAMETER :: ifhome     =  1
      INTEGER (KIND=4), PARAMETER :: ifhome2    =  2
      INTEGER (KIND=4), PARAMETER :: ifinpt_in  =  3
      INTEGER (KIND=4), PARAMETER :: ifout      =  4
      INTEGER (KIND=4), PARAMETER :: iflog      =  5
      INTEGER (KIND=4), PARAMETER :: ifpdfout   =  6
      INTEGER (KIND=4), PARAMETER :: iffldat    =  7
      INTEGER (KIND=4), PARAMETER :: nfiles = 7
      INTEGER (KIND=4) :: ifound(nfiles)
      INTEGER (KIND=4) :: i, ilen
      CHARACTER(LEN=20) :: search(nfiles)
      CHARACTER(LEN=fileln) :: filename(0:nfiles), string
      CHARACTER(LEN=30+fileln-1) :: line

      INTEGER(KIND=4) :: hlen, hlen2
      CHARACTER(LEN=fileln) :: home,home2

!-----  MAIN ACTION

!-----  initialization
      DO i = 1, nfiles
         ifound(i) = 0
      END DO

      filename(0) = ' '


!-----  searched file strings
      search(ifhome    ) = 'HOME                '
      search(ifhome2   ) = 'HOME2               '

      !--- input
      search(ifinpt_in)  = 'input file          '

      !--- output
      search(ifout     ) = 'output file         '
      search(iflog     ) = 'log file            '

      search(ifpdfout  ) = 'PDF output file     '

      !--- chemistry lookup files
      search(iffldat   ) = 'FLAME datafile      '

      IF (imon > 0) WRITE(lout,'(/A)') 'Monitor output in READ_FILENAMES'

!-----  read the contents of the whole file (standard input)
      DO
         READ(linfile,'(A)',ERR=9000,END=1000) line
         IF (imon > 0) WRITE(lout,'(A)') line
         DO i = 1, nfiles
            IF (INDEX(line,search(i)) == 1) THEN
               filename(i) = line(30:30+fileln-1)
               ifound(i) = 1
            END IF
         END DO
      END DO

!-----  entry point after end-of-file was encountered
 1000 CONTINUE

!-----  check if all file names were read
      DO i = 1, nfiles
         IF (ifound(i) == 0) THEN
            WRITE(6,'(/A,A,A)') 'READ_FILENAMES: Warning ', &
              search(i),' not found'
         END IF
      END DO

!-----  add HOME / HOME2 string
      hlen  = LEN_TRIM(filename(ifhome))
      hlen2 = LEN_TRIM(filename(ifhome2))
      IF (imon > 1) THEN
         WRITE(lout,'(A,A)') 'HOME =',filename(ifhome)
         WRITE(lout,'(A,A)') 'HOME2=',filename(ifhome2)
      END IF
      IF (hlen > 0) THEN
         home(1:hlen) = filename(ifhome)(1:hlen)
         IF (imon > 1) WRITE(lout,'(A,I3)') 'hlen = ',hlen
      ELSE
         hlen = 1
         home(1:1) = '.'
         IF (imon > 1) WRITE(lout,'(A)') 'HOME reset to .'
      END IF
      IF (hlen2 > 0) THEN
         home2(1:hlen2) = filename(ifhome2)(1:hlen2)
         IF (imon > 1) WRITE(lout,'(A,I3)') 'hlen2 = ',hlen2
      ELSE
         hlen2 = 1
         home2(1:1) = '.'
         IF (imon > 1) WRITE(lout,'(A)') 'HOME2 reset to .'
      END IF
      IF (imon > 1) THEN
         WRITE(lout,'(A,A)') 'HOME =',home
         WRITE(lout,'(A,A)') 'HOME2=',home2
      END IF

      DO i = 1, nfiles
         IF (filename(i)(1:6) == 'HOME2/') THEN
            ilen = LEN_TRIM(filename(i))-5

            IF (ilen < 1) THEN
               WRITE(6,'(A)') 'Error in FILNAM'
               WRITE(6,'(A)') 'Filename too short'
               WRITE(6,'(A)') search(i)//':'//filename(i)
               STOP
            END IF

            IF (hlen2+ilen > LEN(filename(i))) THEN
               WRITE(6,'(A)') 'Error in FILNAM'
               WRITE(6,'(A)') 'Length of filename too short'
               WRITE(6,'(A)') search(i)//':'//filename(i)
               WRITE(6,'(A,I3)') 'Expected length:', ilen+hlen2
               STOP
            END IF

            string(1:ilen) = filename(i)(6:5+ilen)
            filename(i) = ' '  !  empty the filename string
            filename(i)(1:hlen2+ilen) = home2(1:hlen2)//string(1:ilen)
            IF (imon > 1) THEN
               WRITE(lout,'(/A,I3)') 'ilen  =', ilen
               WRITE(lout,'( A,A)')  'string=', string(1:ilen)
            END IF

         ELSE IF (filename(i)(1:5) == 'HOME/') THEN
            ilen = LEN_TRIM(filename(i))-4

            IF (ilen < 1) THEN
               WRITE(6,'(A)') 'Error in FILNAM'
               WRITE(6,'(A)') 'Filename too short'
               WRITE(6,'(A)') search(i)//':'//filename(i)
               STOP
            END IF

            IF (hlen+ilen > LEN(filename(i))) THEN
               WRITE(6,'(A)') 'Error in FILNAM'
               WRITE(6,'(A)') 'Length of filename too short'
               WRITE(6,'(A)') search(i)//':'//filename(i)
               WRITE(6,'(A,I3)') 'Expected length:', ilen+hlen
               STOP
            END IF

            string(1:ilen) = filename(i)(5:4+ilen)
            filename(i) = ' '
            filename(i)(1:hlen+ilen) = home(1:hlen)//string(1:ilen)
            IF (imon > 1) THEN
               WRITE(lout,'(/A,I3)') 'ilen  =', ilen
               WRITE(lout,'( A,A)')  'string=', string(1:ilen)
            END IF

         END IF
      END DO
  

!-----  substitute file names
      inpt_in  = filename(ifinpt_in )
      output   = filename(ifout     )
      logfile  = filename(iflog     )
      pdfout   = filename(ifpdfout  )
      fldat    = filename(iffldat   )

      IF (imon > 0) THEN
         DO i = 1, nfiles
            WRITE(lout,'(A,A,A)') search(i), ' = ', filename(i)
         END DO
      END IF
      
      RETURN

!-----  entry point after an error occurred in reading the input lines
 9000 CONTINUE
      WRITE(6,'(/A,/A)') 'Error in READ_FILENAMES; line = ', TRIM(line)
      STOP

      END SUBROUTINE read_filenames


!#2
!======================================================================= 
!     subroutine open_infofiles
!
!  function :
!    open PDFD output files depending on logical unit numbers
!
!  includes :
!----------------------------------------------------------------------- 
      SUBROUTINE open_infofiles

      IMPLICIT NONE

!-----  local variables
      INTEGER (KIND=4) :: i

!-----  MAIN ACTION

!-----  open files and write headers
      OPEN(UNIT=loutput, &
           FILE=output, &
           FORM='FORMATTED', &
           ACCESS='SEQUENTIAL')

      IF ((llog /= 6) .AND. (llog /= loutput)) THEN
         OPEN(UNIT=llog, &
              FILE=logfile, &
              FORM='FORMATTED', &
              ACCESS='SEQUENTIAL')
      END IF

      WRITE(loutput, '(72A1)') ('*', i=1, 72)
      WRITE(loutput, '(A)') 'C A L C U L A T I O N   O U T P U T'
      WRITE(loutput, '(72A1)') ('*', i=1, 72)

      WRITE(llog, '(72A1)') ('*', i=1, 72)
      WRITE(llog, '(A)') 'L O G   O U T P U T'
      WRITE(llog, '(72A1)') ('*', i=1, 72)

      END SUBROUTINE open_infofiles


      END MODULE files_module

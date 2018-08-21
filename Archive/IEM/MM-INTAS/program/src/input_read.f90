!#0   subroutine read_input
!#1   subroutine input_particles
!#2   subroutine input_iterations
!#3   subroutine input_mixing
!#4   subroutine input_lookup
!#5   subroutine input_pdf_outp


!#0 
!=======================================================================
!        subroutine read_input
!=======================================================================
!  function :
!    read the input file describing the problem
!
!  note :
!    the input file is divided into blocks and sub-blocks
!-----------------------------------------------------------------------
      SUBROUTINE read_input (linput, inputfile, loutf, &
                             icasestudied, nstep, dt)

!-----  used modules
      USE pdfd_cktool_module, ONLY : cksubs, upcase

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: linput, loutf
      CHARACTER(LEN=*), INTENT(IN) :: inputfile 

!----- output variables
      INTEGER (KIND=4), INTENT(INOUT) :: icasestudied
      INTEGER (KIND=4), INTENT(INOUT) :: nstep
      REAL    (KIND=8), INTENT(INOUT) :: dt

!----- local variables
      INTEGER (KIND=4) :: ind, nsub
      INTEGER (KIND=4) :: iblock
      LOGICAL          :: Lerr, Lexist
      CHARACTER(LEN=100) :: line, sub(2)
      CHARACTER(LEN=8)   :: keywrd

!------  define the set of allowed keyword blocks
      INTEGER (KIND=4), PARAMETER :: nblock = 5
      CHARACTER(LEN=8), PARAMETER :: blknam(nblock) = (/ 'PARTICLE', &
                                                         'LOOKUP  ', &
                                                         'MIXING  ', &
                                                         'ITERATIO', &
                                                         'PDF_OUTP' /)
      INTEGER (KIND=4) :: ifoundblk(nblock)


!-----  MAIN ACTION

      WRITE(loutf,'(/A,/A,/A/)') &
         '*********************', &
         'SUBROUTINE READ_INPUT', &
         '*********************'

!-----  open input file
      INQUIRE(FILE=inputfile, EXIST = Lexist)

      IF (.NOT.Lexist) THEN
         WRITE(loutf,'(/,''Error in subroutine READ_INPUT:'')')
         WRITE(loutf,'(''Could not find keyword file'')')
         WRITE(loutf,'(A)') inputfile
         STOP
      END IF

      OPEN(UNIT=linput, &
           FILE=inputfile, &
           FORM='FORMATTED', &
           ACCESS='SEQUENTIAL', &
           STATUS='OLD')

      ifoundblk = -1  ! (iblock=1,nblock)

      DO iblock = 1, nblock
         REWIND(linput)
         DO
            READ (linput,'(A)',ERR=9000,END=6000) line

            !---  look for comment indicator
            ind = INDEX(line,'!')
            IF (ind == 1) CYCLE
            IF (ind > 1) line(ind:) = ' '
            line = upcase(line,100)

            !---  parse line into substrings
            CALL cksubs (line,loutf,4,sub,nsub,Lerr)
            IF (Lerr) THEN
               WRITE(loutf,'(/A,/A)') &
                  ' Error in subroutine READ_INPUT', &
                  ' Input contains too many substrings, line:'
               WRITE(loutf,'(A,A)') '   ', TRIM(line)
               STOP
            END IF

            !---  no substrings; line is empty; read next line
            IF (nsub == 0) CYCLE


            keywrd = sub(1)(1:8)

            !--- the block has not been found yet.
            IF (ifoundblk(iblock) == -1) THEN

               !---  check for a 'BEGIN'
               IF (keywrd(1:5) == 'BEGIN') THEN           

                  !---  look for the blockname
                  IF (nsub < 2) THEN
                     WRITE(loutf,'(/A)') &
                        ' Undefined beginning of block in READ_INPUT, line:'
                     WRITE(loutf,'(A,A)') '   ', TRIM(line)
                     STOP
                  END IF

                  !---  check whether it is the same name as block iblock
                  keywrd = sub(2)(1:8)
                  IF (keywrd(1:8) == blknam(iblock)(1:8)) THEN
                     ifoundblk(iblock) = 1
                  END IF
                  CYCLE

               END IF
            END IF

            !--- the block has been found
            IF (ifoundblk(iblock) == 1) THEN

               !--- end of block?
               IF (keywrd(1:3) == 'END') EXIT

               SELECT CASE (blknam(iblock))
                  CASE ('PARTICLE')
                     CALL input_particles (loutf, line, sub, nsub)
                  CASE ('LOOKUP  ')
                     CALL input_lookup (loutf, line, sub, nsub, &
                                        icasestudied)
                  CASE ('MIXING  ')
                     CALL input_mixing (loutf, line, sub, nsub)
                  CASE ('ITERATIO')
                     CALL input_iterations (loutf, line, sub, nsub, &
                                            nstep, dt)
                  CASE ('PDF_OUTP')
                     CALL input_pdf_outp (loutf, line, sub, nsub)
               END SELECT

            END IF
         END DO

!--- end of file
6000     CONTINUE

         IF (ifoundblk(iblock) == -1) THEN
            WRITE(*,*) 'Warning: BLOCK ', blknam(iblock), ' not read.'
         END IF

      END DO

      CLOSE(linput)

      RETURN

!-----  error messages
 9000 WRITE(loutf,'(/A)') &
         'File error when reading input data in READ_INPUT'
      STOP

      END SUBROUTINE read_input


!#1
!=======================================================================
!     subroutine INPUT_PARTICLES
!=======================================================================
      SUBROUTINE input_particles(loutf, line, sub, nsub)

!-----  used modules
      USE cmn_gaspdf, ONLY : initscal, mean_mixf,var_mixf,min_mixf,max_mixf, np
      USE cmn_mixing, ONLY : omega

      USE pdfd_cktool_module, ONLY : ckxnum

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: nsub
      CHARACTER(LEN=*), INTENT(IN) :: line, sub(nsub)


!-----  local variables
      INTEGER (KIND=4) :: nval, n
      REAL    (KIND=8) :: rval(3)
      LOGICAL          :: Lerr
      CHARACTER(LEN=8) :: keywrd

!-----  MAIN ACTION

      keywrd = sub(1)(1:8)

!-----  check for single keywords

      !---  number of particles
      IF (keywrd(1:8) == 'INITSCAL') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        initscal = INT(rval(1))
        WRITE(loutf, 5000) keywrd, ' = ', initscal

      !---  mean mixture fraction
      ELSE IF (keywrd(1:8) == 'MEAN_MIX') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        mean_mixf = rval(1)
        WRITE(loutf, 5500) keywrd, ' = ', mean_mixf

      !---  mixture fraction variance
      ELSE IF (keywrd(1:8) == 'VAR_MIXF') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        var_mixf = rval(1)
        WRITE(loutf, 5500) keywrd, ' = ', var_mixf

      !---  1st peak mean mixture fraction
      ELSE IF (keywrd(1:8) == 'MIN_MIXF') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        min_mixf = rval(1)
        WRITE(loutf, 5500) keywrd, ' = ', min_mixf

      !---  2nd peak mean mixture fraction
      ELSE IF (keywrd(1:8) == 'MAX_MIXF') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        max_mixf = rval(1)
        WRITE(loutf, 5500) keywrd, ' = ', max_mixf

      !---  number of particles
      ELSE IF (keywrd(1:2) == 'NP') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        np = INT(rval(1))
        WRITE(loutf, 5000) keywrd, ' = ', np

      !---  turbulence frequency
      ELSE IF (keywrd(1:5) == 'OMEGA') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        omega = rval(1)
        WRITE(loutf, 5500) keywrd, ' = ', omega
 

!-------  illegal keyword
      ELSE
        WRITE(loutf,*) 'Error in INPUT_PARTICLES', &
                       'Unknown keyword ',keywrd,' in line:'
        WRITE(loutf,*) line
        STOP 'INPUT_PARTICLES: unknown keyword'
      END IF

      RETURN

!-----  formats
 5000 FORMAT(4X,A8,A3,I5)
 5500 FORMAT(4X,A8,A3,F8.3)

!-----  error messages
 8000 WRITE(loutf,'(/A,/A,A)') &
         'Error in subroutine CKXNUM', &
         'line  =',TRIM(line)
      DO n = 1, nsub
         WRITE(loutf,'(A,I1,A,A)') 'sub(',n,')=',TRIM(sub(n))
      END DO
      STOP

 8500 WRITE(loutf,'(/A,/A)') &
         ' Error in INPUT_PARTICLES', &
         ' Too many words in line:'
      WRITE(loutf,'(A,A)') '   ', TRIM(line)
      STOP

 8600 WRITE(loutf,'(/A,/A)') &
         ' Error in INPUT_PARTICLES', &
         ' Not enough words in line:'
      WRITE(loutf,'(A,A)') '   ', TRIM(line)
      STOP

      END SUBROUTINE input_particles


!#2
!=======================================================================
!     subroutine INPUT_ITERATIONS
!=======================================================================
      SUBROUTINE input_iterations (loutf, line, sub, nsub, &
     &                             nstep, dt)

!-----  used modules
      USE pdfd_cktool_module, ONLY : ckxnum

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: nsub
      CHARACTER(LEN=*), INTENT(IN) :: line, sub(nsub)

!-----  input/output variables
      INTEGER (KIND=4), INTENT(INOUT) :: nstep
      REAL    (KIND=8), INTENT(INOUT) :: dt

!-----  local variables
      INTEGER (KIND=4) :: nval, n
      REAL    (KIND=8) :: rval(3)
      LOGICAL          :: Lerr
      CHARACTER(LEN=8) :: keywrd


!-----  MAIN ACTION

      keywrd = sub(1)(1:8)

!-----  check for keywords
      IF (keywrd(1:5) == 'NITER') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        nstep = INT(rval(1))
        WRITE(loutf, 5000) keywrd, ' = ', nstep

      ELSE IF (keywrd(1:2) == 'DT') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        dt = rval(1)
        WRITE(loutf, 5100) keywrd, ' = ', dt

!---------  illegal keyword
      ELSE
        WRITE(loutf,*) 'Error in INPUT_ITERATIONS', &
                       'Unknown keyword ',keywrd,' in line:'
        WRITE(loutf,*) line
        STOP 'INPUT_ITERATIONS: unknown keyword'
      END IF

      RETURN

!-------  formats
 5000 FORMAT(4X,A8,A3,I5)
 5100 FORMAT(4X,A8,A3,E12.4)

!-------  error messages
 8000 WRITE(loutf,'(/A,/A,A)') &
         'Error in subroutine CKXNUM', &
         'line  =',TRIM(line)
      DO n = 1, nsub
         WRITE(loutf,'(A,I1,A,A)') 'sub(',n,')=',TRIM(sub(n))
      END DO
      STOP

 8500 WRITE(loutf,'(/A,/A)') &
         ' Error in INPUT_ITERATIONS', &
         ' Too many words in line:'
      WRITE(loutf,'(A,A)') '   ', TRIM(line)
      STOP

      END SUBROUTINE input_iterations


!#3
!=======================================================================
!     subroutine INPUT_MIXING
!=======================================================================
      SUBROUTINE input_mixing(loutf, line, sub, nsub)

!-----  used modules
      USE cmn_mixing, ONLY : mixmod, mixext, cphi

      USE pdfd_cktool_module, ONLY : ckxnum

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: nsub
      CHARACTER(LEN=*), INTENT(IN) :: line, sub(nsub)

!-----  local variables
      INTEGER (KIND=4) :: nval, n
      REAL    (KIND=8) :: rval(3)
      LOGICAL          :: Lerr
      CHARACTER(LEN=8) :: keywrd

!-----  MAIN ACTION

      keywrd = sub(1)(1:8)

!-----  check for single keywords
      IF (keywrd(1:6) == 'MIXMOD') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        mixmod = INT(rval(1))
        WRITE(loutf, 5000) keywrd, ' = ', mixmod

      ELSE IF (keywrd(1:6) == 'MIXEXT') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        mixext = INT(rval(1))
        WRITE(loutf, 5000) keywrd, ' = ', mixext

      ELSE IF (keywrd(1:4) == 'CPHI') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        cphi = rval(1)
        WRITE(loutf, 5500) keywrd, ' = ', cphi

!---------  illegal keyword
      ELSE
        WRITE(loutf,*) 'Error in INPUT_MIXING', &
                       'Unknown keyword ',keywrd,' in line:'
        WRITE(loutf,*) line
        STOP 'INPUT_MIXING: unknown keyword'
      END IF

      RETURN

!-------  formats
 5000 FORMAT(4X,A8,A3,I5)
 5500 FORMAT(4X,A8,A3,F8.3)

!-------  error messages
 8000 WRITE(loutf,'(/A,/A,A)') &
         'Error in subroutine CKXNUM', &
         'line  =',TRIM(line)
      DO n = 1, nsub
         WRITE(loutf,'(A,I1,A,A)') 'sub(',n,')=',TRIM(sub(n))
      END DO
      STOP

 8500 WRITE(loutf,'(/A,/A)') &
         ' Error in INPUT_MIXING', &
         ' Too many words in line:'
      WRITE(loutf,'(A,A)') '   ', TRIM(line)
      STOP

 8600 WRITE(loutf,'(/A,/A)') &
         ' Error in INPUT_MIXING', &
         ' Not enough words in line:'
      WRITE(loutf,'(A,A)') '   ', TRIM(line)
      STOP

      END SUBROUTINE input_mixing


!#4
!=======================================================================
!     subroutine INPUT_LOOKUP
!=======================================================================
      SUBROUTINE input_lookup(loutf, line, sub, nsub, icasestudied)

!-----  used modules
      USE cmn_lookup, ONLY : ilut

      USE pdfd_cktool_module, ONLY : ckxnum

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: nsub
      CHARACTER(LEN=*), INTENT(IN) :: line, sub(nsub)

!-----  input/output variables
      INTEGER (KIND=4), INTENT(INOUT) :: icasestudied

!-----  local variables
      INTEGER (KIND=4) :: nval, n
      REAL    (KIND=8) :: rval(3)
      LOGICAL          :: Lerr
      CHARACTER(LEN=8) :: keywrd


!-----  MAIN ACTION

      keywrd = sub(1)(1:8)

!-----  check for single keywords

      !---  case considered (Freon/air, c2h4/air or ch4/h2)
      IF (keywrd(1:4) == 'CASE') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        icasestudied = INT(rval(1))
        WRITE(loutf, 5000) keywrd, ' = ', icasestudied

      !---  look-up performed every ilut iteration
      ELSE IF (keywrd(1:4) == 'ILUT') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        ilut = INT(rval(1))
        WRITE(loutf, 5000) keywrd, ' = ', ilut


!-------  illegal keyword
      ELSE
        WRITE(loutf,*) 'Error in INPUT_LOOKUP', &
                       'Unknown keyword ',keywrd,' in line:'
        WRITE(loutf,*) line
        STOP 'INPUT_LOOKUP: unknown keyword'
      END IF


      RETURN

!-----  formats
 5000 FORMAT(4X,A8,A3,I5)

!-----  error messages
 8000 WRITE(loutf,'(/A,/A,A)') &
         'Error in subroutine CKXNUM', &
         'line  =',TRIM(line)
      DO n = 1, nsub
         WRITE(loutf,'(A,I1,A,A)') 'sub(',n,')=',TRIM(sub(n))
      END DO
      STOP

 8500 WRITE(loutf,'(/A,/A)') &
         ' Error in INPUT_LOOKUP', &
         ' Too many words in line:'
      WRITE(loutf,'(A,A)') '   ', TRIM(line)
      STOP

 8600 WRITE(loutf,'(/A,/A)') &
         ' Error in INPUT_LOOKUP', &
         ' Not enough words in line:'
      WRITE(loutf,'(A,A)') '   ', TRIM(line)
      STOP

      END SUBROUTINE input_lookup


!#5
!=======================================================================
!     subroutine INPUT_PDF_OUTP
!=======================================================================
      SUBROUTINE input_pdf_outp(loutf, line, sub, nsub)

!-----  used modules
      USE cmn_pdfoutp, ONLY : ipdfop

      USE pdfd_cktool_module, ONLY : ckxnum

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: nsub
      CHARACTER(LEN=*), INTENT(IN) :: line, sub(nsub)

!-----  local variables
      INTEGER (KIND=4) :: nval, n
      REAL    (KIND=8) :: rval(3)
      LOGICAL          :: Lerr
      CHARACTER(LEN=8) :: keywrd

!-----  MAIN ACTION

      keywrd = sub(1)(1:8)

!-----  check for keywords
      IF (keywrd(1:6) == 'IPDFOP') THEN
        IF (nsub > 2) GOTO 8500
        CALL ckxnum(sub(2),1,loutf,nval,rval(1),Lerr)
        IF (Lerr) GOTO 8000
        ipdfop = INT(rval(1))
        WRITE(loutf, 5000) keywrd, ' = ', ipdfop


!---------  illegal keyword
      ELSE
        WRITE(loutf,'(/A,/A,A,A)') &
           ' Error in INPUT_PDF_OUTP', &
           ' Unknown keyword ',TRIM(keywrd),', in line:'
        WRITE(loutf,'(A,A)') '   ', TRIM(line)
        STOP
      END IF

      RETURN

!-------  formats
 5000 FORMAT(4X,A8,A3,I5)

!-------  error messages
 8000 WRITE(loutf,'(/A,/A,A)') &
         'Error in subroutine CKXNUM', &
         'line  =',TRIM(line)
      DO n = 1,nsub
         WRITE(loutf,'(A,I1,A,A)') 'sub(',n,')=',TRIM(sub(n))
      END DO
      STOP

 8500 WRITE(loutf,'(/A,/A)') &
         ' Error in INPUT_PDF_OUTP', &
         ' Too many words in line:'
      WRITE(loutf,'(A,A)') '   ', TRIM(line)
      STOP

      END SUBROUTINE input_pdf_outp

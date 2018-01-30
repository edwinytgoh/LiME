!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! ONLY THE SUBROUTINES USED IN PDFD !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!#1 subroutine ipparr
!#2 function ifirch
!#3 function ilasch
!#4 subroutine ckcomp
!#5 function ipplen
!#6 function upcase
!#7 subroutine cksubs
!#8 subroutine ckxnum
!#9 subroutine upper_case


   MODULE pdfd_cktool_module

   CHARACTER(LEN=1), PARAMETER, PRIVATE :: LCASE(26) = &
      &  (/'a','b','c','d','e','f','g','h','i','j','k','l','m', &
      &    'n','o','p','q','r','s','t','u','v','w','x','y','z'/)

   CHARACTER(LEN=1), PARAMETER, PRIVATE :: UCASE(26) = &
      &  (/'A','B','C','D','E','F','G','H','I','J','K','L','M', &
      &    'N','O','P','Q','R','S','T','U','V','W','X','Y','Z'/)

   CONTAINS


!#1
!***********************************************************************
   SUBROUTINE IPPARR (STRING, ICARD, NEXPEC, RVAL, NFOUND, IERR, LOUT)
!***********************************************************************
!  IPPARR may be used for parsing an input record that contains real
!  values, but was read into a character variable instead of directly
!  into real variables.
!  The following benefits are gained by this approach:
!    - specification of only certain elements of the array is allowed,
!      thus letting the others retain default values
!    - variable numbers of values may be input in a record, up to a
!      specified maximum
!    - control remains with the calling program in case of an input error
!    - diagnostics may be printed by IPPARR to indicate the nature
!      of input errors
!
!  The contents of STRING on input indicate which elements of RVAL
!  are to be changed from their entry values, and values to which
!  they should be changed on exit.  Commas and blanks serve as
!  delimiters, but multiple blanks are treated as a single delimeter.
!  Thus, an input record such as:
!    '   1.,   2,,4.e-5   , ,6.e-6'
!  is interpreted as the following set of instructions by IPGETR:
!
!    (1) set RVAL(1) = 1.0
!    (2) set RVAL(2) = 2.0
!    (3) leave RVAL(3) unchanged
!    (4) set RVAL(4) = 4.0E-05
!    (5) leave RVAL(5) unchanged
!    (6) set RVAL(6) = 6.0E-06
!
!  IPPARR will print diagnostics on the default output device, if desired.
!
!  IPPARR is part of IOPAK, and is written in ANSI FORTRAN 77
!
!  Examples:
!
!    Assume RVAL = (0., 0., 0.) and NEXPEC = 3 on entry:
!
!   input string           RVAL on exit            IERR    NFOUND
!   -------------          ----------------------  ----    ------
!  '  2.34e-3,  3 45.1'    (2.34E-03, 3.0, 45.1)     0       3
!  '2,,3.-5'               (2.0, 0.0, 3.0E-05)       0       3
!  ',1.4,0.028E4'          (0.0, 1.4, 280.0)         0       3
!  '1.0, 2.a4, 3.0'        (1.0, 0.0, 0.0)           1       1
!  '1.0'                   (1.0, 0.0, 0.0)           2       1
!
!  Assume RVAL = (0.,0.,0.,0.) and NEXPEC = -4 on entry:
!
!   input string           RVAL on exit            IERR    NFOUND
!   -------------          ----------------------  ----    ------
!  '1.,2.'                 (1.0, 2.0)                0       2
!  ',,3  4.0'              (0.0, 0.0, 3.0, 4.0)      0       4
!  '1,,3,,5.0'             (0.0, 0.0, 3.0, 0.0)      3       4
!
!  arguments: (I=input,O=output)
!  -----------------------------
!  STRING (I) - the character string to be parsed.
!
!  ICARD  (I) - data statement number, and error processing flag
!         < 0 : no error messages printed
!         = 0 : print error messages, but not ICARD
!         > 0 : print error messages, and ICARD
!
!  NEXPEC (I) - number of real variables expected to be input.  If
!         < 0, the number is unknown, and any number of values
!         between 0 and abs(nexpec) may be input.  (see NFOUND)
!
!  RVAL (I,O) - the real value or values to be modified.  On entry,
!       the values are printed as defaults.  The formal parameter
!       corresponding to RVAL must be dimensioned at least NEXPEC
!       in the calling program if NEXPEC > 1.
!
!  NFOUND (O) - the number of real values represented in STRING,
!         only in the case that there were as many or less than
!         NEXPEC.
!
!  IERR (O) - error flag:
!       = 0 if no errors found
!       = 1 syntax errors or illegal values found
!       = 2 for too few values found (NFOUND < NEXPEC)
!       = 3 for too many values found (NFOUND > NEXPEC)
!
!  LOUT (I) - unit number of output file
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)    :: ICARD, NEXPEC, LOUT
      INTEGER, INTENT(OUT)   :: NFOUND, IERR
      REAL(KIND=8) , INTENT(INOUT) :: RVAL(*)
      CHARACTER(LEN=*), INTENT(IN)  :: STRING

      INTEGER :: NEXP,IE,NC,IBS,IES,NCH
      LOGICAL :: OKINCR
      CHARACTER(LEN=8), PARAMETER :: &
     &   FMT(22) = (/' (E1.0)',' (E2.0)',' (E3.0)',' (E4.0)', &
     &               ' (E5.0)',' (E6.0)',' (E7.0)',' (E8.0)',' (E9.0)', &
     &               '(E10.0)','(E11.0)','(E12.0)','(E13.0)','(E14.0)', &
     &               '(E15.0)','(E16.0)','(E17.0)','(E18.0)','(E19.0)', &
     &               '(E20.0)','(E21.0)','(E22.0)'/)
      CHARACTER(LEN=80) :: ITEMP

      IERR   = 0
      NFOUND = 0
      NEXP = IABS(NEXPEC)
      IE = ILASCH(STRING)
      IF (IE == 0) GOTO 500
      NC = 1

!-----
!  OKINCR is a flag that indicates it's OK to increment NFOUND,
!  the index of the array into which the value should be read.
!  It is set negative when a space follows a real value substring,
!  to keep incrementing from occurring if a comma should be 
!  encountered before the next value.
!------------------------------------
      OKINCR = .TRUE.

!-----  begin overall loop on characters in string
  100 CONTINUE

      IF (STRING(NC:NC) == ',') THEN
        IF (OKINCR) THEN
          NFOUND = NFOUND + 1
        ELSE
          OKINCR = .TRUE.
        END IF
        GOTO 450
      END IF
      IF (STRING(NC:NC) == ' ') GOTO 450

!-----  first good character (non-delimeter) found
!-----  now find last good character
      IBS = NC
      DO
        NC = NC + 1
        IF (NC > IE) EXIT
        SELECT CASE (STRING(NC:NC))
          CASE (' ')
            OKINCR = .FALSE.
            EXIT
          CASE (',')
            OKINCR = .TRUE.
            EXIT
          CASE DEFAULT
            CYCLE
        END SELECT
      END DO

!-----  end of substring found - read value into real array
      NFOUND = NFOUND + 1
      IF (NFOUND > NEXP) THEN
        IERR = 3
        GOTO 500
      END IF

      IES = NC - 1
      NCH = IES - IBS + 1
      ITEMP = ' '
      ITEMP = STRING(IBS:IES)
      READ (ITEMP(1:NCH), FMT(NCH), ERR=400) RVAL(NFOUND)
      GOTO 450
  400 CONTINUE
      WRITE(LOUT,'(A)') STRING(IBS:IES)
      IERR = 1
      GOTO 510
  450 CONTINUE
      NC = NC + 1
      IF (NC <= IE) GOTO 100

  500 CONTINUE
      IF (NEXPEC > 0 .AND. NFOUND < NEXP) IERR = 2
  510 CONTINUE

!-----  normal return
      IF (IERR == 0 .OR. ICARD < 0) RETURN

!-----  error messages
      IF (ICARD /= 0) WRITE (LOUT,'(A,I3)') &
     &  '! Error in data statement number ',ICARD
      IF (IERR == 1) &
     &  WRITE (LOUT,'(A)') 'Syntax error, or illegal value'
      IF (IERR == 2) WRITE(LOUT,'(A,I2,A,I2)') &
     &  ' Too few data items.  Number found = ',NFOUND, &
     &  '  Number expected = ',NEXPEC
      IF (IERR == 3) WRITE(LOUT,'(A,I2)') &
     &  ' Too many data items.  Number expected = ',NEXPEC

   END SUBROUTINE IPPARR


!#2
!***********************************************************************
   FUNCTION IFIRCH (STRING)
!***********************************************************************
!  IFIRCH locates the first non-blank character in a string of
!  arbitrary length.  If no characters are found, IFIRCH is set = 0.
!  When used with the companion routine ILASCH, the length of a string
!  can be determined, and/or a concatenated substring containing the
!  significant characters produced.
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IFIRCH
      INTEGER :: I
      CHARACTER(LEN=*) :: STRING

      IFIRCH = 0

      DO I = 1, LEN(STRING)
        IF (STRING(I:I) /= ' ') THEN
          IFIRCH = I
          EXIT
        END IF
      END DO

   END FUNCTION IFIRCH


!#3
!***********************************************************************
   FUNCTION ILASCH(STRING)
!***********************************************************************
!  ILASCH locates the last non-blank character in a string of
!  arbitrary length.  If no characters are found, ILASCH is set = 0.
!  When used with the companion routine IFIRCH, the length of a string
!  can be determined, and/or a concatenated substring containing the
!  significant characters produced.
!  Note that the FORTRAN intrinsic function LEN returns the length
!  of a character string as declared, rather than as filled.  The
!  declared length includes leading and trailing blanks, and thus is
!  not useful in generating 'significant' substrings.
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: ILASCH
      INTEGER :: I
      CHARACTER(LEN=*) :: STRING

      ILASCH = 0

      DO I = LEN(STRING), 1, -1
        IF (STRING(I:I) /= ' ') THEN
          ILASCH = I
          EXIT
        END IF
      END DO

   END FUNCTION ILASCH


!#4
!***********************************************************************
   SUBROUTINE CKCOMP (IST, IRAY, II, I)
!***********************************************************************
!  Returns the index of an element of a reference character
!  string array which corresponds to a character string;
!  leading and trailing blanks are ignored.
!-----------------------------------------------------------------------
!  INPUT
!    IST   - CHAR*(*) character string.
!    IRAY  - An array of CHAR*(*) character strings;
!            dimension IRAY(*) at least II
!    II    - The length of IRAY.
!-----------------------------------------------------------------------
!  OUTPUT
!    I     - The first integer location in IRAY in which IST
!            corresponds to IRAY(I); if IST is not also an
!            entry in IRAY, I=0.
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: II
      CHARACTER(LEN=*), INTENT(IN) :: IST, IRAY(*)
      INTEGER, INTENT(OUT) :: I

      INTEGER :: N, IS1, IS2, IR1, IR2

      I = 0
      IS1 = IFIRCH(IST)
      IS2 = ILASCH(IST)
      IF (IS2 >= IS1 .AND. IS2 > 0) THEN
        DO N = 1, II
          IR1 = IFIRCH(IRAY(N))
          IR2 = ILASCH(IRAY(N))
          IF (IR2 >= IR1 .AND. IR2 > 0 .AND. &
     &        IST(IS1:IS2) == IRAY(N)(IR1:IR2)) THEN
            I = N
            EXIT
          END IF
        END DO
      END IF

   END SUBROUTINE CKCOMP


!#5
!***********************************************************************
   FUNCTION IPPLEN (LINE)
!***********************************************************************
!  Returns the effective length of a character string, i.e.,
!  the index of the last significant character (NOT a blank!)
!  before an exclamation mark (!) indicating a comment.
!  Trailing blanks are therefore ignored, but leading blanks are retained.
!-----------------------------------------------------------------------
!  INPUT
!     LINE  - A character string.
!-----------------------------------------------------------------------
!  OUTPUT
!     IPPLEN - The effective length of the character string.
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IPPLEN
      INTEGER :: IN
      CHARACTER(LEN=*) :: LINE

      IN = IFIRCH(LINE)
      IF (IN == 0) THEN
        IPPLEN = 0
      ELSE IF (LINE(IN:IN) == '!') THEN
        IPPLEN = 0
      ELSE
        IN = INDEX(LINE, '!')
        IF (IN == 0) THEN
          IPPLEN = ILASCH(LINE)
        ELSE
          IPPLEN = ILASCH(LINE(:IN-1))
        END IF
      END IF

   END FUNCTION IPPLEN


!#6
!***********************************************************************
   FUNCTION UPCASE (STRING, ILEN)
!***********************************************************************

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: STRING
      CHARACTER(LEN=LEN(STRING))   :: UPCASE
      INTEGER, INTENT(IN) :: ILEN

      INTEGER :: JJ, J, N

      UPCASE = ' '
      UPCASE = STRING(:ILEN)
      JJ = MIN(LEN(UPCASE), LEN(STRING), ILEN)
      DO J = 1, JJ
         DO N = 1, 26
            IF (STRING(J:J) == LCASE(N)) UPCASE(J:J) = UCASE(N)
         END DO
      END DO

   END FUNCTION UPCASE


!#7
!***********************************************************************
   SUBROUTINE CKSUBS (LINE, LOUT, NDIM, SUB, NFOUND, KERR)
!***********************************************************************
!  Returns an array of substrings in a character string with blanks
!  as the delimiter
!-----------------------------------------------------------------------
!  input
!     LINE   - CHAR*(*) character string.
!     LOUT   - Output unit for printed diagnostics.
!     NDIM   - Dimension of array SUB(*)*(*)
!-----------------------------------------------------------------------
!  output
!    SUB    - The CHAR*(*) character substrings of LINE.
!             Dimension SUB(*) at least NDIM
!    NFOUND - Number of substrings found in LINE.
!    KERR   - Error flag; dimensioning errors will result in KERR = .TRUE.
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: LOUT, NDIM
      INTEGER, INTENT(OUT) :: NFOUND
      LOGICAL, INTENT(OUT) :: KERR
      CHARACTER(LEN=*), INTENT(IN)  :: LINE
      CHARACTER(LEN=*), INTENT(OUT) :: SUB(*)

      INTEGER :: ILEN, IEND, ISTART, L

      NFOUND = 0
      ILEN = LEN(SUB(1))

      IEND = 0
      KERR = .FALSE.

      OUTERLOOP : DO
        ISTART = IEND + 1
        DO L = ISTART, IPPLEN(LINE)
          IF (LINE(L:L) /= ' ') THEN
            IEND = INDEX(LINE(L:), ' ')
            IF (IEND == 0) THEN
              IEND = IPPLEN(LINE)
            ELSE
              IEND = L + IEND - 1
            END IF
            IF (IEND-L+1 > ILEN) THEN
              WRITE (LOUT,*) ' Error in CKSUBS...substring too long'
              KERR = .TRUE.
            ELSE IF (NFOUND+1 > NDIM) THEN
              WRITE (LOUT,*) ' Error in CKSUBS...NDIM too small'
              KERR = .TRUE.
            ELSE
              NFOUND = NFOUND + 1
              SUB(NFOUND) = LINE(L:IEND)
            END IF
            CYCLE OUTERLOOP
          END IF
        END DO
        EXIT OUTERLOOP
      END DO OUTERLOOP

   END SUBROUTINE CKSUBS


!#8
!***********************************************************************
   SUBROUTINE CKXNUM (LINE, NEXP, LOUT, NVAL, RVAL, KERR)
!***********************************************************************
!  This subroutine is called to parse a character string, LINE,
!  that is composed of several blank-delimited substrings.
!  Each substring is expected to represent a number, which
!  is converted to entries in the array of real numbers, RVAL(*).
!  NEXP is the number of values expected, and NVAL is the
!  number of values found.  This allows format-free input of
!  numerical data.  For example:
!
!     input:  LINE    = " 0.170E+14 0 47780.0"
!             NEXP    = 3, the number of values requested
!             LOUT    = 6, a logical unit number on which to write
!                       diagnostic messages.
!     output: NVAL    = 3, the number of values found
!             RVAL(*) = 1.700E+13, 0.000E+00, 4.778E+04
!             KERR    = .FALSE.
!-----------------------------------------------------------------------
!  input
!     LINE   - A CHAR*(*) character string.
!     NEXP   - Number of real values to be found in character string.
!     LOUT   - Output unit for printed diagnostics.
!-----------------------------------------------------------------------
!  OUTPUT
!     NVAL   - Number of real values found in character string.
!     RVAL   - Array of real values found.
!              Dimension RVAL(*) at least NEXP
!     KERR   - Error flag;  syntax or dimensioning error results
!              in KERR = .TRUE.
!-----------------------------------------------------------------------

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: LINE
      INTEGER, INTENT(IN)  :: LOUT, NEXP
      INTEGER, INTENT(OUT) :: NVAL
      REAL(KIND=8) , INTENT(OUT) :: RVAL(*)
      LOGICAL, INTENT(OUT) :: KERR

      INTEGER :: ILEN,IERR,N
      REAL(KIND=8)  :: RTEMP(80)
      CHARACTER(LEN=80) :: ITEMP

!-----  find comment string (! signifies comment)
      ILEN = IPPLEN(LINE)
      NVAL = 0
      KERR = .FALSE.

      IF (ILEN <= 0) RETURN
      IF (ILEN > 80) THEN
        WRITE (LOUT,*) ' Error in CKXNUM...line length > 80 '
        WRITE (LOUT,'(A)') LINE
        KERR = .TRUE.
        RETURN
      END IF

      ITEMP = LINE(:ILEN)
      IF (NEXP < 0) THEN
        CALL IPPARR (ITEMP, -1, NEXP, RTEMP, NVAL, IERR, LOUT)
      ELSE
        CALL IPPARR (ITEMP, -1, -NEXP, RTEMP, NVAL, IERR, LOUT)
        IF (IERR == 1) THEN
          WRITE(LOUT,*) ' Syntax errors in CKXNUM...'
          WRITE(LOUT,'(A)') LINE
          KERR = .TRUE.
        ELSE IF (NVAL /= NEXP) THEN
          WRITE(LOUT,*) ' Error in CKXNUM...'
          WRITE(LOUT,'(A)') LINE
          KERR = .TRUE.
          WRITE(LOUT,*) NEXP,' values expected;'
          WRITE(LOUT,*) NVAL,' values found.'
        END IF
      END IF
      IF (NVAL <= ABS(NEXP)) THEN
        DO N = 1, NVAL
          RVAL(N) = RTEMP(N)
        END DO
      END IF

   END SUBROUTINE CKXNUM


!#9
!***********************************************************************
      SUBROUTINE UPPER_CASE (string, ilen)
!***********************************************************************
!  This subroutine converts all lowercase letters in the input STRING
!  to uppercase.
!  If ILEN is specified, then only the first ILEN letters are converted.
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER (KIND=4), INTENT(IN), OPTIONAL :: ilen
      CHARACTER(LEN=*), INTENT(INOUT) :: string

      INTEGER (KIND=4) :: jj, j, n

      jj = LEN_TRIM(string)
      IF (PRESENT(ilen)) THEN
         JJ = MIN(jj, ilen)
      END IF

      DO j = 1, jj
         DO n = 1, 26
            IF (string(j:j) == lcase(n)) string(j:j) = ucase(n)
         END DO
      END DO

      END SUBROUTINE UPPER_CASE

   END MODULE pdfd_cktool_module

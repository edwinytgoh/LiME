MODULE splitidnt
!-----------------------------------------------------------------------
!  Identity of f90split utility
!-----------------------------------------------------------------------
      CHARACTER (LEN=*), PARAMETER :: zsccs = &
"@(#)splitidnt.f90	1.4	95/04/04 Michel Olagnon, Phil Garnatz"
      CHARACTER (LEN=*), PARAMETER :: zvers = &
"@(#) splitidnt.f90	V-1.0 95/04/04 Michel Olagnon, Phil Garnatz"
      CHARACTER (LEN=*), PARAMETER :: zusg = &
"( usage: f90split < largefile [ > list_file ] )"
      CHARACTER (LEN=*), PARAMETER :: zhlp  = '( &
&"Fortran 90 utility to split free source form code into"/&
&"as many files as there are procedures. Contained procedures"/&
&"are stored within their container"/&
&"_____________________________________________________________________"/&
&"All rights to this code waived, so that it may be freely distributed"/&
&"as public domain software subject to the condition that these 6 lines"/&
&"are verbatim reproduced. Originally written by Michel Olagnon, from"/&
&"Ifremer, France, who would be pleased to receive your comments and"/&
&"corrections. M. Olagnon (Michel.Olagnon@ifremer.fr) Improved by"/&
&"Phil Garnatz, Cray Research Inc. for makefile generation"/&
&"_____________________________________________________________________"/&
&"                    version 1.0 of 4th April 1995"/&
&"  Split standard input stream, containing source of several fortran90"/&
&"  program units into individual files, each containing a single"/&
&"  program unit, named after it, or main0001.f90-main9999.f90, or"/&
&"  bdta0001.f90-bdta9999.f90. If a file with that name already exists,"/&
&"  it is put in dupl0001.f90-dupl9999.f90."/&
&"  Lists on stdout the use and include dependencies"/&
&"_____________________________________________________________________"/&
&"Note: If you do not like code to start in column 7, remember that,"/&
&"      had Diophantes left a 6 characters margin, then mathematicians"/&
&"      might have spared much efforts on A**N = B**N + C**N ..."/&
&"      My margin is wide to let you put your corrections there."/&
&"____________________________________________________________________")'
!
END MODULE splitidnt

MODULE splitprms
!-----------------------------------------------------------------------
!  Parameters for f90split utility
!-----------------------------------------------------------------------
      CHARACTER (LEN=26), PARAMETER :: zlwc="abcdefghijklmnopqrstuvwxyz"
      CHARACTER (LEN=26), PARAMETER :: zupc="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      CHARACTER (LEN=10), PARAMETER :: zdgt="1234567890"
      CHARACTER (LEN=1) , PARAMETER :: ztab=char(9)
      INTEGER, PARAMETER              :: lend = 3
      CHARACTER (LEN=lend), PARAMETER :: zend = "END"
      INTEGER, PARAMETER              :: lctn = 8
      CHARACTER (LEN=lctn), PARAMETER :: zctn = "CONTAINS"
      INTEGER, PARAMETER              :: lntf = 9
      CHARACTER (LEN=lntf), PARAMETER :: zntf = "INTERFACE"
      INTEGER, PARAMETER              :: lsub = 10
      CHARACTER (LEN=lsub), PARAMETER :: zsub = "SUBROUTINE"
      INTEGER, PARAMETER              :: lpgm = 7
      CHARACTER (LEN=lpgm), PARAMETER :: zpgm = "PROGRAM"
      INTEGER, PARAMETER              :: lmdl = 6
      CHARACTER (LEN=lmdl), PARAMETER :: zmdl = "MODULE"
      INTEGER, PARAMETER              :: lfun = 8
      CHARACTER (LEN=lfun), PARAMETER :: zfun = "FUNCTION"
      INTEGER, PARAMETER              :: lbdt = 9
      CHARACTER (LEN=lbdt), PARAMETER :: zbdt = "BLOCKDATA"
      INTEGER, PARAMETER               :: lbdt1 = 5
      CHARACTER (LEN=lbdt1), PARAMETER :: zbdt1 = "BLOCK"
      INTEGER, PARAMETER               :: lbdt2 = 4
      CHARACTER (LEN=lbdt2), PARAMETER :: zbdt2 = "DATA"
      INTEGER, PARAMETER              :: luse = 3
      CHARACTER (LEN=luse), PARAMETER :: zuse = "USE"
      INTEGER, PARAMETER              :: linc = 7
      CHARACTER (LEN=linc), PARAMETER :: zinc = "INCLUDE"
!
      CHARACTER (LEN=*), PARAMETER :: zbasp = "main0000"
      CHARACTER (LEN=*), PARAMETER :: zbasb = "bdta0000"
      CHARACTER (LEN=*), PARAMETER :: zbasm = "modl0000"
      CHARACTER (LEN=*), PARAMETER :: zbasd = "dupl0000"
      CHARACTER (LEN=*), PARAMETER :: zbask = "dpds0000"
      CHARACTER (LEN=*), PARAMETER :: zfmtn =    "(i4.4)"
      INTEGER, PARAMETER  :: ifmts = 5 ! start pos. in names
      INTEGER, PARAMETER  :: ifmte = 8 ! end  pos. in names
      INTEGER, PARAMETER  :: nnamm = 9999 ! number max in names
      INTEGER, PARAMETER  :: klwc = -1 ! case processing: to lower
      INTEGER, PARAMETER  :: kupc =  1 ! case processing: to upper
      INTEGER, PARAMETER  :: klve =  0 ! case processing: leave as is
      INTEGER, PARAMETER  :: kpgm =  0 ! code for main program
      INTEGER, PARAMETER  :: kbdt =  1 ! code for block data
      INTEGER, PARAMETER  :: ksub =  2 ! code for subroutine
      INTEGER, PARAMETER  :: kfun =  3 ! code for function
      INTEGER, PARAMETER  :: kmdl =  4 ! code for module
      INTEGER, PARAMETER  :: kdup =  5 ! code for duplicate
      INTEGER, PARAMETER  :: kdpd =  6 ! code for dependencies
      INTEGER, PARAMETER  :: kend = -1 ! code for end-of-input
      INTEGER, PARAMETER  :: ktabn = 0 ! assume no tabs
      INTEGER, PARAMETER  :: ktabi = 1 ! accept tabs, no expand
      INTEGER, PARAMETER  :: ktabe = 2 ! expand tabs
      INTEGER, PARAMETER  :: nplam = 3 ! # of plans to expand tabs
      INTEGER, PARAMETER  :: luerr = 0 ! logical unit for stderr
      INTEGER, PARAMETER  :: lutmp = 2 ! logical unit for temp. file
      INTEGER, PARAMETER  :: lufil = 3 ! logical unit for final file
      INTEGER, PARAMETER  :: luinp = 5 ! logical unit for stdin
      INTEGER, PARAMETER  :: ludpd = 7 ! logical unit for depend file
      INTEGER, PARAMETER  :: lnamm = 31 ! max. variable name length
      INTEGER, PARAMETER  :: lfilm = 64 ! max. file name length
      INTEGER, PARAMETER  :: ncntm = 39 ! max. # cont. lines
      INTEGER, PARAMETER  :: linem = 132 ! max. line length
      INTEGER, PARAMETER  :: lsttm = (linem-1)*ncntm+linem
                             ! max. statement length
      INTEGER, PARAMETER  :: ndepm =  100 ! max use/include deps
!-----------------------------------------------------------------------
! The following declaration is technically non-standard in Fortran90:
! (the "max" function is not required to be accepted in a PARAMETER
! statement)  to fix this, I added a contained routine, called at the
! beginning of the main program.
!-----------------------------------------------------------------------
!     INTEGER, PARAMETER, DIMENSION (linem, nplam) :: nxttab =  &
!     RESHAPE (                                                 &
!              (/ MAX( (/ (6+3*((i-6+3)/3), i= 1,linem),        &
!                         (6+2*((i-6+2)/2), i= 1,linem) /),     &
!                      (/ (6, i= 1, 2*linem) /)            ),   &
!                 (/                      (i, i= 1,linem) /) /),&
!               (/ linem, nplam /) )

      INTEGER,            DIMENSION (linem, nplam) :: nxttab =  &
      RESHAPE (                                                 &
               (/  (/ (6+3*((i-6+3)/3), i= 1,linem) /),         &
                   (/ (6+2*((i-6+2)/2), i= 1,linem) /),         &
                   (/               (i, i= 1,linem) /) /),      &
                (/ linem, nplam /) )
CONTAINS
   SUBROUTINE maxnxt
      nxttab(:,1:2) = MAX(6,nxttab(:,1:2))
   END SUBROUTINE maxnxt

END MODULE splitprms


MODULE splitdefs
!-----------------------------------------------------------------------
!  Default settings for f90split utility
!-----------------------------------------------------------------------
USE splitprms
!-----------------------------------------------------------------------
      CHARACTER (LEN=*), PARAMETER :: zsuff = ".f90"
      CHARACTER (LEN=*), PARAMETER :: zsufk = ".mk"
      CHARACTER (LEN=*), PARAMETER :: zsufm = ".mod"
      CHARACTER (LEN=*), PARAMETER :: zsufo = ".o"
      INTEGER  :: ktab =  ktabe
      INTEGER  :: kcas =  klve ! code for case processing
      INTEGER  :: kmkd =  1    ! code for making dependencies
END MODULE splitdefs


MODULE splitcurs
USE splitprms
!-----------------------------------------------------------------------
!  Current status variables in f90split utility
!-----------------------------------------------------------------------
      INTEGER, SAVE  :: nlini =  0 ! Lines input
      INTEGER, SAVE  :: nlins =  0 ! in current sub-unit
      INTEGER, SAVE  :: iplac =  1 ! plan for tab expansion
      INTEGER, SAVE  :: mlins =  0 ! max line length
      INTEGER, SAVE  :: ndep  =  0 ! number of use/includes deps
      INTEGER, SAVE  :: iflina = 0   ! advance line is multiple
      INTEGER, SAVE  :: llina =  -1  ! length of advance stored line
      CHARACTER (LEN=linem) :: zlina ! line in advance
      CHARACTER (LEN=lfilm), DIMENSION (ndepm) :: zdept ! current dependencies
END MODULE splitcurs


PROGRAM f90split
!-----------------------------------------------------------------------
!  Split standard input stream, containing source of several fortran90
!  program units into individual files, each containing a single
!  program unit, named after it, or main0001.f90-main9999.f90, or
!  bdta0001.f90-bdta9999.f90. If a file with that name already exists,
!  it is put in dupl0001.f90-dupl9999.f90.
!-----------------------------------------------------------------------
      USE splitidnt
      USE splitdefs
      USE splitcurs

      CHARACTER (LEN=lfilm) :: zfil, zdpd
      CHARACTER (LEN=lnamm) :: znam
!
      !! WRITE (luerr, "(A)") "This is f90split: " // zvers
      !! WRITE (luerr, "(A)")  zusg
      CALL maxnxt

     body: DO
!
!  Open temporary output file
!
         OPEN (lutmp, STATUS='scratch', IOSTAT=kerr)
         IF (kerr /= 0) THEN
            WRITE (luerr,*) "Unable to open scratch file"
            EXIT body
         END IF
!
!  Open dependencies file
!
         IF (kmkd == 1) THEN
            CALL nxtnam (kdpd, zdpd, kerr)
            IF (kerr /= 0) THEN
               WRITE (luerr,*) "Name space exhausted"
               EXIT body
            END IF
            ! Originally written to ludpd
            ! OPEN (*, FILE=TRIM(zdpd)//zsufk, IOSTAT=kerr)
            ! IF (kerr /= 0) THEN
            !    WRITE (luerr,*) "Unable to open dependencies file"
            !    EXIT body
            ! END IF
            ! WRITE (luerr, "(A, A)") "Writing dependencies to ",   &
            !                         TRIM (zdpd) // zsufk
         END IF
!
         DO
            nlins = 0
            IF (ktab == ktabe) iplac = 1
!
!  Find type and name of unit
!
            CALL fndfst (kunt, znam, kerr)
            IF (kunt == kend) THEN
               WRITE (luerr, *) "Trailing comments removed"
               EXIT body
            END IF
            IF (kerr /= 0) THEN
               EXIT body
            END IF
!
!  Find name for corresponding file
!
            CALL getnam (kunt, znam, zfil, kerr)
            lfil = LEN_TRIM (zfil)
            IF (kerr /= 0) THEN
               WRITE (luerr,*) "Name space exhausted"
               EXIT body
            ELSE
               IF (kunt == kdup) THEN
                  WRITE (luerr, *) TRIM (znam), "->", zfil (1:lfil)
               END IF
            END IF
!
!  Find end of current program unit
!
            IF (kmkd == 1) THEN
               CALL fndiue (kerr)
            END IF
            IF (kerr /= 0) THEN
               WRITE (luerr,*) zfil (1:lfil), " : Missing END statement"
            END IF
            REWIND lutmp
!
!  Copy scratch file to destination
!
            IF (kerr /= 0) THEN
               WRITE (luerr,*) zfil (1:lfil), " : Unable to write"
               EXIT body
            ELSE
               IF (kmkd == 1) THEN
                  lsuf = LEN (zsuff)
                  ! Originally written to ludpd
                  WRITE (*,*) zfil (1:lfil-lsuf) // zsufo // " : "
                  WRITE (*,*) zfil (1:lfil)
                  DO i = 1, ndep
                     WRITE (*,*) " ", TRIM (zdept (i))
                  END DO  
                  IF (zsufm /= zsufo .AND. kunt == kmdl) THEN
                  ! Originally written to ludpd 
                  ! THIS IS FOR THE .MOD FILE , SKIP IT!
                  !WRITE (*,*) zfil (1:lfil-lsuf) // zsufm // " : "
                  !WRITE (*,*) zfil (1:lfil)
                  !DO i = 1, ndep
                  !   WRITE (*,*) " ", TRIM (zdept (i))
                  !END DO
                  END IF
                  ndep = 0
               END IF
               !!!! WRITE (*,*) zfil (1:lfil) ! To write name to stdout
            END IF
!
!  Loop to next program unit
!
            REWIND lutmp
         END DO
!
      END DO body
      CLOSE (lutmp)
END PROGRAM f90split


SUBROUTINE fndfst (kunt, znam, kerr)
!-----------------------------------------------------------------------
!  Read input file, copying it to the scratch file, until the first
!  non-comment statement is found.
!  Analyse this statement, and decide of the type and name of the
!  program unit that starts there.
!-----------------------------------------------------------------------
USE splitprms
USE splitcurs
      INTEGER, INTENT (out) :: kunt           ! type of program unit
      CHARACTER (LEN=*), INTENT (out) :: znam ! name chosen
      INTEGER, INTENT (out) :: kerr           ! error code
!
      CHARACTER (LEN=lsttm) :: zstt
!
!     get next statement
!
      CALL nxtstt (zstt, ksta)

      IF (ksta == 0) THEN
         CALL nlzfst (zstt, kunt, znam)
      ELSE IF (ksta < 0) THEN
         IF (nlins > 0) THEN
            kunt = kend
         ELSE
            kunt = kpgm
         END IF
         kerr = -1
      ELSE
         kunt = kpgm
         znam = ' '
         kerr = ksta
      END IF
END SUBROUTINE fndfst


SUBROUTINE fndiue (kerr)
!-----------------------------------------------------------------------
!  Read input file, copying it to the scratch file, and
!  looking for dependencies, until an END statement is found.
!-----------------------------------------------------------------------
USE splitdefs
USE splitcurs
!
      INTEGER, INTENT (out) :: kerr           ! error code
!
      CHARACTER (LEN=lsttm) :: zstt
      CHARACTER (LEN=lnamm) :: znam
      INTEGER, SAVE         :: jlvl = 0
      INTEGER, SAVE         :: jntf = 0
!
      ifnew = 0
      DO
!
!  Get next statement
!
         CALL nxtstt (zstt, ksta)
         IF (ksta /= 0) THEN
            kerr = 1
            EXIT
         END IF
!
!  Look for USE, INCLUDE, or END of sub-unit or of INTERFACE
!
         CALL nlziue (zstt, klst)
         SELECT CASE (klst)
         CASE (-1)! problem
            kerr = 1
            EXIT
         CASE (0) ! not end of anything
            continue
         CASE (1) ! end of sub-unit
            IF (jlvl <= 0 .AND. jntf == 0) THEN
               kerr = 0
               EXIT
            END IF
            IF (jlvl > 0) jlvl = jlvl - 1
            ifnew = 1
            CYCLE
         CASE (2) ! end of interface
            IF (jntf <= 0) THEN
               WRITE (luerr, *) "END INTERFACE out of place"
            ELSE
               jntf = jntf - 1
            END IF
            ifnew = 0
            CYCLE
         END SELECT
!
!  Look for INTERFACE statement
!
         CALL fndntf (zstt, ifntf)
         IF (ifntf /= 0) THEN
            jntf = jntf + 1
            ifnew = 1
            CYCLE
         END IF
!
!  Look for CONTAINS statement
!
         CALL fndctn (zstt, ifctn)
         IF (ifctn /= 0) THEN
            ifnew = 1
            CYCLE
         END IF
!
!  Look for start of new unit
!
         IF (ifnew /= 0) THEN
            CALL nlzfst (zstt, kunt, znam)
            IF (kunt == ksub .OR. kunt == kfun) THEN
               jlvl = jlvl + 1
               ifnew = 0
               CYCLE
            END IF
         END IF
      END DO
END SUBROUTINE fndiue


SUBROUTINE nxtstt (zstt, ksta)
!-----------------------------------------------------------------------
!  Get (possibly multiple) non-comment statement and extract
!  single statement out of it
!-----------------------------------------------------------------------
USE splitcurs
      CHARACTER (LEN=lsttm), INTENT (out) :: zstt
      INTEGER, INTENT (out)               :: ksta ! status code
!
      CHARACTER (LEN=1)           :: zdlm
      CHARACTER (LEN=lsttm), SAVE :: zmul
      INTEGER, SAVE               :: istt = 0
      INTEGER, SAVE               :: lmul
!
      ksta = 0
      body: DO
         IF (istt == 0) THEN
!
!  Get a (possibly multiple) non-comment statement
!
            CALL reastt (zmul, lmul, kget)
!
            IF (kget /= 0) THEN
               ksta = kget
               istt = 0
               EXIT body
            ELSE
               istt = 1
            END IF
         END IF
!
!  Look for character context
!
         ifchc1 = 0
         iloo   = istt
         lstt   = lmul
         DO
!
!  Outside of character context, truncate at ; if any
!
            IF (ifchc1 == 0) THEN
               ichc0 = SCAN (zmul (iloo:lstt), "'"//'"')
               IF (ichc0 == 0) THEN
                  ismc = INDEX (zmul (iloo:lstt), ';')
               ELSE
                  ismc = INDEX (zmul (iloo:ichc0), ';')
               END IF
               IF (ismc > 0) THEN
                  lstt = iloo + ismc - 2
                  EXIT
               ELSE IF (ichc0 > 0) THEN
                  ifchc1 = 1
                  iloo = iloo + ichc0
                  zdlm = zmul (iloo-1:iloo-1)
               ELSE
                  EXIT
               END IF
            ELSE
!
!  Within character context, look for its termination
!
               ichc1 = SCAN (zmul (iloo:lstt), zdlm)
               IF (ichc1 == 0) THEN
                  EXIT
               ELSE
                  ifchc1 = 0
                  iloo  = iloo + ichc1
               END IF
            END IF
         END DO
!
!  Copy current statement into zstt
!
         zstt = zmul (istt:lstt)
         IF (istt == 1 .AND. lstt == lmul) THEN
            iflina = 0
         ELSE
            iflina = 1
         END IF
         IF (lstt+1 < lmul) THEN
            istt = lstt + VERIFY (zmul (lstt+2:lmul), ' ') + 1
         ELSE
            istt = 0
         END IF
         IF (LEN_TRIM (zstt) > 0) THEN
            EXIT body
         END IF
      END DO body
END SUBROUTINE nxtstt


SUBROUTINE reastt (zmul, lstt, ksta)
!-----------------------------------------------------------------------
!  Read input file, copying it to the scratch file, until a
!  (possibly multiple) non-comment statement is found.
!-----------------------------------------------------------------------
USE splitdefs
USE splitcurs
      CHARACTER (LEN=lsttm), INTENT (out) :: zmul
      INTEGER, INTENT (out)               :: lstt ! istt. length
      INTEGER, INTENT (out)               :: ksta ! status code
!
      CHARACTER (LEN=linem) :: zlin
      CHARACTER (LEN=1)     :: zdlm
!
      lstt  = 0
      ifchc0 = 0
      ifcnt0 = 0
      DO
!
!  Something to write ?
!  Write advance line to scratch file
!
            IF (llina > 0) THEN
               WRITE (lutmp, "(A)", IOSTAT=kwri) zlina (1:llina)
               IF (kwri /= 0) THEN
                  WRITE (luerr,*) "Problem writing scratch file"
                  ksta = 2
                  EXIT
               END IF
            ELSE IF (llina == 0) THEN
               WRITE (lutmp, "()", IOSTAT=kwri)
               IF (kwri /= 0) THEN
                  WRITE (luerr,*) "Problem writing scratch file"
                  ksta = 2
                  EXIT
               END IF
            END IF
!
!  Read a line
!
         READ (luinp, "(A)", IOSTAT=krea) zlin
!
         SELECT CASE (krea)
         CASE (1:)
            ksta = 1
            llina = -1
            WRITE (luerr,*) "Problem reading input"
            EXIT
         CASE (:-1)
            ksta = -1
            llina = -1
            EXIT
         CASE (0)
            ksta = 0
            nlini = nlini + 1
            nlins = nlins + 1
            llin  = LEN_TRIM (zlin)
            llina = llin
            IF (llin <= 0) CYCLE
            zlina (1:llina) = zlin (1:llin)
!
!  process tabs
!
            SELECT CASE (ktab)
            CASE (ktabi)
               CALL rmvtab (zlin, llin)
               mlins = MAX (mlins, llin)
            CASE (ktabn)
               mlins = MAX (mlins, llin)
            CASE (ktabe)
               CALL chktab (zlin, llin)
               CALL rmvtab (zlin, llin)
            END SELECT
!
!  Recognize and skip comments
!
            ifst = VERIFY (zlin (1:llin), ' ')
            IF (ifst == 0) CYCLE
            IF (zlin (ifst:ifst) == '!') CYCLE
!
!  Recognize and skip pre-processing commands
!
            IF (zlin (ifst:ifst) == '$') CYCLE
            IF (zlin (ifst:ifst) == '#') CYCLE
!
!  Do not explore trailing comments if any
!
!  Look for character context
!
            ifchc1 = ifchc0
            iloo = ifst
            lxpl = llin
            DO
!
!  Outside of character context, truncate at ! if any
!
               IF (ifchc1 == 0) THEN
                  ichc0 = SCAN (zlin (iloo:llin), "'"//'"')
                  IF (ichc0 == 0) THEN
                     icmt = INDEX (zlin (iloo:llin), '!')
                  ELSE
                     icmt = INDEX (zlin (iloo:ichc0), '!')
                  END IF
                  IF (icmt > 0) THEN
                     ltmp = iloo + icmt - 2
                     lxpl = LEN_TRIM (zlin (1:ltmp))
                     EXIT
                  ELSE IF (ichc0 > 0) THEN
                     ifchc1 = 1
                     iloo = iloo + ichc0
                     zdlm = zlin (iloo-1:iloo-1)
                  ELSE
                     EXIT
                  END IF
               ELSE
!
!  Within character context, look for its termination
!
                  ichc1 = SCAN (zlin (iloo:llin), zdlm)
                  IF (ichc1 == 0) THEN
                     EXIT
                  ELSE
                     ifchc1 = 0
                     iloo  = iloo + ichc1
                  END IF
               END IF
            END DO
!
!  Look for continuation mark
!
            IF (zlin (lxpl:lxpl) == '&') THEN
               ifcnt1 = 1
               llin   = LEN_TRIM (zlin (1:lxpl-1))
            ELSE
               ifcnt1 = 0
            END IF
!
!  Copy current statement fragment into zmul
!
!  Look for continued mark
!
            IF (zlin (ifst:ifst) == '&') THEN
               ifst = ifst + VERIFY (zlin (ifst+1:llin), ' ')
               IF (ifchc0 == 0) THEN
                  lstt = lstt + 1
                  zmul (lstt:lstt) = ' '
               END IF
            END IF
!
!  Copy
!
            lfrg = llin - ifst + 1
            zmul (lstt+1:lstt+lfrg) = zlin (ifst:llin)
            lstt = lstt + lfrg
            IF (ifcnt1 == 0) EXIT
         END SELECT
!
!  Loop until end of statement
!
         ifcnt0 = ifcnt1
         ifchc0 = ifchc1
      END DO
END SUBROUTINE reastt


SUBROUTINE nlzfst (zstt, kunt, znam)
!-----------------------------------------------------------------------
!  Analyse a statement, and decide of the type (and name) of the
!  program unit that starts there.
!-----------------------------------------------------------------------
USE splitcurs
      CHARACTER (LEN=lsttm), INTENT (in) :: zstt ! the statement
      INTEGER, INTENT (out) :: kunt              ! type of program unit
      CHARACTER (LEN=*), INTENT (out) :: znam    ! name chosen
!
      CHARACTER (LEN=lsttm) :: zsttw
      LOGICAL               :: ifwrk
!
body: DO
!
!  Raise to upper case (No label to be removed)
!
         zsttw = zstt
         CALL raicas (zsttw)
!
!  Look for PROGRAM
!
         lstt = LEN_TRIM (zsttw)
         ipgm = INDEX (zsttw (1:lstt), zpgm)
         IF (ipgm == 1) THEN
            kunt = kpgm
            ikwdf = lpgm
         ELSE
!
!  Look for MODULE
!
            imdl = INDEX (zsttw (1:lstt), zmdl)
            IF (imdl == 1) THEN
               kunt = kmdl
               ikwdf = lmdl
            ELSE
!
!  Look for FUNCTION
!
               ifun = INDEX (zsttw (1:lstt), zfun)
               IF (ifun <= 1) THEN
                  ifwrk = (ifun == 1)
               ELSE
                  ifwrk = (zsttw (ifun-1:ifun-1) == ' ')
               END IF
               IF (ifwrk) THEN
                  kunt = kfun
                  ikwdf = lfun + ifun - 1
               ELSE
!
!  Look for SUBROUTINE
!
                  isub = INDEX (zsttw (1:lstt), zsub)
                  IF (isub <= 1) THEN
                     ifwrk = (isub == 1)
                  ELSE
                     ifwrk = (zsttw (isub-1:isub-1) == ' ')
                  END IF
                  IF (ifwrk) THEN
                     kunt = ksub
                     ikwdf = lsub + isub - 1
                  ELSE
!
!  Look for BLOCK DATA
!
                     ibdt1 = INDEX (zsttw (1:lstt), zbdt1)
                     IF (ibdt1 == 1) THEN
                        ikwdf = lbdt1 &
                              + VERIFY (zsttw (lbdt1+1:lstt), ' ') &
                              - 1
                        IF (ikwdf >= lbdt1) THEN
                           ibdt2 = INDEX (zsttw (ikwdf+1:lstt), zbdt2)
                           IF (ibdt2 == 1) THEN
                              kunt = kbdt
                              ikwdf = ikwdf + lbdt2
                           ELSE
                              kunt = kpgm
                              znam = ' '
                              EXIT body
                           END IF
                        ELSE
                           kunt = kpgm
                           znam = ' '
                           EXIT body
                        END IF
                     ELSE
                        kunt = kpgm
                        znam = ' '
                        EXIT body
                     END IF
                  END IF
               END IF
            END IF
         END IF
!
!  Find name
!
         inams = ikwdf + VERIFY (zsttw (ikwdf+1:lstt), ' ')
         IF (inams < ikwdf+2) THEN
            IF (kunt /= kbdt) kunt = kpgm
            znam = ' '
            EXIT body
         END IF
         iname = inams   &
               + VERIFY (zsttw (inams+1:lstt+1), zupc//zdgt//"_") &
               - 1
         IF (iname < inams) THEN
            IF (kunt /= kbdt) kunt = kpgm
            znam = ' '
            EXIT body
         END IF
         znam = zstt (inams:iname)
         EXIT body
      END DO body
      kwri = 0
      IF (iflina /= 0 .AND. llina < 0) THEN
         WRITE (lutmp, "(A)", IOSTAT=kwri) TRIM (zstt)
      ELSE IF (llina >= 0) THEN
         IF (iflina == 0) THEN
            WRITE (lutmp, "(A)", IOSTAT=kwri) zlina (1:llina)
         ELSE
            WRITE (lutmp, "(A)", IOSTAT=kwri) TRIM (zstt)
         END IF
         llina = -1
      END IF
      IF (kwri /= 0) THEN
         WRITE (luerr,*) "Problem writing scratch file"
      END IF
END SUBROUTINE nlzfst


SUBROUTINE nlziue (zstt, klst)
!-----------------------------------------------------------------------
!  Analyse a statement, and decide if it is use, include, or if current
!  program unit ends there.
!-----------------------------------------------------------------------
USE splitcurs
USE splitdefs
      CHARACTER (LEN=lsttm), INTENT (in) :: zstt ! The statement
      INTEGER, INTENT (out) :: klst              ! result
!
      CHARACTER (LEN=lsttm) :: zsttw
      CHARACTER (LEN=1)     :: zdlm
!
body: DO
         zsttw = zstt
!
!  Remove label and raise to upper case
!
         CALL rmvlbl (zsttw)
         CALL raicas (zsttw)
!
!  Look for first token, to be INCLUDE, USE, or END
!
         itokf = VERIFY (zsttw, zupc) - 1
         lstt = LEN_TRIM (zsttw)
         klst = 0
         IF (itokf == luse) THEN
             IF (zsttw(1:luse) == zuse) THEN
!
!  Look for [space] use name
!
                 itoks = luse + VERIFY (zsttw (luse+1:lstt), ' ')
                 itoke = itoks + VERIFY (zsttw (itoks+1:lstt+1),        &
                                         zupc//zdgt//"_") - 1
                 IF (ndep < ndepm) THEN
                    ndep = ndep + 1
                    zdept (ndep) = zsttw (itoks:itoke) // zsufo
                    CALL lwrcas (zdept (ndep))
                    DO idep = 1, ndep - 1
                      IF (zdept (idep) == zdept (ndep)) THEN
                        ndep = ndep - 1
                        EXIT body
                      END IF
                    END DO
                 END IF
                 EXIT body
             END IF
         END IF

         IF (itokf == linc) THEN
             IF (zsttw(1:linc) == zinc) THEN
!
!  Look for [space] 'include_string' or "include_string"
!
                 itoks = linc + VERIFY (zsttw (linc+1:lstt), ' ')
                 zdlm  = zsttw (itoks:itoks)
                 IF (zdlm /= '"' .AND. zdlm /= "'") THEN
                    EXIT body
                 END IF
                 itoks = itoks + 1
                 itoke = itoks + INDEX (zsttw (itoks+1:lstt), zdlm) - 1
                 IF (itoke == itoks-1) THEN
                     EXIT body    ! no trailing delim found
                 END IF

                 IF (ndep < ndepm) THEN
                    ndep = ndep + 1
                    zdept (ndep) = zsttw (itoks:itoke)
                    DO idep = 1, ndep - 1
                       IF (zdept (idep) == zdept (ndep)) THEN
                          ndep = ndep - 1
                          EXIT body
                       END IF
                    END DO
                 END IF
                 EXIT body
             END IF
         END IF

         IF (itokf < lend) THEN
            klst = 0
            EXIT body
         END IF
         IF (zsttw (1:lend) /= zend) THEN
            klst = 0
            EXIT body
         END IF
!
!  Nothing after END
!
         IF (lstt == lend) THEN
            klst = 1
            EXIT body
         END IF
!
!  Look for [space] unit name
!
         itoks = lend + VERIFY (zsttw (lend+1:lstt), ' ')
         itoke = itoks + INDEX (zsttw (itoks+1:lstt+1), ' ') - 1
         IF (itoke < itoks+2) THEN
            klst = 0
            EXIT body
         END IF
         IF (    (zsttw (itoks:itoke) == zpgm)     &
             .OR.(zsttw (itoks:itoke) == zsub)     &
             .OR.(zsttw (itoks:itoke) == zfun)     &
             .OR.(zsttw (itoks:itoke) == zbdt)     & ! Be laxist
             .OR.(zsttw (itoks:itoke) == zmdl)     ) THEN
            klst = 1
            EXIT body
         ELSE IF (zsttw (itoks:itoke) == zntf) THEN
            klst = 2
            EXIT body
         ELSE IF (zsttw (itoks:itoke) == zbdt1) THEN
            itoks = itoke + VERIFY (zsttw (itoke+1:lstt), ' ')
            IF (itoks < itoke+2) THEN
               klst = 0
               EXIT body
            END IF
            itoke = itoks + INDEX (zsttw (itoks+1:lstt+1), ' ') - 1
            IF (itoke < itoks+2) THEN
               klst = 0
               EXIT body
            END IF
            IF (zsttw (itoks:itoke) == zbdt2) THEN
               klst = 1
               EXIT body
            ELSE
               klst = 0
               EXIT body
            END IF
         ELSE
            klst = 0
            EXIT body
         END IF
      END DO body
      kwri = 0
      IF (iflina /= 0 .AND. llina < 0) THEN
         WRITE (lutmp, "(A)", IOSTAT=kwri) TRIM (zstt)
      ELSE
         IF (iflina == 0) THEN
            WRITE (lutmp, "(A)", IOSTAT=kwri) zlina (1:llina)
         ELSE
            WRITE (lutmp, "(A)", IOSTAT=kwri) TRIM (zstt)
         END IF
         llina = -1
      END IF
      IF (kwri /= 0) THEN
         WRITE (luerr,*) "Problem writing scratch file"
         klst = -1
      END IF
END SUBROUTINE nlziue


SUBROUTINE fndctn (zstt, ifctn)
!-----------------------------------------------------------------------
!  Look for CONTAINS statement
!-----------------------------------------------------------------------
USE splitprms
      CHARACTER (LEN=lsttm), INTENT (in) :: zstt ! The statement
      INTEGER, INTENT (out) :: ifctn
!
      CHARACTER (LEN=lsttm) :: zsttw
!
body: DO
         zsttw = zstt
!
!  Remove label and raise to upper case
!
         CALL rmvlbl (zsttw)
         CALL raicas (zsttw)
!
!  Look for first token, to be CONTAINS
!
         itokf = VERIFY (zsttw, zupc) - 1
         IF (itokf /= lctn) THEN
            ifctn = 0
            EXIT body
         END IF
         IF (zsttw (1:lctn) /= zctn) THEN
            ifctn = 0
            EXIT body
         END IF
!
!  Nothing after CONTAINS
!
         lstt = LEN_TRIM (zsttw)
         IF (lstt == lctn .AND. zsttw (1:lctn) == zctn) THEN
            ifctn = 1
            EXIT body
         ELSE
            ifctn = 0
            EXIT body
         END IF
      END DO body
END SUBROUTINE fndctn


SUBROUTINE fndntf (zstt, ifntf)
!-----------------------------------------------------------------------
!  Look for INTERFACE statement
!-----------------------------------------------------------------------
USE splitprms
      CHARACTER (LEN=lsttm), INTENT (in) :: zstt ! The statement
      INTEGER, INTENT (out) :: ifntf
!
      CHARACTER (LEN=lsttm) :: zsttw
!
body: DO
         zsttw = zstt
!
!  Remove label and raise to upper case
!
         CALL rmvlbl (zsttw)
         CALL raicas (zsttw)
!
!  Look for first token, to be INTERFACE
!
         itokf = VERIFY (zsttw, zupc) - 1
         IF (itokf /= lntf) THEN
            ifntf = 0
            EXIT body
         END IF
         IF (zsttw (1:lntf) /= zntf) THEN
            ifntf = 0
            EXIT body
         END IF
!
!  Nothing after INTERFACE
!
         lstt = LEN_TRIM (zsttw)
         IF (lstt == lntf .AND. zsttw (1:lntf) == zntf) THEN
            ifntf = 1
            EXIT body
         ELSE IF (lstt > lntf .AND. zsttw (1:lntf) == zntf) THEN
            IF (zsttw (lntf+1:lntf+1) == ' ') THEN
               ifntf = 1
               EXIT body
            ELSE
               ifntf = 0
               EXIT body
            END IF
         ELSE
            ifntf = 0
            EXIT body
         END IF
      END DO body
END SUBROUTINE fndntf


SUBROUTINE getnam (kunt, znam, zfil, kerr)
!-----------------------------------------------------------------------
!  Return a file name from the type (and name) of the
!  program unit that is processed.
!-----------------------------------------------------------------------
USE splitdefs
      INTEGER, INTENT (INOUT) :: kunt            ! type of program unit
      CHARACTER (LEN=*), INTENT (in) :: znam     ! name if any
      CHARACTER (LEN=*), INTENT (out) :: zfil    ! file name
      INTEGER, INTENT (out) :: kerr              ! error code
!
      LOGICAL :: ifxst
      CHARACTER (LEN=lnamm) :: znamw
!
!  Change according to desired case
!
      znamw = znam
      lnam = LEN_TRIM (znamw)
      SELECT CASE (kcas)
      CASE (-1)
         CALL lwrcas (znamw (1:lnam))
      CASE (+1)
         CALL raicas (znamw (1:lnam))
      CASE DEFAULT
         CONTINUE
      END SELECT
      IF (lnam > 0) THEN
!
!  Check that name is valid
!
         zfil = znamw (1:lnam) // zsuff
         inquire (file=zfil, exist=ifxst)
         IF (ifxst) THEN
            kunt = kdup
            CALL nxtnam (kunt, znamw, kerr)
            lnam = LEN_TRIM (znamw)
            zfil = znamw (1:lnam) // zsuff
         ELSE
            kerr = 0
         END IF
      ELSE
         CALL nxtnam (kunt, znamw, kerr)
         lnam = LEN_TRIM (znamw)
         zfil = znamw (1:lnam) // zsuff
      END IF
END SUBROUTINE getnam


SUBROUTINE nxtnam (kunt, znam, kerr)
!-----------------------------------------------------------------------
!  Return the next name for the type of the
!  program unit that is processed.
!-----------------------------------------------------------------------
USE splitdefs
      INTEGER, INTENT (in) :: kunt               ! type of program unit
      CHARACTER (LEN=*), INTENT (out) :: znam    ! name
      INTEGER, INTENT (out) :: kerr              ! error code
!
      LOGICAL       :: ifxst
      INTEGER, SAVE :: idpd = 0
      INTEGER, SAVE :: ipgm = 0
      INTEGER, SAVE :: ibdt = 0
      INTEGER, SAVE :: idup = 0
      INTEGER, SAVE :: imdl = 0
      CHARACTER (LEN=lfilm) :: zsuf
!
body: DO
         SELECT CASE (kunt)
         CASE (kdpd)
            idpd = idpd + 1
            inum = idpd
            znam = zbask
            zsuf = zsufk
         CASE (kpgm)
            ipgm = ipgm + 1
            inum = ipgm
            znam = zbasp
            zsuf = zsuff
         CASE (kmdl)
            imdl = imdl + 1
            inum = imdl
            znam = zbasm
            zsuf = zsuff
         CASE (kbdt)
            ibdt = ibdt + 1
            inum = ibdt
            znam = zbasb
            zsuf = zsuff
         CASE (ksub,kfun,kdup)
            idup = idup + 1
            inum = idup
            znam = zbasd
            zsuf = zsuff
         END SELECT
         DO
            IF (inum > nnamm) THEN
               kerr = 1
               EXIT body
            END IF
            WRITE (znam (ifmts:ifmte), zfmtn) inum
            inquire (file=TRIM(znam)//TRIM(zsuf), exist=ifxst)
            IF (ifxst) THEN
               inum = inum + 1
               CYCLE
            ELSE
               EXIT
            END IF
         END DO
         kerr = 0
         EXIT body
      END DO body
END SUBROUTINE nxtnam


SUBROUTINE raicas (zstr)
!-----------------------------------------------------------------------
!  Raise a string to upper case
!-----------------------------------------------------------------------
USE splitprms
      CHARACTER (LEN=*), INTENT (INOUT) :: zstr  ! The string
      LOGICAL :: toggle
      CHARACTER(LEN=1) :: togglechar
!
!   Modified to not do upper case to embedded strings (pgg 21.11.94)
!
      toggle = .TRUE.
      lstr = LEN_TRIM (zstr)
      DO istr = 1, lstr
         IF (toggle) THEN
            IF (zstr (istr:istr) == '"' .OR. zstr (istr:istr) == "'") THEN
               toggle = .NOT. toggle
               togglechar = zstr (istr:istr)
            END IF
            irnk = INDEX (zlwc, zstr (istr:istr))
            IF (irnk > 0) THEN
               zstr (istr:istr) = zupc (irnk:irnk)
            END IF
         ELSE
            IF (zstr (istr:istr) == togglechar) toggle = .NOT. toggle
         END IF
      END DO
END SUBROUTINE raicas


SUBROUTINE lwrcas (zstr)
!-----------------------------------------------------------------------
!  Lower a string to lower case
!-----------------------------------------------------------------------
USE splitprms
      CHARACTER (LEN=*), INTENT (INOUT) :: zstr  ! The string
!
      lstr = LEN_TRIM (zstr)
      DO istr = 1, lstr
         irnk = INDEX (zupc, zstr (istr:istr))
         IF (irnk > 0) THEN
            zstr (istr:istr) = zlwc (irnk:irnk)
         END IF
      END DO
END SUBROUTINE lwrcas


SUBROUTINE rmvlbl (zstt)
!-----------------------------------------------------------------------
!  Remove statement label (Note: Label /= Construct name)
!-----------------------------------------------------------------------
USE splitprms
      CHARACTER (LEN=lsttm), INTENT (INOUT) :: zstt  ! The statement
!
      IF (INDEX (zdgt, zstt (1:1)) > 0) THEN
         istt = VERIFY (zstt, zdgt//' ')
         zstt = zstt (istt:lsttm)
      END IF
END SUBROUTINE rmvlbl


SUBROUTINE rmvtab (zstr, lstr)
!-----------------------------------------------------------------------
!  Remove tabs and replace with spaces
!-----------------------------------------------------------------------
USE splitprms
      CHARACTER (LEN=lstr), INTENT (INOUT) :: zstr  ! The string
      INTEGER, INTENT (INOUT)              :: lstr  ! its trimmed length
!
      lsrc = lstr
      DO
!
!  Search backwards so that trailing tabs eliminated first
!
         lsrc = INDEX (zstr (1:lsrc), ztab, back=.true.)
         IF (lsrc == 0) THEN
            EXIT
         END IF
         zstr (lsrc:lsrc) = ' '
         lsrc = lsrc - 1
      END DO
      lstr = LEN_TRIM (zstr)
END SUBROUTINE rmvtab


SUBROUTINE chktab (zstr, lstr)
!-----------------------------------------------------------------------
!  Verify and possibly update current tab expansion plan
!-----------------------------------------------------------------------
USE splitdefs
USE splitcurs
      CHARACTER (LEN=*), INTENT (INOUT) :: zstr  ! The string
      INTEGER, INTENT (INOUT)           :: lstr  ! its trimmed length
!
      lexp = lstr
!
!  Quick return when possible
!
body: DO
      IF (iplac == nplam) EXIT body
      IF (INDEX (zstr (1:lstr), ztab) == 0) EXIT body
      IF (VERIFY (zstr (1:lstr), ztab//' ') == 0) THEN
         lexp = 0
         lstr = 0
         EXIT body
      END IF
!
!  Loop on expansion plans
!
      DO
         istr = 1
         lexp = 0
         CALL expand
!
!  Check if line fits with current plan
!
         IF (lexp > linem) THEN
            iplac = iplac + 1
            IF (iplac < nplam) CYCLE
            lexp = lstr
         END IF
         EXIT body
      END DO
      END DO body
      mlins = MAX (lexp, mlins)
CONTAINS
   SUBROUTINE expand
!
!  Expand each tab on to next tab mark
!
      DO
         IF (lstr >= istr) THEN
            iwrk = INDEX (zstr (istr:lstr), ztab)
         ELSE
            EXIT
         END IF
         IF (iwrk /= 0) THEN
            lexp = lexp + iwrk - 1
            istr = istr + iwrk
!
!  Expand tab on to next tab mark
!
            iexp = lexp + 1
            lfil = MIN (iexp, linem)
            lexp = nxttab (lfil, iplac)
!
!  Fill-up with spaces
!
         ELSE
            EXIT
         END IF
      END DO
      lexp = lexp + lstr - istr + 1
   END SUBROUTINE expand
END SUBROUTINE chktab

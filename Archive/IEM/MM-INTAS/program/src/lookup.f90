!#1 subroutine lookup

!#1
!=======================================================================
!      subroutine lookup
!=======================================================================
      SUBROUTINE lookup (loutf, llogf, istep, &
                         np, npd, nsv, ndv, &
                         f_kf, f_kdf, f_kvol)

!-----  used modules
      USE cmn_lookup, ONLY : ilut, lu_flam, lu_ildm, lu_isat, ilumthd, &
                             index_depprop, prho

      USE fllookup_module, ONLY : ppropv

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf, llogf
      INTEGER (KIND=4), INTENT(IN) :: istep
      INTEGER (KIND=4), INTENT(IN) :: np, npd, nsv, ndv

!-----  input/output variables
      REAL    (KIND=8), INTENT(INOUT) :: f_kf(npd,nsv)
      REAL    (KIND=8), INTENT(INOUT) :: f_kdf(npd,ndv)

!-----  output variables
      REAL    (KIND=8), INTENT(OUT) :: f_kvol(np)


!-----  MAIN ACTION

      !---  check for look-up method
      SELECT CASE (ilumthd)

      !--------------------------------------------------------------
      !--- FLAME look-up
      !--------------------------------------------------------------
      CASE (lu_flam)

         !---  unconditional lookup every 'ilut' iteration
         IF (istep == 0) THEN  ! this is for the first call in 'init_mc'
            WRITE (llogf,'(2X,A)') 'LOOKUP: calling ppropv'
            CALL ppropv(np, npd, nsv, ndv, &
                        index_depprop, prho, &
                        f_kf, f_kdf, f_kvol)
         ELSE IF (ilut > 0) THEN
            IF (MOD(istep, ilut) == 0) THEN
               CALL ppropv(np, npd, nsv, ndv, &
                           index_depprop, prho, &
                           f_kf, f_kdf, f_kvol)
            END IF
         END IF

      !--------------------------------------------------------------
      !--- ILDM look-up
      !--------------------------------------------------------------
      CASE (lu_ildm)
         STOP 'LOOKUP: ILDM look-up method not implemented here'

      !--------------------------------------------------------------
      !--- ISAT step (integration / look-up)
      !--------------------------------------------------------------
      CASE (lu_isat)
         STOP 'LOOKUP: ISAT look-up method not implemented here'

      !--------------------------------------------------------------
      !--- unknown method
      !--------------------------------------------------------------
      CASE DEFAULT
         WRITE (loutf,*) 'Error in subroutine LOOKUP'
         WRITE (loutf,*) 'Look-up method unknown'
         STOP

      END SELECT


      END SUBROUTINE lookup

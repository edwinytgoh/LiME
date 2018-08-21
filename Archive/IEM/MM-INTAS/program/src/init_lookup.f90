!#1 subroutine init_lookup_size
!#2 subroutine init_lookup


!#1
!=======================================================================
!  subroutine init_lookup_size
!=======================================================================
!  function :
!    depending on the case studied / lookup method used, returns the
!    number of independent scalars (nindv) and dependent scalars (ndepv)
!-----------------------------------------------------------------------
      SUBROUTINE init_lookup_size(icasestudied, nindv, ndepv)

!-----  used modules
      USE cmn_lookup, ONLY : lu_flam, lu_ildm, lu_isat, ilumthd, case_methane
      USE fllookup_module, ONLY : luinit

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: icasestudied

!-----  output variables
      INTEGER (KIND=4), INTENT(OUT) :: nindv, ndepv

!-----  MAIN ACTION

      SELECT CASE (ilumthd)

      CASE (lu_flam)

         !----  initialise scalar variables
         SELECT CASE (icasestudied)
            CASE (0)
               nindv = 1
               ndepv = 0

            CASE (case_methane)
               nindv =  1
               ndepv = 10

         CASE DEFAULT
            STOP 'Case not implemented.'
         END SELECT

         IF (icasestudied > 0) THEN
            CALL luinit
         END IF


      CASE (lu_ildm)
         STOP 'INIT_LOOKUP_SIZE: ILDM look-up method not implemented here'

      CASE (lu_isat)
         STOP 'INIT_LOOKUP_SIZE: ISAT look-up method not implemented here'


      CASE DEFAULT
         STOP 'Look-up method unknown'
      END SELECT


      END SUBROUTINE init_lookup_size


!#2
!=======================================================================
!  subroutine init_lookup
!=======================================================================
!  function :
!    . depending on the case studied / lookup method used, returns the
!      name of independent scalars (name_indprop) and dependent scalars
!      (name_depprop)
!    . build the cross-reference tables 'index_depprop'
!    . set 'prho' (indice of DENSITY in the dependent variable array)
!-----------------------------------------------------------------------
      SUBROUTINE init_lookup (icasestudied, nindv, ndepv, &
                              name_indprop, name_depprop)

!-----  used modules
      USE cmn_lookup, ONLY : index_depprop, prho, &
                             ilumthd, lu_flam, lu_ildm, lu_isat, &
                             case_methane

      USE fllookup_module, ONLY : get_index_depprop

      IMPLICIT NONE

!-----  VARIABLES

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: icasestudied
      INTEGER (KIND=4), INTENT(IN) :: nindv, ndepv

!-----  output variables
      CHARACTER(LEN=*), INTENT(OUT) :: name_indprop(nindv)
      CHARACTER(LEN=*), INTENT(OUT) :: name_depprop(ndepv)

!-----  MAIN ACTION

      SELECT CASE (ilumthd)

      CASE (lu_flam)
     ! Remark: in subroutine 'ppropv' it is assumed that mixture fraction
     !         is the first independent scalar and enthalpy (if used) is
     !         the second one.

         SELECT CASE (icasestudied)
            CASE (0)
               ! Independent variables, nindv = 1
               name_indprop = (/'MEAN_MIXTURE_FRACTION    '/)

            CASE (case_methane)
               ! Independent variables, nindv = 1
               name_indprop = (/'MEAN_MIXTURE_FRACTION    '/)

               ! Dependent variables, ndepv = 10
               name_depprop = (/'RHO             ', &
                                'TEMP            ', &
                                'CPM             ', &
                                'MASS CH4        ', &
                                'MASS O2         ', &
                                'MASS N2         ', &
                                'MASS H2         ', &
                                'MASS H2O        ', &
                                'MASS CO         ', &
                                'MASS CO2        ' /)
               prho = 1

         CASE DEFAULT
            STOP 'Case not implemented.'
         END SELECT

         IF (icasestudied > 0) THEN
            ALLOCATE(index_depprop(ndepv))
            CALL get_index_depprop (ndepv, name_depprop, index_depprop)
         END IF

      CASE (lu_ildm)
         STOP 'INIT_LOOKUP: ILDM look-up method not implemented here'

      CASE (lu_isat)
         STOP 'INIT_LOOKUP: ISAT look-up method not implemented here'

      CASE DEFAULT
         STOP 'Look-up method unknown'
      END SELECT


      END SUBROUTINE init_lookup

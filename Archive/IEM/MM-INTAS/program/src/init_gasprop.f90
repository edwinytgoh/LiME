!#1  subroutine set_property_indices
!#2  subroutine set_property_names


!#1
!=======================================================================
!  subroutine set_property_indices
!=======================================================================
!  Function:
!    Set the indice of the properties of array 'f(particle,property)'
!-----------------------------------------------------------------------
      SUBROUTINE set_property_indices(loutf, nindv, ndepv)

!-----  used modules
      USE cmn_gaspdf, ONLY : nsv, ndv, nkd, kvol, kf, kl, kdf, kdl

      IMPLICIT NONE

!-----  input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf
      INTEGER (KIND=4), INTENT(IN) :: nindv
      INTEGER (KIND=4), INTENT(IN) :: ndepv


!------  MAIN ACTION
      nsv = nindv
      ndv = ndepv

      nkd =  1 + nsv + ndv

      !-----  set property indices
      kf   = 1
      kl   = kf   + nsv - 1
      kvol = kl   + 1
      kdf  = kvol + 1
      kdl  = kdf  + ndv - 1

      !-----  check nkd, array dimension of f(npd,nkd) from param, with
      !       calculated kdl
      IF (nkd /= kdl) THEN
         WRITE(loutf,*) 'SET_PROPERTY_INDICES: nkd=', nkd, '.NE. kdl=', kdl
         IF (nkd < kdl) STOP
      END IF

      END SUBROUTINE set_property_indices


!#2
!=======================================================================
!  subroutine set_property_names
!=======================================================================
!  Function:
!    Initialise 'iprlab': names of the properties of array 'f'
!-----------------------------------------------------------------------
      SUBROUTINE set_property_names(name_indprop, name_depprop)

!-----  used modules
      USE cmn_gaspdf, ONLY : kvol, kf, kl, kdf, kdl, &
                             iprlab, nsv, ndv, nkd

      IMPLICIT NONE

!-----  VARIABLES

!-----  input variables
      CHARACTER(LEN=*), INTENT(IN) :: name_indprop(nsv)
      CHARACTER(LEN=*), INTENT(IN) :: name_depprop(ndv)

!-----  local variables
      INTEGER (KIND=4) :: k

!------  MAIN ACTION
      ALLOCATE(iprlab(nkd)); iprlab = ' '

      iprlab(kvol) = 'SPECIFIC_VOLUME'

      !--  independent scalars
      DO k = kf, kl
         iprlab(k) = name_indprop(k - kf + 1)
      END DO

      !--  dependent scalars
      DO k = kdf, kdl
         iprlab(k) = name_depprop(k - kdf + 1)
      END DO

      END SUBROUTINE set_property_names

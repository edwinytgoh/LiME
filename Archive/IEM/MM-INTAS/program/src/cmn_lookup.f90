
      MODULE cmn_lookup

!======================================================================= 
!
!     LOOK-UP TABLE PROCEDURE
!
!     ilut     : designator table look-up
!                0   - no look-up
!                >=1 - look-up every ilut-th iteration
!     ilumthd  : designates table look-up method
!                1   - FLAME table
!                2   - ILDM ==> not implemented here
!                3   - ISAT ==> not implemented here
!
!-----------------------------------------------------------------------

      INTEGER(KIND=4), SAVE :: ilut = 1           ! set to default value

      INTEGER (KIND=4), PARAMETER :: lu_flam = 1
      INTEGER (KIND=4), PARAMETER :: lu_ildm = 2
      INTEGER (KIND=4), PARAMETER :: lu_isat = 3
      INTEGER (KIND=4), SAVE :: ilumthd = lu_flam ! set to default value

      INTEGER (KIND=4), PARAMETER :: case_methane = 1


!=======================================================================
!
!     CROSS-REFERENCE
!
!     Monte Carlo property kdf + k - 1 corresponds to dependent
!     property index_depprop(k) in look-up table
!
!     index_depprop(prho) corresponds to density in look-up table
!
!-----------------------------------------------------------------------

      INTEGER (KIND=4), ALLOCATABLE, SAVE :: index_depprop(:)  !  (ndv)
      INTEGER (KIND=4), SAVE :: prho = 0


      END MODULE cmn_lookup

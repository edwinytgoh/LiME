
      MODULE cmn_gaspdf

      USE cmn_param, ONLY : nameln


!======================================================================= 
!
!     NUMBER OF PARTICLES
!
!     np     : number of particles
!
!-----------------------------------------------------------------------

      !---  from input file
      INTEGER (KIND=4), SAVE :: np


!======================================================================= 
!
!     NUMBER OF PARTICLE PROPERTIES
!
!     nkd = 1 + nsv + ndv
!
!     "1"    : particle specific volume (kvol)
!     nsv    : number of scalar variables
!     ndv    : number of dependent variables
!     nkd    : maximum number of particle properties
!
!-----------------------------------------------------------------------

      INTEGER (KIND=4), SAVE :: nsv
      INTEGER (KIND=4), SAVE :: ndv
      INTEGER (KIND=4), SAVE :: nkd


!======================================================================= 
!
!     INDICES TO PARTICLE PROPERTIES
!
!     particle properties are contained in the array f(np,nkd).
!     - the first index of f refers to a particle.
!     - the second index of f refers to a particle property.
!
!     particle property indices
!       kvol   : specific volume
!       kf     : first scalar property (often the mixture fraction)
!       kl     : last scalar property
!       kdf    : first dependent property
!       kdl    : last dependent property
!
!-----------------------------------------------------------------------

      INTEGER (KIND=4), SAVE :: kvol, kf, kl, kdf, kdl

 
!======================================================================= 
!
!     PARTICLE PROPERTIES
!
!     f(np,nkd) :
!       particle properties are contained in the array f(np,nkd).
!       - the first index of f refers to a particle.
!       - the second index of f refers to a particle property.
!
!     iprlab :
!       array containing character labels to particle properties
!
!-----------------------------------------------------------------------
      REAL    (KIND=8), ALLOCATABLE, SAVE :: f(:,:)           ! (np,nkd)
      CHARACTER(LEN=nameln), ALLOCATABLE, SAVE :: iprlab(:)   ! (nkd)


!======================================================================= 
!
!     INITIALISATION OF INDEPENDENT SCALARS
!
!     initscal :
!       - 0, one peak (single Dirac PDF)
!       - 1, two equal peaks
!       - 2, beta-function PDF
!
!     parameters:
!         mean_mixf : mean mixture fraction (used if initscal = 0 or 2)
!         var_mixf  : mixture fraction variance (used if initscal = 2)
!         min_mixf  : 1st peak (used if initscal = 1)
!         max_mixf  : 2nd peak (used if initscal = 1)
!
!-----------------------------------------------------------------------

      INTEGER (KIND=4), SAVE :: initscal = 1      ! set to default value
      REAL    (KIND=8), SAVE :: mean_mixf = 0.5d0 ! set to default value
      REAL    (KIND=8), SAVE :: var_mixf  = 0.1d0 ! set to default value
      REAL    (KIND=8), SAVE :: min_mixf  = 0.0d0 ! set to default value
      REAL    (KIND=8), SAVE :: max_mixf  = 1.0d0 ! set to default value

 
!======================================================================= 
!
!     RELATIVE VOLUME REPRESENTED BY PARTICLES
!
!     revol :
!       we assume that each particle "i" has a unit mass: m_i = 1
!       the volume represented by one particle is then the specific
!       volume:
!          v_i = m_i / rho_i = 1 / rho_i = f(i,kvol)
!
!       The total relative volume represented by the particles is:
!
!          relvol = SUM ( f(i,kvol) )
!
!       This volume might change because of mixing.
!
!     revol_init :
!       initial value of relvol
!
!-----------------------------------------------------------------------

      REAL    (KIND=8), SAVE :: relvol
      REAL    (KIND=8), SAVE :: relvol_init


      END MODULE cmn_gaspdf

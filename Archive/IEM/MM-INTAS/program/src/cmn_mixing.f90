
      MODULE cmn_mixing

!======================================================================= 
!
!     MICRO-MIXING PROCEDURE
!
!     mixmod   : designator for mixing model to be used
!                0  - no mixing
!                1  - IEM
!                2  - Coalescence Dispersion (= modified Curl)
!                3  - Mapping closure
!
!     mixext   : designator model determining degree of mixing
!                in case of C/D model (= modified Curl)
!                1  - full mixing
!                2  - mixing uniformly distributed between 0 and 1
!                3  - distribution (b) from Pope (1982, Comb.Sci.Tech.)
!                4  - distribution (c) ""
!                5  - distribution (d) ""
!                6  - distribution (e) ""
!
!     cphi     : constant in IEM model
!
!-----------------------------------------------------------------------

      INTEGER(KIND=4), SAVE :: mixmod = 1         ! set to default value

      INTEGER(KIND=4), PARAMETER :: mx_nomix = 0
      INTEGER(KIND=4), PARAMETER :: mx_iem   = 1
      INTEGER(KIND=4), PARAMETER :: mx_codi  = 2
      INTEGER(KIND=4), PARAMETER :: mx_mcmg  = 3

      INTEGER(KIND=4), SAVE :: mixext = 2         ! set to default value

      REAL   (KIND=8), SAVE :: cphi = 2.d0        ! set to default value

!======================================================================= 
!
!     TURBULENCE FREQUENCY
!
!-----------------------------------------------------------------------

      REAL   (KIND=8), SAVE :: omega


      END MODULE cmn_mixing

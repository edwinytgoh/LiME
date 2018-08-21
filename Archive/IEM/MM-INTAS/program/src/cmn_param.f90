
      MODULE cmn_param

!=======================================================================  
!    
!     STRING LENGTHS
!
!     nameln : length of particle property names
!     fileln : length of filenames
!
!-----------------------------------------------------------------------
      
      INTEGER (KIND=4), PARAMETER :: nameln = 30
      INTEGER (KIND=4), PARAMETER :: fileln = 120


!=======================================================================  
!
!     STANDARD REAL PARAMETERS
!
!-----------------------------------------------------------------------

      REAL    (KIND=8), PARAMETER :: big  = 1.d30
      REAL    (KIND=8), PARAMETER :: small= 1.d-30

      REAL    (KIND=8), PARAMETER :: zero  = 0.d0
      REAL    (KIND=8), PARAMETER :: one   = 1.d0
      REAL    (KIND=8), PARAMETER :: two   = 2.d0
      REAL    (KIND=8), PARAMETER :: three = 3.d0
      REAL    (KIND=8), PARAMETER :: four  = 4.d0
      REAL    (KIND=8), PARAMETER :: ten   = 10.d0
      REAL    (KIND=8), PARAMETER :: d1d2  = 0.5d0
      REAL    (KIND=8), PARAMETER :: d3d2  = 1.5d0
      REAL    (KIND=8), PARAMETER :: d1d3  = 1.d0/3.d0
      REAL    (KIND=8), PARAMETER :: d2d3  = 2.d0/3.d0
      REAL    (KIND=8), PARAMETER :: d4d3  = 4.d0/3.d0
      REAL    (KIND=8), PARAMETER :: d1d4  = 0.25d0
      REAL    (KIND=8), PARAMETER :: d3d4  = 0.75d0
      REAL    (KIND=8), PARAMETER :: half  = d1d2
      REAL    (KIND=8), PARAMETER :: one3rd = d1d3
      REAL    (KIND=8), PARAMETER :: two3rd = d2d3
      REAL    (KIND=8), PARAMETER :: one4th = d1d4
      REAL    (KIND=8), PARAMETER :: pi = 3.14159265358979d0
      REAL    (KIND=8), PARAMETER :: twopi = two*pi


      END MODULE cmn_param

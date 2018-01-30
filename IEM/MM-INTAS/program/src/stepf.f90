!#1 subroutine stepf

!#1
!=======================================================================
!  subroutine stepf
!=======================================================================
!  function :
!    routine to advance scalars for one time step, according to
!    designated model.
!
!  variable mixmod designates model :
!    0 - no mixing
!    1 - IEM
!    2 - Coalescence Dispersion
!    3 - Mapping closure
!
!  input :
!    dt            : timestep
!    omega         : mean turbulent frequency
!    np            : number of particles
!    npd           : first dimension size of array f
!    nsv           : number of independent scalars
!    f_kf(npd,nsv) : particle independant scalars
!
!  output :
!    f_kf(npd,nsv) : particle independant scalars
!
!  calls :
!    iem
!    codi
!    mcmg
!-----------------------------------------------------------------------
      SUBROUTINE stepf (loutf, llogf, &
                        dt, omega, np, npd, nsv, &
                        f_kf)

!-----  used modules
      USE cmn_mixing, ONLY : mixmod, mixext, cphi, &
                             mx_nomix, mx_iem, mx_codi, mx_mcmg

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: loutf, llogf
      INTEGER (KIND=4), INTENT(IN) :: np, npd, nsv

      REAL    (KIND=8), INTENT(IN) :: dt
      REAL    (KIND=8), INTENT(IN) :: omega

!----- input/output variables
      REAL    (KIND=8), INTENT(INOUT) :: f_kf(npd,nsv)

!----- local variables
      INTEGER (KIND=4) :: i

      !--- we consider only 1 cell
      INTEGER (KIND=4), PARAMETER :: ncell = 1
      INTEGER (KIND=4) :: nf(ncell), nl(ncell)
      REAL    (KIND=8) :: omean_cell(ncell)
      REAL    (KIND=8) :: f_kwt(np)


!----- MAIN ACTION

      !--- We consider only ONE cell
      nf(1) = 1               ! first particle in the cell
      nl(1) = np              ! last particle in the cell
      omean_cell(1) = omega   ! mean turbulence frequency in the cell
      DO i = 1, np
         f_kwt(i) = 1.d0      ! all particles have the same mass
      END DO
 

!----- choice between models
      SELECT CASE (mixmod)

      CASE (mx_nomix)
         RETURN

      CASE (mx_iem)
         WRITE(llogf,'(2X,A)') 'STEPF: calling iem'
         CALL iem (dt, cphi, &
                   ncell, &
                   omean_cell, &
                   np, npd, nsv, nf, nl, &
                   f_kwt, f_kf)

      CASE (mx_codi)
         WRITE(llogf,'(2X,A)') 'STEPF: calling codi'
         CALL codi (loutf, &
                    mixext, dt, cphi, &
                    ncell, &
                    omean_cell, &
                    np, npd, nsv, nf, nl, &
                    f_kwt, f_kf)

      CASE (mx_mcmg)
         WRITE(llogf,'(2X,A)') 'STEPF: calling mcmg'
         CALL mcmg(loutf, &
                   dt, cphi, &
                   ncell, &
                   omean_cell, &
                   np, npd, nsv, nf, nl, &
                   f_kwt, f_kf)

      CASE DEFAULT
         WRITE(loutf, *) 'Illegal value for mixmod in routine stepf', mixmod
         STOP

      END SELECT


      END SUBROUTINE stepf

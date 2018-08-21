!#1 subroutine iem

!#1
!=======================================================================
!  subroutine iem
!=======================================================================
!  function :
!    routine to advance scalars for one time step.
!    physical model : evolution towards the mean value
!
!  input :
!    dt            : timestep
!    cphi          : constant in IEM model
!    ncell         : number of cells in the domain
!    omean_cell(l) : mean turbulent frequency in cell l
!    np            : number of particles
!    npd           : first dimension size of array f
!    nsv           : number of independent scalars
!    nf(l)         : first particle in cell l
!    nl(l)         : last particle in cell l
!    f_kwt(np)     : particle weights
!    f_kf(npd,nsv) : particle independent scalars
!
!  output :
!    f_kf(npd,nsv) : particle independent scalars
!
!  includes:
!    cmn_param
!-----------------------------------------------------------------------
      SUBROUTINE iem (dt, cphi, &
                      ncell, &
                      omean_cell, &
                      np, npd, nsv, nf, nl, &
                      f_kwt, f_kf)

!-----  used modules
      USE cmn_param, ONLY : half

      IMPLICIT NONE

!----- input variables
      INTEGER (KIND=4), INTENT(IN) :: ncell
      INTEGER (KIND=4), INTENT(IN) :: np, npd, nsv
      INTEGER (KIND=4), INTENT(IN) :: nf(ncell)
      INTEGER (KIND=4), INTENT(IN) :: nl(ncell)
      REAL    (KIND=8), INTENT(IN) :: dt, cphi
      REAL    (KIND=8), INTENT(IN) :: omean_cell(ncell)
      REAL    (KIND=8), INTENT(IN) :: f_kwt(np)

!----- input/output variables
      REAL    (KIND=8), INTENT(INOUT) :: f_kf(npd,nsv)

!----- local variables
      INTEGER (KIND=4) :: i, k, lc
      REAL    (KIND=8) :: omean, decay, fmean
      REAL    (KIND=8) :: sumwt, sumf


!----- MAIN ACTION


      !--- loop over cells 
      DO lc = 1, ncell

         !--- calculate total weight in cell
         sumwt = 0.d0
         DO i = nf(lc), nl(lc)
            sumwt = sumwt + f_kwt(i)
         END DO
         IF (sumwt <= 0.d0) CYCLE

         !--- calculate mean omega in this cell 
         omean = omean_cell(lc)

         !--- calculate mean decay rate in this cell
         decay = DEXP(-half * cphi * omean * dt)

!-------- QUESTION: is the factor 2.d0 correct?
!         ANSWER  : yes, the variance of the scalar decays as
!         exp(-cphi * omean * dt), so the distance from the
!         mean decays as exp(-cphi/2 * omean * dt)


         !--- loop over scalars
         DO k = 1, nsv

            !--- calculate mean value of scalar in this cell
            sumf = 0.d0

            !--- sum scalar over particles in this cell
            DO i = nf(lc), nl(lc)
               sumf = sumf + f_kwt(i) * f_kf(i,k)
            END DO

            fmean = sumf / sumwt

            !--- second step : calculate step in f
            !    (physical model : decay to the cell mean)

            !--- loop over particles in cell
            DO i = nf(lc), nl(lc)
               f_kf(i,k) = fmean + (f_kf(i,k) - fmean) * decay
            END DO
         END DO

      END DO


      END SUBROUTINE iem

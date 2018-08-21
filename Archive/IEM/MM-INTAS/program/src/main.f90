
!======================================================================!
!                                                                      !
!  Main program for MM-INTAS                                           !
!                                                                      !
!                                                                      !
!  Author: Bertrand Naud                                               !
!                                                                      !
!                                                                      !
!  Contact: Prof. Dirk Roekaerts   /   dirkr@ws.tn.tudelft.nl          !
!                                                                      !
!           Delft University of Technology                             !
!           Thermal and Fluids Sciences Section                        !
!           Lorentzweg 1                                               !
!           2628 CJ  Delft                                             !
!           THE NETHERLANDS                                            !
!                                                                      !
!                                                                      !
!                                                                      !
!  1. INTRODUCTION                                                     !
!                                                                      !
!  The computer program MM-INTAS was first provided to the groups      !
!  contributing to a INTAS project (www.intas.be), grant 00-353.       !
!  In this project, groups from Delft (The Netherlands), Moscow        !
!  (Russia), Minsk (White Russia), Zaragoza (Spain), Marseille and     !
!  Rouen (France) are involved. The project concerns the development   !
!  of advanced models for turbulent mixing.                            !
!                                                                      !
!  This computer program was provided as a tool in order to test       !
!  micro-mixing models in the context of Lagrangian modeling. It is    !
!  meant to be modified for one's specific needs. It was written from  !
!  the Monte-Carlo PDF code "PDFD" (Delft University of Technology;    !
!  authors: D. Roekaerts, P.A. Nooren, H.A. Wouters, B. Naud; based on !
!  the routine "pdf2ds" by S.B. Pope). Only a small part of the        !
!  original code remains here.                                         !
!                                                                      !
!  Three mixing models are available here: IEM/LMSE model, modified    !
!  Curl's model and mapping closure model. The subroutines available   !
!  for those three models correspond to the implementations of Nooren  !
!  and Wouters. For more information on those models, see:             !
!       [1] Nooren, P.A. (1998) "Stochastic modeling of turbulent      !
!           natural-gas flames".  PhD thesis, Delft University of      !
!           Technology.                                                !
!       [2] Wouters, H.A. (1998) "Lagrangian models for turbulent      !
!           reacting flows".  PhD thesis, Delft University of          !
!           Technology.                                                !
!       [3] Nooren, P.A., H.A. Wouters, T.W.J. Peeters, D. Roekaerts,  !
!           U. Maas and D. Schmidt (1997). Monte Carlo PDF modeling of !
!           a turbulent natural-gas diffusion flame.                   !
!           Combust. Theory Modelling 1, 79-96.                        !
!  Note that although all particles have the same computational weight !
!  in this code, we kept the general implementation of the models for  !
!  non-uniform particle weights.                                       !
!                                                                      !
!  B. Naud developed this computer code such that it should be         !
!  (hopefully) easy to use, modify and extend.                         !
!                                                                      !
!                                                                      !
!  2. DESCRIPTION OF THE COMPUTER PROGRAM                              !
!                                                                      !
!  We consider the evolution of the composition PDF in a small domain  !
!  where a frozen homogeneous turbulence exists. This turbulent field  !
!  is characterized by a turbulence frequency "omega":                 !
!     omega = epsilon / k                                              !
!     with: epsilon, dissipation of turbulent kinetic energy           !
!            and  k, turbulent kinetic energy                          !
!  The user specifies the number 'nstep' and the duration 'dt' of      !
!  time-steps to be performed.                                         !
!                                                                      !
!  The composition PDF is represented by 'np' stochastic particles:    !
!  each particle has a composition and a density. All particles have   !
!  the same mass which remains constant during the calculation. The    !
!  volume of particles (and therefore the volume of the domain         !
!  considered) may change.                                             !
!                                                                      !
!  The composition PDF (i.e. the evolution of the composition of each  !
!  of the 'np' particles) evolves due to mixing and possibly because   !
!  of chemical reaction.                                               !
!                                                                      !
!  We consider 'nsv' scalar variables from which dependent scalars can !
!  be obtained via state equations ('ndv' dependent variables are      !
!  considered).                                                        !
!                                                                      !
!  The information is stored in an array f(np,nkd), where np is the    !
!  number of particles considered and nkd=nsv+ndv is the number of     !
!  scalar properties carried by each particle.                         !
!  (actually nkd=nsv+ndv+1: the extra property being the inverse of    !
!   density)                                                           !
!                                                                      !
!  ------------------------------------------------------------------  !
! | RESTRICTION: In the current version of the code, only one scalar | !
! | variable is considered (nsv=1): mixture fraction. This conserved | !
! | scalar only evolves due to mixing.                               | !
!  ------------------------------------------------------------------  !
!                                                                      !
!                                                                      !
!  3. WHAT CAN BE DONE WITH THE CURRENT VERSION                        !
!                                                                      !
!  3.1. "Test case" possibilities.                                     !
!                                                                      !
!    \A. The program can be used to test mixing models for a constant  !
!  density flow by considering mixture fraction only (nsv=1 / ndv=0).  !
!                                                                      !
!    \B. Two look-up tables are provided, giving properties of a       !
!  methane/air mixture (at 298K) as a function of mixture fraction,    !
!  corresponding to two situations: 1. pure mixing, and 2. chemical    !
!  equilibrium.                                                        !
!  Different species mass fractions, temperature and density can be    !
!  retrieved from those tables as function of mixture fraction         !
!  (nsv=1 / ndv>0).                                                    !
!  The program can be used to test mixing models in those variable     !
!  density flow cases and the evolution of the shape of e.g. the       !
!  temperature PDF induced by the mixing model applied on mixture      !
!  fraction can be studied.                                            !
!                                                                      !
!  3.2. Initial conditions.                                            !
!                                                                      !
!  The initial distribution for the mixing scalar (mixture fraction)   !
!  can be chosen as:                                                   !
!                                                                      !
!    \A. one peak (useless for mixing studies)                         !
!                                                                      !
!    \B. two equal peaks                                               !
!                                                                      !
!    \C. beta-function PDF corresponding to a given mean and variance. !
!        (note that in this case, it is not sure that the mean scalar  !
!         value obtained from the particle ensemble exactly corresponds!
!         to the specified mean value, due to statistical error)       !
!                                                                      !
!  Case \B is a standard case in order to test mixing models.          !
!  Case \C can be interesting in order to evaluate the performance of  !
!  mixing models in "real" simulation situations: how will a given     !
!  mixing model influence the original beta-shape of the composition   !
!  PDF in comparison with other mixing models.                         !
!                                                                      !
!  3.3 Outputs of the program                                          !
!                                                                      !
!  The output of the program is simply the list of the scalar values   !
!  of all particles at different time steps (1 file per time step).    !
!  The user can extract the wanted information (PDFs, scatter plots,   !
!  ...) from those files.                                              !
!                                                                      !
!                                                                      !
!  4. POSSIBLE EXTENSIONS                                              !
!                                                                      !
!  Some straightforward extensions of this computer code should be     !
!  possible:                                                           !
!    - more than one mixing scalar                                     !
!    - different initial conditions for the scalar PDF                 !
!    - new mixing models                                               !
!    - non-uniform particle weights                                    !
!    - "processed" outputs                                             !
!    - ...                                                             !
!                                                                      !
!======================================================================!


!=======================================================================
      PROGRAM MM_INTAS
!=======================================================================

      USE cmn_param, ONLY : nameln, fileln

      IMPLICIT NONE

      !---  Output files, Input file
      INTEGER (KIND=4) :: loutf, llogf, linpt
      CHARACTER(LEN=fileln) :: inpt

      !--- Case studied
      INTEGER (KIND=4) :: icasestudied

      INTEGER (KIND=4) :: nindv, ndepv
      CHARACTER(LEN=nameln), ALLOCATABLE :: name_indprop(:)   ! (nindv)
      CHARACTER(LEN=nameln), ALLOCATABLE :: name_depprop(:)   ! (ndepv)

      !---  Time step/iterations from input file
      !  nstep : user-specified number of timesteps to be performed
      !  dt    : user-specified time-step
      INTEGER (KIND=4) :: nstep
      REAL    (KIND=8) :: dt


!======  MAIN ACTION

!----------------------------------------------------------------------- 
!  CHAPTER 1 : INITIALISATION
!----------------------------------------------------------------------- 

      CALL init_files (loutf, llogf, linpt, inpt)

   !---------------------------
   !  INPUTS OF THE PROBLEM
   !---------------------------

      icasestudied = - 1
      nstep        = - 1
      dt           = - 1.d0

      WRITE(*,*) '----------------'
      WRITE(*,*) 'Call input'
      CALL read_input (linpt, inpt, loutf, icasestudied, nstep, dt)
      WRITE(*,*) 'Exit input'
      WRITE(*,*) '----------------'
      CALL check_input (loutf, icasestudied, nstep, dt)


   !-----------------------------------------------
   !  INITIALISE PROPERTIES, INITIALISE LOOK-UP
   !-----------------------------------------------

      !---- Number of properties
      !---  depending on the case studied, set the number of dependent
      !     scalars and independent sacalars
      CALL init_lookup_size (icasestudied, nindv, ndepv)

      !---  set the number of particle properties and the indices of
      !     those properties (kvol, kf, ...) in the array f
      CALL set_property_indices (loutf, nindv, ndepv)

      !---- Names of properties
      !---  depending on the case studied, set the names of dependent
      !     scalars and independent scalars,
      !     and build the cross reference table (index_depprop) when
      !     a look-up table is used
      ALLOCATE(name_indprop(nindv))
      ALLOCATE(name_depprop(ndepv))
      CALL init_lookup (icasestudied, nindv, ndepv, &
                        name_indprop, name_depprop)

      !---  set the names of particle properties (array iprlab)
      CALL set_property_names (name_indprop, name_depprop)


!----------------------------------------------------------------------- 
!  CHAPTER 2 : STARTING VALUES
!----------------------------------------------------------------------- 

      WRITE(*,'(a)')
      WRITE(*,'(a)')
      WRITE(*,'(a)')
      WRITE(*,'(a)') 'STARTING VALUES'


   !------------------------------------------------------------
   !  INITIALISE THE MONTE CARLO PROCEDURE
   !   (the number of particles was read from the input file)
   !   initialise particle property arrays
   !------------------------------------------------------------
      CALL init_mc (loutf, llogf, icasestudied, nindv)


!----------------------------------------------------------------------- 
!  CHAPTER 3 : MAIN LOOP
!----------------------------------------------------------------------- 
      WRITE(*,'(a)')
      WRITE(*,'(a)')
      WRITE(*,'(a)')
      WRITE(*,'(a)') 'MAIN LOOP'

      CALL mainloop (loutf, llogf, icasestudied, nstep, dt)


!----------------------------------------------------------------------- 
!  CHAPTER 4  EXIT
!----------------------------------------------------------------------- 

!------  write tails and close files: output and loginf
      CALL close_infofiles (loutf, llogf)


!=======================================================================
      END PROGRAM MM_INTAS
!=======================================================================

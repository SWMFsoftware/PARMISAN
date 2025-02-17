  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModParticle
   use PT_ModConst
   use PT_ModFieldline, ONLY: get_random_shock_location, setup_lagrangian_grid, &
      compute_weight, compute_timestep, update_lagrangian_grid, nS, nT, &
      rMin, rMax, tMax, nDim

   use PT_ModSolver, ONLY: rk2_sde, milstein
   implicit none
   SAVE

   real, allocatable :: Particle_V(:)
   integer :: nParticlePerProc

   integer, parameter :: LagrCoord_    = 1,  &
                         Momentum_     = 2,  &                      
                         Time_         = 3,  &
                         TimeStep_     = 4,  &
                         S_            = 5,  &
                         R_            = 6,  &
                         Energy_       = 7,  &
                         Weight_       = 8,  &
                         SOld_         = 9,  &
                         ROld_         = 10, &
                         LagrCoordOld_ = 11, &
                         MomentumOld_  = 12, &
                         nVar          = 12

   real         :: E0             ! = 10.0*keV
   integer      :: nSplitLev, nSplitMax      ! defaults: 40, 80
   real         :: eSplitLevelMin ! = 1.d0*MeV      ! energy of first split level
   real         :: eSplitLevelMax ! = 20000.d0*MeV  ! energy of last split level
   real, allocatable :: eSplitLev_I(:)
contains
   !============================================================================
   subroutine read_param(NameCommand)

      use ModReadParam, ONLY: read_var
      use ModUtilities, ONLY: CON_stop

      character(len=*), intent(in):: NameCommand ! From PARAM.in
      character(len=*), parameter:: NameSub = 'read_param'
      !--------------------------------------------------------------------------
      select case(NameCommand)
      case('#PARTICLE')
         call read_var('nParticlePerProc', nParticlePerProc )
         call read_var('initialEnergy', E0)
         call read_var('nSplitLev', nSplitLev)
         call read_var('nSplitMax', nSplitMax)
         call read_var('eSplitLevelMin', eSplitLevelMin)
         call read_var('eSplitLevelMax', eSplitLevelMax)

         ! Convert to erg
         E0 = E0*ckeV
         eSplitLevelMax = eSplitLevelMax*cMeV
         eSplitLevelMin = eSplitLevelMin*cMeV
         
      case default
         call CON_stop(NameSub//' Unknown command '//NameCommand)
      end select

   end subroutine read_param
   !============================================================================
   subroutine initialize_particle()

      if(.not.allocated(Particle_V)) allocate(Particle_V(1:nVar)) ! modparticle init function?

      ! Get particle initial time, radial distance, and distance along fieldline
      ! Samples from 1/r^2 distribution in [rMin, rMax]    
      call get_random_shock_location(Particle_V(Time_), Particle_V(R_), Particle_V(S_))
      
      ! Get initial momentum
      Particle_V(Energy_) = E0
      call energy_to_momentum(E0, Particle_V(Momentum_))
      
      ! Get initial lagrangian coordinate and set up lagrangian grid
      call setup_lagrangian_grid(Particle_V(Time_), Particle_V(S_), Particle_V(LagrCoord_))

      ! Get initial timestep
      call compute_timestep(Particle_V(Time_), Particle_V(LagrCoord_), & 
                            Particle_V(Momentum_), Particle_V(TimeStep_))

      ! Get initial factor (B/ds)
      call compute_weight(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(Weight_))

      ! Save particle's current position as previous position
      call save_state()

   end subroutine initialize_particle
   !============================================================================
   subroutine init_split_particle(Time, S, Momentum, R)

      real, intent(in) :: Time, S, Momentum, R
      
      Particle_V(Time_) = Time
      Particle_V(S_) = S
      Particle_V(R_) = R
      Particle_V(Momentum_) = Momentum

      call momentum_to_energy(Momentum, Particle_V(Energy_))

      ! Get initial lagrangian coordinate and set up lagrangian grid
      call setup_lagrangian_grid(Particle_V(Time_), Particle_V(S_), Particle_V(LagrCoord_))
      ! Get initial timestep
      call compute_timestep(Particle_V(Time_), Particle_V(LagrCoord_), & 
                            Particle_V(Momentum_), Particle_V(TimeStep_))

      call compute_weight(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(Weight_))

      call save_state

   end subroutine init_split_particle
   !============================================================================
   subroutine save_state()

      Particle_V(SOld_) = Particle_V(S_)
      Particle_V(ROld_) = Particle_V(R_)
      Particle_V(LagrCoordOld_) = Particle_V(LagrCoord_)
      Particle_V(MomentumOld_) = Particle_V(Momentum_)

   end subroutine save_state
   !============================================================================
   subroutine advance_particle()

      real :: XOld(nDim), XNew(nDim)
      ! Save particle's current position as previous position
      call save_state

      ! Solve SDE - updates Lagrangian coordinate and momentum
      ! Equation advances p^3/3 not p
      XOld(LagrCoord_) = Particle_V(LagrCoord_)
      XOld(Momentum_)  = Particle_V(Momentum_)**3.0 / 3.0
      
      ! Timestep is calculated after moving particle below
      call rk2_sde(XOld, Particle_V(Time_), Particle_V(TimeStep_), XNew)
      
      ! Milstein method allows for calculation of timestep and weight during 
      ! timestep.
      !call milstein(XOld, Particle_V(Time_), Particle_V(TimeStep_), Particle_V(Weight_), XNew)

      Particle_V(LagrCoord_) = XNew(LagrCoord_)
      Particle_V(Momentum_) = (3.0*XNew(Momentum_))**(1.0/3.0)

      call momentum_to_energy(Particle_V(Momentum_), Particle_V(Energy_))

      ! make this more robust ?
      if(Particle_V(LagrCoord_) < 1.0 .or. Particle_V(LagrCoord_) > nS - 1.0) return

      ! Update distance along field line (S) and radial distance (R)
      call update_lagrangian_grid(Particle_V(LagrCoord_), Particle_V(TimeStep_), &
                                  Particle_V(S_), Particle_V(R_))

      ! Update particle time then calculate new timestep
      Particle_V(Time_) = Particle_V(Time_) + Particle_V(TimeStep_)
      call compute_timestep(Particle_V(Time_), Particle_V(LagrCoord_), & 
                            Particle_V(Momentum_), Particle_V(TimeStep_))

      ! Update particle weight
      call compute_weight(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(Weight_))

   end subroutine advance_particle
   !============================================================================
   subroutine init_split_grid()
   use PT_ModPlot, ONLY: OutputDir, splitFile
   use PT_ModProc, ONLY: iProc
   real :: dLogEsplit
   ! Loop variable:
   integer :: iLev
   !--------------------------------------------------------------------------

   if(.not.allocated(eSplitLev_I)) allocate(eSplitLev_I(1:nSplitLev+1))

   eSplitLev_I(1)           = eSplitLevelMin
   eSplitLev_I(1+nSplitLev) = eSplitLevelMax
   dLogEsplit = log10(eSplitLevelMax/eSplitLevelMin)/real(nSplitLev)

   do iLev = 2, nSplitLev
      eSplitLev_I(iLev) = eSplitLev_I(iLev-1)*10**dLogEsplit
   end do

   ! Write energy splitting grid to file
   if(iProc==0)then
      open(45,file=OutputDir//splitFile,status='unknown')
      do iLev = 1, nSplitLev+1
         write(45,*)iLev, eSplitLev_I(iLev)/cMeV
      end do
      close(45)
   end if

   end subroutine init_split_grid
   !============================================================================
   subroutine check_boundary_conditions(IsOutside)
      logical, intent(out) :: IsOutside

      IsOutside = .False.
      if(Particle_V(R_) < rMin .or. Particle_V(R_) > rMax) IsOutside = .True.
      if(Particle_V(LagrCoord_) < 1.0 .or. Particle_V(LagrCoord_) > nS-1.0) IsOutside = .True.
      if((Particle_V(Time_) + Particle_V(TimeStep_)) > tMax) IsOutside = .True.

   end subroutine check_boundary_conditions
   !============================================================================
   subroutine get_random_initial_energy(Energy)
      real, intent(out) :: Energy
      real :: RandUniform
      real :: Emin, Emax
      real, parameter :: cThreeHalf = 1.5
      ! outputs random energy in [keV]
      ! f(E) is proportional to v^-5 or E^(-5/2)
      ! Random selection is done in energy to get better statistics
      ! at higher energies

      ! energy limits [keV]
      Emin = 10.0
      Emax = 100.0

      call random_number(RandUniform)

      Energy = (1.0 - RandUniform) * Emin**(-cThreeHalf) + &
               RandUniform * Emax**(-cThreeHalf)
      Energy = Energy**(-1.0/cThreeHalf)

   end subroutine get_random_initial_energy
   !=============================================================================
   subroutine energy_to_momentum(Energy, Momentum)
      real, intent(in) :: Energy  ! erg 
      real, intent(out) :: Momentum ! g*cm/s

      Momentum = sqrt(2.0*cProtonMass*Energy + (Energy/cLightSpeed)**2)

   end subroutine energy_to_momentum
   !=============================================================================
   subroutine momentum_to_energy(Momentum, Energy)
      real, intent(in) :: Momentum   ! g*cm/s 
      real, intent(out) :: Energy  ! erg

      Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - &
               cProtonRestEnergy

   end subroutine momentum_to_energy
   !=============================================================================
   subroutine debug_particle(Message)
      character(len = *), intent(in) :: Message

      write(*,*) '#', repeat('-', 25), '#'
      write(*,*) Message
      write(*, '(A, F14.2, A)') '  Time: ', Particle_V(Time_), ' s' 
      write(*, '(A, F10.2, A)') '  Timestep: ', Particle_V(TimeStep_), ' s'
      write(*, '(A, F10.2, A)') '  Energy: ', Particle_V(Energy_)/ckeV, ' keV'
      write(*, '(A, F16.2, A)') '  R: ', Particle_V(R_)/cRsun, ' Rs'
      write(*, '(A, F16.2, A)') '  S: ', Particle_V(S_)/cRsun, ' Rs'
      write(*, '(A, F11.2, A)') '  LagrCoord: ', Particle_V(LagrCoord_)
      write(*,*) '#', repeat('-', 25), '#'

   end subroutine
   !=============================================================================
   end module PT_ModParticle
   !=============================================================================

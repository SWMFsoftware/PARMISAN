  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModParticle
   use PT_ModConst   
   use ModKind
   ! use PT_ModFieldline, ONLY: get_random_shock_location, compute_weight, &
   !                            get_particle_location, nS, rMin, rMax, tMax, nDim
   use PT_ModTestFieldline, ONLY: get_random_shock_location, compute_weight, &
                                  get_particle_location, nS, rMin, rMax, tMax, nDim

   use PT_ModSolver, ONLY: euler_sde, milstein_sde
   
   implicit none
   SAVE

   real, allocatable :: Particle_V(:,:)
   integer :: nParticlePerProc

   integer, parameter :: LagrCoord_    = 1,  &
                         Momentum_     = 2,  &                      
                         Time_         = 3,  &
                         TimeStep_     = 4,  &
                         R_            = 5,  &
                         Energy_       = 6,  &
                         Weight_       = 7,  &
                         ROld_         = 8,  &
                         LagrCoordOld_ = 9,  &
                         MomentumOld_  = 10, &
                         nVar          = 10

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

         ! Convert to joules
         E0 = E0*ckeV
         eSplitLevelMax = eSplitLevelMax*cMeV
         eSplitLevelMin = eSplitLevelMin*cMeV
         
      case default
         call CON_stop(NameSub//' Unknown command '//NameCommand)
      end select

   end subroutine read_param
   !============================================================================
   subroutine initialize_particles()
      integer :: i
      
      if(.not.allocated(Particle_V)) allocate(Particle_V(1:nParticlePerProc, 1:nVar)) ! modparticle init function?

      ! Get particle initial time, radial distance, and distance along fieldline    
      do i = 1, nParticlePerProc
         call get_random_shock_location(Particle_V(i, Time_), Particle_V(i, LagrCoord_), Particle_V(i, R_))
         Particle_V(i, Energy_) = E0
         call energy_to_momentum(Particle_V(i, Energy_), Particle_V(i, Momentum_))
         ! Get weight: T*rho
         ! TODO: check to see if rho is proton only, if not, need number density of protons
         call compute_weight(Particle_V(i, Time_), Particle_V(i, LagrCoord_), Particle_V(i, Weight_))
         call save_state(i)
      end do

   end subroutine initialize_particles
   !============================================================================
   ! subroutine init_split_particle(Time, LagrCoord, Momentum, splitLev)

         ! real, intent(in) :: Time, LagrCoord, Momentum
         ! integer, intent(in) :: splitLev

         ! Particle_V(Time_) = Time
         ! Particle_V(LagrCoord_) = LagrCoord
         ! Particle_V(Momentum_) = Momentum

         ! call momentum_to_energy(Momentum, Particle_V(Energy_))

         ! call get_particle_location(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(R_))

         ! call compute_weight(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(Weight_))
         ! Particle_V(Weight_) = Particle_V(Weight_)*(2.d0**(-real(splitLev-1)))

         ! call save_state()

   ! end subroutine init_split_particle
   !============================================================================
   subroutine save_state(iParticle)
      integer, intent(in) :: iParticle

      Particle_V(iParticle, ROld_) = Particle_V(iParticle, R_)
      Particle_V(iParticle, LagrCoordOld_) = Particle_V(iParticle, LagrCoord_)
      Particle_V(iParticle, MomentumOld_) = Particle_V(iParticle, Momentum_)

   end subroutine save_state
   !============================================================================
   subroutine advance_particle(iParticle)
      integer, intent(in) :: iParticle

      ! Vectors solved by SDE
      real :: XOld(nDim), XNew(nDim)
      
      ! Save particle's current position as previous position
      call save_state(iParticle)

      ! Equation advances lagrcoord and p^3/3 not p
      XOld(LagrCoord_) = Particle_V(iParticle, LagrCoord_)
      XOld(Momentum_)  = Particle_V(iParticle, Momentum_)**3.0 / 3.0
      
      ! Advance psuedo-particle one time step
      call euler_sde(XOld, Particle_V(iParticle, Time_), Particle_V(iParticle, TimeStep_), XNew)
      ! call milstein_sde(XOld, Particle_V(Time_), Particle_V(TimeStep_), XNew)

      ! Update Lagrangian coordinate and momentum
      Particle_V(iParticle, LagrCoord_) = XNew(LagrCoord_)
      Particle_V(iParticle, Momentum_) = (3.0*XNew(Momentum_))**(1.0/3.0)
      call momentum_to_energy(Particle_V(iParticle, Momentum_), Particle_V(iParticle, Energy_))

      ! Update particle time
      Particle_V(iParticle, Time_) = Particle_V(iParticle, Time_) + Particle_V(iParticle, TimeStep_)

      ! reflection at inner boundary
      if(Particle_V(iParticle, LagrCoord_) < 2.0) Particle_V(iParticle, LagrCoord_) = 4.0 - Particle_V(iParticle, LagrCoord_)

      ! Update radial distance (R)
      call get_particle_location(Particle_V(iParticle, Time_), Particle_V(iParticle, LagrCoord_), Particle_V(iParticle, R_))

   end subroutine advance_particle
   !============================================================================
   subroutine init_split_grid()
      
      use PT_ModProc, ONLY: iProc
      
      character(len=*), parameter :: OutputDir = 'PT/IO2/'
      character(len=*), parameter :: SplitFile = 'split_levels.dat'

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
   subroutine check_boundary_conditions(iParticle, IsOutside)
      integer, intent(in) :: iParticle
      logical, intent(out) :: IsOutside

      IsOutside = .False.
      if(Particle_V(iParticle, R_) < rMin .or. Particle_V(iParticle, R_) > rMax) IsOutside = .True.
      if(Particle_V(iParticle, LagrCoord_) < 2.0 .or. Particle_V(iParticle, LagrCoord_) > nS-2.0) IsOutside = .True.
      if((Particle_V(iParticle, Time_) + Particle_V(iParticle, TimeStep_)) > tMax) IsOutside = .True.

   end subroutine check_boundary_conditions
   !============================================================================
   subroutine get_random_initial_energy(Energy)
      real, intent(out) :: Energy
      real :: RandUniform
      real :: Emin, Emax
      real, parameter :: cThreeHalf = 1.5
      ! outputs random energy in [J]
      ! f(E) is proportional to v^-5 or E^(-5/2)
      ! Random selection is done in energy to get better statistics
      ! at higher energies

      ! energy limits [keV]
      Emin = 10.0
      Emax = 100.0

      call random_number(RandUniform)

      Energy = (1.0 - RandUniform) * Emin**(-cThreeHalf) + &
               RandUniform * Emax**(-cThreeHalf)
      Energy = Energy**(-1.0/cThreeHalf) * ckeV

   end subroutine get_random_initial_energy
   !=============================================================================
   subroutine energy_to_momentum(Energy, Momentum)
      real, intent(in) :: Energy  ! Joules 
      real, intent(out) :: Momentum ! kg*m/s

      Momentum = sqrt(2.0*cProtonMass*Energy + (Energy/cLightSpeed)**2)

   end subroutine energy_to_momentum
   !=============================================================================
   subroutine momentum_to_energy(Momentum, Energy)
      real, intent(in) :: Momentum   ! kg*m/s 
      real, intent(out) :: Energy  ! joules

      Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - &
               cProtonRestEnergy

   end subroutine momentum_to_energy
   !=============================================================================
   end module PT_ModParticle
   !=============================================================================

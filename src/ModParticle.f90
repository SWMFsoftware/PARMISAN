  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModParticle
   use PT_ModConst, only: cProtonRestEnergy, cProtonMass, cLightSpeed, &
                          cMu, cMeV, ckeV
   use PT_ModSize, ONLY: nDim
   use PT_ModTime, ONLY: PTTime

   use PT_ModSolver, ONLY: euler_sde
   
   implicit none
   SAVE

   real, allocatable :: Particle_IV(:,:)
   integer :: nInject
   integer, allocatable :: nParticleOnLine(:)

   integer, parameter :: LagrCoord_    = 1,  &
                         Momentum_     = 2,  &                      
                         Time_         = 3,  &
                         TimeStep_     = 4,  &
                         R_            = 5,  &
                         Energy_       = 6,  &
                         Weight_       = 7,  &
                         ParentIndex_  = 8,  &
                         NumChildren_  = 9,  &
                         SplitLevel_   = 10, &
                         nSplitVar     = 10, & ! number of variables if splitting
                         nVar          = 7     ! number of variables if not splitting

   real         :: E0             ! = 10.0*keV
   integer      :: nSplitLev, nSplitMax      ! defaults: 40, 80
   real         :: eSplitLevelMin ! = 1.d0*MeV      ! energy of first split level
   real         :: eSplitLevelMax ! = 20000.d0*MeV  ! energy of last split level
   real, allocatable :: eSplitLev_I(:)
   logical :: UseSplit = .false.

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
         call read_var('nInject', nInject)
         call read_var('initialEnergy', E0)
         call read_var('UseSplit', UseSplit)
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
   subroutine init
      use PT_ModGrid, only: nLine
      use PT_ModSize, only: nVertexMax
      integer :: MaxTotalParticles

      MaxTotalParticles = nVertexMax*nInject
      if(UseSplit) then
         MaxTotalParticles = MaxTotalParticles + MaxTotalParticles*nSplitMax
         allocate(Particle_IV(1:MaxTotalParticles, 1:nSplitVar))
         call init_split_grid
      else
         allocate(Particle_IV(1:MaxTotalParticles, 1:nVar))
      end if

      allocate(nParticleOnLine(1:nLine))
      nParticleOnLine = 0
     
   end subroutine init
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
   subroutine inject_particles(iLine, Time, LagrCoord)
      use PT_ModFieldline, only: get_particle_location, compute_weight
      use PT_ModDistribution, only: increase_total_weight
      use PT_ModUnit, only: kinetic_energy_to_momentum
      
      integer, intent(in) :: iLine
      real, intent(in) :: Time, LagrCoord
      integer :: i, iStart, iEnd


      iStart = nParticleOnLine(iLine) + 1
      iEnd = nParticleOnLine(iLine) + nInject

      do i = iStart, iEnd
         ! Set particle initial time, radial distance, and lagrangian coord on fieldline   
         Particle_IV(i, Time_) = Time
         Particle_IV(i, LagrCoord_) = LagrCoord
         call get_particle_location(Time, Particle_IV(i, LagrCoord_), Particle_IV(i, R_))
         ! Set particle energy/momentum
         Particle_IV(i, Energy_) = E0
         Particle_IV(i, Momentum_) = kinetic_energy_to_momentum(Particle_IV(i, Energy_))

         ! RIGHT NOW WEIGHT IS SET TO 1
         call compute_weight(Particle_IV(i, Time_), Particle_IV(i, LagrCoord_), Particle_IV(i, Weight_))
         call increase_total_weight(Particle_IV(i, Weight_))

         ! particle splitting variables
         ! parent index = index of first parent (to track and limit maximum number of splits)
         ! numchildren = number of child  particles (to track and limit maximum number of splits)
         ! splitlevel = number of splits removed from first parent (to track next split)
         ! if(UseSplit) then
         !    Particle_IV(i, ParentIndex_) = i
         !    Particle_IV(i, NumChildren_) = 0
         !    Particle_IV(i, SplitLevel_) = 1
         ! end if
         nParticleOnLine(iLine) = nParticleOnLine(iLine) + 1
      end do
      
   end subroutine inject_particles
   !============================================================================
   subroutine check_split(iParticle, DoSplit)
      integer, intent(in) :: iParticle
      logical, intent(out) :: DoSplit

      integer :: SplitLevel, ParentNumChildren

      DoSplit = .false.
      ! current energy threshold index of splitting
      SplitLevel = int(Particle_IV(iParticle, SplitLevel_))
      ! total number of children from first parent
      ParentNumChildren = int(Particle_IV(int(Particle_IV(iParticle, ParentIndex_)), &
                                         NumChildren_))
      
      ! if max split energy threshold has not yet been reached
      ! if particle crosses next energy threshold
      ! if max child particles per first parent has not yet been reached
      if(SplitLevel.lt.nSplitLev.and. &
         Particle_IV(iParticle, Energy_).gt.eSplitLev_I(SplitLevel) &
         .and.ParentNumChildren.le.nSplitMax) DoSplit = .true.

   end subroutine check_split
   !============================================================================
   subroutine split_particle(iParticle, iLine)

      ! index of particle being split, total number of current particles in simulation
      integer, intent(in) :: iParticle, iLine
      integer :: ParentIndex
      integer :: SplitInd

      ! index of newly split particle
      SplitInd = nParticleOnLine(iLine) + 1
      ! index of first parent
      ParentIndex = int(Particle_IV(iParticle, ParentIndex_))
      
      ! increase split level of current particle
      Particle_IV(iParticle, SplitLevel_) = Particle_IV(iParticle, SplitLevel_) + 1
      ! increase total number of children from first parent
      Particle_IV(ParentIndex, NumChildren_) = Particle_IV(ParentIndex, NumChildren_) + 1

      ! Copy particle information
      Particle_IV(SplitInd, :) = Particle_IV(iParticle, :)

      ! adjust weights of split particles
      ! conserves total weight
      Particle_IV(SplitInd, Weight_) = Particle_IV(SplitInd, Weight_) * 0.5
      Particle_IV(iParticle, Weight_) = Particle_IV(iParticle, Weight_) * 0.5
      
      ! Increase number of particles on fieldline
      nParticleOnLine(iLine) = nParticleOnLine(iLine) + 1

   end subroutine split_particle
   !============================================================================
   subroutine advance_particles(iLine, TimeLimit, BinTime)

      use PT_ModDistribution, ONLY: bin_particle, TimeWindow
      use PT_ModFieldline, ONLY: check_boundary_conditions

      integer, intent(in) :: iLine
      real, intent(in) :: TimeLimit, BinTime
      integer :: iParticle, nLagr, iShock, iShockOld, nProgress, iProgress
      logical :: DoSplit, IsOutside

      ! Loop over particles on this line
      PARTICLE: do iParticle = 1, nParticleOnLine(iLine)
         ! particle time loop    
         do while(Particle_IV(iParticle, Time_).lt.TimeLimit) 
            ! move particle one timestep
            call advance_particle(iParticle, TimeLimit)

            ! check if particle left spatial boundaries
            call check_boundary_conditions(Particle_IV(iParticle, Time_),      &
                                           Particle_IV(iParticle, LagrCoord_), &
                                           Particle_IV(iParticle, R_),         &
                                           IsOutside)
            
            ! bin particle in space/time/momentum at end of time step
            if(Particle_IV(iParticle, Time_).ge.(BinTime-TimeWindow)) then
               call bin_particle(Particle_IV(iParticle, LagrCoord_), &
                                 Particle_IV(iParticle, Energy_),    &
                                 Particle_IV(iParticle, Weight_))
            end if

            ! if particle leaves spatial boundary set time to large number
            if(IsOutside) then
               Particle_IV(iParticle, Time_) = huge(0.0)
            end if

            ! particle splitting:
            ! total weight is conserved
            ! if(UseSplit) then
            !    call check_split(iParticle, DoSplit)
            !    if(DoSplit) then          
            !       call split_particle(iParticle, iLine)
            !    end if
            ! end if

         end do ! particle time loop
         
      end do PARTICLE 

   end subroutine advance_particles
   !============================================================================
   subroutine advance_particle(iParticle, tStepMax)
      use PT_ModFieldline, only : get_particle_location
      use PT_ModUnit, only: momentum_to_kinetic_energy
      integer, intent(in) :: iParticle
      real, intent(in) :: tStepMax

      ! Vectors solved by SDE
      real :: XOld(nDim), XNew(nDim)

      ! Equation advances lagrcoord and p^3/3 not p
      XOld(LagrCoord_) = Particle_IV(iParticle, LagrCoord_)
      XOld(Momentum_)  = Particle_IV(iParticle, Momentum_)**3.0 / 3.0

      ! Advance psuedo-particle one time step
      call euler_sde(XOld, tStepMax, Particle_IV(iParticle, Time_), &
                     Particle_IV(iParticle, TimeStep_), XNew)
                     
      ! Update Lagrangian coordinate and momentum
      Particle_IV(iParticle, LagrCoord_) = XNew(LagrCoord_)
      Particle_IV(iParticle, Momentum_) = (3.0*XNew(Momentum_))**(1.0/3.0)
      Particle_IV(iParticle, Energy_) = momentum_to_kinetic_energy(Particle_IV(iParticle, Momentum_))

      ! Update particle time
      Particle_IV(iParticle, Time_) = Particle_IV(iParticle, Time_) + Particle_IV(iParticle, TimeStep_)
      ! Update radial distance (R)
      call get_particle_location(Particle_IV(iParticle, Time_), &
                                 Particle_IV(iParticle, LagrCoord_), &
                                 Particle_IV(iParticle, R_))


   end subroutine advance_particle
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
   end module PT_ModParticle
   !=============================================================================

program dsa_1D_spherical

   ! Solves the Parker transport equation using stochastic integration method
   ! 1D lagrangian frame along single IMF line
   ! solves acceleration and transport
   !
   ! Particle splitting is used
   ! Uses either Euler-Maruyama or Milstein method 
   ! Particles are injected uniformly in r
   ! Particles can be injected at constant energy (default) or from 1/v^5 distribution
   use ModKind
   use ModMpi
   use ModReadParam, ONLY: read_file, read_init
   use ModUtilities, ONLY: remove_file, touch_file
 
   use PT_ModTime,   ONLY: iIter, init_time => init, PTTime
  
   use PT_ModMain,   ONLY: &
         IsLastRead, IsStandAlone,          &
         TimeMax, nIterMax, CpuTimeMax,     &
         PT_read_param => read_param,       &
         PT_check      => check,            &
         PT_initialize => initialize,       &
         PT_run        => run,              &
         PT_finalize   => finalize
   
   use PT_ModGrid,   ONLY: init_stand_alone, init_grid => init
   
   use PT_ModConst, ONLY: cRsun, ckeV, cMeV
   
   use PT_ModPlot, ONLY : set_bins, bin_particle, save_output, increase_total_weight
   use PT_ModProc, ONLY: iProc, nProc, iComm, iError
   use PT_ModParticle

   ! use PT_ModFieldline, ONLY: read_fieldline, tMax
   use PT_ModTestFieldline, ONLY: read_fieldline, tMax

   use PT_ModRandom, ONLY: init_random_seed, read_seed_file, &
                           save_seed, UseInputSeedFile
      
   implicit none

   integer      :: iSession = 1
   real(Real8_) :: CpuTimeStart, CpuTime
   logical      :: IsFirstSession = .true.

   logical :: IsOutside, DoSplit
   ! particle loop variable
   integer :: iParticle, nParticle
   
   ! particle-splitting variables
   integer :: lev, iSplit, iSplitCounter
   
   ! to be replaced by SWMF coupled variables
   real :: SimTime = 0.0, tStepMax
   real :: SimTimeStep = 120.0

   character(len = 20) :: iProcStr ! for testing

   ! mpi initialization routines
   !----------------------------------------------------------------------------
   call MPI_INIT( iError )
   
   ! Assign communicator
   iComm = MPI_COMM_WORLD
   call MPI_COMM_RANK( iComm, iProc, iError )
   call MPI_COMM_SIZE( iComm, nProc, iError )

   ! for testing purposes
   write(iProcStr, *) iProc
   iProcStr = adjustl(iProcStr)

   ! Mark the run as a stand alone
   IsStandAlone = .true.

   ! Read PARAM.in file. Provide default restart file for #RESTART
   call read_file('PARAM.in', iComm)

   SESSIONLOOP: do
      call read_init('  ', iSessionIn=iSession)
   
      if(iProc==0)&
         write(*,*)'----- Starting Session ',iSession,' ------'
      ! Set and check input parameters for this session
      call PT_read_param  ! Identical to SP_set_param('READ')
      call PT_check       ! Similar to SP_set_param('CHECK'), but see init_time

      ! Use input seed file or randomize PRNG seed
      if(UseInputSeedFile) then
         call read_seed_file()
      else
         call init_random_seed()
      end if

      call save_seed() ! saves seed to file

      ! read in field line data
      ! remove when coupled?
      ! sets tMin, tMax, rMin, rMax
      call read_fieldline()
      
      if(IsFirstSession)then
         ! Time execution (timing parameters set by SP_read_param)
         ! call init_grid        ! Similar to SP_set_param('GRID')
         call init_stand_alone ! Distinctions from SWMF (CON_bline) version
         call PT_initialize    ! Similar to SP_init_session, NO origin points
         call init_time        ! StartTime, StartTimeJulian from StartTime_I
         ! DataInputTime=0, PTTime=0 unless set in PARAM.in or restart.H
         ! In SWMF, PTTime is set in PT_set_param('CHECK'), DataInput time is
         ! set either in coupling or in reading the MHD data from file.
         ! definititions, simulation paramters (cgs units)

         call init_split_grid
   
         ! set up bin matrix
         call set_bins ! move to PT initialize?

      end if
      
      ! Initialize time which is used to check CPU time
      CpuTimeStart = MPI_WTIME()

      ! init all particles
      call initialize_particles()

      ! sum total weight 
      ! old function from previous loop structure 
      ! can be improved
      do iParticle = 1, nParticlePerProc
         call increase_total_weight(Particle_V(iParticle, Weight_))
      end do

      ! initial number of particles
      nParticle = nParticlePerProc

      ! SWMF time loop 
      do while(SimTime.le.tMax)  
         ! upper bound of current time step
         tStepMax = SimTime + SimTimeStep
         
         ! update fieldline data goes here?

         ! particle loop
         do iParticle = 1, nParticle
            ! particle time loop
            do while(Particle_V(iParticle, Time_).le.tStepMax) 
               
               call advance_particle(iParticle)

               ! bin particle in space/time/momentum
               call bin_particle(Particle_V(iParticle, LagrCoord_), &
                                 Particle_V(iParticle, Time_), &
                                 Particle_V(iParticle, Energy_), &
                                 Particle_V(iParticle, Weight_))
               ! check if particle left spatial boundaries
               call check_boundary_conditions(iParticle, IsOutside)
               ! if particle leaves spatial boundary set time to large number
               if(IsOutside) Particle_V(iParticle, Time_) = tMax*2.0

               ! particle splitting:
               ! total weight is conserved
               if(UseSplit) then
                  call check_split(iParticle, DoSplit)
                  if(DoSplit) then          
                     call split_particle(iParticle, nParticle)
                     nParticle = nParticle + 1
                  end if
               end if

            end do ! particle time loop

            ! Missing CPU time check - what loop to put this in?
         end do ! particle loop
         SimTime = tStepMax
      end do ! time loop

      if(IsLastRead)EXIT SESSIONLOOP

      if(iProc==0) &
            write(*,*)'----- End of Session   ',iSession,' ------'
            
      iSession       = iSession + 1
      IsFirstSession = .false.
   end do SESSIONLOOP

   ! finialize mpi routine
   call save_output
   call MPI_finalize(iError)
   
end program dsa_1D_spherical
!==============================================================================


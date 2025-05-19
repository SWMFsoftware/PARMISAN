program dsa_1D_spherical

   ! solves the Parker transport equation using stochastic integration method
   ! 1D lagrangian frame along single IMF line
   ! solves acceleration and transport
   !
   ! particle splitting is used
   !
   ! Can solve equivalent Ito or Stratonovich using RK2 (or Milstein for Ito) approach
   !
   ! Particles are injected uniformly in r
   !
   ! Particles can be injected at constant energy (default) or from 1/v^5 distribution
   !
   ! this is an mpi version
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

   use PT_ModFieldline, ONLY: read_fieldline
   ! use PT_ModTestFieldline, ONLY: read_fieldline

   use PT_ModRandom, ONLY: init_random_seed, read_seed_file, &
                           save_seed, UseInputSeedFile
      
   implicit none

   integer      :: iSession = 1
   real(Real8_) :: CpuTimeStart, CpuTime
   logical      :: IsFirstSession = .true.
   logical      :: UseSplit = .false.

   logical :: IsOutside

   integer :: n, lev, iSplit, iSplitCounter
   integer, allocatable :: lev_save_split(:)
   real :: weight
   real, allocatable ::  t_save_split(:), p_save_split(:), lagr_save_split(:)

   character(len = 20) :: iProcStr
   integer :: counter = 0
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

      if(UseInputSeedFile) then
         call read_seed_file()
      else
         call init_random_seed()
      end if

      call save_seed() ! saves seed to file

      ! read in field line data
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

         if(.not.allocated(lev_save_split)) allocate(lev_save_split(1:nSplitMax))
         if(.not.allocated(t_save_split)) allocate(t_save_split(1:nSplitMax))
         if(.not.allocated(lagr_save_split)) allocate(lagr_save_split(1:nSplitMax))
         if(.not.allocated(p_save_split)) allocate(p_save_split(1:nSplitMax))

         call init_split_grid
   
         ! set up bin matrix
         call set_bins ! move to PT initialize?

      end if
      
      ! Initialize time which is used to check CPU time
      CpuTimeStart = MPI_WTIME()

      ! particle loop
      do n = 1, nParticlePerProc

         iSplit = 0
         SPLITLOOP: do
            if(iSplit.eq.0) then

               call initialize_particle
               ! if(n.eq.1) call save_particle_trajectory(iProc, n)
               lev = 1
               iSplitCounter = 0

            else
               iSplitCounter = iSplitCounter + 1

               if(iSplitCounter > iSplit) EXIT SPLITLOOP

               call init_split_particle(t_save_split(iSplitCounter), &
                                        lagr_save_split(iSplitCounter), &
                                        p_save_split(iSplitCounter), lev)
               
            end if
            ! time loop
            call increase_total_weight(Particle_V(Weight_))
            ! call save_initial_state(iProc)
            TIMELOOP: do
               call advance_particle()

               ! if(n.eq.1) call save_particle_trajectory(iProc, n)

               ! Checks particle position, lagrange coordinate, and time
               call check_boundary_conditions(IsOutside)               
               if(IsOutside) then 
                  ! call debug_particle('')
                  ! call save_final_state(iProc)
                  EXIT TIMELOOP
               end if

               call bin_particle(Particle_V(LagrCoord_), Particle_V(Time_), &
                                 Particle_V(Energy_), Particle_V(Weight_))

               ! calculate cpu time
               CpuTime =  mpi_wtime() - CpuTimeStart
               if (CpuTime > CpuTimeMax)then
                  write(*,*) 'CPU Max Time Reached: iProc, N_completed: ', iProc, n-1
                  call save_output
                  call MPI_finalize(iError)
                  stop
               end if
               
               ! particle splitting
               if(UseSplit)then
                  if(Particle_V(Energy_) > eSplitLev_I(lev).and.lev<nSplitLev &
                        .and.iSplit <= nSplitMax)then
                     lev=lev+1
                     iSplit=iSplit+1
                     t_save_split(iSplit) = Particle_V(Time_)
                     p_save_split(iSplit) = Particle_V(Momentum_)
                     lagr_save_split(iSplit) = Particle_V(LagrCoord_)
                     lev_save_split(iSplit) = lev
                  end if
               end if
            end do TIMELOOP
           
            if( lev==1 ) EXIT SPLITLOOP ! no particles were split, exit
            ! need to go back and do all split particles
         end do SPLITLOOP

      end do ! PARTICLE LOOP 

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


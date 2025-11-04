!  Copyright (C) 2002 Regents of the University of Michigan
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModMain

  use PT_ModProc,       ONLY: iProc
  use PT_ModReadMhData, ONLY: DoReadMhData
  use ModUtilities,     ONLY: CON_stop

  implicit none

  SAVE
  private ! except

  ! Indicator of stand alone mode
  logical:: IsStandAlone = .false.

  ! Stopping conditions. These variables are only used in stand alone mode.
  real    :: TimeMax     = -1.0, CpuTimeMax = -1.0
  integer :: nIterMax    = -1
  logical :: UseStopFile = .true.
  logical :: IsLastRead  = .false.

  ! Logicals for actions
  !----------------------
  ! run test analytical fieldline
  logical :: DoRunTest = .false.
  ! run the component
  logical :: DoRun = .true.
  ! restart the run
  logical :: DoRestart = .false.
  ! Methods and variables from ModReadMhData
  public  :: DoReadMhData

  ! Methods and variables from this module
  public:: read_param, initialize, finalize, run, check, DoRestart,          &
       IsLastRead, UseStopFile, CpuTimeMax, TimeMax, nIterMax, IsStandAlone, &
       DoRunTest
contains
   !============================================================================
   subroutine read_param
      use PT_ModGrid,          ONLY: read_param_grid       => read_param
      use PT_ModOriginPoints,  ONLY: read_param_origin     => read_param
      use PT_ModReadMHData,    ONLY: read_param_mhdata     => read_param
      use PT_ModTime,          ONLY: read_param_time       => read_param
      use PT_ModDistribution,  ONLY: read_param_distribution => read_param
      use PT_ModRandom,        ONLY: read_param_random     => read_param
      use PT_ModParticle,      ONLY: read_param_particle   => read_param
      use PT_ModSolver,        ONLY: read_param_solver     => read_param
      use PT_ModShock,         ONLY: read_param_shock      => read_param
      use PT_ModFieldline,     ONLY: read_param_fieldline

      ! Read input parameters for PT component
      use ModReadParam, ONLY: &
            read_var, read_line, read_command, i_session_read, read_echo_set

      ! aux variables
      integer:: nParticleCheck, nLonCheck, nLatCheck
      logical:: DoEcho

      ! The name of the command
      character (len=100) :: NameCommand
      ! Read the corresponding section of input file
      character(len=*), parameter:: NameSub = 'read_param'
      !--------------------------------------------------------------------------
      do
         if(.not.read_line() ) then
            IsLastRead = .true.
            EXIT
         end if
         if(.not.read_command(NameCommand)) CYCLE
         select case(NameCommand)
         case('#RESTART')
            call read_var('DoRestart', DoRestart)
            ! read parameters for each module
         case('#ORIGIN')
            if(IsStandAlone) CYCLE
            call read_param_origin
         case('#COORDSYSTEM', '#COORDINATESYSTEM', '#TESTPOS', &
               '#CHECKGRIDSIZE', '#GRIDNODE', "#BOUNDARY")
            ! Currently we do not need '#DOSMOOTH'
            if(i_session_read() /= 1) CYCLE
            call read_param_grid(NameCommand)
         case('#READMHDATA','#MHDATA')
            call read_param_mhdata(NameCommand)
         case('#TRACESHOCK', '#IDENTIFYSHOCK')
            call read_param_shock(NameCommand)
         case('#DORUN')
            call read_var('DoRun', DoRun)
         case('#TEST')
            call check_stand_alone
            call read_var('DoRunTest', DoRunTest)
            !if(DoRunTest) call read_param_test(NameCommand)
         case('#END')
            call check_stand_alone
            IsLastRead=.true.
            EXIT
         case('#RUN')
            call check_stand_alone
            IsLastRead=.false.
            EXIT
         case('#STOP')
            call check_stand_alone
            call read_var('nIterMax', nIterMax)
            call read_var('TimeMax' , TimeMax)
         case('#CPUTIMEMAX')
            call check_stand_alone
            call read_var('CpuTimeMax', CpuTimeMax)
         case('#CHECKSTOPFILE')
            call check_stand_alone
            call read_var('UseStopFile',UseStopFile)
         case('#ECHO')
            call check_stand_alone
            call read_var('DoEcho', DoEcho)
            if(iProc==0)call read_echo_set(DoEcho)
         case("#STARTTIME",'#NSTEP','#TIMESIMULATION')
            if(i_session_read() /= 1)CYCLE
            call read_param_time(NameCommand)
         case("#SETREALTIME")
            if(i_session_read() /= 1)CYCLE
            call check_stand_alone
            call read_param_time(NameCommand)
         case("#TIMEACCURATE")
            call check_stand_alone
            call read_param_time(NameCommand)
         case("#DISTRIBUTION")
            call read_param_distribution(NameCommand)
         case("#PARTICLE")
            call read_param_particle(NameCommand)
         case("#INPUTSEED")
            call read_param_random(NameCommand)
         case("#SDE")
            call read_param_solver(NameCommand)
         case("#DIFFUSION")
            call read_param_fieldline(NameCommand)
         case default
            call CON_stop(NameSub//': Unknown command '//NameCommand)
         end select
      end do
   contains
      !==========================================================================
      subroutine check_stand_alone

         ! certain options are only available for stand alone mode;
         ! check whether the mode is used and stop the code if it's no the case
         !------------------------------------------------------------------------
         if(IsStandAlone)RETURN
         call CON_stop(NameSub//': command '//trim(NameCommand)//&
            ' is only allowed in stand alone mode, correct PARAM.in')
      end subroutine check_stand_alone
      !==========================================================================
   end subroutine read_param
   !============================================================================
   subroutine initialize
      use PT_ModReadMhData,   ONLY: init_mhdata       => init
      use PT_ModRandom,       ONLY: init_seed         => init
      use PT_ModDistribution, ONLY: init_distribution => init
      use PT_ModParticle,     ONLY: init_particle     => init
      use PT_ModShock,        ONLY: init_shock        => init
      use PT_ModPlot,         ONLY: save_bin_arrays
      
      ! initialize the model
      character(len=*), parameter:: NameSub = 'initialize'
      !--------------------------------------------------------------------------
      if(iProc.eq.0)then
         write(*,'(a)')'PT: initializing'
      end if

      ! Reads in first data file if reading in files, else nothing
      call init_mhdata
      call init_shock ! Init shock arrays

      ! Sets random seeds for all processors
      call init_seed
      ! Initializes phase space bins and counts array
      call init_distribution
      call save_bin_arrays
      ! Init particle array and split grid
      call init_particle
      
      if(iProc.eq.0)then
         write(*,'(a)')'PT: initialized'
      end if

   end subroutine initialize
   !============================================================================
   subroutine finalize
      use PT_ModReadMhData, ONLY: finalize_mhdata    => finalize
      ! finalize the model
      ! if(IsStandAlone)call stand_alone_final_restart
      character(len=*), parameter:: NameSub = 'finalize'
      !--------------------------------------------------------------------------

      call finalize_mhdata
   end subroutine finalize
   !============================================================================
   subroutine run

      use PT_ModGrid,          ONLY: get_other_state_var, copy_old_state, Used_B
      use PT_ModShock,         ONLY: get_shock_location, get_dLogRho, &
                                     DoTraceShock, set_initial_shock
      use PT_ModReadMhData,    ONLY: read_mh_data
      use PT_ModTime,          ONLY: PTTime, DataInputTime, iIter
      use PT_ModAdvance,       ONLY: advance
      use PT_ModPlot,          only: save_analytic_solution
      
      ! advance the solution in time
      logical, save :: IsFirstCall = .true.
      real :: Dt ! time increment in the current call

      ! write the initial background state to the output file
      !--------------------------------------------------------------------------
      if(IsFirstCall)then
         ! recompute the derived components of state vector, e.g.
         ! magnitude of magnetic field and velocity etc. Smooth if needed.
         if(.not.DoReadMhData) call get_other_state_var
         
         if(DoRunTest) then 
            call set_initial_shock
            call save_analytic_solution
         end if
         IsFirstCall = .false.
      end if

      ! May need to read background data from files
      if(DoReadMhData)then
         ! copy some variables from the previous time step
         call copy_old_state
         ! Read the background data from file
         call read_mh_data()
         ! Read from file: MHData_VIB(0:nMHData,::) for the time moment
         ! DataInputTime
      end if

      ! recompute the derived components of state vector, e.g.
      ! magnitude of magnetic field and velocity etc
      call get_other_state_var

      ! if no new background data loaded, don't advance in time
      if(DataInputTime <= PTTime) RETURN
      
      if(DoTraceShock) then
         call get_dLogRho
         call get_shock_location
      end if

      ! run the model from PTTime to min(DataInputTime, TimeMax)
      if(iProc.eq.0) write(*,*) 'Advancing from: ', PTTime, min(DataInputTime, TimeMax)
      call advance(min(DataInputTime, TimeMax))

      ! update time & iteration counters
      iIter = iIter + 1
      Dt = min(DataInputTime, TimeMax) - PTTime
      PTTime = PTTime + Dt
      
      ! call save_plot_all

      ! save restart if needed
      ! call check_save_restart(Dt)
   end subroutine run
   !============================================================================
   subroutine check

      use ModUtilities, ONLY: make_dir
      ! Make output directory

      character(len=*), parameter:: NameSub = 'check'
      !--------------------------------------------------------------------------
   end subroutine check
   !============================================================================
end module PT_ModMain
!==============================================================================

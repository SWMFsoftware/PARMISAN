!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
program MITTENS
   use ModKind
   use PT_ModProc,   ONLY: iProc, nProc, iComm
   use ModUtilities, ONLY: remove_file, touch_file, CON_Stop
   use PT_ModTime,   ONLY: iIter, init_time  => init
   use PT_ModMain,   ONLY: &
         IsLastRead, IsStandAlone, DoRunTest, &
         TimeMax, nIterMax,                   &
         PT_read_param => read_param,         &
         PT_check      => check,              &
         PT_initialize => initialize,         &
         PT_run        => run,                &
         PT_finalize   => finalize           
   use PT_ModGrid,   ONLY: init_stand_alone, init_grid=>init
   use ModReadParam, ONLY: read_file, read_init
   use ModMpi
   ! use PT_CreateSyntheticData, only: create_files

   implicit none

   integer      :: iError
   integer      :: iSession = 1
   real(Real8_) :: CpuTimeStart, IterStart, IterStop
   logical      :: IsFirstSession = .true.

   ! Initialization of MPI/parallel message passing.
   !----------------------------------------------------------------------------
   call MPI_INIT(iError)
   iComm=MPI_COMM_WORLD
   call MPI_COMM_RANK(iComm, iProc, iError)
   call MPI_COMM_SIZE(iComm, nProc, iError)

   ! Initialize time which is used to check CPU time
   CpuTimeStart = MPI_WTIME()

   ! Delete MITTENS.SUCCESS and MITTENS.STOP files if found
   if(iProc==0)then
      call remove_file('MITTENS.SUCCESS')
      call remove_file('MITTENS.STOP')
   end if

   ! Mark the run as a stand alone
   IsStandAlone = .true.

   ! Read PARAM.in file. Provide default restart file for #RESTART
   call read_file('PARAM.in',iComm)

   SESSIONLOOP: do
      call read_init('  ', iSessionIn=iSession)

      if(iProc==0)&
            write(*,*)'----- Starting Session ',iSession,' ------'

      ! Set and check input parameters for this session
      call PT_read_param  ! Identical to SP_set_param('READ')
      call PT_check       ! Similar to SP_set_param('CHECK'), but see init_time

      if(IsFirstSession)then
         ! Time execution (timing parameters set by SP_read_param)
         call init_grid        ! Similar to SP_set_param('GRID')
         call init_stand_alone ! Distinctions from SWMF (CON_bline) version
         call PT_initialize    ! Similar to SP_init_session, NO origin points
         call init_time        ! StartTime, StartTimeJulian from StartTime_I
      end if
      TIMELOOP: do
         if(stop_condition_true()) EXIT TIMELOOP
         if(is_time_to_stop()) EXIT SESSIONLOOP
         IterStart =  MPI_WTIME()
         call PT_run
         if(iProc.eq.0) write(*,*) "iteration took: ", MPI_WTIME() - IterStart
         ! call show_progress
      end do TIMELOOP

      if(IsLastRead)EXIT SESSIONLOOP
      if(iProc==0) &
            write(*,*)'----- End of Session   ',iSession,' ------'
      iSession       = iSession + 1
      IsFirstSession = .false.

   end do SESSIONLOOP
   if(iProc==0)then
      write(*,'(a)')'    -----------------------------'
      write(*,'(a)')'    Finished Numerical Simulation'
      write(*,'(a)')'    -----------------------------'
   end if

   ! Finish writing to log file
   call PT_finalize

   ! Touch MITTENS.SUCCESS
   if(iProc==0) call touch_file('MITTENS.SUCCESS')
   
   ! Finalize MPI
   call MPI_Finalize(iError)

contains
  !============================================================================
  function stop_condition_true() result(IsStopCondition)
    use PT_ModMain, ONLY: nIterMax
    use PT_ModTime, ONLY: PTTime
    logical :: IsStopCondition
    !--------------------------------------------------------------------------
    IsStopCondition = .false.

    if(nIterMax >= 0  .and. iIter >= nIterMax) IsStopCondition = .true.
    if(TimeMax >  0.0 .and. PTTime >= TimeMax) IsStopCondition = .true.

  end function stop_condition_true
  !============================================================================
  function is_time_to_stop() result(IsTimeToStop)
    use PT_ModMain, ONLY: CpuTimeMax, UseStopFile
    logical :: IsTimeToStop
    !--------------------------------------------------------------------------
    IsTimeToStop = .false.

    if(iProc==0)then
       if(CpuTimeMax > 0.0 .and. MPI_WTIME()-CpuTimeStart >= CpuTimeMax)then
          write(*,*)'CPU time exceeded:',CpuTimeMax,MPI_WTIME()-CpuTimeStart
          IsTimeToStop=.true.
       end if
       if(.not.IsTimeToStop .and. UseStopFile) then
          inquire(file='MITTENS.STOP',exist=IsTimeToStop)
          if (IsTimeToStop) &
               write(*,*)'MITTENS.STOP file exists: received stop signal'
       end if
    end if
    if(nProc==1) RETURN
    call MPI_BCAST(IsTimeToStop,1,MPI_LOGICAL,0,iComm,iError)

  end function is_time_to_stop
  !============================================================================
end program MITTENS
!==============================================================================


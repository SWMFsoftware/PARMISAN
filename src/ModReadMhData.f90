!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModReadMhData

  ! This module contains methods for reading input MH data

  use PT_ModGrid,    ONLY: iblock_to_lon_lat, get_other_state_var,   &
       nMHData, nLine, Z_, Used_B, FootPoint_VB, nVertex_B, MHData_VIB, &
       LagrID_, MinLagr, MaxLagr
  use PT_ModTime,    ONLY: PTTime, DataInputTime
  use ModPlotFile,   ONLY: read_plot_file
  use ModUtilities,  ONLY: fix_dir_name, open_file, close_file, CON_stop
  use ModIoUnit,     ONLY: io_unit_new

  implicit none

  SAVE

  private ! except

  ! Public members
  public:: init         ! Initialize module variables
  public:: read_param   ! Read module variables
  public:: read_mh_data ! Read MH_data from files
  public:: finalize     ! Finalize module variables DoReadMhData
  ! If the folliwing logical is true, read MH_data from files
  logical, public :: DoReadMhData
  ! the input directory
  character(len=100)         :: NameInputDir=""
  ! the name with list of file tags
  character(len=100)         :: NameTagFile=""
  ! the input file name base
  character(len=4)           :: NameFileExtension
  character(len=20)          :: TypeMhDataFile

  ! IO unit for file with list of tags
  integer :: iIOTag
  integer :: iLonFile = 0, iLatFile = 0

contains
   !============================================================================
   subroutine read_param(NameCommand)

      use ModReadParam, ONLY: read_var
      ! set parameters of input files with background data
      ! integer :: nFileRead
      character (len=*), intent(in):: NameCommand ! From PARAM.in
      character(len=*), parameter:: NameSub = 'read_param'
      !--------------------------------------------------------------------------
      select case(NameCommand)
      case('#READMHDATA')
         ! determine whether to read the MHD data
         call read_var('DoReadMhData', DoReadMhData)

         if(.not. DoReadMhData)&
            RETURN
         ! the input directory
         call read_var('NameInputDir', NameInputDir)
         call fix_dir_name(NameInputDir) ! adds "/" if not present

         call read_var('iLonFile', iLonFile)
         call read_var('iLatFile', iLatFile)
      case('#MHDATA')
         ! type of data files
         call read_var('TypeFile', TypeMhDataFile)
         ! the format of output file must be set
         select case(trim(TypeMhDataFile))
         case('tec')
            NameFileExtension='.dat'
         case('idl','ascii')
            TypeMhDataFile = 'ascii'
            NameFileExtension='.out'
         case('real4','real8')
            NameFileExtension='.out'
         case default
            call CON_stop(NameSub//': input format was not set in PARAM.in')
         end select
         ! number of input files
         ! call read_var('nFileRead', nFileRead)
         ! name of the file with the list of tags
         call read_var('NameTagFile', NameTagFile)
      end select

   end subroutine read_param
   !============================================================================
   subroutine init
      ! initialize by setting the time and interation index of input files
      integer:: iTag
      character(len=50):: StringAux
      character(len=*), parameter:: NameSub = 'init'
      !--------------------------------------------------------------------------
      if(.not.DoReadMhData) RETURN

      ! open the file with the list of tags
      iIOTag = io_unit_new()
      call open_file(iUnitIn=iIOTag, &
         file=trim(NameInputDir)//trim(NameTagFile), status='old')
      ! if nTag > 0, need to skip nTag lines
      ! if(nTag>0)then
      !   do iTag = 1, nTag-1
      !      read(iIOTag,'(a)') StringAux
      !   end do
      ! end if

      ! read the first input file
      call read_mh_data(DoOffsetIn = .false.)
      call get_other_state_var
      PTTime = DataInputTime

   end subroutine init
  !============================================================================
  subroutine finalize

    ! close currentl opend files
    !--------------------------------------------------------------------------
    if(DoReadMhData) call close_file(iUnitIn=iIOTag)
  end subroutine finalize
  !============================================================================
  subroutine read_mh_data(DoOffsetIn)

   ! use SP_ModPlot,    ONLY: NameMHData
   character(len=*), parameter :: NameMHData = "MH_data"
   logical, optional, intent(in ):: DoOffsetIn
   ! read 1D MH data, which are produced by write_mh_1d n ModWrite
   ! separate file is read for each field line, name format is
   ! (usually)MH_data_<iLon>_<iLat>_t<ddhhmmss>_n<iIter>.{out/dat}
   ! name of the input file
   character(len=100):: NameFile
   ! loop variables
   integer:: iLine
   ! indexes of corresponding node, latitude and longitude
   integer:: iLat, iLon
   ! size of the offset to apply compared to the previous state
   integer:: iOffset
   integer :: startIndex = 1
   ! local value of DoOffset
   logical:: DoOffset
   ! auxilary variables to apply positive offset for appended particles
   ! auxilary parameter index
   integer, parameter:: RShock_ = Z_ + 2
   integer, parameter:: StartTime_  = RShock_ + 1
   integer, parameter:: StartJulian_= StartTime_ + 1
   ! additional parameters of lines
   real:: Param_I(LagrID_:StartJulian_)
   ! data input time before reading new data file
   real:: DataInputTimeOld
   ! timetag
   character(len=50):: StringTag
   ! status of reading .lst file
   integer :: ioStat, SleepCounter = 0
   integer :: TimeToWait, TimeToQuit

   ! check whether need to apply offset, default is .true.

   character(len=*), parameter:: NameSub = 'read_mh_data'
   !--------------------------------------------------------------------------
   if(present(DoOffsetIn))then
      DoOffset = DoOffsetIn
   else
      DoOffset = .true.
   end if

   ! get the tag for files
   ! read(iIOTag,'(a)') StringTag

   ! For Feb 2026 Artemis real-time demonstration
   ! Continually check that .lst file is updated and wait until it is
   ! Terminate if .lst is not updated after 2 minutes
   TimeToWait = 3   ! seconds
   TimeToQuit = 120 ! seconds
   do 
      read(iIOtag, '(a)', iostat = ioStat) StringTag
      ! sometimes the file returns empty string and the next read is successful
      ! if the string is empty, the fieldline will be marked as unused
      if(StringTag.eq.'') read(iIOtag, '(a)', iostat = ioStat) StringTag
      if(ioStat.lt.0) then
         backspace(iIOTag)
         SleepCounter = SleepCounter + TimeToWait
         call sleep(TimeToWait)
         if(SleepCounter.gt.TimeToQuit) &
            call CON_Stop(NameSub//': .lst file not updated for 120 seconds.')
         cycle
      else if(ioStat.gt.0) then
         call CON_Stop(NameSub//': error reading .lst file.')
      else
         print *, 'Successful: ', StringTag
         exit
      end if

   end do
   ! ---------------------------------------------

   ! save the current data input time
   DataInputTimeOld = DataInputTime
   ! read the data
   line:do iLine = 1, nLine
      
      if(.not.Used_B(iLine))then
         nVertex_B(iLine) = 0
         MinLagr(iLine) = 1
         MaxLagr(iLine) = 1
         CYCLE line
      end if
      
      call iblock_to_lon_lat(iLine, iLon, iLat)

      ! set the file name
      write(NameFile,'(a,i3.3,a,i3.3,a)') &
         trim(NameInputDir)//NameMHData//'_',iLonFile,&
         '_',iLatFile, '_'//trim(StringTag)//NameFileExtension

      inquire(file=NameFile,exist=Used_B(iLine))

      if(.not.Used_B(iLine))then
         write(*,'(a)')NameSub//': the file '//NameFile//' is not found!'
         write(*,'(a)')NameSub//': the line marked as unused'
         nVertex_B(iLine) = 0
         MinLagr(iLine) = 1
         MaxLagr(iLine) = 1
         CYCLE line
      end if
      ! read the header first

      call read_plot_file(NameFile          ,&
         TypeFileIn = TypeMhDataFile      ,&
         TimeOut    = DataInputTime       ,&
         n1out      = nVertex_B(iLine)    ,&
         ParamOut_I = Param_I(LagrID_:StartJulian_))
      ! find offset in data between new and old states
      if(DoOffset)then
         ! check consistency: time counter MUST advance
         if(DataInputTimeOld >= DataInputTime)then
            call CON_stop(NameSub//&
               ': time counter didnt advance when reading mh data file '//&
               'with tag '//trim(StringTag)//&
               '; the tag may be repeated in '//trim(NameTagFile))
         end if
         ! amount of the offset is determined from difference
         ! in LagrID_
         iOffset = nint(FootPoint_VB(LagrID_,iLine) - Param_I(LagrID_))
      else
         iOffset = 0
      end if
      ! Parameters
      FootPoint_VB(LagrID_:Z_,iLine) = Param_I(LagrID_:Z_)
      ! read MH data
      call read_plot_file(NameFile           ,&
         TypeFileIn = TypeMhDataFile       ,&
         Coord1Out_I= MHData_VIB(LagrID_   ,&
         1:nVertex_B(iLine),iLine)         ,&
         VarOut_VI  = MHData_VIB(1:nMHData ,&
         1:nVertex_B(iLine),iLine))
         
      ! Shift data such that index == lagr coordinate
      
      ! sometimes negative lagr coords appear
      ! ignoring those for now - consult igor
      do while(MhData_VIB(LagrID_, startIndex, iLine).lt.1)
         startIndex = startIndex + 1
      end do
      
      MinLagr(iLine) = MhData_VIB(LagrID_, startIndex, iLine)
      MaxLagr(iLine) = MhData_VIB(LagrID_, nVertex_B(iLine), iLine)
      MHData_VIB(:, MinLagr(iLine):MaxLagr(iLine), iLine) = &
         MHData_VIB(:, startIndex:nVertex_B(iLine), iLine)

      ! zero out lagr coord with no data
      MHData_VIB(:, 1:MinLagr(iLine)-1, iLine) = 0.0

    end do line

  end subroutine read_mh_data
  !============================================================================
end module PT_ModReadMhData
!==============================================================================

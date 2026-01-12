!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModShock

  ! This module contains subroutines for determining the shock location,
  ! and steepening the density and magnetic field strength at shock front.
  use PT_ModGrid,   ONLY: nLine, nLineAll, Used_B, MinLagr, MaxLagr, &
       LagrID_, X_, Z_, State_VIB, MhData_VIB, &
       iShock_IB, NoShock_, Shock_, ShockOld_, check_line_ishock, &
       ShockUp_, ShockDown_, ShockUpOld_, ShockDownOld_

  use PT_ModSize,   ONLY: nVertexMax
  use ModUtilities, ONLY: CON_stop

  implicit none

  SAVE

  PRIVATE ! except

  ! Public members:
  public:: read_param           ! Read parameters
  public:: init                 ! Initialize arrays on the grid
  public:: get_dLogRho          ! calculate temporal change in density 
  public:: get_shock_location   ! finds shock location on all lines
!   public:: steepen_shock        ! steepen the density profile at the shock
  public:: set_initial_shock

  ! If the shock wave is traced, the advance algorithms are modified
  logical, public :: DoTraceShock = .true.
  ! divergence of velocity \vec{U}: for determining the shock locations
  real, public, allocatable :: dLogRho_II(:, :), dLogRhoOld_II(:, :)

  ! Shock algorithm parameters:
  real,    public :: dLogRhoThreshold = 0.0001 ! Empirical value
  integer, public :: nShockWidth = 10, nSearchMax = 50

  ! Parameters for the shock coordinates
  integer, public, parameter :: nShockVar = 18, &
       ShockID_    = 1, & ! Shock index
       XShock_     = 2, & ! Shock X coordinates
       YShock_     = 3, & ! Shock Y coordinates
       ZShock_     = 4, & ! Shock Z coordinates
       RShock_     = 5, & ! Shock radial distance
       RhoShock_   = 6, & ! Shock number density
       TShock_     = 7, & ! Shock temperature
       UxShock_    = 8, & ! Shock Ux
       UyShock_    = 9, & ! Shock Uy
       UzShock_    =10, & ! Shock Uz
       BxShock_    =11, & ! Shock Bx
       ByShock_    =12, & ! Shock By
       BzShock_    =13, & ! Shock Bz
       Wave1Shock_ =14, & ! Shock Wave1
       Wave2Shock_ =15, & ! Shock Wave2
       LonShock_   =16, & ! Shock longitude
       LatShock_   =17, & ! Shock latitude
       CompRatio_  =18    ! Compression ratio
  real, public, allocatable :: StateShock_VIB(:,:)
  logical, public :: DoSaveStateShock = .false.

  ! Shock variable names
  character(len=10), public, parameter:: NameVarShock_V(ShockID_:CompRatio_) &
       = ['ShockID   ', &
       'XShock    ', &
       'YShock    ', &
       'ZShock    ', &
       'RShock    ', &
       'nShock    ', &
       'TShock    ', &
       'UxShock   ', &
       'UyShock   ', &
       'UzShock   ', &
       'BxShock   ', &
       'ByShock   ', &
       'BzShock   ', &
       'Wave1Shock', &
       'Wave2Shock', &
       'LonShock  ', &
       'LatShock  ', &
       'CompRatio ']

  ! Unit for all the shock variables: Length is in the unit of Rsun
  character(len=6), public :: NameVarShockUnit_V(ShockID_:CompRatio_) = [&
       'none  ', &
       'RSun  ', &
       'RSun  ', &
       'RSun  ', &
       'RSun  ', &
       'amu/m3', &
       'kev   ', &
       'm/s   ', &
       'm/s   ', &
       'm/s   ', &
       'T     ', &
       'T     ', &
       'T     ', &
       'J/m3  ', &
       'J/m3  ', &
       'Deg   ', &
       'Deg   ', &
       'none  ']

  ! Shock skeleton for visualization
  real, public, allocatable :: XyzShockEffUnit_DG(:, :)
  ! Arrays to construct a triangular mesh on a sphere
  logical, public :: IsShockTriMade = .false., DoOutputShock = .false.
  integer, public :: nShockTriMesh, lidShockTri, ridShockTri
  integer, public, allocatable :: &
       iListShock_I(:), iPointerShock_I(:), iEndShock_I(:)
contains
     !============================================================================
     subroutine read_param(NameCommand)

          use ModReadParam, ONLY: read_var
          character(len=*), intent(in):: NameCommand ! From PARAM.in
          character(len=*), parameter:: NameSub = 'read_param'
          !--------------------------------------------------------------------------
          select case(NameCommand)
          case('#TRACESHOCK')
               call read_var('DoTraceShock', DoTraceShock)
          case('#IDENTIFYSHOCK')
               call read_var('nShockWidth', nShockWidth)
               call read_var('nSearchMax', nSearchMax)
               call read_var('dLogRhoThreshold', dLogRhoThreshold)
               call read_var('DoOutputShock', DoOutputShock)
          case default
               call CON_stop(NameSub//': Unknown command '//NameCommand)
          end select

     end subroutine read_param
     !============================================================================
     subroutine init

          ! initialize arrays related to the shock
          use ModUtilities, ONLY: check_allocate
          use PT_ModProc,   ONLY: iError
          character(len=*), parameter:: NameSub = 'init'
          integer :: iLine
          !--------------------------------------------------------------------------
          do iLine = 1, nLine
               iShock_IB(:, iLine) = MinLagr(iLine)
          end do

          if(allocated(dLogRho_II)) deallocate(dLogRho_II)
          allocate(dLogRho_II(1:nVertexMax, 1:nLine)) ! divU
          call check_allocate(iError, 'dLogRho_II')
          dLogRho_II = 0.0

          if(allocated(dLogRhoOld_II)) deallocate(dLogRhoOld_II)
          allocate(dLogRhoOld_II(1:nVertexMax, 1:nLine)) ! divU
          call check_allocate(iError, 'dLogRhoOld_II')
          dLogRhoOld_II = 0.0

          ! allocate and assign values only when we save states for the shock
          if(DoSaveStateShock) then
               if(allocated(StateShock_VIB)) deallocate(StateShock_VIB)
               allocate(StateShock_VIB(XShock_:CompRatio_, 1:nLine)) ! States
               call check_allocate(iError, 'StateShock_VIB')
               StateShock_VIB = -1.0

               if(allocated(XyzShockEffUnit_DG)) deallocate(XyzShockEffUnit_DG)
               allocate(XyzShockEffUnit_DG(XShock_:LatShock_, 1:nLineAll+2)) ! Effective
               call check_allocate(iError, 'XyzShockEffUnit_DG')
               XyzShockEffUnit_DG = -1.0
          end if

     end subroutine init
     !============================================================================
     subroutine get_dLogRho

          ! divU = dLogRho_I for time-accurate run
          use PT_ModGrid, ONLY: Rho_, RhoOld_, MinLagrOld
          use PT_ModTime, ONLY: PTTime, DataInputTime

          ! Loop variables
          integer :: iLine, iMin, iMax
          !--------------------------------------------------------------------------

          do iLine = 1, nLine
          ! go line by line and get divU if active
          if(.not.Used_B(iLine)) then
               iShock_IB(Shock_,iLine) = NoShock_
               CYCLE
          end if
          iMin = MinLagr(iLine)
          iMax = MaxLagr(iLine)
          
          dLogRhoOld_II(:, iLine) = 0.0
          dLogRhoOld_II(iMin:iMax, iLine) = dLogRho_II(iMin:iMax, iLine)

          dLogRho_II(:, iLine) = 0.0
          dLogRho_II(iMin:iMax, iLine) = &
               log(MhData_VIB(Rho_, iMin:iMax,iLine) / &
                   State_VIB(RhoOld_, iMin:iMax,iLine)) / &
               (DataInputTime - PTTime)
          ! if new particles were introduced - avoid NaNs
          if(MinLagr(iLine).lt.MinLagrOld(iLine)) &
               dLogRho_II(iMin:MinLagrOld(iLine), iLine) = 0.0
          end do

     end subroutine get_dLogRho
     !============================================================================
     subroutine get_shock_location

          ! find location of a shock wave on a given line (line)
          ! shock front is assumed to be location of max log(Rho/RhoOld)
          use PT_ModConst,       ONLY: cRadToDeg
          use ModCoordTransform, ONLY: xyz_to_rlonlat
          use PT_ModGrid,        ONLY: R_, Rho_, Wave2_, iLineAll0, ROld_
          ! Do not search too close to the Sun
          real, parameter :: RShockMin = 1.2  ! *RSun
          integer         :: iShockMin
          ! Do not search too close to the heliosphere boundary
          integer :: iShockMax
          ! Misc
          integer :: iShockCandidate
          ! Loop variables
          integer :: iLine, iEnd, iShockForward, i

          character(len=*), parameter:: NameSub = 'get_shock_location'
          !--------------------------------------------------------------------------

          do iLine = 1, nLine
               ! go line by line and get the shock location if active
               if(.not.Used_B(iLine)) CYCLE
               
               ! Number of the active particles on the line
               iEnd = MaxLagr(iLine)
               
               iShockMin = iShock_IB(ShockOld_, iLine)
               
               ! iShockMax = iShockMax = iEnd - nShockMargin - 1
               iShockMax = iEnd - 1
 
               ! get the forward grid index for iShockCandidate
               if (any(State_VIB(R_,iShockMin:iShockMax,iLine) > RShockMin .and. &
                    dLogRho_II(iShockMin:iShockMax, iLine) > dLogRhoThreshold)) then
                    iShockForward = maxloc( &
                         dLogRho_II(iShockMin:iShockMax, iLine), DIM=1, MASK= &
                         State_VIB(R_,iShockMin:iShockMax,iLine) > RShockMin .and. &
                         dLogRho_II(iShockMin:iShockMax, iLine) > dLogRhoThreshold, BACK = .true.)
               else
                    iShockForward = 0
               end if
               iShockCandidate = iShockMin - 1 + iShockForward

               ! if shock has moved forward
               if(iShockCandidate >= iShockMin) then
                    iShock_IB(Shock_, iLine) = iShockCandidate
                    
                    ! Determine upstream extent of shock
                    i = 0
                    do while(dLogRho_II(iShockCandidate + i, iLine).gt.dLogRhoThreshold.and. &
                             (iShockCandidate+i).lt.iEnd.and.i.le.nSearchMax)
                         i = i + 1
                    end do
                    iShock_IB(ShockUp_, iLine) = iShockCandidate + i
                    
                    ! Determine downstream extent of shock
                    i = 0
                    do while(dLogRho_II(iShockCandidate - i, iLine).gt.dLogRhoThreshold.and. &
                            (iShockCandidate-i).gt.1.and.i.le.nSearchMax) 
                         i = i + 1
                    end do
                    iShock_IB(ShockDown_, iLine) = iShockCandidate - i   
               else
                    ! no change to shock location or width
                    iShock_IB(Shock_, iLine) = iShockMin
                    iShock_IB(ShockUp_, iLine) = iShock_IB(ShockUpOld_, iLine)
                    iShock_IB(ShockDown_, iLine) = iShock_IB(ShockDownOld_, iLine)
               end if

               ! check_line_ishock: update Used_B(iLine)
               call check_line_ishock(iLine)
               if(.not.Used_B(iLine)) CYCLE
               ! calculate values only when we save states for the shock
               if(DoSaveStateShock) then
                    ! get the coordinates
                    StateShock_VIB(XShock_:ZShock_, iLine) = &
                         MHData_VIB(X_:Z_, iShockCandidate, iLine)
                    call xyz_to_rlonlat(StateShock_VIB(XShock_:ZShock_, iLine), &
                         StateShock_VIB(RShock_, iLine), &
                         StateShock_VIB(LonShock_, iLine), &
                         StateShock_VIB(LatShock_, iLine))
                    if(StateShock_VIB(RShock_, iLine) == 0.0) then
                         write(*,*) "On the field line, iLineAll=", iLineAll0+iLine
                         call CON_Stop(NameSub//": Error of the shock location (R=0.0).")
                    end if

                    ! convert units for angles
                    StateShock_VIB([LonShock_,LatShock_], iLine) = &
                         StateShock_VIB([LonShock_,LatShock_], iLine) * cRadToDeg
                    ! get MHD VARs
                    StateShock_VIB(RhoShock_:Wave2Shock_, iLine) = &
                         MHData_VIB(Rho_:Wave2_, iShockCandidate, iLine)
                    ! also get compression ratio at shock surface
                    StateShock_VIB(CompRatio_, iLine) = &
                         maxval(MHData_VIB(Rho_, &
                         iShockCandidate-nShockWidth+1:iShockCandidate+1, iLine), &
                         MASK=dLogRho_II(iShockCandidate-nShockWidth+1: &
                         iShockCandidate+1, iLine) > dLogRhoThreshold)/ & ! post shock
                         minval(MHData_VIB(Rho_, iShockCandidate+1: &
                         iShockCandidate+nShockWidth, iLine), &
                         MASK=dLogRho_II(iShockCandidate+1:iShockCandidate+nShockWidth, &
                         iLine) > dLogRhoThreshold .and. &
                         MHData_VIB(Rho_, iShockCandidate+1: &
                         iShockCandidate+nShockWidth, iLine) > 0.0)        ! pre shock
               end if
          end do

     end subroutine get_shock_location
     !============================================================================
     ! subroutine steepen_shock(iLine, nX, iShock, BSi_I, dLogRhoIn_I)

     !      ! change the density profile near the shock front
     !      ! so it becomes steeper for the current line
     !      use PT_ModConst,   ONLY: cTiny, cRsun
     !      use PT_ModGrid, ONLY: D_

     !      ! INPUTs
     !      integer, intent(in) :: iLine, iShock ! Indices of line and shock front
     !      integer, intent(in) :: nX            ! Number of meshes along s_L axis
     !      real, intent(inout) :: BSi_I(nX)     ! Magnetic field strength
     !      real, optional, intent(inout) :: dLogRhoIn_I(nX) ! for time-accurate run
     !      ! Local VARs
     !      real :: DsSi_I(1:nX-1), dLogRho_I(nX)
     !      real :: dLogRhoExcess_I(iShock-nShockWidth:iShock+nShockWidth-1)
     !      real :: dLogRhoExcessSum
     !      !--------------------------------------------------------------------------
     !      DsSi_I(1:nX-1) = State_VIB(D_, 1:nX-1, iLine)*cRsun
          
     !      ! get dLogRho_I if given; otherwise = -divU
     !      if(present(dLogRhoIn_I)) then
     !           dLogRho_I = dLogRhoIn_I
     !      else
     !           ! steady state: we do not have dlogrho/dt since dt=0 so we keep divU
     !           dLogRho_I = -divU_II(1:nX, iLine)
     !      end if

     !      ! find the excess of dLogRho within the shock compared
     !      ! to background averaged over length
     !      dLogRhoExcess_I = max(0.5*( &
     !           dLogRho_I(iShock-nShockWidth:iShock+nShockWidth-1) + &
     !           dLogRho_I(iShock-nShockWidth+1:iShock+nShockWidth)) - &
     !           dLogRhoThreshold, 0.0)

     !      ! a jump (dLogRhoExcess>0) in velocity accross the shock wave * \Delta t
     !      dLogRhoExcessSum = sum(dLogRhoExcess_I* &
     !           DsSi_I(iShock-nShockWidth:iShock+nShockWidth-1))

     !      ! check for zero excess
     !      !if(abs(dLogRhoExcessSum) <= cTiny) RETURN
          
     !      ! nullify excess within the smoothed shock
     !      dLogRho_I(iShock-nShockWidth:iShock+nShockWidth) = min( &
     !           dLogRhoThreshold, dLogRho_I(iShock-nShockWidth:iShock+nShockWidth))
          
     !      ! ... and concentrate it at the shock front, applying the whole jump
     !      ! in the velocity at a single grid point
     !      dLogRho_I(iShock) = dLogRhoThreshold + &
     !           dLogRhoExcessSum/DsSi_I(iShock)
     !      ! dLogRho_I(iShock-nShockWidth:iShock+nShockWidth) = dLogRho_I(iShock-nShockWidth:iShock+nShockWidth)
     !      ! dLogRho_I(iShock-nShockWidth:iShock+nShockWidth) = dLogRhoThreshold + &
     !      !      dLogRhoExcessSum/DsSi_I(iShock-nShockWidth:iShock+nShockWidth)
          
     !      ! also, sharpen the magnetic field magnitude
     !      ! post shock part
     !      BSi_I(iShock-nShockWidth+1:iShock+1) = &
     !           maxval(BSi_I(iShock+1-nShockWidth:iShock+1))
     !      ! pre shock part
     !      BSi_I(iShock+1:iShock+nShockWidth  ) = &
     !           minval(BSi_I(iShock+1:iShock+nShockWidth))

     !      ! update dLogRhoIn (if given) and divU
     !      if(present(dLogRhoIn_I)) then
     !           dLogRhoIn_I = dLogRho_I
     !           divU_II(1:nX, iLine) = -dLogRhoIn_I
     !      else
     !           ! steady state: we do not have dlogrho/dt since dt=0 so we keep dLogRho
     !           divU_II(1:nX, iLine) = -dLogRho_I
     !      end if

     ! end subroutine steepen_shock
     !============================================================================
     subroutine set_initial_shock
          use PT_ModGrid, only: U_, D_

          integer :: iLine, iEnd, i, iShock
          real :: divU(nVertexMax)


          do iLine = 1, nLine
               divU = 0.0

               ! D_ and U_ are not defined at MaxLagr
               iEnd = MaxLagr(iLine) - 1
               divU(MinLagr(iLine):iEnd-1) = (State_VIB(U_, MinLagr(iLine)+1:iEnd, iLine) - &
                                              State_VIB(U_, MinLagr(iLine):iEnd-1, iLine)) / &
                                              State_VIB(D_, MinLagr(iLine):iEnd-1, iLine)

               iShock = minloc(divU(MinLagr(iLine):iEnd-1), DIM = 1) + MinLagr(iLine) - 1
               iShock_IB(Shock_, iLine) = iShock
              
               i = 0
               do while(divU(iShock + i).lt.-dLogRhoThreshold.and.(iShock+i).lt.iEnd.and.i.le.nSearchMax)
                    i = i + 1
               end do
               iShock_IB(ShockUp_, iLine) = iShock + 1
               
               i = 0
               do while(divU(iShock - i).lt.-dLogRhoThreshold.and.(iShock-i).gt.MinLagr(iLine).and.i.le.nSearchMax) 
                    i = i + 1
               end do
               iShock_IB(ShockDown_, iLine) = iShock - i   
               write(*,*) 'INITIAL: ', iShock, iShock_IB(ShockUp_, iLine), iShock_IB(ShockDown_, iLine)

          end do

     end subroutine set_initial_shock
     !============================================================================
end module PT_ModShock

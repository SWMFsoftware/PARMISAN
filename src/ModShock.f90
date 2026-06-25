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
  public:: sharpen_shock        ! steepen the density profile at the shock
  public:: get_test_shock
  public:: scale_shock

  ! If the shock wave is traced, the advance algorithms are modified
  logical, public :: DoScaleShock = .false., DoSharpenShock = .false.
  ! divergence of velocity \vec{U}: for determining the shock locations
  real, public, allocatable :: dLogRho_II(:, :), dLogRhoOld_II(:, :)

  ! Shock algorithm parameters:
  real,    public :: dLogRhoThreshold = 0.0001 ! Empirical value
  integer, public :: nSearchMax = 50

  ! Parameters for the shock coordinates
  integer, public, parameter :: nShockVar = 2, &
     RShock_     = 1, &
     ShockSpeed_ = 2

  real, public, allocatable :: StateShock_II(:,:)
  logical, public :: DoSaveStateShock = .false.

  ! Shock variable names
!   character(len=10), public, parameter:: NameVarShock_V(ShockID_:CompRatio_) &
!        = ['ShockID   ', &
!        'XShock    ', &
!        'YShock    ', &
!        'ZShock    ', &
!        'RShock    ', &
!        'nShock    ', &
!        'TShock    ', &
!        'UxShock   ', &
!        'UyShock   ', &
!        'UzShock   ', &
!        'BxShock   ', &
!        'ByShock   ', &
!        'BzShock   ', &
!        'Wave1Shock', &
!        'Wave2Shock', &
!        'LonShock  ', &
!        'LatShock  ', &
!        'CompRatio ']

  ! Unit for all the shock variables: Length is in the unit of Rsun
!   character(len=6), public :: NameVarShockUnit_V(ShockID_:CompRatio_) = [&
!        'none  ', &
!        'RSun  ', &
!        'RSun  ', &
!        'RSun  ', &
!        'RSun  ', &
!        'amu/m3', &
!        'kev   ', &
!        'm/s   ', &
!        'm/s   ', &
!        'm/s   ', &
!        'T     ', &
!        'T     ', &
!        'T     ', &
!        'J/m3  ', &
!        'J/m3  ', &
!        'Deg   ', &
!        'Deg   ', &
!        'none  ']

contains
     !============================================================================
     subroutine read_param(NameCommand)

          use ModReadParam, ONLY: read_var
          character(len=*), intent(in):: NameCommand ! From PARAM.in
          character(len=*), parameter:: NameSub = 'read_param'
          !--------------------------------------------------------------------------
          select case(NameCommand)
          case('#SHOCK')
               call read_var('DoScaleShock', DoScaleShock)
               call read_var('DoSharpenShock', DoSharpenShock)
               call read_var('nSearchMax', nSearchMax)
               call read_var('dLogRhoThreshold', dLogRhoThreshold)

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


          if(allocated(StateShock_II)) deallocate(StateShock_II)
          allocate(StateShock_II(1:nShockVar, 1:nLine))
          call check_allocate(iError, 'StateShock_II')
          StateShock_II = 0.0


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
     subroutine scale_shock
          use PT_ModGrid, only: get_density_gradient, Rho_, R_, D_, S_, SOld_, ROld_
          use PT_ModProc, only: iProc
          use PT_ModTime, only: DataInputTime, PTTime

          integer :: iLine, iMin, iMax
          integer :: iShockDown, iShock, iShockUp
          real    :: GradientInRho, CompressionRatio
          real :: ShockSpeed

          do iLine = 1, nLine
               iMin = MinLagr(iLine)
               iMax = MaxLagr(iLine)
               if(.not.Used_B(iLine)) CYCLE
               ! If shock does not move - do not scale
               if(iShock_IB(Shock_, iLine).eq.iShock_IB(ShockOld_, iLine)) CYCLE

               iShock = iShock_IB(Shock_, iLine)
               iShockDown = iShock_IB(ShockDown_, iLine)
               iShockUp = iShock_IB(ShockUp_, iLine)

               ! Calculate compression ratio. Density is normalized by local gradient calculated at 
               ! the start of the simulation
               call get_density_gradient(iLine, State_VIB(R_, iShock, iLine), GradientInRho)
               CompressionRatio = MhData_VIB(Rho_, iShockDown, iLine) / &
                                  MhData_VIB(Rho_, iShockUp, iLine)
               CompressionRatio = CompressionRatio * (State_VIB(R_, iShockDown, iLine) / &
                                                      State_VIB(R_, iShockUp, iLine))**GradientInRho
             
               ! limit compression ratio
               if(CompressionRatio.lt.1.0.or.CompressionRatio.gt.4.0.and.iProc.eq.0) &
                    write(*,*) 'Limiting CompressionRatio: ', CompressionRatio
               CompressionRatio = min(4.0, max(1.0, CompressionRatio))

               ! dLogRho_II(iShock_IB(ShockDown_, iLine):iShock_IB(ShockUp_, iLine), iLine) = &
               !      dLogRho_II(iShock_IB(ShockDown_, iLine):iShock_IB(ShockUp_, iLine), iLine) * &
               !      (DataInputTime - PTTime) / (iShock - iShock_IB(ShockOld_, iLine))

               dLogRho_II(iShock_IB(ShockDown_, iLine):iShock_IB(ShockUp_, iLine), iLine) = &
                    dLogRho_II(iShock_IB(ShockDown_, iLine):iShock_IB(ShockUp_, iLine), iLine) * &
                    log(CompressionRatio) / sum(dLogRho_II(iShock_IB(ShockDown_, iLine):iShock_IB(ShockUp_, iLine), iLine))

               
          end do
     end subroutine scale_shock
     !============================================================================        
     subroutine sharpen_shock
          use PT_ModGrid, only: get_density_gradient, Rho_, R_, B_, dB_, D_
          use PT_ModProc, only: iProc

          integer :: iLine
          real    :: GradientInRho, dLogRhoSum

          do iLine = 1, nLine

               if(.not.Used_B(iLine)) CYCLE
               ! If shock does not move - do not sharpen
               if(iShock_IB(Shock_, iLine).eq.iShock_IB(ShockOld_, iLine)) CYCLE

               dLogRhoSum = sum(dLogRho_II(iShock_IB(ShockDown_, iLine):iShock_IB(ShockUp_, iLine), iLine))

               dLogRho_II(iShock_IB(ShockDown_, iLine):iShock_IB(ShockUp_, iLine), iLine) = dLogRhoThreshold
               dLogRho_II(iShock_IB(Shock_, iLine), iLine) = dLogRhoSum * 2.0 / 3.0
               dLogRho_II(iShock_IB(Shock_, iLine) + 1, iLine) = dLogRhoSum / 3.0

               ! Create step-functions of dS, B, dB at shock location
               ! State_VIB(D_, iShock_IB(ShockDown_, iLine):iShock_IB(Shock_, iLine), iLine) = &
               !      State_VIB(D_, iShock_IB(ShockDown_, iLine), iLine)
               ! State_VIB(D_, iShock_IB(Shock_, iLine)+2:iShock_IB(ShockUp_, iLine), iLine) = &
               !      State_VIB(D_, iShock_IB(ShockUp_, iLine), iLine)
               ! State_VIB(D_, iShock_IB(Shock_, iLine)+1, iLine) = 0.5 * &
               !      (State_VIB(D_, iShock_IB(ShockDown_, iLine), iLine) + State_VIB(D_, iShock_IB(ShockUp_, iLine), iLine))

               ! State_VIB(B_, iShock_IB(ShockDown_, iLine):iShock_IB(Shock_, iLine), iLine) = &
               !      State_VIB(B_, iShock_IB(ShockDown_, iLine), iLine)
               ! State_VIB(B_, iShock_IB(Shock_, iLine)+2:iShock_IB(ShockUp_, iLine), iLine) = &
               !      State_VIB(B_, iShock_IB(ShockUp_, iLine), iLine)
               ! State_VIB(B_, iShock_IB(Shock_, iLine)+1, iLine) = 0.5 * &
               !      (State_VIB(B_, iShock_IB(ShockDown_, iLine), iLine) + State_VIB(B_, iShock_IB(ShockUp_, iLine), iLine))

               ! State_VIB(dB_, iShock_IB(ShockDown_, iLine):iShock_IB(Shock_, iLine), iLine) = &
               !      State_VIB(dB_, iShock_IB(ShockDown_, iLine), iLine)
               ! State_VIB(dB_, iShock_IB(Shock_, iLine)+2:iShock_IB(ShockUp_, iLine), iLine) = &
               !      State_VIB(dB_, iShock_IB(ShockUp_, iLine), iLine)
               ! State_VIB(dB_, iShock_IB(Shock_, iLine)+1, iLine) = 0.5 * &
               !      (State_VIB(dB_, iShock_IB(ShockDown_, iLine), iLine) + State_VIB(dB_, iShock_IB(ShockUp_, iLine), iLine))
               
               ! Set new width of shock
               iShock_IB(ShockDown_, iLine) = iShock_IB(Shock_, iLine) - 2
               iShock_IB(ShockUp_, iLine) = iShock_IB(Shock_, iLine) + 3
          end do

     end subroutine sharpen_shock
     !============================================================================
     subroutine get_shock_location

          ! find location of a shock wave on a given line (line)
          ! shock front is assumed to be location of max log(Rho/RhoOld)
          use PT_ModConst,       ONLY: cRadToDeg, cRsun
          use ModCoordTransform, ONLY: xyz_to_rlonlat
          use PT_ModGrid,        ONLY: R_, Rho_, Wave2_, iLineAll0, ROld_
          use PT_ModTime,        ONLY: PTTime, DataInputTime

          ! Do not search too close to the Sun
          real, parameter :: RShockMin = 1.2  ! *RSun
          integer         :: iShockMin
          ! Do not search too close to the heliosphere boundary
          integer :: iShockMax
          ! Misc
          integer :: iShockCandidate
          ! Loop variables
          integer :: iLine, iEnd, iShockForward, i
          real :: Rshock, RshockOld

          character(len=*), parameter:: NameSub = 'get_shock_location'
          !--------------------------------------------------------------------------

          do iLine = 1, nLine
               ! go line by line and get the shock location if active
               if(.not.Used_B(iLine)) CYCLE
               
               ! Number of the active particles on the line
               iEnd = MaxLagr(iLine)

               iShockMin = iShock_IB(ShockOld_, iLine)

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
                    iShock_IB(ShockUp_, iLine) = iShockCandidate + i + 1
                    
                    ! Determine downstream extent of shock
                    i = 0
                    do while(dLogRho_II(iShockCandidate - i, iLine).gt.dLogRhoThreshold.and. &
                            (iShockCandidate-i).gt.MinLagr(iLine).and.i.le.nSearchMax) 
                         i = i + 1
                    end do
                    iShock_IB(ShockDown_, iLine) = iShockCandidate - i - 1
               else
                    ! no change to shock location - shock does not exist
                    ! Set shock location to old shock location
                    ! Set width equal to 0
                    iShock_IB(Shock_, iLine)     = iShock_IB(ShockOld_, iLine)
                    iShock_IB(ShockUp_, iLine)   = iShock_IB(ShockOld_, iLine)
                    iShock_IB(ShockDown_, iLine) = iShock_IB(ShockOld_, iLine)
               end if
               
               ! check_line_ishock: update Used_B(iLine)
               call check_line_ishock(iLine)
               if(.not.Used_B(iLine)) CYCLE
          end do

     end subroutine get_shock_location
     !============================================================================
     subroutine get_test_shock
          use PT_ModGrid, only: U_, D_, DensityGradient

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

               ! slighly dubious here, calculating divU in this manner for the first timestep isn't perfect
               ! location of shock from divU is actually one lagr coordinate ahead - investigate!
               iShock_IB(Shock_, iLine) = iShock - 1

               i = 0
               do while(divU(iShock + i).lt.-dLogRhoThreshold.and.(iShock+i).lt.iEnd.and.i.le.nSearchMax)
                    i = i + 1
               end do
               iShock_IB(ShockUp_, iLine) = iShock + i 
               
               i = 0
               do while(divU(iShock - i).lt.-dLogRhoThreshold.and.(iShock-i).gt.MinLagr(iLine).and.i.le.nSearchMax) 
                    i = i + 1
               end do
               iShock_IB(ShockDown_, iLine) = iShock - i - 2 

               ! Set density gradient to zero. Need to do this because the shock exists at t = 0, 
               ! skewing the calculation
               DensityGradient(:, iLine) = 0.0

          end do

     end subroutine get_test_shock
     !============================================================================
end module PT_ModShock

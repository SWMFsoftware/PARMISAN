!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModGrid

  ! Multi-line grid, D.Borovikov & I.Sokolov, Dec,17, 2017.
  ! Dec.23 2017: exclude fluxes from the state vector.
  ! Dec.25 2017: standard init and read_param
  ! Dec.25 2017: rename nVarRead=>nMhData, add NoShock_ param.

#ifdef OPENACC
   use ModUtilities, ONLY: norm2
#endif
   use ModUtilities, ONLY: CON_stop
   use PT_ModSize,  ONLY: nVertexMax, nDim
   use PT_ModProc,  ONLY: iProc
   use PT_ModConst, ONLY: cTwoPi, cPi
   use PT_ModConst, ONLY: cMu, cRsun, cProtonMass
   use PT_ModTime, ONLY: PTTime, DataInputTime

   implicit none
   SAVE

   PRIVATE ! except
   ! Public members:
   public:: read_param          ! read parameters related to grid
   public:: init                ! Initialize arrays on the grid
   public:: init_stand_alone    ! Initialize arrays on the grid
   public:: copy_old_state      ! save old arrays before getting new ones
   public:: get_other_state_var ! Auxiliary components of state vector
   public:: check_line_ishock   ! check if iShock exceeds iEnd of the field line

   ! Coordinate system and geometry
   character(len=3), public :: TypeCoordSystem = 'HGR'
   !
   ! Grid info
   ! Angular grid at origin surface
   integer, public :: nLon  = 1
   integer, public :: nLat  = 1
   integer:: iLatOffset = 0     ! Offset of line index along Lat grid
   integer:: iLonOffset = 0     ! Offset of line index along Lon grid

   ! Total number of magnetic field lines on all PEs
   ! (just a product of nLat * nLon)
   integer, public :: nLineAll = 1

   ! We do MPI on field lines:
   ! All nodes are enumerated. Last node number on the previous proc = iProc-1,
   ! equals (iProc*nLineAll)/nProc. Store this:
   integer, public :: iLineAll0

   ! The nodes on a given PE have node numbers ranging from iLineAll0 +1 to
   ! iNodeLast =((iProc + 1)*nLineAll)/nProc. The iLine index to enumerate
   ! lines on a given proc ranges from 1 to iNodeLast.
   ! nLine = nNodeLast - iLineAll0 is the number of
   ! lines (blocks) on this processor. For iLine=1:nLine
   ! iLineAll = iLineAll0+1:iNodeLast
   integer, public :: nLine
   
   ! Number of particles (vertexes, Lagrangian meshes) per line (line):
   integer, public, pointer :: nVertex_B(:), nVertex_BOld(:)

   ! Min and max indices to be used on fieldlines
   ! Index == Lagrangian Coordinate
   integer, public, pointer :: MinLagr(:), MaxLagr(:), &
                               MinLagrOld(:), MaxLagrOld(:)

   ! Function converting line number to lon-lat location of the line
   public :: iBlock_to_lon_lat

   ! Array for current and present location of shock wave
   integer, public, parameter:: nShockParam = 6,  &
         Shock_   = 1, & ! Current location of a shock wave
         ShockOld_= 2, &   ! Old location of a shock wave
         ShockUp_ = 3, &   ! Index of upstream extent of shock
         ShockDown_ = 4, & ! Index of downstream extent of shock
         ShockUpOld_ = 5, &
         ShockDownOld_ = 6

   integer, public, allocatable:: iShock_IB(:,:)
   integer, public, parameter:: NoShock_ = 1
   !
   ! Information about the magnetic field line foot point:
   ! the Lagrangian (0) and Cartesian (1:3) coordinates, and
   integer, public, parameter :: & ! init length of segment 1-2:
         Length_ = 4               ! control appending  new particles
   real, public, pointer :: FootPoint_VB(:,:)
   
   ! Logical to mark unusable lines
   logical, public, pointer :: Used_B(:)

   ! MHD state vector;
   ! 1st index - identification of variable (LagrID_:Wave2_)
   ! 2nd index - particle index along the field line
   ! 3rd index - local line number
   real, public, pointer     :: MhData_VIB(:,:,:)
   
   ! Aux state vector;
   ! 1st index - identification of variable (D_:BOld_)
   ! 2nd index - particle index along the field line
   ! 3rd index - local line number
   real, public, pointer     :: State_VIB(:,:,:)

   ! Boundary conditions for SDE - move to mod particle or create mod_bc
   real, public :: Rmax, Rmin

   ! Number of variables in the state vector and the identifications
   integer, public, parameter :: nMhData = 13, nVar = 27, &
         !
         LagrID_     = 0, & ! Lagrangian id           ^saved/   ^set to 0
         X_          = 1, & !                         |read in  |in copy_
         Y_          = 2, & ! Cartesian coordinates   |restart  |old_stat
         Z_          = 3, & !                         v/        |saved to
         Rho_        = 4, & ! Background plasma density         |mhd1
         T_          = 5, & ! Background temperature            |
         Ux_         = 6, & !                                   |may be
         Uy_         = 7, & ! Background plasma bulk velocity   |read from
         Uz_         = 8, & !                                   |mhd1
         Bx_         = 9, & !                                   |or
         By_         =10, & ! Background magnetic field         |received
         Bz_         =11, & !                                   |from
         Wave1_      =12, & !                                   |coupler
         Wave2_      =13, & !-Alfven wave turbulence            v
         !-
         R_          =14, & ! Heliocentric distance          ^derived from
         D_          =15, & ! Distance to the next particle  |MHdata in
         S_          =16, & ! Distance from the foot point   |get_other_
         U_          =17, & ! Plasma speed along field line  |state_var
         B_          =18, & ! Magnitude of magnetic field    v
         dB_         =19, & ! total alfven wave amplitude
         RhoOld_     =20, & ! Background plasma density      ! copy_
         UOld_       =21, & ! Background plasma bulk speed   ! old_
         BOld_       =22, & ! Magnitude of magnetic field    ! state
         SOld_       =23, & 
         DOld_       =24, &
         dBOld_      =25, &
         ROld_       =26, &
         LagrOld_    =27

   ! variable names
   character(len=10), public, parameter:: NameVar_V(LagrID_:nVar)&
      = ['LagrID    ', &
         'X         ', &
         'Y         ', &
         'Z         ', &
         'Rho       ', &
         'T         ', &
         'Ux        ', &
         'Uy        ', &
         'Uz        ', &
         'Bx        ', &
         'By        ', &
         'Bz        ', &
         'Wave1     ', &
         'Wave2     ', &
         'R         ', &
         'D         ', &
         'S         ', &
         'U         ', &
         'B         ', &
         'dB        ', &
         'RhoOld    ', &
         'UOld      ', &
         'BOld      ', &
         'SOld      ', &
         'DOld      ', &
         'dBOld     ', &
         'ROld      ', &
         'LagrOld   ' ]

   logical:: DoInit = .true.

contains
   !============================================================================
   subroutine read_param(NameCommand)
      use ModReadParam, ONLY: read_var
      character(len=*), intent(in):: NameCommand ! From PARAM.in
      ! Misc
      integer :: nParticleCheck, nLonCheck, nLatCheck
      character(len=*), parameter:: NameSub = 'read_param'
      !--------------------------------------------------------------------------
      select case(NameCommand)
      case('#CHECKGRIDSIZE')
         call read_var('nVertexMax',nParticleCheck)
         call read_var('nLon',     nLonCheck)
         call read_var('nLat',     nLatCheck)
         if(iProc==0.and.any([nLon,     nLat] /= [nLonCheck,nLatCheck])) &
            write(*,'(a,2I5)') 'nLon,nLat are reset to ', nLonCheck, nLatCheck
         nLon = nLonCheck
         nLat = nLatCheck
         nLineAll = nLon*nLat
         if(nParticleCheck > nVertexMax)then
            if(iProc==0)write(*,*)&
               'nVertexMax is too small, use ./Config.pl -g=',nParticleCheck
            call CON_stop('Code stopped')
         end if
      case('#COORDSYSTEM','#COORDINATESYSTEM')
         call read_var('TypeCoordSystem', TypeCoordSystem, &
            IsUpperCase=.true.)
      case('#GRIDNODE')
         call read_var('nLat',  nLat)
         call read_var('nLon',  nLon)
         nLineAll = nLat * nLon
      case('#BOUNDARY')
         call read_var('Rmin', Rmin)
         call read_var('Rmax', Rmax)

      case default
         call CON_stop(NameSub//' Unknown command '//NameCommand)
      end select

   end subroutine read_param
   !============================================================================
   subroutine init

      ! allocate the grid used in this model
      use ModUtilities,      ONLY: check_allocate
      use PT_ModProc,        ONLY: nProc

      integer:: iError
      integer:: iNodeLast

      character(len=*), parameter:: NameSub = 'init'
      !--------------------------------------------------------------------------
      if(.not.DoInit)RETURN
      DoInit = .false.

      ! distribute nodes between processors
      if(nLineAll >= nProc) then
         iLineAll0 = ( iProc   *nLineAll)/nProc
         iNodeLast = ((iProc+1)*nLineAll)/nProc
         nLine = iNodeLast-iLineAll0
      else
         ! there is one processor for each field line: we keep
         ! iProc = 0~nLineAll-1 working and others for the last line
         ! we also send the warning message for this over-request
         nLine = 1
         if(iProc < nLineAll) then
            iLineAll0 = iProc
         else
            ! Some work/trial has been done, but just partially. One can refer
            ! to the code version (915'th commit) on August 22, 2024.
            iLineAll0 = nLineAll-1
            write(*,*) "Here we keep iProc's >", nLineAll, 'on the last line.'
         end if
      end if

      ! check consistency
      if(nLat <= 0 .or. nLon <= 0)&
         call CON_stop(NameSub//': Origin surface grid is invalid')

      ! allocate data and grid containers
      allocate(iShock_IB(nShockParam, nLine), stat=iError)
      call check_allocate(iError, NameSub//'iShock_IB')
      iShock_IB = NoShock_

   end subroutine init
   !============================================================================
   subroutine init_stand_alone

      ! allocate the grid used in this model
      use ModUtilities,      ONLY: check_allocate
      integer :: iVertex, iError
      character(len=*), parameter:: NameSub = 'init_stand_alone'
      !--------------------------------------------------------------------------
      ! Allocate here if stand alone
      allocate(MhData_VIB(LagrID_:nMhData, 1:nVertexMax, nLine))
      MhData_VIB(1:nMhData,:,:) = 0.0

      ! reset lagrangian ids
      do iVertex = 1, nVertexMax
         MhData_VIB(LagrID_, iVertex, 1:nLine) = real(iVertex)
      end do

      ! Allocate auxiliary State vector
      allocate(State_VIB(nMhData+1:nVar, 1:nVertexMax, nLine))
      State_VIB = -1

      allocate(nVertex_B(nLine))
      allocate(nVertex_BOld(nLine))
      nVertex_B = 0
      nVertex_BOld = 0
      
      allocate(MinLagrOld(nLine))
      allocate(MaxLagrOld(nLine))
      allocate(MinLagr(nLine))
      allocate(MaxLagr(nLine))
      MinLagrOld = 0
      MinLagr = 0
      MaxLagrOld = 0
      MaxLagr = 0

      allocate(FootPoint_VB(LagrID_:Length_, nLine))
      FootPoint_VB = -1

      allocate(Used_B(nLine))
      Used_B = .true.

   end subroutine init_stand_alone
   !============================================================================
   subroutine iblock_to_lon_lat(iBlockIn, iLonOut, iLatOut)

      ! return angular grid's indexes corresponding to this line
      integer, intent(in) :: iBlockIn
      integer, intent(out):: iLonOut
      integer, intent(out):: iLatOut

      integer :: iLineAll
      !--------------------------------------------------------------------------
      !
      ! Get node number from line number
      !
      iLineAll = iBlockIn + iLineAll0
      iLatOut = 1 + (iLineAll - 1)/nLon
      iLonOut = iLineAll - nLon*(iLatOut - 1)

   end subroutine iblock_to_lon_lat
   !============================================================================
   subroutine copy_old_state

      ! copy current state to old state for all field lines
      integer:: i1, i2, iLine
      !--------------------------------------------------------------------------
      do iLine = 1, nLine
         
         if(.not.Used_B(iLine))CYCLE
         
         i1 = MinLagr(iLine)
         i2 = MaxLagr(iLine)
         
         iShock_IB(ShockOld_, iLine) = iShock_IB(Shock_, iLine)
         iShock_IB(ShockDownOld_, iLine) = iShock_IB(ShockDown_, iLine)
         iShock_IB(ShockUpOld_, iLine)   = iShock_IB(ShockUp_, iLine)

         State_VIB(RhoOld_, i1:i2, iLine) = MhData_VIB(Rho_, i1:i2, iLine)
         State_VIB(UOld_, i1:i2, iLine) = State_VIB(U_, i1:i2, iLine)
         State_VIB(BOld_, i1:i2, iLine) = State_VIB(B_, i1:i2, iLine)
         State_VIB(SOld_, i1:i2, iLine) = State_VIB(S_, i1:i2, iLine)
         State_VIB(DOld_, i1:i2, iLine) = State_VIB(D_, i1:i2, iLine)
         State_VIB(dBOld_, i1:i2, iLine) = State_VIB(dB_, i1:i2, iLine)
         State_VIB(ROld_, i1:i2, iLine) = State_VIB(R_, i1:i2, iLine)
         State_VIB(LagrOld_, i1:i2, iLine) = MhData_VIB(LagrID_, i1:i2, iLine)
         
         nVertex_BOld(iLine) = nVertex_B(iLine)
         MinLagrOld(iLine) = MinLagr(iLine)
         MaxLagrOld(iLine) = MaxLagr(iLine)

         ! reset variables read from file or received via coupler
         MhData_VIB(1:nMhData, i1:i2, iLine) = 0.0

      end do

   end subroutine copy_old_state
   !============================================================================
   subroutine get_other_state_var

      integer:: iLine, iVertex, iEnd
      integer:: iAux1, iAux2
      real   :: XyzAux1_D(x_:z_), XyzAux2_D(x_:z_)
      character(len=*), parameter:: NameSub = 'get_other_state_var'
      !--------------------------------------------------------------------------
      do iLine = 1, nLine
         if(.not.Used_B(iLine))CYCLE
         iEnd = MaxLagr(iLine)

         do iVertex = MinLagr(iLine), MaxLagr(iLine)
            ! Heliocentric Distance
            State_VIB(R_, iVertex, iLine) = &
               norm2(MhData_VIB(X_:Z_, iVertex, iLine))
            ! magnetic field
            State_VIB(B_,iVertex, iLine) = &
               norm2(MhData_VIB(Bx_:Bz_,iVertex,iLine))
            ! total wave amplitude
            State_VIB(dB_, iVertex, iLine) = cMu*(MhData_VIB(Wave1_, iVertex, iLine) + &
                                                  MhData_VIB(Wave2_, iVertex, iLine))
            ! Velocity vector along fieldline
            State_VIB(U_,iVertex,iLine) = sum(MhData_VIB(Ux_:Uz_,iVertex,iLine) * &
                     MHData_VIB(Bx_:Bz_,iVertex,iLine)) / State_VIB(B_,iVertex,iLine)
            ! distance between lagrangian points
            if(iVertex /= MaxLagr(iLine))then
               State_VIB(D_, iVertex, iLine) = norm2(&
                  MhData_VIB(X_:Z_, iVertex  , iLine) - &
                  MhData_VIB(X_:Z_, iVertex+1, iLine))
            end if
            ! distance from the beginning of the line
            if(iVertex == MinLagr(iLine))then
               State_VIB(S_, iVertex, iLine) = 0.0
            else
               State_VIB(S_, iVertex, iLine) = &
                  State_VIB(S_, iVertex-1, iLine) + &
                  State_VIB(D_, iVertex-1, iLine)
            end if

         end do
      end do

   end subroutine get_other_state_var
   !============================================================================
   subroutine check_line_ishock(iLine)

      ! check if the shock front is beyond the last coordinate of this field line
      ! in case that we accidentally have fewer particles on this field line
      integer, intent(in) :: iLine            ! index of line
      !--------------------------------------------------------------------------
      if(iShock_IB(Shock_, iLine) > MaxLagr(iLine)) then
         write(*,*) 'Shock in front of last coordinate', iShock_IB(Shock_, iLine), MaxLagr(iLine)
         Used_B(iLine) = .false.
      end if
   end subroutine check_line_ishock
   !============================================================================
end module PT_ModGrid
!==============================================================================

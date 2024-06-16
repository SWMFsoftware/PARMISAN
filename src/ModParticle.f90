  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModParticle
  use PT_ModConst
  implicit none
  PRIVATE ! Except
  integer, public :: nParticle = 0, nParticleMax = 1000000
  ! Particle arrays:
  ! 1. Named indexes:
  integer, parameter :: Coord_ = 1, p_ = 2, Weight_ = 3
  ! 2. Phase coords of particles:
  real, allocatable  :: Particle_VI(:,:)
  ! 3. To sort out the unused particles, left from the computational domain
  integer, allocatable :: iUsedParticle_I(:)
  integer :: nUsedParticle = 0
  ! 4. To add split particles:
  integer, allocatable :: iLev_I(:)
  integer :: nSplitParticle = 0
  ! Public members:
  public :: init_particles, complete_particle_loop
  ! Random numbers:
  integer :: Seed_I(630)
  ! Splitting:
  integer, public :: nSplitLev = 40
  real :: eSplitLevelMin = 1.d0*MeV      ! energy of first split level
  real :: eSplitLevelMax = 20000.d0*MeV  ! energy of last split level
  real, public :: eSplitLev_I(100)
contains
  !============================================================================
  subroutine init_particles
    use PT_ModPlot, ONLY: NamePath, name_sl
    use PT_ModProc, ONLY: iProc
    real :: dLogEsplit, xi
    ! Loop variable:
    integer :: iLev
    !--------------------------------------------------------------------------

    if(allocated(Particle_VI))RETURN
    ! Init random numbers
    Seed_I(630) = 123 + iProc
    call random_seed( put=Seed_I )
    call random_number( xi )

    allocate(Particle_VI(Coord_:Weight_,1:nParticleMax))
    allocate(iUsedParticle_I(1:nParticleMax))
    allocate(iLev_I(1:nParticleMax))
    nSplitParticle = 0; nUsedParticle = 0
    nParticle = 0
    ! Splitting:
    eSplitLev_I(1)           = eSplitLevelMin
    eSplitLev_I(1+nSplitLev) = eSplitLevelMax
    dLogEsplit = log(eSplitLev_I(1+nSplitLev)/eSplitLev_I(1))/real(nSplitLev)
    do iLev = 2, nSplitLev
       eSplitLev_I(iLev) = eSplitLev_I(iLev-1)*exp(dLogEsplit)
    end do
    if(iProc==0)then
       open(45,file=NamePath//name_sl,status='unknown')
       do iLev = 1, nSplitLev+1
          write(45,*)iLev, eSplitLev_I(iLev)/MeV
       end do
       close(45)
    end if
  end subroutine init_particles
  !============================================================================
  subroutine complete_particle_loop

    !--------------------------------------------------------------------------
    if(nParticle==0)then
       nSplitParticle = 0; nUsedParticle = 0
       RETURN
    end if
    if(nUsedParticle==0.and.nSplitParticle==0)then
       nParticle = 0
       RETURN
    end if
    if(nUsedParticle /= nParticle)then
       ! Wash out unused particles
       Particle_VI(Coord_:Weight_,1:nUsedParticle) = &
            Particle_VI(Coord_:Weight_,iUsedParticle_I(1:nUsedParticle))
       iLev_I(1:nUsedParticle) = iLev_I(iUsedParticle_I(1:nUsedParticle))
       if(nSplitParticle > 0)then
          ! Shift split particles to the left
          Particle_VI(Coord_:Weight_,&
               nUsedParticle+1:nUsedParticle+nSplitParticle) = &
               Particle_VI(Coord_:Weight_,&
               nParticle+1:nParticle+nSplitParticle)
          iLev_I(nUsedParticle+1:nUsedParticle+nSplitParticle) = &
               iLev_I(nParticle+1:nParticle+nSplitParticle)
       end if
       nParticle = nUsedParticle + nSplitParticle
       nSplitParticle = 0; nUsedParticle = 0
    else
       nParticle = nParticle + nSplitParticle
       nSplitParticle = 0; nUsedParticle = 0
    end if
  end subroutine complete_particle_loop
  !============================================================================
end module PT_ModParticle
!==============================================================================

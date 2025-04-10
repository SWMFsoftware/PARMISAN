!=====================================================================!
! Numerical Test of 1D Isotropic Parker Transport in Lagrangian Frame (SDE)
! Uses numerical shock used in Sokolov 2023
!        Modified slightly - no 0.5 offset between i, s

! Location of key parameters to change:
! SDE Type, shock width, Dxx - module Fieldline
! timestep - module Particle - subroutine initialize_particle
! number of particles to run (per core) - program TestParmisan
!=====================================================================!
module Fieldline
   use ModKind
   implicit none
   SAVE

   real, parameter :: tMin = 0.0, tMax = 10000.0
   real, parameter :: sMin = 0.0, sMax = 10000.0
   real, parameter :: Dxx = 10.0

   integer, parameter :: nDim = 2
   integer, parameter :: LagrCoord_ = 1, Momentum_ = 2
   real, parameter :: StratoFactor = 0.0 ! 1 if Ito, 0 if Stratonovich
   real, parameter :: dShock = 1.0 ! originally 0.5. how far in sL to calculate spacing between grid points
                                   ! this effectively changes the width of the shock = 2*dShock + 1
                                   ! 1 = acceleration region (constant)

contains
   !=====================================================================!
   subroutine get_distance_along_fieldline(Time, LagrCoord, S)
      real, intent(in) :: Time, LagrCoord
      real, intent(out) :: S
      real :: ShockTime, DeltaT

      call get_shock_arrival(LagrCoord, ShockTime)
      DeltaT = Time - ShockTime

      ! Instantaneous acceleration
      ! Change dSL to 1d-10 for true discontinuity
      ! --------------------------------
      ! if(DeltaT < 0) then                      ! stationary before shock arrival
      !    S = LagrCoord - 0.5
      ! else
      !    S = LagrCoord - 0.5 + 0.75 * DeltaT
      !    ! S = LagrCoord - 0.5 + 0.5 * DeltaT   ! constant speed after shock arrival
      ! end if

      ! Igor's method
      ! --------------------------------
      if(DeltaT < 0) then 
         S = LagrCoord 
      else if(DeltaT > 1) then
         S = LagrCoord + 0.375 + 0.75*(DeltaT - 1)
      else
         S = LagrCoord + 0.375*DeltaT**2
      end if
      
   end subroutine get_distance_along_fieldline
   !=====================================================================!
   subroutine get_ds(Time, LagrCoord, dS)
      real, intent(in) :: Time, LagrCoord
      real, intent(out) :: dS

      real :: Sup, Sdown, ShockTime, DeltaT

      call get_distance_along_fieldline(Time, LagrCoord - dShock, Sdown)
      call get_distance_along_fieldline(Time, LagrCoord + dShock, Sup)
      dS = (Sup - Sdown) / (2 * dShock)

   end subroutine get_ds
   !=====================================================================!
   subroutine get_rho(Time, LagrCoord, rho)
      real, intent(in) :: Time, LagrCoord
      real, intent(out) :: rho
      real :: dS

      call get_ds(Time, LagrCoord, dS)
      rho = 1.0 / dS

   end subroutine get_rho
   !=====================================================================!
   subroutine get_shock_arrival(LagrCoord, ShockTime)
      real, intent(in) :: LagrCoord
      real, intent(out) :: ShockTime

      ShockTime = LagrCoord + dShock

   end subroutine get_shock_arrival
   !=====================================================================!
   subroutine get_sde_coeffs(X_I, Time, Timestep, Wk, DriftCoeff, DiffCoeff)
         
      real, intent(in) :: X_I(nDim), Time, Timestep, Wk
      real, intent(out) :: DriftCoeff(nDim), DiffCoeff(nDim)

      real :: rho, rhoFuture, dRhodTau, dS, dSdown, dSup
      ! real :: newLagrCoord

      call get_ds(Time, X_I(LagrCoord_), dS)
      rho = 1.0 / dS

      ! found that this derivative converges around 0.01
      call get_rho(Time + 0.001, X_I(LagrCoord_), rhoFuture)
      dRhodTau = (rhoFuture - rho) / 0.001

      if(StratoFactor.eq.1) then ! ito

         call get_ds(Time, X_I(LagrCoord_) - dShock, dSdown)
         call get_ds(Time, X_I(LagrCoord_) + dShock, dSup) 
         DriftCoeff(LagrCoord_) = Dxx * (1.0/dSup -  1.0/dSdown) / (dS * 2 * dShock)
      else ! stratonovich
         DriftCoeff(LagrCoord_) = 0.0
      end if

      DiffCoeff(LagrCoord_) = sqrt(2.0 * Dxx) / dS

      DriftCoeff(Momentum_) = X_I(Momentum_) * dRhodTau / rho
      DiffCoeff(Momentum_) = 0.0

   end subroutine get_sde_coeffs
   !=====================================================================!
   subroutine get_sde_coeffs_milstein(X_I, Time, Timestep, DriftCoeff, DiffCoeff, dDiffdX)
      real, intent(in) :: X_I(nDim), Time, Timestep
      real, intent(out) :: DriftCoeff(nDim), DiffCoeff(nDim), dDiffdX(nDim)
      real :: LagrCoord, Momentum, dS, dSup, dSdown
      real :: rho, rhoFuture, dRhodTau

      LagrCoord = X_I(LagrCoord_)
      Momentum  = (3.0*X_I(Momentum_))**(1.0/3.0)

      call get_ds(Time, X_I(LagrCoord_), dS)
      call get_ds(Time, X_I(LagrCoord_) - dShock, dSdown)
      call get_ds(Time, X_I(LagrCoord_) + dShock, dSup)

      call get_rho(Time, X_I(LagrCoord_), rho)
      call get_rho(Time + 0.001, X_I(LagrCoord_), rhoFuture)

      dRhodTau = (rhoFuture - rho) / 0.001
      
      ! ==============================
      ! calculate sde coefficients 
      DriftCoeff(LagrCoord_) = Dxx * ( 1.0/dSup -  1.0/dSdown) / (2 * dShock * dS)
      DriftCoeff(Momentum_) = X_I(Momentum_) * dRhodTau / rho

      DiffCoeff(LagrCoord_) = sqrt(2.0 * Dxx) / dS
      DiffCoeff(Momentum_) = 0.0

      dDiffdX(LagrCoord_) = sqrt(2.0) * Dxx * (1.0/dSup**2 -  1.0/dSdown**2) / (2 * dShock)
      dDiffdX(Momentum_) = 0.0

   end subroutine get_sde_coeffs_milstein
   !=====================================================================!
   subroutine get_random_shock_location(Time, LagrCoord, S)
      real, intent(out) :: Time, LagrCoord, S
      real :: RandomUniform
      
      call random_number(RandomUniform)

      LagrCoord = sMin + (sMax - sMin) * RandomUniform
      call get_shock_arrival(LagrCoord, Time)

      call get_distance_along_fieldline(Time, LagrCoord, S)

   end subroutine get_random_shock_location
   !=====================================================================!
   subroutine compute_weight(Time, LagrCoord, Weight)
      real, intent(in) :: Time, LagrCoord
      real, intent(out) :: Weight
      real :: dS

      ! F = dS / B * f
      call get_ds(Time, LagrCoord, dS)
      Weight = dS

   end subroutine compute_weight
   !=====================================================================!
end module Fieldline
!==============================================================================!
module Solver
   use ModKind
   use Fieldline, ONLY: get_sde_coeffs, nDim, StratoFactor, get_sde_coeffs_milstein

   use PT_ModConst, ONLY: cFourPi
   implicit none
   SAVE
   real :: Wk, Sk

contains
   !==============================================================================
   subroutine milstein(X_I, Time, Timestep, Xnew_I)
       real, intent(in) :: X_I(nDim), Time, Timestep
       real, intent(out) :: Xnew_I(nDim)

       real :: DriftCoeff(nDim), DiffCoeff(nDim), dDiffdX(nDim)
       real :: RandomNormal

       call get_sde_coeffs_milstein(X_I, Time, Timestep, DriftCoeff, DiffCoeff, dDiffdX)
       
       call get_random_normal(RandomNormal)
       Wk = sqrt(Timestep)*RandomNormal

       Xnew_I = X_I + DriftCoeff * Timestep + DiffCoeff * Wk + &
                dDiffdX * (Wk**2 - Timestep)

   end subroutine milstein
   !==============================================================================
   subroutine get_random_normal(RandNormal1)
      ! returns random number sampled from normal distribution
      ! with mean = 0 and std = 1

      ! Box-Muller transformation limits the random variable to rn < ~6 
      ! physically limiting the size of the random diffusive process

      real, intent(out) :: RandNormal1 !, RandNormal2
      real :: RandUniform1, RandUniform2

      ! uniform random numbers over [0,1)
      call random_number(RandUniform1)
      call random_number(RandUniform2)

      ! redistribute to (0, 1] to avoid 0
      RandUniform1 = 1 - RandUniform1
      RandUniform2 = 1 - RandUniform2
      
      ! Box-Muller transformation
      ! two independent random variable with standard normal distribution
      ! only need one for 1-D version
      
      RandNormal1 = sqrt(-2*log(RandUniform1))*cos(0.5*cFourPi*RandUniform2)
      ! RandNormal2 = sqrt(-2*log(RandUniform1))*sin(0.5*fourpi*RandUniform2)

   end subroutine get_random_normal
   !==============================================================================
   subroutine rk2_sde(X_I, Time, Timestep, Xnew_I)
       ! Runge Kutta 2 Scheme for SDEs
       ! Roberts 2012  "Modify the improved Euler scheme..."
       real, intent(in) :: X_I(nDim), Time, Timestep
       real, intent(out) :: Xnew_I(nDim)

       real :: RandomNormal, RandomUniform
       real :: K1_I(nDim), K2_I(nDim) ! runge-kutta coefficients

       ! the same Wk and Sk needs to be used for both K1, K2
       call get_random_normal(RandomNormal)
       call random_number(RandomUniform)
       
       Sk = StratoFactor * sign(sqrt(Timestep), RandomUniform-0.5)
       Wk = sqrt(Timestep)*RandomNormal - Sk

       call calculate_k(X_I, Time, TimeStep, K1_I)

       Wk = Wk + 2*Sk

       call calculate_k(X_I + K1_I, Time + Timestep, Timestep, K2_I)

       Xnew_I = X_I + 0.5 * (K1_I + K2_I)

   end subroutine rk2_sde
   !==============================================================================
   subroutine calculate_k(X_I, Time, Timestep, K_I)
       ! Calculate K1 and K2 for second order RK scheme for SDE
       real, intent(in) :: X_I(nDim), Time, Timestep
       real, intent(out) :: K_I(nDim)

       integer :: i
       real(8) :: DriftCoeff(nDim), DiffCoeff(nDim) !, pDiff

       call get_sde_coeffs(X_I, Time, Timestep, Wk, DriftCoeff, DiffCoeff)

       do i = 1, nDim
           K_I(i) = Timestep * DriftCoeff(i) + Wk * DiffCoeff(i)
       end do

   end subroutine calculate_k
   !==============================================================================
end module Solver
!==============================================================================!
module Particle
   use PT_ModConst
   use Fieldline, ONLY: tMin, tMax, sMin, sMax, nDim, &
    get_distance_along_fieldline, compute_weight, get_random_shock_location

   use Solver, ONLY: rk2_sde, milstein

   implicit none
   SAVE

   real, allocatable :: Particle_V(:)

   integer, parameter :: LagrCoord_    = 1,  &
                         Momentum_     = 2,  &                      
                         Time_         = 3,  &
                         TimeStep_     = 4,  &
                         S_            = 5,  &
                         Energy_       = 6,  &
                         Weight_       = 7,  &
                         SOld_         = 8,  &
                         LagrCoordOld_ = 9,  &
                         MomentumOld_  = 10, &
                         nVar          = 10

contains
   !============================================================================
   subroutine initialize_particle
      
      if(.not.allocated(Particle_V)) allocate(Particle_V(1:nVar))

      call get_random_shock_location(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(S_))

      ! Get initial momentum
      Particle_V(Energy_) = 10.0 * ckeV
      call energy_to_momentum(Particle_V(Energy_), Particle_V(Momentum_))

      Particle_V(TimeStep_) = 0.01

      ! Get initial factor (1/ds)
      call compute_weight(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(Weight_))

      ! Save particle's current position as previous position
      call save_state()

   end subroutine initialize_particle
   !============================================================================
   subroutine advance_particle()

      real :: XOld(nDim), XNew(nDim)
      ! Save particle's current position as previous position
      call save_state

      ! Solve SDE - updates Lagrangian coordinate and momentum
      ! Equation advances p^3/3 not p
      XOld(LagrCoord_) = Particle_V(LagrCoord_)
      XOld(Momentum_)  = Particle_V(Momentum_)**3.0 / 3.0

      call rk2_sde(XOld, Particle_V(Time_), Particle_V(TimeStep_), XNew)
      
      ! Milstein method allows for calculation of timestep and weight during 
      ! timestep.
      ! call milstein(XOld, Particle_V(Time_), Particle_V(TimeStep_), XNew)

      Particle_V(LagrCoord_) = XNew(LagrCoord_)
      Particle_V(Momentum_) = (3.0*XNew(Momentum_))**(1.0/3.0)

      ! Update particle time then calculate new timestep
      Particle_V(Time_) = Particle_V(Time_) + Particle_V(TimeStep_)
      call get_distance_along_fieldline(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(S_))

      call momentum_to_energy(Particle_V(Momentum_), Particle_V(Energy_))

      ! Update particle weight
      call compute_weight(Particle_V(Time_), Particle_V(LagrCoord_), Particle_V(Weight_))

   end subroutine advance_particle
   !============================================================================
   subroutine save_state()

      Particle_V(SOld_) = Particle_V(S_)
      Particle_V(LagrCoordOld_) = Particle_V(LagrCoord_)
      Particle_V(MomentumOld_) = Particle_V(Momentum_)

   end subroutine save_state
   !============================================================================
   subroutine energy_to_momentum(Energy, Momentum)
      real, intent(in) :: Energy  ! Joules 
      real, intent(out) :: Momentum ! kg*m/s

      Momentum = sqrt(2.0*cProtonMass*Energy + (Energy/cLightSpeed)**2)

   end subroutine energy_to_momentum
   !=============================================================================
   subroutine momentum_to_energy(Momentum, Energy)
      real, intent(in) :: Momentum   ! kg*m/s 
      real, intent(out) :: Energy  ! joules

      Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - &
               cProtonRestEnergy

   end subroutine momentum_to_energy
   !=============================================================================
   subroutine check_boundary_conditions(isInside)
      logical, intent(out) :: isInside

      isInside = .true.
      if(Particle_V(S_) < sMin .or. Particle_V(S_) > sMax) isInside = .false.
      if((Particle_V(Time_) + Particle_V(TimeStep_)) > tMax) isInside = .false.

   end subroutine check_boundary_conditions
   !============================================================================
   subroutine debug_particle(Message)
      character(len = *), intent(in) :: Message

      write(*,*) '#', repeat('-', 25), '#'
      write(*,*) Message
      write(*, '(A, F14.4, A)') '  Time: ', Particle_V(Time_), ' s' 
      write(*, '(A, F10.4, A)') '  Timestep: ', Particle_V(TimeStep_), ' s'
      write(*, '(A, F16.4, A)') '  S(previous): ', Particle_V(SOld_)
      write(*, '(A, F16.4, A)') '  S(now): ', Particle_V(S_)
      write(*, '(A, F11.4, A)') '  LagrCoord(previous): ', Particle_V(LagrCoordOld_)
      write(*, '(A, F11.4, A)') '  LagrCoord(now): ', Particle_V(LagrCoord_)
      write(*, '(A, F10.4, A)') '  Energy: ', Particle_V(Energy_)/ckeV, ' keV'
      write(*, '(A, E10.4, A)') '  Momentum: ', Particle_V(Momentum_), ' kg*m/s'
      write(*,*) '#', repeat('-', 25), '#'

   end subroutine debug_particle
   !=============================================================================
   subroutine save_particle(FileName, FileUnit)
      character(len = *), intent(in) :: FileName
      integer, intent(in) :: FileUnit

      open(FileUnit, file = 'PT/IO2/'//FileName, status = 'unknown', position = 'append')
      write(FileUnit, *)  Particle_V(Time_), Particle_V(LagrCoord_), &
                          Particle_V(S_), Particle_V(Energy_)/ckeV
      close(FileUnit)

   end subroutine save_particle
   !=============================================================================
end module Particle
!==============================================================================!
module Plot

   use PT_ModConst
   use PT_ModProc, ONLY : iProc, iComm, iError
   use Fieldline, ONLY: tMin, tMax, sMin, sMax

   implicit none
   SAVE 
   
   character(len=*), parameter :: OutputDir = 'PT/IO2/'
   character(len=*), parameter :: TimeBinFile = 'time.dat'
   character(len=*), parameter :: EnergyBinFile = 'energy.dat'
   character(len=*), parameter :: rCrossFile = 'lagrCoord.dat'

   ! time and energy are bins: particle fell within t_i and t_i+1
   ! so timebins and energybins(nbins + 1)
   ! lagrcoord is crossing so lagrBins(nLagrbins)
   integer :: nTimeBins = 100, nEnergyBins = 100, nLagrbins = 100
   real    :: tBinMin = tMin, tBinMax = tMax
   real    :: eBinMin = 10.0, eBinMax = 100000.0
   real    :: iBinMin = sMin, iBinMax = sMax

   ! integer :: nTimeBins = 2, nEnergyBins = 200, nLagrbins = 500
   ! real    :: tBinMin = 0.0, tBinMax = 6000.0
   ! real    :: eBinMin = 10.0, eBinMax = 100000
   ! real    :: iBinMin = 1.0, iBinMax = 10000.0

   real, allocatable :: EnergyBins(:), TimeBins(:), LagrBins(:)
   real, allocatable :: Counts(:,:,:)

contains
   !============================================================================
   subroutine set_output_arrays

      integer :: iBin
      real :: dT, dLogE, dLagr
      real :: EnergyGrid(nEnergyBins), TimeGrid(nTimeBins)

      eBinMin = eBinMin * ckeV
      eBinMax = eBinMax * ckeV

      allocate(EnergyBins(nEnergyBins + 1))
      allocate(TimeBins(nTimeBins + 1))

      allocate(LagrBins(nLagrbins))
      allocate(Counts(nTimeBins, nEnergyBins, nLagrbins))
      Counts = 0.0

      dT = (tBinMax - tBinMin) / (nTimeBins)
      do iBin = 1, nTimeBins+1
         TimeBins(iBin) = tBinMin + dT * (iBin-1) 
      end do

      dLogE = (log10(eBinMax) - log10(eBinMin)) / (nEnergyBins)
      do iBin = 1, nEnergyBins+1
         EnergyBins(iBin) = 10.0**(log10(eBinMin)+ dLogE*(iBin-1))
      end do
  
      dLagr = (iBinMax - iBinMin) / (nLagrbins-1)
      do iBin = 1, nLagrbins
         LagrBins(iBin) = iBinMin + dLagr * (iBin - 1)
      end do
  
      TimeGrid = 0.5 * (TimeBins(1:nTimeBins) + TimeBins(2:nTimeBins+1))
      EnergyGrid = 10**(0.5 * (log10(EnergyBins(1:nEnergyBins)) + log10(EnergyBins(2:nEnergyBins+1))))

      if(iProc==0) then
         open(801,file=OutputDir//EnergyBinFile,status='unknown')
         write(801,'(1000e15.6)')EnergyGrid/ckeV
         close(801)
         open(801,file=OutputDir//TimeBinFile,status='unknown')
         write(801,'(1000e15.6)')TimeGrid
         close(801)
         open(801,file=OutputDir//rCrossFile,status='unknown')
         write(801,'(1000e15.6)')LagrBins
         close(801)
      end if

   end subroutine set_output_arrays
   !============================================================================
   subroutine bin_particle(LagrCoordOld, LagrCoord, Time, Energy, Weight)
      real, intent(in) :: LagrCoordOld, LagrCoord, Time, Energy, Weight
      
      logical :: DidCross(nLagrbins)
      integer :: iT, iE

      DidCross = (LagrCoord - LagrBins) * (LagrCoordOld - LagrBins) < 0.0
      
      if(.not.any(DidCross)) return

      if(Time.ge.TimeBins(nTimeBins+1)) return
      if(Energy.ge.EnergyBins(nEnergyBins+1)) return

      if(Time.lt.TimeBins(1)) return
      if(Energy.lt.EnergyBins(1)) return

      ! iT = minloc(abs(TimeBins - Time), DIM = 1)
      ! iE = minloc(abs(EnergyBins - Energy), DIM = 1)

      iT = minloc(Time - TimeBins, mask = (Time-TimeBins > 0), dim = 1)
      iE = minloc(Energy - EnergyBins, mask = (Energy-EnergyBins > 0), dim = 1)

      where(DidCross) Counts(iT, iE, :) = Counts(iT, iE, :) + Weight

   end subroutine bin_particle
   !============================================================================
   subroutine save_counts
      
      use ModMpi, ONLY: MPI_reduce_real_array, MPI_SUM
      
      integer :: iTimeBin, iCross
      character(len = 50) :: CrossingFile

      integer :: totalBins

      totalBins = nTimeBins * nEnergyBins * nLagrbins
      call MPI_reduce_real_array(Counts, totalBins, MPI_SUM, 0, iComm, iError)
 
      if(iProc/=0) RETURN
      
      do iCross = 1, nLagrbins
         
         write(CrossingFile, '(A13, I0)') 'fluxes_cross_', iCross
         CrossingFile = adjustl(CrossingFile)
         open(901,file=OutputDir//trim(CrossingFile),&
               status='unknown',action="READWRITE")
         do iTimeBin = 1, nTimeBins
            write(901,'(1000e15.6)')Counts(iTimeBin,:,iCross)
         end do
         close(901)
      end do
      
   end subroutine save_counts
   !============================================================================
end module Plot
!==============================================================================!
program TestParmisan

   use ModKind
   use ModMpi
   use ModUtilities, ONLY: CON_stop

   use PT_ModProc, ONLY: iProc, nProc, iComm, iError

   use Particle
   use Solver
   use Plot

   implicit none
   
   character(len = 20) :: iProcStr, nParticleStr
   character(len = 30) :: ParticleFileName
   logical :: isInside

   integer, allocatable :: Seed_I(:)
   integer :: nSeed, fileUnit, istat

   integer :: n, nParticlePerProc = 1, nSaveParticle = 1
   !---------------------------------------------------------- !
   call MPI_INIT( iError )
   
   ! Assign communicator
   iComm = MPI_COMM_WORLD
   call MPI_COMM_RANK( iComm, iProc, iError )
   call MPI_COMM_SIZE( iComm, nProc, iError )

   ! for testing purposes
   write(iProcStr, *) iProc
   iProcStr = adjustl(iProcStr)

   ! -------------------------------------------------------------------- !
   ! get size of seed (compiler-dependent)
   call random_seed(size = nSeed)
   allocate(Seed_I(nSeed))

   ! use dev/urandom to set seed for each processor
   open(newunit=fileUnit, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
       read(fileUnit) Seed_I
       close(fileUnit)
       call random_seed(put = Seed_I)
   else
       call CON_stop('/dev/urandom not found')
   end if
   ! -------------------------------------------------------------------- !

   call set_output_arrays

   ! ------------------- Run particles to save -------------------- !
   do n = 1, nSaveParticle
      write(nParticleStr, *) n
      nParticleStr = adjustl(nParticleStr)
      ParticleFileName = 'particle_'//trim(nParticleStr)//'_iProc_'//trim(iProcStr)//'.dat'
      
      call initialize_particle
      call save_particle(ParticleFileName, iProc)
      isInside = .true.

      do while(isInside)
         
         call advance_particle
         call save_particle(ParticleFileName, iProc)
         call bin_particle(Particle_V(LagrCoordOld_), Particle_V(LagrCoord_), &
                           Particle_V(Time_), Particle_V(Energy_), Particle_V(Weight_))

         call check_boundary_conditions(isInside)
      end do
   end do

   ! ------------------ Run remainder of particles ---------------------- !
   do n = nSaveParticle + 1, nParticlePerProc

      call initialize_particle
      isInside = .true.

      do while(isInside)
         
         call advance_particle
         call bin_particle(Particle_V(LagrCoordOld_), Particle_V(LagrCoord_), &
                           Particle_V(Time_), Particle_V(Energy_), Particle_V(Weight_))

         call check_boundary_conditions(isInside)
      end do
   end do

   call save_counts
   call MPI_finalize(iError)  
end program TestParmisan
!==============================================================================!
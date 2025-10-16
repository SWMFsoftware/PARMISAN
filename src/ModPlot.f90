  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModPlot
  
  use PT_ModConst, only: ckeV
  use PT_ModProc, ONLY : iProc, iComm, iError
  use PT_ModUnit, ONLY : kinetic_energy_to_momentum
  
  implicit none
  SAVE
  
  character(len=*), parameter :: OutputDir = 'PT/IO2/'
  character(len=*), parameter :: SplitFile = 'split_levels.dat'
  character(len=*), parameter :: TimeFile =  'time.dat'
  character(len=*), parameter :: EnergyBinFile = 'energy_bin.dat'
  character(len=*), parameter :: LagrBinFile = 'lagr_bin.dat'
  character(len=*), parameter :: RFile = 'R_bin.dat'
  
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: CON_stop

    character(len=*), intent(in):: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#PLOT')
      ! Maybe something will go here one day...
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine save_bin_arrays
    use ModUtilities, only: touch_file
    use PT_ModDistribution, only: EnergyBin_I, LagrBin_I

    open(801,file=OutputDir//EnergyBinFile,status='unknown')
    write(801,'(3000e15.6)')EnergyBin_I/ckeV
    close(801)
    
    ! Create time file - will be appended with the output time each timestep
    call touch_file(OutputDir//TimeFile)
    
    open(801,file=OutputDir//LagrBinFile,status='unknown')
    write(801,'(3000e15.6)')LagrBin_I
    close(801)

    ! Create file that will hold radial distance - time dependent
    call touch_file(OutputDir//RFile)

  end subroutine save_bin_arrays
  !============================================================================
  subroutine save_position_bins(Time)

    use PT_ModDistribution, only: LagrBin_I, nLagrBins
    use PT_ModFieldline, only: get_particle_location

    real, intent(in) :: Time
    real :: RBins_I(nLagrBins+1)
  
    integer :: iBin

    do iBin = 1, nLagrBins+1
      call get_particle_location(Time, LagrBin_I(iBin), RBins_I(iBin))
    end do

    open(911, file=OutputDir//RFile, position="append", action="write")
      write(911, '(4000e15.6)') RBins_I
    close(911)

  end subroutine save_position_bins
  !============================================================================
  subroutine save_distribution_function(Iter, Time)
    use PT_ModDistribution, only: Counts_II, nEnergyBins, nLagrBins, TotalWeight, &
                                  calculate_distribution_function
    use ModMpi, ONLY: MPI_reduce_real_array, MPI_SUM, MPI_reduce_real_scalar
    real, intent(in) :: Time
    integer, intent(in) :: Iter
    integer :: iTime, iLagr
    character(len = 50) :: outputFile

    ! Save fluxes stored at different radial distances
    !--------------------------------------------------------------------------
    call MPI_reduce_real_array(Counts_II, nEnergyBins*nLagrBins, MPI_SUM, 0,&
         iComm, iError)
    call MPI_reduce_real_scalar(TotalWeight, MPI_SUM, 0, iComm, iError)

    ! Save data
    ! Calculate conversion from F = ds/B*f 
    ! Divide counts by total weight
    if(iProc/=0) then
      ! Reset Counts for next timestep
      Counts_II = 0.0
    else
      ! Convert counts to distribution function - variable name does not change
      call calculate_distribution_function(Time)

      ! append time of output to time file
      open(901, file=OutputDir//TimeFile, position="append", action="write")
        write(901, *) Time
      close(901)

      call save_position_bins(Time)
      ! create output distribution function file name
      write(outputFile, '(A15, I0)') 'distfunc_iter_', int(Time)
      outputFile = adjustl(outputFile)

      ! output distribution function (nEnergy x nLagrCoord)
      open(902, file=OutputDir//trim(outputFile), status='unknown', action="READWRITE")
      do iLagr = 1, nLagrBins
        write(902,'(3000e15.6)') Counts_II(:, iLagr) 
      end do
      close(902)

      ! Reset Counts for next timestep
      Counts_II = 0.0
    end if

  end subroutine save_distribution_function
  !============================================================================
  subroutine save_analytic_solution()
    ! currently assumes: 
    !     constant Dxx along entire fieldline
    use PT_ModGrid, only: State_VIB, U_, nVertex_B
    use PT_ModFieldline, only: DxxConst
    use PT_ModDistribution, only: EnergyBin_I, nEnergyBins
    use PT_ModUnit, only: kinetic_energy_to_momentum
    use PT_ModParticle, only: E0

    real :: UpstreamU, DownstreamU, PowerLaw, P0, Pnorm
    real, allocatable :: fSteadyState(:), Tacc(:)
    integer :: i

    if(iProc.ne.0) return

    allocate(fSteadyState(1:nEnergyBins+1))
    allocate(Tacc(1:nEnergyBins+1))

    ! Shock moves at constant lagrangian speed - shift to shock frame (-1.0)
    UpstreamU = State_VIB(U_, nVertex_B(1)-5, 1) - 1.0
    DownstreamU = State_VIB(U_, 5, 1) - 1.0
    PowerLaw = -3.0 * UpstreamU / (UpstreamU - DownstreamU)
    P0 = kinetic_energy_to_momentum(E0)

    do i = 1, nEnergyBins+1
      Pnorm = kinetic_energy_to_momentum(EnergyBin_I(i)) / P0
      fSteadyState(i) = Pnorm ** PowerLaw
      Tacc(i) = 3.0 * DxxConst * (UpstreamU**-1.0 + DownstreamU**-1.0) &
                / (UpstreamU - DownstreamU) * log(Pnorm)
    end do

    open(801,file=OutputDir//'steadystate_f.dat',status='unknown')
    write(801,'(1000e15.6)') fSteadyState
    close(801)

    open(801,file=OutputDir//'acceleration_time.dat',status='unknown')
    write(801,'(1000e15.6)') Tacc
    close(801)


  end subroutine save_analytic_solution
  !============================================================================
end module PT_ModPlot
!==============================================================================

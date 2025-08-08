  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModPlot
  
  use ModKind
  use PT_ModConst
  use PT_ModProc, ONLY : iProc, iComm, iError
  ! use PT_ModFieldline, ONLY : compute_conversion_factor
  use PT_ModTestFieldline, ONLY : compute_conversion_factor
  use PT_ModParticle, ONLY : energy_to_momentum

  implicit none
  SAVE
  
  character(len=*), parameter :: OutputDir = 'PT/IO2/'
  character(len=*), parameter :: SplitFile = 'split_levels.dat'
  character(len=*), parameter :: TimeFile =  'time.dat'
  character(len=*), parameter :: EnergyBinFile = 'energy_bin.dat'
  character(len=*), parameter :: LagrBinFile = 'lagr_bin.dat'
  
  integer :: nSaveTimes, nEnergyBins, nLagrBins
  real :: tSaveMin, tSaveMax, eBinMin, eBinMax, lagrBinMin, lagrBinMax

  real :: TotalWeight

  ! Time bin
  real, parameter :: timeWindow = 5.0
  ! Array for storing the counts
  real, allocatable :: Counts_III(:,:,:)
  
  ! Bins
  real, allocatable :: Ebin_I(:), LagrBin_I(:), Time_I(:)
  
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

       call read_var('nLagrBins', nLagrBins)
       call read_var('lagrBinMin', lagrBinMin)
       call read_var('lagrBinMax', lagrBinMax)
       call read_var('nSaveTimes', nSaveTimes)
       call read_var('tSaveMin', tSaveMin)
       call read_var('tSaveMax', tSaveMax)
       call read_var('nEnergyBins', nEnergyBins)
       call read_var('eBinMin', eBinMin)
       call read_var('eBinMax', eBinMax)

      ! Convert to seconds
      tSaveMin = tSaveMin * 3600.
      tSaveMax = tSaveMax * 3600.

      ! Convert to keV
      eBinMin = eBinMin * ckeV
      eBinMax = eBinMax * ckeV

    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine set_bins()

    integer :: iBin
    real :: dLogE, dT, dR

    allocate(Counts_III(nSaveTimes, nEnergyBins, nLagrBins)); Counts_III = 0.0
    allocate(Time_I(nSaveTimes+1))
    allocate(Ebin_I(nEnergyBins+1))
    allocate(LagrBin_I(nLagrBins+1))

    dT = (tSaveMax - tSaveMin) / (nSaveTimes)
    do iBin = 1, nSaveTimes+1
      Time_I(iBin) = tSaveMin + dT * (iBin-1) 
    end do

    dLogE = (log10(eBinMax) - log10(eBinMin)) / (nEnergyBins)
    do iBin = 1, nEnergyBins+1
       Ebin_I(iBin) = 10.0**(log10(eBinMin)+ dLogE*(iBin-1))
    end do

    dR = (lagrBinMax - lagrBinMin) / (nLagrBins)
    do iBin = 1, nLagrBins+1
      LagrBin_I(iBin) = lagrBinMin + dR * (iBin - 1)
    end do

    if(iProc==0) then
       open(801,file=OutputDir//EnergyBinFile,status='unknown')
       write(801,'(1000e15.6)')Ebin_I/ckeV
       close(801)
       open(801,file=OutputDir//TimeFile,status='unknown')
       write(801,'(1000e15.6)')Time_I
       close(801)
       open(801,file=OutputDir//LagrBinFile,status='unknown')
       write(801,'(1000e15.6)')LagrBin_I
       close(801)
    end if

  end subroutine set_bins
  !============================================================================
  subroutine increase_total_weight(Weight)
    real, intent(in) :: Weight
    TotalWeight = TotalWeight + Weight
  end subroutine increase_total_weight
  !============================================================================
  subroutine bin_particle(LagrCoord, Time, Energy, Weight)
    ! Bin in time, energy, and location
    real, intent(in) :: LagrCoord, Time, Energy, Weight
    integer :: iL, iE, iT, i

    if(LagrCoord.lt.LagrBin_I(1).or.LagrCoord.ge.LagrBin_I(nLagrBins+1)) return
    if(Time.lt.Time_I(1).or.Time.ge.Time_I(nSaveTimes+1)) return
    if(Energy.lt.Ebin_I(1).or.Energy.ge.Ebin_I(nEnergyBins+1)) return

    !--------------------------------------------------------------------------
    iT = minloc(Time - Time_I, mask = (Time - Time_I > 0), dim = 1) 
    iL = minloc(LagrCoord - LagrBin_I, mask = (LagrCoord - LagrBin_I > 0), dim = 1)
    iE = minloc(Energy - Ebin_I, mask = (Energy - Ebin_I > 0), dim = 1)
    Counts_III(iT, iE, iL) = Counts_III(iT, iE, iL) + Weight
  end subroutine bin_particle
  !============================================================================
  subroutine calculate_distribution_function
    ! Convert counts (sum of weights) in bin to distribution function
    ! divide by total weight, divide by bin widths, multiply by conversion factor
    ! Bins are momentum, time, and lagrcoord
    ! Conversion factor is F = dS/B * f. Solved for F, want f.
    
    real :: dP, p1, p2
    real :: dLagr, factor, conversion, dTime
    integer :: iE, iT, iL

    do iE = 1, nEnergyBins
      call energy_to_momentum(Ebin_I(iE), p1)
      call energy_to_momentum(Ebin_I(iE+1), p2)
      dP = p2**3.0 / 3.0 - p1**3.0 / 3.0

      do iT = 1, nSaveTimes
        dTime = Time_I(iT+1) - Time_I(iT)
        do iL = 1, nLagrBins
          dLagr = LagrBin_I(iL+1) - LagrBin_I(iL)
          factor = TotalWeight * dLagr * dP * dTime
          call compute_conversion_factor(Time_I(iT)+0.5*dTime, LagrBin_I(iL)+0.5*dLagr, conversion)
          Counts_III(iT, iE, iL) = Counts_III(iT, iE, iL) * conversion / factor
        end do
      end do

    end do

  end subroutine
  !============================================================================
  subroutine save_output

    use ModMpi, ONLY: MPI_reduce_real_array, MPI_SUM, MPI_reduce_real_scalar
    integer :: iTime, iLagr
    character(len = 50) :: outputFile

    ! Save fluxes stored at different radial distances
    !--------------------------------------------------------------------------
    call MPI_reduce_real_array(Counts_III, nSaveTimes*nEnergyBins*nLagrBins, MPI_SUM, 0,&
         iComm, iError)
    call MPI_reduce_real_scalar(TotalWeight, MPI_SUM, 0, iComm, iError)

    ! Save data
    ! Calculate conversion from F = ds/B*f 
    ! Divide counts by total weight
    if(iProc/=0) RETURN

    write(*,*) 'TotalWeight = ', TotalWeight
    call calculate_distribution_function

    do iTime = 1, nSaveTimes
      
      write(outputFile, '(A15, I0)') 'distfunc_iTime_', iTime
      outputFile = adjustl(outputFile)

      open(901,file=OutputDir//trim(outputFile),status='unknown',action="READWRITE")

      do iLagr = 1, nLagrBins
        write(901,'(2000e15.6)') Counts_III(iTime, :, iLagr) 
      end do

      close(901)

    end do

  end subroutine save_output
!============================================================================
end module PT_ModPlot
!==============================================================================

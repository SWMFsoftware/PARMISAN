  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModPlot
  
  use PT_ModConst
  use PT_ModProc, ONLY : iProc, iComm, iError
  use PT_ModFieldline, ONLY : compute_conversion_factor, get_lagr_coord
  ! use PT_ModTestFieldline, ONLY : compute_conversion_factor, get_lagr_coord
  use PT_ModParticle, ONLY : energy_to_momentum

  implicit none
  SAVE
  
  character(len=*), parameter :: OutputDir = 'PT/IO2/'
  character(len=*), parameter :: SplitFile = 'split_levels.dat'
  character(len=*), parameter :: TimeFile =  'time.dat'
  character(len=*), parameter :: EnergyBinFile = 'energy_bin.dat'
  character(len=*), parameter :: DistBinFile = 'r_bin.dat'
  
  integer :: NumTimes, nEnergyBins, NumLocations
  real :: tBinMin, tBinMax, eBinMin, eBinMax, rBinMin, rBinMax
  real :: LagrBinSize, TimeBinSize
  real :: TotalWeight

  ! Array for storing the counts
  real, allocatable :: Counts_III(:,:,:)
  
  ! Bins
  real, allocatable :: eBin_I(:), rBin_I(:), tBin_I(:), lagrBin_I(:,:)
  
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

       call read_var('NumLocations', NumLocations)
       call read_var('rBinMin', rBinMin)
       call read_var('rBinMax', rBinMax)
       call read_var('LagrBinSize', LagrBinSize)
       call read_var('NumTimes', NumTimes)
       call read_var('tBinMin', tBinMin)
       call read_var('tBinMax', tBinMax)
       call read_var('TimeBinSize', TimeBinSize)
       call read_var('nEnergyBins', nEnergyBins)
       call read_var('eBinMin', eBinMin)
       call read_var('eBinMax', eBinMax)

      ! Convert to seconds
      tBinMin = tBinMin * 3600.
      tBinMax = tBinMax * 3600.

      ! Convert to meters
      rBinMin = rBinMin * cRsun
      rBinMax = rBinMax * cRsun

      ! Convert to keV
      eBinMin = eBinMin * ckeV
      eBinMax = eBinMax * ckeV

    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine set_bins()

    integer :: iE, iR, iT
    real :: dLogE, dT, dR

    allocate(Counts_III(NumTimes, nEnergyBins, NumLocations)); Counts_III = 0.0
    allocate(tBin_I(NumTimes))
    allocate(eBin_I(nEnergyBins+1))
    allocate(rBin_I(NumLocations))
    allocate(lagrBin_I(NumTimes, NumLocations))

    ! centers of time and spatial bins
    ! lagr bins are function of time
    ! width of time, lagr bins are from param.in and used in bin_particle
    dT = (tBinMax - tBinMin) / (NumTimes-1)
    dR = (rBinMax - rBinMin) / (NumLocations-1)
    do iT = 1, NumTimes
      tBin_I(iT) = tBinMin + dT * (iT-1)
      do iR = 1, NumLocations
        rBin_I(iR) = rBinMin + dR * (iR - 1)
        call get_lagr_coord(tBin_I(iT), rBin_I(iR), lagrBin_I(iT, iR))
      end do 
    end do

    ! edges of energy bins
    dLogE = (log10(eBinMax) - log10(eBinMin)) / (nEnergyBins)
    do iE = 1, nEnergyBins+1
       eBin_I(iE) = 10.0**(log10(eBinMin)+ dLogE*(iE-1))
    end do

    if(iProc==0) then
       open(801,file=OutputDir//EnergyBinFile,status='unknown')
       write(801,'(1000e15.6)')eBin_I/ckeV
       close(801)
       open(801,file=OutputDir//TimeFile,status='unknown')
       write(801,'(1000e15.6)')tBin_I
       close(801)
       open(801,file=OutputDir//DistBinFile,status='unknown')
       write(801,'(1000e15.6)')rBin_I/cRsun
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
    ! Time and energy, to determine the bin to put the contribution:
    real, intent(in) :: LagrCoord, Time, Energy, Weight
    integer :: iL, iE, iT, i

    !--------------------------------------------------------------------------
    if(Energy.lt.eBin_I(1).or.Energy.ge.eBin_I(nEnergyBins+1)) return
    if(Time.lt.tBin_I(1).or.Time.ge.tBin_I(NumTimes)) return

    !--------------------------------------------------------------------------
    iT = -1
    do i = 1, NumTimes
        if((Time.le.tBin_I(i) + 0.5*TimeBinSize).and.(Time.gt.tBin_I(i) - 0.5*TimeBinSize)) iT = i
    end do
 
    if(iT.eq.-1) return

    !--------------------------------------------------------------------------
    iL = -1
    do i = 1, NumLocations
      if((LagrCoord.le.lagrBin_I(iT, i) + 0.5*LagrBinSize).and.(LagrCoord.gt.lagrBin_I(iT, i) - 0.5*LagrBinSize)) iL = i
    end do

    if(iL.eq.-1) return

    !--------------------------------------------------------------------------
    iE = minloc(Energy - eBin_I, mask = ((Energy - eBin_I).ge.0), dim = 1)
    Counts_III(iT, iE, iL) = Counts_III(iT, iE, iL) + Weight

  end subroutine bin_particle
  !============================================================================
  subroutine calculate_distribution_function
    ! Convert counts (sum of weights) in bin to distribution function
    ! divide by total weight, divide by bin widths, multiply by conversion factor
    ! Bins are momentum and lagrcoord taken at a snapshot in time
    ! Conversion factor is F = dS/B * f. Solved for F, want f.
    real :: dP, p1, p2
    real :: dLagr, factor, conversion
    integer :: iE, iT, iL

    do iE = 1, nEnergyBins
      call energy_to_momentum(eBin_I(iE), p1)
      call energy_to_momentum(eBin_I(iE+1), p2)
      dP = p2**3.0 / 3.0 - p1**3.0 / 3.0

      do iT = 1, NumTimes

        do iL = 1, NumLocations
          factor = TotalWeight * LagrBinSize * dP * TimeBinSize
          call compute_conversion_factor(tBin_I(iT), lagrBin_I(iT, iL) - 0.5 * LagrBinSize, &
                                         lagrBin_I(iT, iL) + 0.5*LagrBinSize, conversion)
          ! TODO: Missing multiplication by timestep - not sure how to do this with adaptive stepping
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
    call MPI_reduce_real_array(Counts_III, NumTimes*nEnergyBins*NumLocations, MPI_SUM, 0,&
         iComm, iError)
    call MPI_reduce_real_scalar(TotalWeight, MPI_SUM, 0, iComm, iError)

    ! Save data
    ! Calculate conversion from F = ds/B*f 
    ! Divide counts by total weight
    !  TODO: divide by bin widths here so output is f
    if(iProc/=0) RETURN

    write(*,*) 'TotalWeight = ', TotalWeight
    call calculate_distribution_function

    do iTime = 1, NumTimes
      
      write(outputFile, '(A15, I0)') 'distfunc_iTime_', iTime
      outputFile = adjustl(outputFile)

      open(901,file=OutputDir//trim(outputFile),status='unknown',action="READWRITE")

      do iLagr = 1, NumLocations
        write(901,'(1000e15.6)') Counts_III(iTime, :, iLagr)
      end do

      close(901)

    end do

  end subroutine save_output
!============================================================================
end module PT_ModPlot
!==============================================================================

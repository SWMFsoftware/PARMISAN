!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModDistribution
    use ModKind
    use PT_ModConst
    use PT_ModProc, ONLY : iProc, iComm, iError

    implicit none
    SAVE

    integer :: nEnergyBins, nLagrBins
    real :: eBinMin, eBinMax, lagrBinMin, lagrBinMax

    real :: TimeWindow
    real :: TotalWeight

    ! Array for storing the counts
    real, allocatable :: Counts_II(:,:)

    ! Bins
    real, allocatable :: EnergyBin_I(:), LagrBin_I(:)

contains
!============================================================================
    subroutine read_param(NameCommand)

        use ModReadParam, ONLY: read_var
        use ModUtilities, ONLY: CON_stop

        character(len=*), intent(in):: NameCommand ! From PARAM.in
        character(len=*), parameter:: NameSub = 'read_param'
        !--------------------------------------------------------------------------
        select case(NameCommand)
        case('#DISTRIBUTION')

        call read_var('nLagrBins', nLagrBins)
        call read_var('lagrBinMin', lagrBinMin)
        call read_var('lagrBinMax', lagrBinMax)
        call read_var('nEnergyBins', nEnergyBins)
        call read_var('eBinMin', eBinMin)
        call read_var('eBinMax', eBinMax)
        call read_var('TimeWindow', TimeWindow)

        ! Convert to joules
        eBinMin = eBinMin * ckeV
        eBinMax = eBinMax * ckeV

        case default
        call CON_stop(NameSub//' Unknown command '//NameCommand)
        end select
    end subroutine read_param
!============================================================================
    subroutine init

        integer :: iBin
        real :: dLogE, dL

        allocate(Counts_II(nEnergyBins, nLagrBins)); Counts_II = 0.0
        allocate(EnergyBin_I(nEnergyBins+1))
        allocate(LagrBin_I(nLagrBins+1))

        dLogE = (log10(eBinMax) - log10(eBinMin)) / (nEnergyBins)
        do iBin = 1, nEnergyBins+1
        EnergyBin_I(iBin) = 10.0**(log10(eBinMin)+ dLogE*(iBin-1))
        end do

        dL = (lagrBinMax - lagrBinMin) / (nLagrBins)
        do iBin = 1, nLagrBins+1
        LagrBin_I(iBin) = lagrBinMin + dL * (iBin - 1)
        end do

        TotalWeight = 0.0

    end subroutine init
!============================================================================
    subroutine increase_total_weight(Weight)    
        real, intent(in) :: Weight
        TotalWeight = TotalWeight + Weight
    end subroutine increase_total_weight
!============================================================================
    subroutine bin_particle(LagrCoord, Energy, Weight)
        ! Bin in time, energy, and location
        real, intent(in) :: LagrCoord, Energy, Weight
        integer :: iL, iE, i

        if(LagrCoord.lt.LagrBin_I(1).or.LagrCoord.ge.LagrBin_I(nLagrBins+1)) return
        if(Energy.lt.EnergyBin_I(1).or.Energy.ge.EnergyBin_I(nEnergyBins+1)) return
        !--------------------------------------------------------------------------
        iL = minloc(LagrCoord - LagrBin_I, mask = (LagrCoord - LagrBin_I > 0), dim = 1)
        iE = minloc(Energy - EnergyBin_I, mask = (Energy - EnergyBin_I > 0), dim = 1)

        Counts_II(iE, iL) = Counts_II(iE, iL) + Weight
    end subroutine bin_particle
!============================================================================
    subroutine calculate_distribution_function(Time)
        ! Convert counts (sum of weights) in bin to distribution function
        ! divide by total weight, divide by bin widths, multiply by conversion factor
        ! Bins are momentum cubed over 3, time, and lagrcoord
        ! Conversion factor is F = dS/B * f. Solved for F, want f.

        use PT_ModFieldline, only: compute_conversion_factor
        use PT_ModUnit, only: kinetic_energy_to_momentum

        real, intent(in) :: Time
        real :: dP, p1, p2
        real :: dLagr, factor, conversion, dTime
        integer :: iE, iT, iL

        do iE = 1, nEnergyBins
            p1 = kinetic_energy_to_momentum(EnergyBin_I(iE))
            p2 = kinetic_energy_to_momentum(EnergyBin_I(iE+1))
            dP = (p2**3.0 / 3.0) - (p1**3.0 / 3.0)

            do iL = 1, nLagrBins
                dLagr = LagrBin_I(iL+1) - LagrBin_I(iL)
                factor = dLagr * dP * TimeWindow * TotalWeight
                call compute_conversion_factor(Time - 0.5*TimeWindow,    &
                                               LagrBin_I(iL)+0.5*dLagr,  &
                                               conversion)
                Counts_II(iE, iL) = Counts_II(iE, iL) * conversion / factor
            end do

        end do
    end subroutine
!============================================================================
end module PT_ModDistribution
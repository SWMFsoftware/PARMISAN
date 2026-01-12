!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModDistribution
    use ModKind
    use PT_ModConst
    use PT_ModProc, ONLY : iProc, iComm, iError
    use PT_ModUnit, only: kinetic_energy_to_momentum

    implicit none
    SAVE

    integer :: nEnergyBins, nLagrBins
    integer :: InjEnergyIndex
    real :: eBinMin, eBinMax
    real :: TimeWindow, TotalWeight, InjEnergy, InjCoeff 
    ! Array for storing the counts
    real, allocatable :: Counts_II(:,:), Counts_Outside_I(:)
    ! Bins
    real, allocatable :: EnergyBin_I(:), LagrBin_I(:)
    ! How many total injected real protons 
    ! integral of f(p_inj) over phase space)
    real :: IntegralOverfInj

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
        call read_var('nEnergyBins', nEnergyBins)
        call read_var('eBinMin', eBinMin)
        call read_var('eBinMax', eBinMax)
        call read_var('InjEnergy', InjEnergy)
        call read_var('InjCoeff', InjCoeff)
        call read_var('TimeWindow', TimeWindow)

        ! Convert to joules
        eBinMin = eBinMin * ckeV
        eBinMax = eBinMax * ckeV
        InjEnergy = InjEnergy * ckeV

        case default
        call CON_stop(NameSub//' Unknown command '//NameCommand)
        end select
    end subroutine read_param
!============================================================================
    subroutine init

        integer :: iBin
        real :: dLogE, dL

        ! allocate bin and solution arrays
        allocate(Counts_II(nEnergyBins, nLagrBins)); Counts_II = 0.0
        allocate(Counts_Outside_I(nEnergyBins));     Counts_Outside_I = 0.0
        allocate(EnergyBin_I(nEnergyBins+1))
        allocate(LagrBin_I(nLagrBins+1))

        ! energy bins equally space in log10
        dLogE = (log10(eBinMax) - log10(eBinMin)) / (nEnergyBins)
        EnergyBin_I(1) = eBinMin
        do iBin = 2, nEnergyBins+1
            EnergyBin_I(iBin) = 10.0**(log10(eBinMin)+ dLogE*(iBin-1))
        end do

        ! index corresponding to bin of injected distribution function
        InjEnergyIndex = minloc(InjEnergy - EnergyBin_I, &
                                mask = (InjEnergy - EnergyBin_I >= 0), dim = 1)

        ! initialize total weight and particle number integral
        TotalWeight = 0.0
        IntegralOverfInj = 0.0

    end subroutine init
!============================================================================
    subroutine set_lagr_bins(iLine)
        ! Set the lagrangian phase space bins for this timestep
        ! Bins the entire fieldline
        use PT_ModGrid, only: MinLagr, MaxLagr
        
        integer, intent(in) :: iLine
        
        integer :: iBin
        real :: dL

        dL = real(MaxLagr(iline) - 1 - MinLagr(iLine)) / real(nLagrBins)
        do iBin = 1, nLagrBins+1
            LagrBin_I(iBin) = real(MinLagr(iLine)) + dL * (iBin - 1)
        end do

    end subroutine set_lagr_bins
!============================================================================
    subroutine update_integral_over_finj(Time, LagrInject)
        
        ! Equation 124 of Sokolov et al., 2023 "High resolution finite...".
        ! Integrate injection distribution function over phase space.
        ! Integration over momentum is only over inj momentum bin.
        ! dSl assumed to be 1. 
        ! Factor of 4pi*Psi does not need to be 
        !   included as it will cancel during normalization.

        ! TODO: This will need to be iLine-dependent

        use PT_ModFieldline, only: calc_thermal_energy, calc_dSOverB
        use PT_ModUnit, only: kinetic_energy_to_momentum
        real, intent(in) :: Time, LagrInject

        real :: pInj, fPinj, p1, p2, dP 
        real :: dSOverB, ThermalEnergy

        ! this is kb*rho*T
        call calc_thermal_energy(Time, LagrInject, ThermalEnergy)
        call calc_dSOverB(Time, LagrInject, dSOverB)
        
        ! Inject particles with p**-5 distribution commonly seen in 
        ! suprathermal tail
        pInj = kinetic_energy_to_momentum(InjEnergy)
        fPinj = InjCoeff * ThermalEnergy / (cPi * pInj**5.0)
        ! call save_finject(Time, LagrInject, fPinj)

        ! Momentum bin
        p1 = kinetic_energy_to_momentum(EnergyBin_I(InjEnergyIndex))
        p2 = kinetic_energy_to_momentum(EnergyBin_I(InjEnergyIndex+1))
        dP = (p2**3.0 - p1**3.0) / 3.0

        IntegralOverfInj = IntegralOverfInj + fPinj * dP * dSOverB

    end subroutine update_integral_over_finj
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

        ! if particle is outside bounds of phase space bins - return
        if(LagrCoord.lt.LagrBin_I(1).or.LagrCoord.ge.LagrBin_I(nLagrBins+1)) return
        if(Energy.lt.EnergyBin_I(1).or.Energy.ge.EnergyBin_I(nEnergyBins+1)) return
        !--------------------------------------------------------------------------
        ! index of bins particle is in
        iL = minloc(LagrCoord - LagrBin_I, mask = (LagrCoord - LagrBin_I > 0), dim = 1)
        iE = minloc(Energy - EnergyBin_I, mask = (Energy - EnergyBin_I > 0), dim = 1)
        ! print *, LagrCoord, Energy
        ! Increase the count in this bin by the weight of particle
        Counts_II(iE, iL) = Counts_II(iE, iL) + Weight

    end subroutine bin_particle
!============================================================================
    subroutine calculate_distribution_function(Time)
        ! Conversion of SDE solution --> PDE solution (F = dS/B * f)
        ! PDF of SDE trajectories in phase space = solution of FP equation
        ! PDF of SDE = P(traj in bin) divided by phase space bin volume
        ! P(traj in bin) = sum of weights in bin / totalWeight
        ! Renormalization of the distribution function is necessary
        ! Renormalize by ensuring the integral over phase space of
        ! f_solution = integral over phase space of f_inject
        ! TODO: need to account for escaped particles!

        use PT_ModFieldline, only: calc_dSOverB
        use PT_ModUnit, only: kinetic_energy_to_momentum

        real, intent(in) :: Time
        real :: dP, dLagr, p1, p2
        real :: BinVolume, IntegralOverF, dSOverB
        integer :: iE, iL

        IntegralOverF = 0.0
        do iL = 1, nLagrBins
            dLagr = LagrBin_I(iL+1) - LagrBin_I(iL)
            do iE = 1, nEnergyBins

                call calc_dSOverB(Time, LagrBin_I(iL) + 0.5 * dLagr, &
                                  dSOverB)

                p1 = kinetic_energy_to_momentum(EnergyBin_I(iE))
                p2 = kinetic_energy_to_momentum(EnergyBin_I(iE+1))
                dP = (p2**3.0 / 3.0) - (p1**3.0 / 3.0)
                BinVolume = dLagr * dP

                ! probability / phase space bin volume
                Counts_II(iE, iL) = Counts_II(iE, iL) / &
                                    (BinVolume * TotalWeight)

                ! integrate over bin and add to total integral
                IntegralOverF = IntegralOverF + Counts_II(iE, iL) * dP * dLagr
                ! Multiple by (B/dS) to convert from F -> f
                Counts_II(iE, iL) = Counts_II(iE, iL) / dSOverB
            end do
        end do 

        ! renormalize such that the integral over phase space is conserved
        Counts_II = Counts_II * IntegralOverfInj / IntegralOverF

    end subroutine
!============================================================================
end module PT_ModDistribution
module PT_ModFieldline
    use ModKind
    use ModMpi
    use ModUtilities, ONLY: find_cell
    use PT_ModRandom, ONLY: get_random_normal

    use PT_ModConst

    implicit none
    
    SAVE
    ! number of dimensions solved by SDE, (passed to ModParticle and ModSolver !)
    ! should probably not be here, need to put in ModMain, ModInitialize, or something higher up
    integer, parameter :: nDim = 2
    integer, parameter :: LagrCoord_ = 1, Momentum_ = 2

    integer :: nT, nS ! number of time and spatial (lagr) grid points
    real :: tMin, tMax, rMin, rMax, shockMaxR, rMinInject, rMaxInject ! extreme values
    real :: TimeStepFactor, MaxTimeStep
    real :: StratoFactor ! 1 if Ito, 0 if Stratonovich

    ! Arrays are data read in from field line files: size (nT, nLagr)
    ! Already in lagrangian coordinates!
    ! TODO: Change array names to match SWMF naming convention
    real, allocatable :: time_Array(:), S_Array(:,:), rho_Array(:,:), &
                         B_Array(:,:), U_Array(:,:), R_Array(:,:), &
                         shockLocation(:), AlfvenTurb_Array(:,:), &
                         Temp_Array(:, :)

contains

   !============================================================================
    subroutine read_param(NameCommand)

        use ModReadParam, ONLY: read_var
        use ModUtilities, ONLY: CON_stop

        character(len=*), intent(in):: NameCommand ! From PARAM.in
        character(len=*), parameter:: NameSub = 'read_param'
        character(len=1) :: Scheme
        !--------------------------------------------------------------------------
        select case(NameCommand)
        case("#SCHEME")
            call read_var('Scheme', Scheme)
        case('#BC')
            call read_var('rMin', rMin )
            call read_var('rMax', rMax)
            call read_var('tMax', tMax)
            call read_var('rMinInject', rMinInject)
            call read_var('rMaxInject', rMaxInject)
            call read_var('TimeStepFactor', TimeStepFactor)
            call read_var('MaxTimeStep', MaxTimeStep)

            ! Convert to m
            rMin = rMin * cRsun
            rMax = rMax * cRsun
            rMinInject = rMinInject * cRsun
            rMaxInject = rMaxInject * cRsun

        case default
            call CON_stop(NameSub//' Unknown command '//NameCommand)
        end select

        select case(Scheme)
        case('I')
            StratoFactor = 1.0
        case('S')
            StratoFactor = 0.0
        case default
            write(*,*) 'Unknown scheme chosen. Use "I" or "S". Using Stratonovich'
            StratoFactor = 0.0
        end select

    end subroutine read_param
    !============================================================================
    subroutine read_fieldline
        use PT_ModProc, ONLY: iProc
        character(len=*), parameter :: FileLocation = 'PT/Param/'
        integer :: i, row, col
        
        ! ====================================
        ! read in time array
        open(8, file = FileLocation//'Lagr_time.dat', status = 'old')
        read(8,*) nT, nS

        allocate(time_Array(nT))
        
        do i = 1, nT
            read(8,*) time_Array(i)
        end do
        close(8)

        tMin = time_Array(1)

        if(tMax.gt.(time_Array(nT)-2*MaxTimeStep)) then
            tMax = time_Array(nT)-2*MaxTimeStep
           if(iProc.eq.0) write(*,*) 'tMax too large, setting tMax = ', tMax, ' s'
        end if

        ! ====================================
        ! read in data along fieldline at each time step
        allocate(S_Array(nT, nS))
        open(8, file = FileLocation//'Lagr_S.dat', status = 'old')
        read(8, *) ((S_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(rho_Array(nT, nS))
        open(8, file = FileLocation//'Lagr_rho.dat', status = 'old')
        read(8, *) ((rho_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(B_Array(nT, nS))
        open(8, file = FileLocation//'Lagr_B.dat', status = 'old')
        read(8, *) ((B_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(U_Array(nT, nS))
        open(8, file = FileLocation//'Lagr_U.dat', status = 'old')
        read(8, *) ((U_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(R_Array(nT, nS))
        open(8, file = FileLocation//'Lagr_R.dat', status = 'old')
        read(8, *) ((R_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(Temp_Array(nT, nS))
        open(8, file = FileLocation//'Lagr_T.dat', status = 'old')
        read(8, *) ((Temp_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(AlfvenTurb_Array(nT, nS))
        open(8, file = FileLocation//'Lagr_AlfvenTurb.dat', status = 'old')
        read(8, *) ((AlfvenTurb_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        ! ====================================
        ! convert to SI
        S_Array = S_Array * cRsun  ! [m]
        R_Array = R_Array * cRsun  ! [m]
        ! U = [m/s], rho = [kg/m^3]
        ! B = [T]
        ! AlfenTurb = [T**2] dB**2 = mu0 * (wave1 + wave2)
        ! T = [keV] ! convert to K?

        if(maxval(R_Array).lt.rMax) then
            rMax = maxVal(R_Array)
            if(iProc.eq.0) write(*, *) 'rMax too large, setting rMax = ', rMax/cRsun
        end if

        if(minval(R_Array).gt.rMin) then
            rMin = minVal(R_Array)
            if(iProc.eq.0) write(*, *) 'rMin too small, setting rMin = ', rMin/cRsun
        end if
        
        ! ====================================
        ! read in shock location
        allocate(shockLocation(nT))
        open(8, file = FileLocation//'shockLoc.dat', status = 'old')
        do i = 1, nT
            read(8,*) shockLocation(i)
        end do
        close(8)
        shockLocation = shockLocation * cRsun
        shockMaxR = shockLocation(nT)

    end subroutine read_fieldline
    !============================================================================
    subroutine calculate_ds(Time, LagrCoord, dS)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: dS
        real :: sUp, sDown

        call get_array_value(Time, LagrCoord + 0.5, S_Array, sUp)
        call get_array_value(Time, LagrCoord - 0.5, S_Array, sDown)
        if(sUp.lt.sDown) write(*,*) 'sup < sdown: ', Time, LagrCoord
        dS = sUp - sDown

    end subroutine calculate_ds
    !============================================================================
    subroutine calc_drhodt(Time, LagrCoord, Timestep, dRhodT)
        real, intent(in) :: Time, LagrCoord, Timestep
        real, intent(out) :: dRhodT
        real :: RhoNow, RhoPast

        call get_array_value(Time, LagrCoord, rho_Array, RhoNow)
        call get_array_value(Time - Timestep, LagrCoord, rho_Array, RhoPast)

        dRhodT = (RhoNow - RhoPast) / (Timestep)

    end subroutine calc_drhodt
    !============================================================================
    subroutine get_sde_coeffs(X_I, Time, Timestep, DriftCoeff, DiffCoeff)
        
        real, intent(in) :: X_I(nDim), Time, Timestep
        real, intent(out) :: DriftCoeff(nDim), DiffCoeff(nDim)
    
        integer :: iS
        real :: LagrCoord, Momentum
        real :: S, dS, sUp, sDown, iSfrac, Sfuture, Spast
        real :: B, Bup, Bdown, dBdS 
        real :: Dxx, DxxUp, DxxDown, dDxxdS
        real :: rho, dRhodT, RhoFuture, RhoPast
        real :: dtmax

        !! TODO: How to sync nDim and indices with ModParticle

        LagrCoord = X_I(LagrCoord_)
        Momentum  = (3.0*X_I(Momentum_))**(1.0/3.0)
   
        ! ==============================
        ! get variables at particle location
        call get_array_value(Time, LagrCoord, B_Array, B)
        call get_array_value(Time, LagrCoord, rho_Array, rho)
        call get_diffusion_coeff(Time, LagrCoord, Momentum, Dxx)

        ! ==============================
        ! Calculate drho/dt - central difference with dt = 1.0 second
        call calc_drhodt(Time, LagrCoord, 1.0, dRhodT)
        ! ==============================
        ! calculate spatial derivatives at i + 0.5, i - 0.5
        call calculate_ds(Time, LagrCoord, dS)

        call get_array_value(Time, LagrCoord + 0.5, B_Array, Bup)
        call get_diffusion_coeff(Time, LagrCoord + 0.5, Momentum, DxxUp)

        call get_array_value(Time, LagrCoord - 0.5, B_Array, Bdown)
        call get_diffusion_coeff(Time, LagrCoord - 0.5, Momentum, DxxDown)

        dBdS = (Bup - Bdown) / dS
        dDxxdS = (DxxUp - DxxDown) / dS

        ! ==============================
        ! calculate sde coefficients
        ! (0.5 + 0.5*StratoFactor) takes into account difference in drift
        ! coefficient for Ito (S = 1) and Stratonovich (S = 0)
        DriftCoeff(LagrCoord_) = (0.5 + 0.5*StratoFactor) * dDxxdS / dS - Dxx * dBdS / (B*dS)
        DiffCoeff(LagrCoord_) = sqrt(2.0 * Dxx) / dS

        DriftCoeff(Momentum_) = X_I(Momentum_) * dRhodT / rho
        DiffCoeff(Momentum_) = 0.0

    end subroutine get_sde_coeffs
    !============================================================================
    ! subroutine get_sde_coeffs_milstein(X_I, Timestep, Weight, DriftCoeff, DiffCoeff, dDiffdX)
        
        !     real, intent(in) :: X_I(nDim)
        !     real, intent(out) :: Timestep, Weight, DriftCoeff(nDim), DiffCoeff(nDim), dDiffdX(nDim)
        
        !     integer :: iS
        !     real :: LagrCoord, Momentum
        !     real :: S, dS, sUp, sDown, dSup, dSdown
        !     real :: B, Bup, Bdown
        !     real :: aDerivative, cDerivative 
        !     real :: Dxx, DxxUp, DxxDown
        !     real :: rho, dRhodT
        !     real :: RhoPast, RhoFuture

        !     !! TODO: How to sync nDim and indices with ModParticle

        !     LagrCoord = X_I(LagrCoord_)
        !     Momentum  = (3.0*X_I(Momentum_))**(1.0/3.0)
        !     ! ==============================
        !     ! particle location at [Time, LagrCoord]
        !     ! Time can be different from LagrTime because of RK2 scheme
        !     call get_distance_along_fieldline(LagrCoord, LagrTime, S)

        !     ! ==============================
        !     ! get variables at particle location
        !     call get_array_value(S, LagrTime, B_Array, B)
        !     call get_array_value(S, LagrTime, rho_Array, rho)
        !     call get_diffusion_coeff(S, LagrTime, Momentum, Dxx)
        !     call calculate_ds(LagrCoord, LagrTime, dS)

        !     ! ==============================
        !     ! Calculate drho/dt - forward difference
        !     ! call calc_drhodt(LagrCoord, LagrTime, 1.0, dRhodT)

        !     call get_array_value(S, LagrTime-Timestep, rho_array, RhoPast)
        !     call get_array_value(S, LagrTime+Timestep, rho_Array, RhoFuture)
        !     dRhodT = (RhoFuture - RhoPast) / (2.0)


        !     ! ==============================
        !     ! calculate quantities at i + 0.5, i - 0.5
        !     call get_distance_along_fieldline(LagrCoord + 0.5, LagrTime, sUp)
        !     call calculate_ds(LagrCoord + 0.5, LagrTime, dSup)
        !     call get_array_value(sUp, LagrTime, B_Array, Bup)
        !     call get_diffusion_coeff(sUp, LagrTime, Momentum, DxxUp)

        !     call get_distance_along_fieldline(LagrCoord - 0.5, LagrTime, sDown)
        !     call calculate_ds(LagrCoord - 0.5, LagrTime, dSdown)
        !     call get_array_value(sDown, LagrTime, B_Array, Bdown)
        !     call get_diffusion_coeff(sDown, LagrTime, Momentum, DxxDown)

        !     aDerivative = DxxUp / (Bup * dSup) - DxxDown / (Bdown * dSdown)
        !     cDerivative = DxxUp / dSup**2 - DxxDown / dSdown**2
        !     ! ==============================
        !     ! calculate sde coefficients 
        !     DriftCoeff(LagrCoord_) = B * aDerivative / dS
        !     DriftCoeff(Momentum_) = X_I(Momentum_) * dRhodT / rho

        !     DiffCoeff(LagrCoord_) = sqrt(2.0 * Dxx) / dS
        !     DiffCoeff(Momentum_) = 0.0

        !     dDiffdX(LagrCoord_) = 0.5 * cDerivative
        !     dDiffdX(Momentum_) = 0.0

        !     Timestep = TimeStepFactor * (DiffCoeff(LagrCoord_)/DriftCoeff(LagrCoord_))**2
        !     ! Set maximum allowed timestep, 10 seconds
        !     TimeStep = min(TimeStep, MaxTimeStep)

        !     Weight = B / dS

    ! end subroutine get_sde_coeffs_milstein
    !============================================================================
    subroutine get_diffusion_coeff(Time, LagrCoord, Momentum, Dxx)
        real, intent(in) :: Time, LagrCoord, Momentum
        real, intent(out) :: Dxx

        real :: R, B, dB
        ! real :: Rshock, DeltaR
        ! integer :: iTime
        ! real :: iTimeFrac
        
        ! get heliocentric distance
        call get_array_value(Time, LagrCoord, R_Array, R)
        ! write(*,*) 'R'
        call get_array_value(Time, LagrCoord, B_Array, B)
        ! write(*,*) 'B'
        call get_array_value(Time, LagrCoord, AlfvenTurb_Array, dB)
        ! write(*,*) 'dB'
        call calc_mhd_dxx(R, B, dB, Momentum, Dxx)

        ! call get_empirical_dxx(R, Momentum, Dxx)
        ! get shock location
        ! call find_cell(1, nT, Time, iTime, iTimeFrac, time_Array)
        ! Rshock = shockLocation(iTime) + iTimeFrac * (shockLocation(iTime+1) - shockLocation(iTime))
        
        ! DeltaR = R - Rshock

        ! if(DeltaR.gt.DxxShockDistance) then
        !     call get_empirical_dxx(R, Momentum, Dxx)
        ! else

        !     call get_array_value(S, Time, B_Array, B)
        !     call get_array_value(S, Time, AlfvenTurb_Array, dB)
        !     call get_mhd_dxx(R, B, dB, Momentum, Dxx)
        ! end if

    end subroutine get_diffusion_coeff
    !============================================================================
    subroutine calc_mhd_dxx(R, B, dB, Momentum, Dxx)
        real, intent(in) :: R, Momentum, B, dB
        real, intent(out) :: Dxx

        real :: ConstantFactor = 81.0 / (7.0*cPi) * (0.5/cPi)**(2.0/3.0)
        real :: Velocity, MeanFreePath, Lmax, BTotal

        Velocity = Momentum / cProtonMass
        Velocity = sqrt(Velocity**2 / (1 + (Velocity/cLightSpeed)**2))

        Lmax = 0.4*R
        Btotal = sqrt(B**2 + dB)
        MeanFreePath = ConstantFactor * BTotal**2 * &
                       (Momentum*Lmax**2/(BTotal * cElectronCharge))**(1.0/3.0) / dB

        Dxx = MeanFreePath * Velocity / 3.0
        Dxx = max(Dxx, 1.0d4 * cRsun)

    end subroutine calc_mhd_dxx
    !============================================================================
    subroutine get_empirical_dxx(R, Momentum, Dxx)
        real, intent(in) :: R, Momentum
        real, intent(out) :: Dxx
        real :: Dxx0, Energy

        ! Uses formula from Chen et al., 2024 for Dxx
        ! Hard-coded numbers are taken straight from Chen et al., 2024 Equation 4
        Dxx0 = 5.16d14 ! [m^2/s]
        Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - cProtonRestEnergy

        Dxx = Dxx0 * (R / cAU)**(1.17) * (Energy / ckeV)**(0.71)

        ! Add uncertainity in PSP measurements?? Only slowed code down and didn't change results
        ! Needs more testing. 
        ! Added in uncertainity by sampling from normal distributions with given errors

        ! call get_random_normal(RandNormal)
        ! Dxx0 = 5.16d14 + 1.22d14*RandNormal! cm^2 / s

        ! ! Cap at 3-sigma value, otherwise this can go negative
        ! Dxx0 = max(Dxx0, 1.5d14)
        
        ! Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - cProtonRestEnergy

        ! call get_random_normal(RandNormal)
        ! rExponent = 1.17 + 0.08*RandNormal

        ! Dxx = Dxx0 * (R / cAU)**(rExponent) * (Energy / ckeV)**(0.71)

    end subroutine get_empirical_dxx
    !============================================================================
    subroutine get_array_value(Time, LagrCoord, Array, InterpValue)
        ! get value of array at (t, i)
        ! uses bilinear interpolation
        real, intent(in) :: Time, LagrCoord, Array(:,:)
        real, intent(out) :: InterpValue

        integer :: iTime, iLagr        ! (lower) index of desired t, i
        real :: iTimeFrac, iLagrFrac ! fractional indices
        real :: f1, f2 ! intermediate values, interpolated in lagr but not t

        if(Time.gt.tMax.or.Time.lt.tMin) write(*,*) 'Time = ', Time
        if(LagrCoord.gt.nS.or.LagrCoord.lt.1) write(*,*) 'LagrCoord = ', LagrCoord
        call find_cell(1, nT, Time, iTime, iTimeFrac, time_Array)

        iLagr = floor(LagrCoord)
        iLagrFrac = LagrCoord - iLagr

        f1 = (1-iLagrFrac) * Array(iTime, iLagr) + iLagrFrac * Array(iTime, iLagr+1)
        f2 = (1-iLagrFrac) * Array(iTime+1, iLagr) + iLagrFrac * Array(iTime+1, iLagr+1)
        InterpValue = (1-iTimeFrac) * f1 + iTimeFrac*f2
    
    end subroutine get_array_value
    !============================================================================
    subroutine get_random_shock_location(Time, LagrCoord, R)
        real, intent(out) :: Time, LagrCoord, R

        real :: RandomUniform
        real :: iTimeFrac, iSFrac1, iSFrac2
        real :: func1, func2 

        integer :: iTime, iS1, iS2
        logical :: doExtrapolate = .true.
 
        ! uniform random numbers over [0,1)
        call random_number(RandomUniform)

        ! inverse transform sampling
        ! desired pdf is ~1/r^2 from [rMinInject, rMaxInject]
        ! equation for r is F(F^-1) = ru
        ! where F is the CDF of the desired normalized pdf
        ! R = rMaxInject*rMinInject / (rMaxInject*(1-RandomUniform) + rMinInject*RandomUniform)

        ! Uniform sampling in R
        R = rMinInject + (rMaxInject - rMinInject) * RandomUniform
     
        ! Find time when shock is at this radial distance
        call find_cell(1, nT, R, iTime, iTimeFrac, shockLocation)
        Time = time_Array(iTime) + iTimeFrac * (time_Array(iTime+1) - time_Array(iTime))

        ! Get lagrangian coordinate at this time
        call get_lagr_coord(Time, R, LagrCoord)

    end subroutine get_random_shock_location
    !=============================================================================
    subroutine get_lagr_coord(Time, R, LagrCoord)
        real, intent(in) :: Time, R
        real, intent(out) :: LagrCoord

        real :: iTimeFrac, iSFrac1, iSFrac2
        real :: f1, f2

        integer :: iTime, iS1, iS2
        logical :: doExtrapolate = .true.

        real :: deltaT, tempR(nS)
        
        ! Find time when shock is at this radial distance
        call find_cell(1, nT, Time, iTime, iTimeFrac, time_Array)

        deltaT = Time - time_Array(iTime)
        tempR = R_Array(iTime, :) + U_Array(iTime, :) * deltaT

        call find_cell(1, nS, R, iS1, iSFrac1, tempR, doExtrapolate)
        LagrCoord = iS1 + iSFrac1

        ! ! Get lagrangian coordinate at this time
        ! call find_cell(1, nS, R, iS1, iSFrac1, R_Array(iTime, :), doExtrapolate)
        ! call find_cell(1, nS, R, iS2, iSFrac2, R_Array(iTime+1, :), doExtrapolate)    
        ! write(*,*) 'iS1, iS2: ', iS1, iS2
        ! f1 = R - (1-iTimeFrac)*R_Array(iTime, iS1) - iTimeFrac * R_Array(iTime+1, iS1)
        ! if(f1.gt.0) then
        !     f2 = (1-iTimeFrac) * (R_Array(iTime, iS1+1) - R_Array(iTime, iS1)) + &
        !           iTimeFrac * (R_Array(iTime+1, iS1+1) - R_Array(iTime+1, iS1))
        !     LagrCoord = iS1 + f1/f2
        ! else
        !     f1 = R - (1-iTimeFrac)*R_Array(iTime, iS2) - iTimeFrac * R_Array(iTime+1, iS2)
        !     ! if(f1.lt.0) write(*,*) iTimeFrac, iSFrac1, iSFrac2
        !     f2 = (1-iTimeFrac) * (R_Array(iTime, iS2+1) - R_Array(iTime, iS2)) + &
        !          iTimeFrac * (R_Array(iTime+1, iS2+1) - R_Array(iTime+1, iS2))
        !     LagrCoord = iS2 + f1/f2
        ! endif 

    end subroutine get_lagr_coord
    !=============================================================================
    subroutine compute_timestep(Time, LagrCoord, Momentum, TimeStep)
        
        real, intent(in) :: Time, LagrCoord, Momentum
        real, intent(out) :: TimeStep

        real :: S, dS, Dxx, B, dBdS, dDxxdS
        real :: Bup, Bdown, DxxUp, DxxDown

        ! Values at particle location
        call calculate_ds(Time, LagrCoord, dS)
        call get_array_value(Time, LagrCoord, B_Array, B)
        call get_diffusion_coeff(Time, LagrCoord, Momentum, Dxx)

        ! Values upstream of particle
        call get_array_value(Time, LagrCoord + 0.5, B_Array, Bup)
        call get_diffusion_coeff(Time, LagrCoord + 0.5, Momentum, DxxUp)

        ! Values downstream of particle
        call get_array_value(Time, LagrCoord - 0.5, B_Array, Bdown)
        call get_diffusion_coeff(Time, LagrCoord - 0.5, Momentum, DxxDown)
 
        ! Calculate derivatives
        dBdS = (Bup - Bdown) / dS
        dDxxdS = (DxxUp - DxxDown) / dS

        ! (diffusion step size > drift step size) --> (dt << function) --> (dt = TimeStepFactor*function)
        TimeStep = TimeStepFactor * 2.0 * Dxx / &
                   ((0.5 + 0.5*StratoFactor)*dDxxdS - Dxx*dBdS/B)**2

        ! Set maximum allowed timestep - from PARAM
        TimeStep = min(TimeStep, MaxTimeStep)

    end subroutine compute_timestep
    !=============================================================================
    subroutine compute_weight(Time, LagrCoord, Weight)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: Weight
        real :: T, rho

        call get_array_value(Time, LagrCoord, Temp_Array, T)
        call get_array_value(Time, LagrCoord, rho_Array, rho)
    
        Weight = T * rho
        Weight = 1.0
    end subroutine compute_weight
    !=============================================================================
    subroutine get_particle_location(Time, LagrCoord, R)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: R

        call get_array_value(Time, LagrCoord, R_Array, R)

    end subroutine get_particle_location
    !=============================================================================
    subroutine compute_conversion_factor(Time, LagrCoord1, LagrCoord2, ConversionFactor)
        real, intent(in) :: Time, LagrCoord1, LagrCoord2
        real, intent(out) :: ConversionFactor
        real :: dS1, dS2, B1, B2

        call get_array_value(Time, LagrCoord1, B_Array, B1)
        call get_array_value(Time, LagrCoord2, B_Array, B2)

        call calculate_ds(Time, LagrCoord1, dS1)
        call calculate_ds(Time, LagrCoord2, dS2)

        ConversionFactor = 0.5 * (B1/dS1 + B2/dS2)

    end subroutine compute_conversion_factor
    !=============================================================================
end module PT_ModFieldline      
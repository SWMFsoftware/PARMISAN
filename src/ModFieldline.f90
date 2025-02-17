module PT_ModFieldline
    use ModKind
    use ModMpi
    use ModUtilities, ONLY: find_cell 

    use PT_ModConst

    implicit none
    
    SAVE
    ! number of dimensions solved by SDE, (passed to ModParticle and ModSolver !)
    ! should probably not be here, need to put in ModMain, ModInitialize, or something higher up
    integer, parameter :: nDim = 2
    integer, parameter :: LagrCoord_ = 1, Momentum_ = 2

    integer :: nT, nS ! number of time and spatial grid points
    real :: tMin, tMax, rMin, rMax, shockMaxR, rMinInject, rMaxInject ! extreme values
    real :: TimeStepFactor, MaxTimeStep
    real :: StratoFactor ! 1 if Ito, 0 if Stratonovich

    ! Arrays are data read in from field line files: size (nT, nS)
    ! TODO: Change array names to match SWMF naming convention
    real, allocatable :: time_Array(:), S_Array(:,:), rho_Array(:,:), &
                         B_Array(:,:), U_Array(:,:), R_Array(:,:), &
                         shockLocation(:)

    ! lagrangian arrays: change with each time step
    ! indicate the S, U values at LagrCoord at the current time step
    ! TODO: Change array names to match SWMF naming convention
    real, allocatable :: S_lagrangian(:), U_lagrangian(:)
    real :: LagrTime

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

            ! Convert to cm
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
        open(8, file = FileLocation//'time.dat', status = 'old')
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
        open(8, file = FileLocation//'S.dat', status = 'old')
        read(8, *) ((S_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(rho_Array(nT, nS))
        open(8, file = FileLocation//'rho.dat', status = 'old')
        read(8, *) ((rho_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(B_Array(nT, nS))
        open(8, file = FileLocation//'B.dat', status = 'old')
        read(8, *) ((B_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(U_Array(nT, nS))
        open(8, file = FileLocation//'U.dat', status = 'old')
        read(8, *) ((U_Array(row, col), col=1,nS), row=1,nT)
        close(8)

        allocate(R_Array(nT, nS))
        open(8, file = FileLocation//'R.dat', status = 'old')
        read(8, *) ((R_Array(row, col), col=1,nS), row=1,nT)
        close(8)
        
        ! ====================================
        ! convert to correct units [cgs]
        S_Array = S_Array * cRsun  ! [cm]
        U_Array = U_Array * 1.0e5 ! [cm/s]
        R_Array = R_Array * cRsun  ! [cm]
        ! rho = [g/cm**3]
        ! B = [G]

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
    subroutine setup_lagrangian_grid(t0, s0, LagrCoord)
        ! given an initial time and position of particle (t0, s0)
        ! create lagrangian grid at that time and interpolate s,U onto it
        ! return (sL): the lagrangian coordinate of the particle

        real, intent(in) :: t0, s0       ! initial time and position of particle
        real, intent(out) :: LagrCoord   ! initial lagrangian coordinate
        integer :: iTime, iS, i          ! indices 
        real :: iTimeFrac, iSfrac        ! fractional indices
        real :: sMin, sMax, dS           ! variables needed to create s grid

        if(.not.allocated(S_lagrangian)) allocate(S_lagrangian(1:nS))
        if(.not.allocated(U_lagrangian)) allocate(U_lagrangian(1:nS))

        ! t0 is betweeen time_Array[iTime, iTime + 1]
        call find_cell(1, nT, t0, iTime, iTimeFrac, time_Array)

        sMin = max(S_Array(iTime, 1), S_Array(iTime+1, 1))
        sMax = min(S_Array(iTime, nS), S_Array(iTime+1, nS))
      
        dS = (sMax - sMin) / (nS-1)

        do i = 1, nS 
            S_lagrangian(i) = sMin + real(i-1) * dS
            call get_array_value(S_lagrangian(i), t0, U_Array, U_lagrangian(i))
        end do

        ! find particle position along lagrangian grid and store it's lagrangian coordinate in sL
        call find_cell(1, nS, s0, iS, iSfrac, S_lagrangian)
        LagrCoord = real(iS) + iSfrac
        LagrTime = t0

    end subroutine setup_lagrangian_grid
    !============================================================================
    ! subroutine check_invariant(LagrCoord)
    !     real, intent(in) :: LagrCoord
        
    !     integer :: iS, i
    !     real :: Sup, Sdown, dS, rho, B, S, invariant, Time
        
    !     open(18, file= 'PT/IO2/invariant.dt', status = 'unknown')
    !     do i = 1, 153700
    !         Time = real(i)
    !         call get_distance_along_fieldline(LagrCoord, Time, S)
    !         call get_array_value(S, Time, B_Array, B)
    !         call get_array_value(S, Time, rho_Array, rho)

    !         call get_distance_along_fieldline(LagrCoord+0.5, Time, Sup)
    !         call get_distance_along_fieldline(LagrCoord-0.5, Time, Sdown)
    !         dS = Sup - Sdown

    !         invariant = rho*dS/B
    !         write(18,*) LagrCoord, Time, invariant
    !         call update_lagrangian_grid(LagrCoord, 1.0, rho, B)

    !     end do
    !     close(18)
    ! end subroutine check_invariant
    !============================================================================
    subroutine get_distance_along_fieldline(LagrCoord, Time, S)

        use ModUtilities, ONLY: CON_stop
        ! given lagrangian coordinate: compute current distance along fieldline
        real, intent(in) :: LagrCoord, Time
        real, intent(out) :: S

        integer :: iS
        real :: iSfrac, dS, Timestep, SLagrTemp(nS)
        character(len=*), parameter :: NameSub = 'get_distance_along_fieldline'

        ! Uses current velocity to backpropagate if Time < LagrTime
        ! Not exactly correct, so warn if timestep is too large
        Timestep = Time - LagrTime
        if(abs(Timestep).gt.60) write(*,*) Namesub//' Timestep > 60 s'

        iS = floor(LagrCoord)
        iSfrac = LagrCoord - iS
        SLagrTemp = S_lagrangian + U_lagrangian*Timestep

        dS = SLagrTemp(iS+1) - SLagrTemp(iS)
        S = SLagrTemp(iS) + iSfrac * dS

        ! dS = S_lagrangian(iS+1) - S_lagrangian(iS)
        ! S = S_lagrangian(iS) + iSfrac * dS

    end subroutine get_distance_along_fieldline
    !============================================================================
    subroutine calculate_ds(LagrCoord, Time, dS)
        real, intent(in) :: LagrCoord, Time
        real, intent(out) :: dS
        real :: sUp, sDown

        call get_distance_along_fieldline(LagrCoord + 0.5, Time, sUp)
        call get_distance_along_fieldline(LagrCoord - 0.5, Time, sDown)
        dS = sUp - sDown

    end subroutine calculate_ds
    !============================================================================
    subroutine calc_drhodt(LagrCoord, Time, Timestep, dRhodT)
        real, intent(in) :: LagrCoord, Time, Timestep
        real, intent(out) :: dRhodT
        real :: Sfuture, Spast, RhoFuture, RhoPast

        call get_distance_along_fieldline(LagrCoord, Time + Timestep, Sfuture)
        call get_distance_along_fieldline(LagrCoord, Time - Timestep, Spast)

        call get_array_value(Sfuture, Time + Timestep, rho_Array, RhoFuture)
        call get_array_value(Spast, Time - Timestep, rho_Array, RhoPast)

        dRhodT = (RhoFuture - RhoPast) / (2.0 * Timestep)

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
        ! particle location at [LagrCoord, Time]
        ! Time can be different from LagrTime because of RK2 scheme
        call get_distance_along_fieldline(LagrCoord, Time, S)

        ! ==============================
        ! get variables at particle location
        call get_array_value(S, Time, B_Array, B)
        call get_array_value(S, Time, rho_Array, rho)
        call get_diffusion_coeff(S, Time, Momentum, Dxx)

        ! ==============================
        ! Calculate drho/dt - central difference with dt = 1.0 second
        !call calc_drhodt(LagrCoord, Time, 1.0, dRhodT)

        call get_array_value(S, Time-Timestep, rho_array, RhoPast)
        call get_array_value(S, Time+Timestep, rho_Array, RhoFuture)

        dRhodT = (RhoFuture - RhoPast) / (2.0 * Timestep)

        ! ==============================
        ! calculate spatial derivatives at i + 0.5, i - 0.5
        call get_distance_along_fieldline(LagrCoord + 0.5, Time, sUp)
        call get_array_value(sUp, Time, B_Array, Bup)
        call get_diffusion_coeff(sUp, Time, Momentum, DxxUp)

        call get_distance_along_fieldline(LagrCoord - 0.5, Time, sDown)
        call get_array_value(sDown, Time, B_Array, Bdown)
        call get_diffusion_coeff(sDown, Time, Momentum, DxxDown)

        dS = sUp - sDown
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
    subroutine get_sde_coeffs_milstein(X_I, Timestep, Weight, DriftCoeff, DiffCoeff, dDiffdX)
        
        real, intent(in) :: X_I(nDim)
        real, intent(out) :: Timestep, Weight, DriftCoeff(nDim), DiffCoeff(nDim), dDiffdX(nDim)
    
        integer :: iS
        real :: LagrCoord, Momentum
        real :: S, dS, sUp, sDown, dSup, dSdown
        real :: B, Bup, Bdown
        real :: aDerivative, cDerivative 
        real :: Dxx, DxxUp, DxxDown
        real :: rho, dRhodT

        real :: RhoPast, RhoFuture, dRhodT2

        !! TODO: How to sync nDim and indices with ModParticle

        LagrCoord = X_I(LagrCoord_)
        Momentum  = (3.0*X_I(Momentum_))**(1.0/3.0)
        ! ==============================
        ! particle location at [LagrCoord, Time]
        ! Time can be different from LagrTime because of RK2 scheme
        call get_distance_along_fieldline(LagrCoord, LagrTime, S)

        ! ==============================
        ! get variables at particle location
        call get_array_value(S, LagrTime, B_Array, B)
        call get_array_value(S, LagrTime, rho_Array, rho)
        call get_diffusion_coeff(S, LagrTime, Momentum, Dxx)
        call calculate_ds(LagrCoord, LagrTime, dS)

        ! ==============================
        ! Calculate drho/dt - forward difference
        call calc_drhodt(LagrCoord, LagrTime, 1.0, dRhodT)

        
        call get_array_value(S, LagrTime-1.0, rho_array, RhoPast)
        call get_array_value(S, LagrTime+1.0, rho_Array, RhoFuture)

        dRhodT2 = (RhoFuture - RhoPast) / (2.0)
        write(*,*) 'drhodtau, drhodt: ', dRhodT, dRhodT2, dRhodT / dRhodT2

        ! ==============================
        ! calculate quantities at i + 0.5, i - 0.5
        call get_distance_along_fieldline(LagrCoord + 0.5, LagrTime, sUp)
        call calculate_ds(LagrCoord + 0.5, LagrTime, dSup)
        call get_array_value(sUp, LagrTime, B_Array, Bup)
        call get_diffusion_coeff(sUp, LagrTime, Momentum, DxxUp)

        call get_distance_along_fieldline(LagrCoord - 0.5, LagrTime, sDown)
        call calculate_ds(LagrCoord - 0.5, LagrTime, dSdown)
        call get_array_value(sDown, LagrTime, B_Array, Bdown)
        call get_diffusion_coeff(sDown, LagrTime, Momentum, DxxDown)

        aDerivative = DxxUp / (Bup * dSup) - DxxDown / (Bdown * dSdown)
        cDerivative = DxxUp / dSup**2 - DxxDown / dSdown**2
        ! ==============================
        ! calculate sde coefficients 
        DriftCoeff(LagrCoord_) = B * aDerivative / dS
        DriftCoeff(Momentum_) = X_I(Momentum_) * dRhodT / rho

        DiffCoeff(LagrCoord_) = sqrt(2.0 * Dxx) / dS
        DiffCoeff(Momentum_) = 0.0

        dDiffdX(LagrCoord_) = 0.5 * cDerivative
        dDiffdX(Momentum_) = 0.0

        Timestep = TimeStepFactor * (DiffCoeff(LagrCoord_)/DriftCoeff(LagrCoord_))**2
        ! Set maximum allowed timestep, 10 seconds
        TimeStep = min(TimeStep, MaxTimeStep)

        Weight = B / dS

    end subroutine get_sde_coeffs_milstein
    !============================================================================
    subroutine get_diffusion_coeff(S, Time, Momentum, Dxx)
        real, intent(in) :: S, Time, Momentum
        real, intent(out) :: Dxx

        real :: Dxx0, R, Energy
        
        Dxx0 = 5.16d18 ! cm^2 / s
        
        call get_array_value(S, Time, R_Array, R)
        
        Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - cProtonRestEnergy
        
        Dxx = Dxx0 * (R / cAU)**(1.17) * (Energy / ckeV)**(0.71)

    end subroutine get_diffusion_coeff
    !============================================================================
    subroutine update_lagrangian_grid(LagrCoord, Timestep, S, R)
        ! move lagrangian grid from t --> t + dt
        ! also given sL outputs s
        use PT_ModProc, ONLY: iProc, nProc

        real, intent(in) :: Timestep, LagrCoord
        real, intent(out) :: S, R
        integer :: i
    
        ! move grid
        S_lagrangian = S_lagrangian + U_lagrangian*Timestep
        LagrTime = LagrTime + Timestep

        ! calculate new position
        call get_distance_along_fieldline(LagrCoord, LagrTime, S)

        ! update velocity grid
        do i = 1, nS           
            call get_array_value(S_lagrangian(i), LagrTime, U_Array, U_lagrangian(i))
        end do
        
        call get_array_value(S, LagrTime, R_Array, R)


    end subroutine update_lagrangian_grid
    !============================================================================
    subroutine get_array_value(S, Time, Array, InterpValue)
        ! get value of array at (t, s)
        ! uses bilinear interpolation
        real, intent(in) :: S, Time, Array(:,:)
        real, intent(out) :: InterpValue

        integer :: iTime, iS1, iS2          ! (lower) index of desired s, t
        real :: iTimeFrac, iSFrac1, iSFrac2 ! fractional indices
        real :: f1, f2 ! intermediate values, interpolated in s but not t
        logical :: doExtrapolate = .true.
        ! find six indices (2 time, 2 position at each time)
        ! value at (t) lies in time between iTime and iTime+1
        ! value of (s) at iTime lies between iS1 and iS1+1
        ! value of (s) at iTime+1 lies between iS2 and iS2+1
        
        ! doExtrapolate is set to true. may need to revisit this
        ! since the code stops if the particle leaves the domain,
        ! it doesn't matter if grid elements leave the domain
        ! another way would be to shrink the array dynamically as the grid moves out of the domain
        call find_cell(1, nT, Time, iTime, iTimeFrac, time_Array)
        call find_cell(1, nS, S, iS1, iSFrac1, S_Array(iTime, :), doExtrapolate)
        call find_cell(1, nS, S, iS2, iSFrac2, S_Array(iTime+1, :), doExtrapolate)

        f1 = (1-iSFrac1) * Array(iTime, iS1) + iSFrac1 * Array(iTime, iS1+1)
        f2 = (1-iSFrac2) * Array(iTime+1, iS2) + iSFrac2 * Array(iTime+1, iS2+1)

        InterpValue = (1-iTimeFrac) * f1 + iTimeFrac*f2
    
    end subroutine get_array_value
    !============================================================================
    subroutine get_random_shock_location(Time, R, S)
        real, intent(out) :: Time, R, S

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
        R = rMaxInject*rMinInject / (rMaxInject*(1-RandomUniform) + rMinInject*RandomUniform)

        ! Find s-value (distance along field line) that corresponds
        ! to when the shock is at radial distance: r at time: t
        call find_cell(1, nT, R, iTime, iTimeFrac, shockLocation)

        Time = time_Array(iTime) + iTimeFrac * (time_Array(iTime+1) - time_Array(iTime))
        call find_cell(1, nS, R, iS1, iSFrac1, R_Array(iTime, :), doExtrapolate)
        call find_cell(1, nS, R, iS2, iSFrac2, R_Array(iTime+1, :), doExtrapolate)
        
        func1 = (1-iSFrac1) * S_Array(iTime, iS1) + iSFrac1 * S_Array(iTime, iS1+1)
        func2 = (1-iSFrac2) * S_Array(iTime+1, iS2) + iSFrac2 * S_Array(iTime+1, iS2+1)

        S = (1-iTimeFrac) * func1 + iTimeFrac*func2

    end subroutine get_random_shock_location
    !=============================================================================
    subroutine compute_timestep(Time, LagrCoord, Momentum, TimeStep)
        
        real, intent(in) :: Time, LagrCoord, Momentum
        real, intent(out) :: TimeStep

        real :: S, dS, Dxx, B, dBdS, dDxxdS
        real :: Sup, Sdown
        real :: Bup, Bdown, DxxUp, DxxDown

        ! Values at particle location
        call get_distance_along_fieldline(LagrCoord, Time, S)
        call get_array_value(S, Time, B_Array, B)
        call get_diffusion_coeff(S, Time, Momentum, Dxx)

        ! Values upstream of particle
        call get_distance_along_fieldline(LagrCoord + 0.5, Time, Sup)
        call get_array_value(Sup, Time, B_Array, Bup)
        call get_diffusion_coeff(Sup, Time, Momentum, DxxUp)

        ! Values downstream of particle
        call get_distance_along_fieldline(LagrCoord - 0.5, Time, Sdown)
        call get_array_value(Sdown, Time, B_Array, Bdown)
        call get_diffusion_coeff(Sdown, Time, Momentum, DxxDown)

        ! Calculate derivatives
        dS = Sup - Sdown
        dBdS = (Bup - Bdown) / dS
        dDxxdS = (DxxUp - DxxDown) / dS

        ! (diffusion > drift) --> (dt << function) --> (dt = TimeStepFactor*function)
        TimeStep = TimeStepFactor * 2.0 * Dxx / &
                   ((0.5 + 0.5*StratoFactor)*dDxxdS - Dxx*dBdS/B)**2

        ! Set maximum allowed timestep, 10 seconds
        TimeStep = min(TimeStep, MaxTimeStep)

    end subroutine compute_timestep
    !=============================================================================
    subroutine compute_weight(Time, LagrCoord, Weight)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: Weight
        real :: dS, Bfield, S, Sup, Sdown
        ! F = f*dS/B

        call get_distance_along_fieldline(LagrCoord, Time, S)
        call get_distance_along_fieldline(LagrCoord + 0.5, Time, Sup)
        call get_distance_along_fieldline(LagrCoord - 0.5, Time, Sdown)
       
        call get_array_value(S, Time, B_Array, Bfield)
       
        dS = Sup - Sdown
        Weight = Bfield/dS
    end subroutine compute_weight
    !=============================================================================
end module PT_ModFieldline
!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModFieldline

    use PT_ModGrid
    use PT_ModTime, ONLY: PTTime, DataInputTime
    use PT_ModConst, ONLY: cPi, ckeV, cAU, cRsun, cProtonMass, &
                           cElectronCharge, cLightSpeed, cProtonRestEnergy

    use PT_ModSize,  ONLY: nVertexMax
    use PT_ModShock, ONLY: dLogRho_II, dLogRhoOld_II, dLogRhoThreshold
    
    implicit none
    save

    integer :: iLine = 1

    logical :: UseConstantDiffusion = .false.
    real    :: DxxConst = 10.0
    real    :: DxxFactor = 1.0

    real, allocatable :: PreviousState(:, :), CurrentState(:,:), &
                         MhdState1(:, :), MhdState2(:,:)

    integer, parameter :: nStateVar    = 8, nStateAdvect = 7, &
                          dLogRho_     = 1, &
                          BState_      = 2, &
                          dBState_     = 3, &
                          RhoState_    = 4, &
                          UState_      = 5, &
                          dSState_     = 6, &
                          TState_      = 7, &
                          RState_      = 8
                          
    real    :: PreviousTime, CurrentTime
    integer :: iShockNew, iShock1, iShock2
    integer :: iShock1Up, iShock2Up, iShock1Down, iShock2Down
    integer :: WidthUp, WidthDown
    real    :: dLogRhoLimit
contains
!============================================================================
    subroutine read_param_fieldline(NameCommand)

        use ModReadParam, ONLY: read_var
        use ModUtilities, ONLY: CON_stop

        character(len=*), intent(in):: NameCommand ! From PARAM.in
        character(len=*), parameter:: NameSub = 'read_param'
        !--------------------------------------------------------------------------
        select case(NameCommand)
        case('#DIFFUSION')
            call read_var('DxxFactor', DxxFactor)
            call read_var('UseConstantDiffusion', UseConstantDiffusion)
            if(UseConstantDiffusion) &
                call read_var('DxxConst', DxxConst)

        case default
            call CON_stop(NameSub//' Unknown command '//NameCommand)
        end select

    end subroutine read_param_fieldline
    !============================================================================
    subroutine set_fieldline(iLineIn)
        integer, intent(in) :: iLineIn
        integer :: iVertex, iLagr

        iLine = iLineIn

        if(allocated(PreviousState)) deallocate(PreviousState)
        if(allocated(CurrentState)) deallocate(CurrentState)
        if(allocated(MhdState1)) deallocate(MhdState1)
        if(allocated(MhdState2)) deallocate(MhdState2)

        allocate(PreviousState(1:nStateVar, 1:nVertexMax))
        allocate(CurrentState(1:nStateVar, 1:nVertexMax))
        allocate(MhdState1(1:nStateVar, 1:nVertexMax))
        allocate(MhdState2(1:nStateVar, 1:nVertexMax))

        MhdState1(RhoState_, :) = State_VIB(RhoOld_, :, iLine)
        MhdState1(BState_, :) = State_VIB(BOld_, :, iLine)
        MhdState1(dBState_, :) = State_VIB(dBOld_, :, iLine)
        MhdState1(dSState_, :) = State_VIB(DOld_, :, iLine)
        MhdState1(UState_, :) = State_VIB(UOld_, :, iLine)
        MhdState1(RState_, :) = State_VIB(ROld_, :, iLine)
        MhdState1(TState_, :) = State_VIB(TOld_, :, iLine)
        MhdState1(dLogRho_, :) = dLogRhoOld_II(:, iLine)

        MhdState2(RhoState_, :) = MhData_VIB(Rho_, :, iLine)
        MhdState2(BState_, :) = State_VIB(B_, :, iLine)
        MhdState2(dBState_, :) = State_VIB(dB_, :, iLine)
        MhdState2(dSState_, :) = State_VIB(D_, :, iLine)
        MhdState2(UState_, :) = State_VIB(U_, :, iLine)
        MhdState2(RState_, :) = State_VIB(R_, :, iLine)
        MhdState2(TState_, :) = MhData_VIB(T_, :, iLine)
        MhdState2(dLogRho_, :) = dLogRho_II(:, iLine)

        PreviousState = MhdState1
        PreviousTime = PTTime

        CurrentTime = PreviousTime
        CurrentState = PreviousState
    
        iShock1 = iShock_IB(ShockOld_, iLine)
        iShock1Up = iShock_IB(ShockUpOld_, iLine)
        iShock1Down = iShock_IB(ShockDownOld_, iLine)
        
        iShock2 = iShock_IB(Shock_, iLine)
        iShock2Up = iShock_IB(ShockUp_, iLine)
        iShock2Down = iShock_IB(ShockDown_, iLine)

        dLogRhoLimit = min(minval(MhdState2(dLogRho_, :)), &
                           minval(MhdState1(dLogRho_, :)))
    end subroutine
    !============================================================================
    subroutine advect_fieldline(Alpha, iShockNewIn, NewTime)
        real, intent(in) :: Alpha, NewTime
        integer, intent(in) :: iShockNewIn

        integer :: iVar, iLagr
        integer :: iShockNewUp, iShockNewDown

        real :: Slope, Ratio
        real :: dLogRhoMax

        iShockNew = iShockNewIn
        PreviousState = CurrentState
        PreviousTime = CurrentTime
        CurrentState = 0.0
        CurrentTime = NewTime

        ! Interpolate width of shock in time
        WidthUp   = floor((1-Alpha) * (iShock1Up - iShock1) + &
                            Alpha   * (iShock2Up - iShock2))
        WidthDown = floor((1-Alpha) * (iShock1 - iShock1Down) + &
                            Alpha   * (iShock2 - iShock2Down))
        ! check bounds
        WidthDown = min(iShockNew - MinLagr(iLine), WidthDown)
        WidthUp   = min(WidthUp, MaxLagr(iLine) - iShock2)

        ! !  no shock
        if(WidthDown.eq.0.and.WidthUp.eq.0) then
            CurrentState = MhdState2
            return
        end if

        ! New indices of shock extent
        iShockNewUp = iShockNew + WidthUp
        iShockNewDown = iShockNew - WidthDown

        ! Take log of B, dB, and rho before interpolation - helps at inner boundary
        ! where(MhdState1(BState_:RhoState_, :).ne.0) &
        !     MhdState1(BState_:RhoState_, :) = log10(MhdState1(BState_:RhoState_, :))
        ! where(MhdState2(BState_:RhoState_, :).ne.0) &
        !     MhdState2(BState_:RhoState_, :) = log10(MhdState2(BState_:RhoState_, :))
    
        ! interpolate fieldline in time
        CurrentState = (1 - Alpha) * MhdState1 + Alpha * MhdState2

        ! Do not advect shock during time period when shock first appears on fieldline
        ! Undo log and return time interpolated fieldline
        if((iShock1-iShock1Down).eq.0) then 
            where(MhdState1(BState_:RhoState_, :).ne.0) &
                MhdState1(BState_:RhoState_, :) = 10.0**MhdState1(BState_:RhoState_, :)
            where(MhdState2(BState_:RhoState_, :).ne.0) &
                MhdState2(BState_:RhoState_, :) = 10.0**MhdState2(BState_:RhoState_, :)
            where(CurrentState(BState_:RhoState_, :).ne.0) &
                CurrentState(BState_:RhoState_, :) = 10.0**CurrentState(BState_:RhoState_, :)           
            return
        end if

        ! if(iShock1-WidthDown.lt.MinLagr(iLine)) then
        !     ! Case when shock first appears
        !     CurrentState(1:nStateAdvect, iShockNewDown:iShockNewUp) = &
        !         (1-alpha) * MhdState1(1:nStateAdvect, iShockNewDown:iShockNewUp) + &
        !         alpha * MhdState2(1:nStateAdvect, iShock2-WidthDown:iShock2+WidthUp)
        ! else
        
        !Insert advected shock
        CurrentState(1:nStateAdvect, iShockNewDown:iShockNewUp) = &
            (1-alpha) * MhdState1(1:nStateAdvect, iShock1-WidthDown:iShock1+WidthUp) + &
            alpha   * MhdState2(1:nStateAdvect, iShock2-WidthDown:iShock2+WidthUp)

        ! end if

        ! interpolate over old mhd shock regions
        do iVar = 1, nStateAdvect
            if(iShockNewDown.gt.iShock1Down) then
                Slope = (CurrentState(iVar, iShockNewDown)-CurrentState(iVar, iShock1Down)) / &
                    (iShockNewDown - iShock1Down)
                do iLagr = iShock1Down,  iShockNewDown
                    CurrentState(iVar, iLagr) = CurrentState(iVar, iShock1Down) + &
                        Slope * (iLagr - iShock1Down)
                end do
            end if
            if(iShock2Up.gt.(iShockNew+WidthUp+1)) then
                Slope = (CurrentState(iVar, iShock2Up)-CurrentState(iVar, iShockNewUp-1)) / &
                    (iShock2Up - iShockNew-WidthUp+1)
                do iLagr = iShockNewUp-1, iShock2Up
                    CurrentState(iVar, iLagr) = CurrentState(iVar, iShockNewUp-1) + &
                        Slope * (iLagr - iShockNew-WidthUp+1)
                end do
            end if
        end do

        ! Undo log
        ! where(MhdState1(BState_:RhoState_, :).ne.0) &
        !     MhdState1(BState_:RhoState_, :) = 10.0**MhdState1(BState_:RhoState_, :)
        ! where(MhdState2(BState_:RhoState_, :).ne.0) &
        !     MhdState2(BState_:RhoState_, :) = 10.0**MhdState2(BState_:RhoState_, :)
        ! where(CurrentState(BState_:RhoState_, :).ne.0) &
        !     CurrentState(BState_:RhoState_, :) = 10.0**CurrentState(BState_:RhoState_, :)

        ! SHOCK SHARPENING ALGORITHM GOES HERE
        ! Currently increase dLogRho by the maximum dLogRho calculated from the advected rho
        ! ------------------------------------------------------------------- !
        dLogRhoMax = maxval(log(CurrentState(RhoState_, :)/PreviousState(RhoState_, :))/ (CurrentTime-PreviousTime))
        where(CurrentState(dLogRho_, iShockNewDown:iShockNewUp).gt.0) &
            CurrentState(dLogRho_, iShockNewDown:iShockNewUp) = &
                CurrentState(dLogRho_, iShockNewDown:iShockNewUp) * &
                dLogRhoMax / maxval(CurrentState(dLogRho_, iShockNewDown:iShockNewUp))
        ! ! ------------------------------------------------------------------- !

        ! energy loss restriction - arbitrary
        ! dLogRhoLimit = -0.0001
        ! dLogRhoLimit = 0.0
        where(PreviousState(dLogRho_, :).lt.dLogRhoLimit) PreviousState(dLogRho_, :) = dLogRhoLimit
        where(CurrentState(dLogRho_, :).lt.dLogRhoLimit)  CurrentState(dLogRho_, :)  = dLogRhoLimit
        ! ------------------------------------------------------------------- !
    end subroutine advect_fieldline
    !============================================================================
    subroutine interpolate_statevar(Time, LagrCoord, Var, InterpValue)
        real, intent(in) :: Time, LagrCoord
        integer, intent(in) :: Var
        real, intent(out) :: InterpValue
        
        integer :: iLagr
        real :: iLagrFrac, TimeFrac, f1, f2
        
        ! bilinear interpolation in time and space
        TimeFrac = (Time - PreviousTime) / (CurrentTime - PreviousTime)

        ! get integer and fractional parts of lagrangian coordinate
        iLagr = floor(LagrCoord)
        iLagrFrac = LagrCoord - iLagr

        ! if either are less than 1 return -1.0
        ! should only occur after particle timestep and particle location is
        ! calculated from new lagrangian coordinate
        ! boundary condition check occurs immediately after
        if(iLagr.lt.MinLagr(iLine).or.iLagr+1.gt.MaxLagr(iLine)) then
            InterpValue = -1.0    
        else
            ! spatial interpolation at old time
            f1 = (1-iLagrFrac) * PreviousState(Var, iLagr) + & 
                iLagrFrac * PreviousState(Var, iLagr+1)

            ! spatial interplation at next time
            f2 = (1-iLagrFrac) * CurrentState(Var, iLagr) + & 
                iLagrFrac * CurrentState(Var, iLagr+1)
            ! interpolation in time
            InterpValue = (1 - TimeFrac) * f1 + TimeFrac * f2
        end if
    end subroutine interpolate_statevar
    !============================================================================
    subroutine calculate_dLogRho(LagrCoord, dLogRhodTau)
        real, intent(in) :: LagrCoord
        real, intent(out) :: dLogRhodTau

        integer :: iLagr
        real :: iLagrFrac, TimeFrac, f1, f2
        real :: RhoCurrent, RhoPrevious

        iLagr = floor(LagrCoord)
        iLagrFrac = LagrCoord - iLagr

        RhoCurrent = CurrentState(RhoState_, iLagr) * (1-iLagrFrac) + &
                     CurrentState(RhoState_, iLagr+1) * iLagrFrac
        RhoPrevious = PreviousState(RhoState_, iLagr) * (1-iLagrFrac) + &
                     PreviousState(RhoState_, iLagr+1) * iLagrFrac                    
        
        dLogRhodTau = log(RhoCurrent/RhoPrevious) / (CurrentTime - PreviousTime)
        ! dLogRhodTau = max(dLogRhodTau, dLogRhoLimit)
   
    end subroutine calculate_dLogRho
    !============================================================================
    subroutine get_dxx(Time, LagrCoord, Momentum, Dxx)
        real, intent(in) :: Time, LagrCoord, Momentum
        real, intent(out) :: Dxx

        real :: R, B, dB
        real :: R0 = 20.0

        if(UseConstantDiffusion) then
            ! constant Rs**2 is to offset the unit conversion
            ! of dS in the SDE coefficients
            ! only for analytical test!
            Dxx = DxxConst * cRsun**2
            return
        end if

        call interpolate_statevar(Time, LagrCoord, RState_, R)

        ! upstream of shock - use empirical PSP values from Chen et al 2024
        if(R.gt.CurrentState(RState_, iShockNew+WidthUp)) then
            ! call get_psp_dxx(R, Momentum, Dxx)

            call get_upstream_dxx(R, Momentum, Dxx)

        ! downstream of shock - use MHD turbulence
        else
            call interpolate_statevar(Time, LagrCoord, BState_, B)
            call interpolate_statevar(Time, LagrCoord, dBState_, dB)
            call get_mhd_dxx(R, B, dB, Momentum, Dxx)
            
            ! Within shock region, reduce Dxx by factor set in PARAM
            ! Do not let Dxx increase. 
            ! R0 = R value where DxxFactor = 1
            if(R.gt.CurrentState(RState_, iShockNew-WidthDown)) &
                Dxx = Dxx / max(1.0, (DxxFactor / R**(log(DxxFactor)/log(R0))))
        end if

    end subroutine get_dxx
    !============================================================================ 
    subroutine get_mhd_dxx(R, B, dB, Momentum, Dxx)
        real, intent(in) :: R, B, dB, Momentum
        real, intent(out) :: Dxx
        real :: ConstantFactor = 81.0 / (7.0*cPi) * (0.5/cPi)**(2.0/3.0)
        real :: Velocity, Lmax, Btotal, MeanFreePath
        
        ! relativistic velocity)
        Velocity = Momentum / sqrt(cProtonMass**2.0 + &
                   (Momentum/cLightSpeed)**2.0)

        ! Calculate mean free path
        Lmax = 0.4*R*cRsun

        ! dB limiter
        ! if(sqrt(dB)/B.gt.1) B = sqrt(dB)

        Btotal = sqrt(B**2 + dB)

        MeanFreePath = ConstantFactor * Btotal**2 * &
                        (Momentum*Lmax**2/(Btotal * cElectronCharge))**(1.0/3.0) / dB
        ! Calculate Dxx 
        Dxx = MeanFreePath * Velocity / 3.0
        Dxx = max(Dxx, 1.0d4 * cRsun)

    end subroutine get_mhd_dxx
    !============================================================================ 
    subroutine get_psp_dxx(R, Momentum, Dxx)
        real, intent(in) :: R, Momentum
        real, intent(out) :: Dxx

        real :: Dxx0, Energy
        ! Chen et al., 2024
        Dxx0 = 5.16d14 ! [m^2/s]
        Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - cProtonRestEnergy

        Dxx = Dxx0 * (R * cRsun / cAU )**(1.17) * (Energy / ckeV)**(0.71)

    end subroutine get_psp_dxx
    !============================================================================   
    subroutine get_upstream_dxx(R, Momentum, Dxx)
        real, intent(in) :: R, Momentum
        real, intent(out) :: Dxx

        real :: Lambda0 = 0.3
        real :: Energy, fac1, fac2, GeV

        Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - cProtonRestEnergy
        GeV = 1e6 * ckeV

        fac1 = Energy * (Energy + 2.0 * cProtonRestEnergy) / GeV**2.0
        fac2 = Energy * (Energy + 2.0 * cProtonRestEnergy) / (Energy + cProtonRestEnergy)**2.0

        Dxx = Lambda0 * cLightSpeed * R * cRsun * fac1 ** (1.0/6.0) * sqrt(fac2) / (3.0)
    
    end subroutine get_upstream_dxx
    !============================================================================   
    subroutine get_sde_coeffs_euler(X_I, Time, TimeStep, DriftCoeff, DiffCoeff)
        use PT_ModSize, only: nDim
        use PT_ModProc, only: iProc

        real, intent(in) :: X_I(nDim), Time
        real, intent(out) :: Timestep, DriftCoeff(nDim), DiffCoeff(nDim)

        real :: Momentum, B, Dxx, dS, ShockWidth
        real :: Bup, DxxUp, dSup
        real :: Bdown, DxxDown, dSdown
        real :: Rho1, RhoOld1, dLogRhodTau

        integer :: iVertex
        ! Need to figure out where to put X_I index variables - hardcoded 2 = Momentum_
        Momentum = (3.0*X_I(2))**(1.0/3.0)

        ! Need to figure out where to put X_I index variables - hardcoded 1 = LagrCoord_
        call calculate_dLogRho(X_I(1), dLogRhodTau)
        ! call interpolate_statevar(Time, X_I(1), dLogRho_, dLogRhodTau)

        ! get values at particle current location
        call interpolate_statevar(Time, X_I(1), BState_, B)
        call interpolate_statevar(Time, X_I(1), dSState_, dS)
        call get_dxx(Time, X_I(1), Momentum, Dxx)
        
        ! get values one lagrangian coordinate upstream of particle location
        call interpolate_statevar(Time, X_I(1) + 1.0, BState_, Bup)
        call interpolate_statevar(Time, X_I(1) + 1.0, dSState_, dSup)
        call get_dxx(Time, X_I(1) + 1.0, Momentum, DxxUp)

        ! convert from [Rsun] to [m]
        dS = dS * cRsun
        dSup = dSup * cRsun

        ! lagr coord sde coefficients
        DriftCoeff(1) = B * ((DxxUp / (Bup * dSup)) - (Dxx / (B * dS))) / dS
        DiffCoeff(1) = sqrt(2.0 * Dxx) / dS

        ! momentum sde coefficients
        DriftCoeff(2) = X_I(2) * dLogRhodTau
        DiffCoeff(2) = 0.0
                
        ! pri
        ! calculate timestep based on coefficients
        ! diffusion >> drift
        ! Maximum spatial step size is less than shockwidth
        ShockWidth = real(WidthUp + WidthDown)
        if(ShockWidth.eq.0) ShockWidth = 1e9
        if(abs(X_I(1) - real(iShockNew)).gt.ShockWidth) then
            Timestep = DiffCoeff(1)**2.0 / DriftCoeff(1)**2.0
        else
            Timestep = min(DiffCoeff(1)**2.0 / DriftCoeff(1)**2.0, &
                           ShockWidth / abs(DriftCoeff(1)),        &
                           ShockWidth**2.0 / DiffCoeff(1)**2.0)
        end if

    end subroutine get_sde_coeffs_euler
    !============================================================================
    subroutine calc_thermal_energy(Time, LagrCoord, ThermalEnergy)
        ! Calculate thermal energy density: kb*mp*np*Tp
        ! Temp [keV] Rho [amu/m^3]
        ! Output: ThermalEnergy [J kg / m^3]
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: ThermalEnergy

        real :: Temp, Rho

        call interpolate_statevar(Time, LagrCoord, TState_, Temp)
        call interpolate_statevar(Time, LagrCoord, RhoState_, Rho)

        Temp = Temp * ckeV
        ThermalEnergy = Temp * Rho * cProtonMass

    end subroutine calc_thermal_energy
    !============================================================================
    subroutine calc_weight(Time, LagrCoord, Weight)
        ! Statistical weight of each particle: rho*T*dS/B
        real, intent(in)  :: Time, LagrCoord
        real, intent(out) :: Weight

        real :: ThermalEnergy, dSOverB

        call calc_thermal_energy(Time, LagrCoord, ThermalEnergy)
        call calc_dSOverB(Time, LagrCoord, dSOverB)

        Weight = ThermalEnergy * dSOverB
    end subroutine calc_weight
    !============================================================================
    subroutine get_particle_location(Time, LagrCoord, R)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: R
        integer :: iLagr, iVertex, iVertexOld
        real :: iLagrFrac, TimeFrac

        call interpolate_statevar(Time, LagrCoord, RState_, R)

    end subroutine get_particle_location
    !============================================================================
    subroutine calc_dSOverB(Time, LagrCoord, dSOverB)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: dSOverB
        
        real :: dS, B

        call interpolate_statevar(Time, LagrCoord, dSState_, dS)
        call interpolate_statevar(Time, LagrCoord, BState_, B)
        dSOverB = dS * cRsun / B

    end subroutine calc_dSOverB
    !============================================================================
    subroutine check_boundary_conditions(Time, LagrCoord, R, IsOutside)
        real, intent(in) :: Time
        real, intent(inout) :: LagrCoord, R
        logical, intent(out) :: IsOutside 
        real :: Rtemp
        IsOutside = .false.
        
        ! reflecting boundary condition - inner spatial boundary
        if(LagrCoord.lt.MinLagr(iLine)) then
             LagrCoord = MinLagr(iLine) + 1.0
             call interpolate_statevar(Time, LagrCoord, RState_, R)
        end if

        ! absorbing boundary conditions - outer spatial boundary
        if(LagrCoord.lt.MinLagr(iLine)) IsOutside = .true.
        ! absorbing boundary conditions - outer spatial boundary
        if(LagrCoord.gt.MaxLagr(iLine)) IsOutside = .true.
        
    end subroutine check_boundary_conditions
    !============================================================================
    subroutine save_fieldline_data(iProgress, NextTimeStep)
        integer, intent(in) :: iProgress
        real, intent(in) :: NextTimeStep
        character(len = 20) :: OutString
        integer :: iVar

        write(OutString, '(I6.6, A, I3.3)') int(DataInputTime), '_', iProgress
        OutString = adjustl(OutString)

        ! Current shock location, shock location at end of coupling time,
        ! shock location at previous coupling time
        open(102, file = 'PT/IO2/'//'shockloc_'//trim(OutString)//'.dat')
        write(102, *) PTTime, &
                    DataInputTime, PreviousTime, CurrentTime, &
                    iShockNew, iShock1, iShock2, MinLagr(iLine), MaxLagr(iLine), &
                    WidthUp, WidthDown, iShock1Down, &
                    iShock1Up, iShock2Down, &
                    iShock2Up
        close(102)
        ! ------------------------------------------------------------------- 
        open(101, file = 'PT/IO2/'//'shockdata_'//trim(OutString)//'.dat')
        
        do iVar = 1, nStateVar
            write(101, '(11000e15.6)') MhdState1(iVar, MinLagr(iLine):MaxLagr(iLine))
            write(101, '(11000e15.6)') MhdState2(iVar, MinLagr(iLine):MaxLagr(iLine))
            write(101, '(11000e15.6)') PreviousState(iVar, MinLagr(iLine):MaxLagr(iLine))
            write(101, '(11000e15.6)') CurrentState(iVar, MinLagr(iLine):MaxLagr(iLine))
        end do
        close(101)

    end subroutine save_fieldline_data
    !============================================================================
end module PT_ModFieldline
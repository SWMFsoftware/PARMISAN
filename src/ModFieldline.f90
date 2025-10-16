!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModFieldline

    use PT_ModGrid
    use PT_ModTime, ONLY: PTTime, DataInputTime
    use PT_ModConst, ONLY: cPi, ckeV, cAU, cRsun, cProtonMass, &
                           cElectronCharge, cLightSpeed, cProtonRestEnergy
    use PT_ModSize,  ONLY: nVertexMax
    
    implicit none
    save

    integer :: iLine = 1

    logical :: UseConstantDiffusion = .false.
    real    :: DxxConst = 10.0

    real, allocatable :: PreviousState(:, :), CurrentState(:,:), &
                         MhdState1(:, :), MhdState2(:,:)

    integer, parameter :: nStateVar    = 6, nStateAdvect = 5, &
                          RhoState_    = 1, &
                          BState_      = 2, &
                          dBState_     = 3, &
                          dSState_     = 4, &
                          UState_      = 5, &
                          RState_      = 6
                          
    real    :: PreviousTime, CurrentTime
    integer :: iShock1, iShock2
    integer :: WidthUp, WidthDown

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
            call read_var('UseConstantDiffusion', UseConstantDiffusion)
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

        MhdState2(RhoState_, :) = MhData_VIB(Rho_, :, iLine)
        MhdState2(BState_, :) = State_VIB(B_, :, iLine)
        MhdState2(dBState_, :) = State_VIB(dB_, :, iLine)
        MhdState2(dSState_, :) = State_VIB(D_, :, iLine)
        MhdState2(UState_, :) = State_VIB(U_, :, iLine)
        MhdState2(RState_, :) = State_VIB(R_, :, iLine)

        PreviousState = MhdState1
        PreviousTime = PTTime

        CurrentTime = PreviousTime
        CurrentState = PreviousState
    
        iShock1 = iShock_IB(ShockOld_, iLine)
        iShock2 = iShock_IB(Shock_, iLine)

    end subroutine
    !============================================================================
    subroutine advect_fieldline(Alpha, iShockNew, NewTime)
        real, intent(in) :: Alpha, NewTime
        integer, intent(in) :: iShockNew

        integer :: iVar, iLagr
        integer :: iShockNewUp, iShockNewDown

        real :: Slope, Ratio

        PreviousState = CurrentState
        PreviousTime = CurrentTime
        CurrentState = 0.0
        CurrentTime = NewTime

        ! Interpolate widths in time
        WidthUp   = floor((1-Alpha) * (iShock_IB(ShockUpOld_, iLine) - iShock1) + &
                            Alpha   * (iShock_IB(ShockUp_, iLine) - iShock2))
        WidthDown = floor((1-Alpha) * (iShock1 - iShock_IB(ShockDownOld_, iLine)) + &
                            Alpha   * (iShock2 - iShock_IB(ShockDown_, iLine)))
        ! check bounds
        WidthDown = min(iShock1 - MinLagr(iLine), WidthDown)
        WidthUp   = min(WidthUp, MaxLagr(iLine) - iShock2)

        ! write(*,*) iShock_IB(ShockDownOld_, iLine), iShock1, iShock_IB(ShockUpOld_, iLine)
        ! write(*,*) iShock_IB(ShockDown_, iLine), iShock2, iShock_IB(ShockUp_, iLine)
        ! write(*,*) iShockNew, WidthDown, WidthUp
        ! write(*,*) '# --------------------------------------- #'
    
        ! interpolate fieldline in time
        CurrentState = (1 - Alpha) * MhdState1 + Alpha * MhdState2

        !Insert advected shock
        CurrentState(1:nStateAdvect, iShockNew-WidthDown:iShockNew+WidthUp) = &
            (1-alpha) * MhdState1(1:nStateAdvect, iShock1-WidthDown:iShock1+WidthUp) + &
              alpha   * MhdState2(1:nStateAdvect, iShock2-WidthDown:iShock2+WidthUp)

        ! Need to smooth over previous shock regions and ensure smooth interpolation
        ! with inserted shock!

        ! do iVar = 1, nStateAdvect
            
        !     Ratio = CurrentState(iVar, iShockNew-WidthDown-1) / &
        !             CurrentState(iVar, iShockNew-WidthDown) 
        !     CurrentState(iVar, iShockNew-WidthDown:iShockNew+WidthUp) = &
        !         CurrentState(iVar, iShockNew-WidthDown:iShockNew+WidthUp) * Ratio
            ! Ratio =  CurrentState(iVar,iShockNew+WidthUp+1) / &
            !          CurrentState(iVar, iShockNew+WidthUp) 
            ! CurrentState(iVar, iShockNew:iShockNew+WidthUp) = &
            !     CurrentState(iVar, iShockNew:iShockNew+WidthUp) * Ratio
            ! do iLagr = iShockNew-2, iShockNew+2
            !     CurrentState(iVar, iLagr) = sum(CurrentState(iVar,iLagr-1:iLagr+1)) / 3.0
            ! end do

        ! end do

        ! SHOCK SHARPENING ALGORITHM GOES HERE
    
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
        ! being calculated from new lagrangian coordinate
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
        
    end subroutine calculate_dLogRho
    !============================================================================
    subroutine get_dxx(Time, LagrCoord, Momentum, Dxx)
        real, intent(in) :: Time, LagrCoord, Momentum
        real, intent(out) :: Dxx

        real :: R, ShockR, B, dB

        if(UseConstantDiffusion) then
            ! constant Rs**2 is to offset the unit conversion
            ! of dS in the SDE coefficients
            ! only for analytical test!
            Dxx = DxxConst * cRsun**2
            return
        end if

        call interpolate_statevar(Time, LagrCoord, RState_, R)
        call interpolate_statevar(Time, LagrCoord, BState_, B)
        call interpolate_statevar(Time, LagrCoord, dBState_, dB)
        call get_mhd_dxx(R, B, dB, Momentum, Dxx)

        ! not quite right - shock location not moving during sharpening...
        ! ShockR = CurrentState(RState_, iShock2)
        ! ShockR + 0.5?
        ! if(R.gt.(ShockR+0.5)) then
        !     call get_psp_dxx(R, Momentum, Dxx)
        ! end if

    end subroutine get_dxx
    !============================================================================ 
    subroutine get_mhd_dxx(R, B, dB, Momentum, Dxx)
        real, intent(in) :: R, B, dB, Momentum
        real, intent(out) :: Dxx
        real :: ConstantFactor = 81.0 / (7.0*cPi) * (0.5/cPi)**(2.0/3.0)
        real :: Velocity, Lmax, Btotal, MeanFreePath
        
        ! relativistic velocity
        Velocity = Momentum / cProtonMass
        Velocity = sqrt(Velocity**2 / (1 + (Velocity/cLightSpeed)**2))

        ! Calculate mean free path
        ! TODO: add maximum dB/B - shouldn't be more than 1
        Lmax = 0.4*R*cRsun
        Btotal = sqrt(B**2 + dB)

        ! if((Btotal/sqrt(dB)).lt.1.0) Btotal = sqrt(dB)
        ! if((Btotal**2/dB).gt.1) write(*,*) 'B/dB > 1'
        MeanFreePath = ConstantFactor * Btotal**2 * &
                        (Momentum*Lmax**2/(B * cElectronCharge))**(1.0/3.0) / dB

        ! Calculate Dxx 
        Dxx = MeanFreePath * Velocity / 3.0
        Dxx = max(Dxx, 1.0d4 * cRsun)

    end subroutine get_mhd_dxx
    !============================================================================ 
    subroutine get_psp_dxx(R, Momentum, Dxx)
        real, intent(in) :: R, Momentum
        real, intent(out) :: Dxx

        real :: Dxx0, Energy

        Dxx0 = 5.16d14 ! [m^2/s]
        Energy = sqrt((Momentum*cLightSpeed)**2 + cProtonRestEnergy**2) - cProtonRestEnergy

        Dxx = Dxx0 * (R * cRsun / cAU )**(1.17) * (Energy / ckeV)**(0.71)

    end subroutine get_psp_dxx
    !============================================================================   
    subroutine get_sde_coeffs_euler(X_I, Time, TimeStep, DriftCoeff, DiffCoeff)
        use PT_ModSize, only: nDim

        real, intent(in) :: X_I(nDim), Time
        real, intent(out) :: Timestep, DriftCoeff(nDim), DiffCoeff(nDim)

        real :: Momentum, B, Dxx, dS
        real :: Bup, DxxUp, dSup
        real :: Bdown, DxxDown, dSdown
        real :: Rho1, RhoOld1, dLogRhodTau

        integer :: iVertex
        ! Need to figure out where to put X_I index variables - hardcoded 2 = Momentum_
        Momentum = (3.0*X_I(2))**(1.0/3.0)

        ! Need to figure out where to put X_I index variables - hardcoded 1 = LagrCoord_
        call calculate_dLogRho(X_I(1), dLogRhodTau)
        
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

        ! calculate timestep based on coefficients
        ! diffusion >> drift
        Timestep = DiffCoeff(1)**2.0 / DriftCoeff(1)**2.0
      
    end subroutine get_sde_coeffs_euler
    !============================================================================
    subroutine compute_weight(Time, LagrCoord, Weight)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: Weight

        real :: Temp, Rho
        ! will need to change to Weight = T*rho
        Weight = 1.0
    end subroutine compute_weight
    !============================================================================
    subroutine get_particle_location(Time, LagrCoord, R)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: R
        integer :: iLagr, iVertex, iVertexOld
        real :: iLagrFrac, TimeFrac

        call interpolate_statevar(Time, LagrCoord, RState_, R)

    end subroutine get_particle_location
    !============================================================================
    subroutine compute_conversion_factor(Time, LagrCoord, ConversionFactor)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: ConversionFactor
        
        real :: B, dS

        call interpolate_statevar(Time, LagrCoord, dSState_, dS)
        call interpolate_statevar(Time, LagrCoord, BState_, B)
        ! use B or Bsteepen????
        ConversionFactor = B / dS
    end subroutine compute_conversion_factor
    !============================================================================
    subroutine check_boundary_conditions(Time, LagrCoord, R, IsOutside)
        real, intent(in) :: Time
        real, intent(inout) :: LagrCoord, R
        logical, intent(out) :: IsOutside 

        IsOutside = .false.

        ! reflecting boundary condition - inner spatial boundary
        ! if(LagrCoord.lt.MinLagr(iLine)) then
        !     LagrCoord = MinLagr(iLine) + 1.0
        !     call interpolate_statevar(Time, LagrCoord, RState_, R)
        ! end if

        ! absorbing boundary conditions - outer spatial boundary
        if(LagrCoord.lt.MinLagr(iLine)) IsOutside = .true.
        ! absorbing boundary conditions - outer spatial boundary
        if(LagrCoord.gt.MaxLagr(iLine)) IsOutside = .true.
      
    end subroutine check_boundary_conditions
    !============================================================================
    subroutine save_fieldline_data(iShockNew, NextTimeStep)
        integer, intent(in) :: iShockNew
        real, intent(in) :: NextTimeStep
        character(len = 20) :: OutString

        write(OutString, *) int(NextTimeStep)
        OutString = adjustl(OutString)

        ! Current shock location, shock location at end of coupling time,
        ! shock location at previous coupling time
        open(102, file = 'PT/IO2/'//'shockloc_'//trim(OutString)//'.dat')
        write(102, *) PTTime, &
                    DataInputTime, PreviousTime, CurrentTime, &
                    iShockNew, iShock1, iShock2, MinLagr(iLine), MaxLagr(iLine), &
                    WidthUp, WidthDown, iShock_IB(ShockDownOld_, iLine), &
                    iShock_IB(ShockUpOld_, iLine), iShock_IB(ShockDown_, iLine), &
                    iShock_IB(ShockUp_, iLine)
        close(102)
        ! ------------------------------------------------------------------- 
        open(101, file = 'PT/IO2/'//'shockdata_'//trim(OutString)//'.dat')
        
        write(101, '(7000e15.6)') MhdState1(RhoState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') MhdState2(RhoState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') PreviousState(RhoState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') CurrentState(RhoState_, MinLagr(iLine):MaxLagr(iLine))
        
        write(101, '(7000e15.6)') MhdState1(dSState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') MhdState2(dSState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') PreviousState(dSState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') CurrentState(dSState_, MinLagr(iLine):MaxLagr(iLine))

        write(101, '(7000e15.6)') MhdState1(dBState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') MhdState2(dBState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') PreviousState(dBState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') CurrentState(dBState_, MinLagr(iLine):MaxLagr(iLine))

        write(101, '(7000e15.6)') MhdState1(RState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') MhdState2(RState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') PreviousState(RState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') CurrentState(RState_, MinLagr(iLine):MaxLagr(iLine))

        write(101, '(7000e15.6)') MhdState1(UState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') MhdState2(UState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') PreviousState(UState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') CurrentState(UState_, MinLagr(iLine):MaxLagr(iLine))

        write(101, '(7000e15.6)') MhdState1(BState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') MhdState2(BState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') PreviousState(BState_, MinLagr(iLine):MaxLagr(iLine))
        write(101, '(7000e15.6)') CurrentState(BState_, MinLagr(iLine):MaxLagr(iLine))

        close(101)

    end subroutine save_fieldline_data
    !============================================================================
end module PT_ModFieldline
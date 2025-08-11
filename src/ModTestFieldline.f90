module PT_ModTestFieldline
    ! Written by Alex Shane
    ! Numerical shock test fieldline
    ! Simulates 1D shock with constant upstream/downstream profiles
    ! Magnetic field assumed to be constant
    ! Shock moves at unit speed, upstream flow is stationary, downstream flow can be varied

    use ModKind
    use ModMpi
    use PT_ModConst
    use PT_ModProc, only : iProc

    implicit none

    SAVE

    character(len=*), parameter :: OutputDir = 'PT/IO2/'
    real, parameter :: tMin = 0.0, tMax = 3600.0    ! seconds
    real, parameter :: rMin = -100.0, rMax = 3700.0 ! lagrcoord bounds (r = lagrcoord)
    real, parameter :: Dxx0 = 10.0  ! either constant or Dxx at injection momentum
 
    integer, parameter :: nS = 10000 ! number of indices of fieldline (stand-in because ModParticle needs this)
    integer, parameter :: nDim = 2 ! number of dimensions of solution vector
    integer, parameter :: LagrCoord_ = 1, Momentum_ = 2 ! indices of solution vector
    real, parameter :: dShock = 4.0 ! how far in sL to calculate spacing between grid points
                                    ! this effectively changes the width of the shock = 2*dShock + 1
                                    ! 1 = acceleration region (constant)
    real, parameter :: velDownstream = 0.75 ! downstream speed (U2) 1-U2 = 1 / compressionRatio
                                            ! assumes upstream flow is 0 or -1 in shock frame

contains
    !=====================================================================!
    subroutine read_fieldline
        ! No data needs to read in for test
        ! This will output generic Pnorm array, steady-state distribution function,
        ! and the theoretical acceleration time.
        ! This is only true for Dxx = constant!
        integer :: i
        integer :: NumPNorm= 200

        real :: PnormMin = 0.0, PnormMax = 4.0 ! log10
        real :: dPNorm, PowerLaw, U1, U2
        real, allocatable :: Pnorm(:), f(:), Tacc(:)


        if(iProc.eq.0) then
            write(*,*) "Running numerical test"
            
            allocate(Pnorm(1:NumPNorm))
            allocate(f(1:NumPNorm))
            allocate(Tacc(1:NumPNorm))
            
            U1 = -1.0 ! upstream speed in shock frame
            U2 = velDownstream - 1.0 ! downstream speed in shock frame
            PowerLaw = -3.0 * U1 / (U1 - U2)
            dPNorm = (PnormMax - PnormMin) / (NumPNorm-1)
            do i = 1, NumPNorm
                Pnorm(i) = 10**(PnormMin + dPNorm * (i-1) )
                f(i) = Pnorm(i) ** PowerLaw
                Tacc(i) = 3.0 * Dxx0 * (U1**-1.0 + U2**-1.0) / (U1 - U2) * log(Pnorm(i))
            end do

            open(801,file=OutputDir//'pnorm.dat',status='unknown')
            write(801,'(1000e15.6)') Pnorm
            close(801)
            open(801,file=OutputDir//'steadystate_f.dat',status='unknown')
            write(801,'(1000e15.6)') f
            close(801)
            open(801,file=OutputDir//'acceleration_time.dat',status='unknown')
            write(801,'(1000e15.6)') Tacc
            close(801)

        end if
    end subroutine read_fieldline
    !=====================================================================!
    subroutine get_particle_location(Time, LagrCoord, S)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: S
        real :: ShockTime, DeltaT

        call get_shock_arrival(LagrCoord, ShockTime)
        DeltaT = Time - ShockTime

        ! From Sokolov 2023
        ! Modified to allow compression ratio variation
        if(DeltaT < 0) then 
            S = LagrCoord 
        else if(DeltaT > 1) then
            S = LagrCoord + 0.5*velDownstream + velDownstream*(DeltaT - 1)
        else
            S = LagrCoord + 0.5 * velDownstream * DeltaT**2
        end if
        
    end subroutine get_particle_location
    !=====================================================================!
    subroutine get_lagr_coord(Time, R, LagrCoord)
        real, intent(in) :: Time, R
        real, intent(out) :: LagrCoord
        ! R = Lagr
        LagrCoord = R / cRsun
    end subroutine get_lagr_coord
    !=====================================================================!
    subroutine get_ds(Time, LagrCoord, dS)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: dS
  
        real :: Sup, Sdown, ShockTime, DeltaT
  
        call get_particle_location(Time, LagrCoord - dShock, Sdown)
        call get_particle_location(Time, LagrCoord + dShock, Sup)
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
    subroutine get_dxx(Time, Momentum, Dxx)
        real, intent(in) :: Time, Momentum
        real, intent(out) :: Dxx
        real :: injMomentum, injEnergy
        
        ! injEnergy = 1d-15
        ! injMomentum = sqrt(2.0 * cProtonMass * injEnergy + (injEnergy/cLightSpeed)**2.0)
        ! if(Time.gt.0.0) then
        !     Dxx = Dxx0 * (Momentum / injMomentum)**(1.3333)
        ! else
        !     Dxx = Dxx0
        ! end if
        Dxx = Dxx0

    end subroutine get_dxx
    !=====================================================================!
    subroutine get_shock_arrival(LagrCoord, ShockTime)
        real, intent(in) :: LagrCoord
        real, intent(out) :: ShockTime
  
        ShockTime = LagrCoord + dShock
  
    end subroutine get_shock_arrival
    !=====================================================================!
    subroutine get_sde_coeffs_euler(X_I, Time, Timestep, DriftCoeff, DiffCoeff)
        real, intent(in) :: X_I(nDim), Time
        real, intent(out) :: Timestep, DriftCoeff(nDim), DiffCoeff(nDim)

        real :: Momentum, Dxx, dS, dSup, dSdown
        real :: rho, rhoFuture, rhoPast, dRhodTau

        Momentum  = (3.0*X_I(Momentum_))**(1.0/3.0)
        call get_dxx(Time, Momentum, Dxx)

        call get_ds(Time, X_I(LagrCoord_), dS)
        rho = 1.0 / dS

        ! found this converges around dt = 0.001
        call get_rho(Time + 0.001, X_I(LagrCoord_), rhoFuture)
        call get_rho(Time - 0.001, X_I(LagrCoord_), rhoPast)
        dRhodTau = (rhoFuture - rhoPast) / 0.002
        
        call get_ds(Time, X_I(LagrCoord_) + 0.5, dSup)
        call get_ds(Time, X_I(LagrCoord_) - 0.5, dSdown)

        DriftCoeff(LagrCoord_) = Dxx/dS * (dSup**(-1.0) - dSdown**(-1.0))
        DiffCoeff(LagrCoord_) = sqrt(2.0 * Dxx) / dS

        DriftCoeff(Momentum_) = X_I(Momentum_) * dRhodTau / rho
        DiffCoeff(Momentum_) = 0.0

        Timestep = 0.01 * DiffCoeff(LagrCoord_) / DriftCoeff(LagrCoord_)**2.0
        Timestep = min(Timestep, 0.01)

    end subroutine get_sde_coeffs_euler   
    !=====================================================================!
    subroutine get_sde_coeffs_milstein(X_I, Time, Timestep, DriftCoeff, DiffCoeff, dDiffdX)
        real, intent(in) :: X_I(nDim), Time
        real, intent(out) :: Timestep, DriftCoeff(nDim), DiffCoeff(nDim), dDiffdX(nDim)

        real :: Momentum, Dxx, dS, dSup, dSdown
        real :: rho, rhoFuture, rhoPast, dRhodTau

        Momentum  = (3.0*X_I(Momentum_))**(1.0/3.0)
        call get_dxx(Time, Momentum, Dxx)

        call get_ds(Time, X_I(LagrCoord_), dS)
        rho = 1.0 / dS

        ! found this converges around dt = 0.001
        call get_rho(Time + 0.001, X_I(LagrCoord_), rhoFuture)
        call get_rho(Time - 0.001, X_I(LagrCoord_), rhoPast)
        dRhodTau = (rhoFuture - rhoPast) / 0.002
        
        call get_ds(Time, X_I(LagrCoord_) + 0.5, dSup)
        call get_ds(Time, X_I(LagrCoord_) - 0.5, dSdown)

        DriftCoeff(LagrCoord_) = Dxx/dS * (dSup**(-1.0) - dSdown**(-1.0))
        DiffCoeff(LagrCoord_) = sqrt(2.0 * Dxx) / dS
        dDiffdX(LagrCoord_) = Dxx * (dSup**(-2.0) - dSdown**(-2.0))

        DriftCoeff(Momentum_) = X_I(Momentum_) * dRhodTau / rho
        DiffCoeff(Momentum_) = 0.0
        dDiffdX(Momentum_) = 0.0

        Timestep = 0.01 * DiffCoeff(LagrCoord_) / DriftCoeff(LagrCoord_)**2.0
        Timestep = min(Timestep, 0.01)

    end subroutine get_sde_coeffs_milstein
    !=====================================================================!
    subroutine get_random_shock_location(Time, LagrCoord, S)
        real, intent(out) :: Time, LagrCoord, S
        real :: RandomUniform
        
        call random_number(RandomUniform)
  
        LagrCoord = tMin + (tMax - tMin) * RandomUniform
        ! call get_shock_arrival(LagrCoord, Time)
        Time = LagrCoord - 2.0 * dShock ! inject particles upstream
        call get_particle_location(Time, LagrCoord, S)
  
    end subroutine get_random_shock_location
    !=====================================================================!
    subroutine compute_weight(Time, LagrCoord, Weight)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: Weight
        real :: dS

        Weight = 1.0

    end subroutine compute_weight
    !=====================================================================! 
    subroutine compute_conversion_factor(Time, LagrCoord, ConversionFactor)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: ConversionFactor
        real :: dS
        call get_ds(Time, LagrCoord, dS)
        ConversionFactor = 1.0 / dS
    end subroutine compute_conversion_factor
    !=====================================================================!
end module PT_ModTestFieldline
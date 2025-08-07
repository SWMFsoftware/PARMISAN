module PT_ModTestFieldline
    use ModKind
    use ModMpi
    use PT_ModConst

    implicit none

    SAVE

    real, parameter :: tMin = 0.0, tMax = 3600.0
    real, parameter :: rMin = -100.0, rMax = 3700.0
    real, parameter :: Dxx0 = 10.0
 
    integer, parameter :: nS = 10000
    integer, parameter :: nDim = 2
    integer, parameter :: LagrCoord_ = 1, Momentum_ = 2
    real, parameter :: StratoFactor = 0.0 ! 1 if Ito, 0 if Stratonovich
    real, parameter :: dShock = 4.0 ! originally 0.5. how far in sL to calculate spacing between grid points
                                    ! this effectively changes the width of the shock = 2*dShock + 1
                                    ! 1 = acceleration region (constant)
    real, parameter :: velDownstream = 0.75 ! downstream speed (U2) 1-U2 = 1 / compressionRatio

contains
    !=====================================================================!
    subroutine read_fieldline
        write(*,*) "Running numerical test"
    end subroutine read_fieldline
    !=====================================================================!
    subroutine get_particle_location(Time, LagrCoord, S)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: S
        real :: ShockTime, DeltaT

        call get_shock_arrival(LagrCoord, ShockTime)
        DeltaT = Time - ShockTime

        ! Instantaneous acceleration
        ! Change dSL to 1d-10 for true discontinuity
        ! --------------------------------
        ! if(DeltaT < 0) then                      ! stationary before shock arrival
        !    S = LagrCoord - 0.5
        ! else
        !    S = LagrCoord - 0.5 + 0.75 * DeltaT
        !    ! S = LagrCoord - 0.5 + 0.5 * DeltaT   ! constant speed after shock arrival
        ! end if

        ! Igor's method
        ! --------------------------------
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
    subroutine get_sde_coeffs_rk2(X_I, Time, Timestep, DriftCoeff, DiffCoeff)
         
        real, intent(in) :: X_I(nDim), Time, Timestep
        real, intent(out) :: DriftCoeff(nDim), DiffCoeff(nDim)
  
        real :: rho, rhoFuture, rhoPast, dRhodTau
        real :: dS, dSdown, dSup, Momentum, Dxx

        Momentum = (3.0 * X_I(Momentum_))**(1.0/3.0)
        call get_dxx(Time, Momentum, Dxx)

        call get_ds(Time, X_I(LagrCoord_), dS)
        rho = 1.0 / dS
  
        ! found that this derivative converges around 0.01
        call get_rho(Time + 0.001, X_I(LagrCoord_), rhoFuture)
        call get_rho(Time - 0.001, X_I(LagrCoord_), rhoPast)
        dRhodTau = (rhoFuture - rhoPast) / 0.002
  
        if(StratoFactor.eq.1) then ! ito
  
           call get_ds(Time, X_I(LagrCoord_) - dShock, dSdown)
           call get_ds(Time, X_I(LagrCoord_) + dShock, dSup) 
           DriftCoeff(LagrCoord_) = Dxx * (1.0/dSup -  1.0/dSdown) / (dS * 2 * dShock)
        else ! stratonovich
           DriftCoeff(LagrCoord_) = 0.0
        end if
  
        DiffCoeff(LagrCoord_) = sqrt(2.0 * Dxx) / dS
  
        DriftCoeff(Momentum_) = X_I(Momentum_) * dRhodTau / rho
        DiffCoeff(Momentum_) = 0.0
  
    end subroutine get_sde_coeffs_rk2
    !=====================================================================!
    subroutine get_random_shock_location(Time, LagrCoord, S)
        real, intent(out) :: Time, LagrCoord, S
        real :: RandomUniform
        
        call random_number(RandomUniform)
  
        LagrCoord = tMin + (tMax - tMin) * RandomUniform
        call get_shock_arrival(LagrCoord, Time)
        Time = LagrCoord - 2.0 * dShock
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
    subroutine compute_timestep(Time, LagrCoord, Momentum, Timestep)
        real, intent(in) :: Time, LagrCoord, Momentum
        real, intent(out) :: Timestep
        real :: MaxTimeStep, Dxx
        
        call get_dxx(Time, Momentum, Dxx)

        Timestep = 0.01
        ! if(Time.gt.1000.0) then
        !     MaxTimeStep = 0.5 * (0.5*dShock + 0.25)**2 / Dxx
        !     Timestep = MaxTimeStep / 10.
        ! else
        !     Timestep = 0.01
        ! end if

    end subroutine compute_timestep
    !=====================================================================! 
    subroutine compute_conversion_factor(Time, LagrCoord, ConversionFactor)
        real, intent(in) :: Time, LagrCoord
        real, intent(out) :: ConversionFactor
        real :: dS
        call get_ds(Time, LagrCoord, dS)
        ConversionFactor = 1.0 / dS
    end subroutine compute_conversion_factor
    !=============================================================================
    subroutine save_location_properties(iProc, Time, LagrCoord, Momentum)
        integer, intent(in) :: iProc
        real, intent(in) :: Time, LagrCoord, Momentum
        integer :: FileUID

        character(len = 20) :: iProcStr, nStr

        real :: rho, Dxx, dS

        call get_ds(Time, LagrCoord, dS)
        call get_dxx(Time, Momentum, Dxx)
        call get_rho(Time, LagrCoord, rho)

        FileUID = iProc
        write(iProcStr, *) iProc
        iProcStr = adjustl(iProcStr)

        open(FileUID, file = 'PT/IO2/location_'//trim(iProcStr)//'.dat', &
            status = 'unknown', position = 'append')
        write(FileUID, *) Time, LagrCoord, Momentum, rho, Dxx, dS
        close(FileUID)

    end subroutine save_location_properties
    !=====================================================================! 
end module PT_ModTestFieldline
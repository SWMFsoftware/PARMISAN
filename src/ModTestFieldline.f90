module PT_ModTestFieldline
    use ModKind
    use ModMpi
    use PT_ModConst

    implicit none

    SAVE

    real, parameter :: tMin = 0.0, tMax = 10000.0
    real, parameter :: rMin = -5000.0, rMax = 15000.0
    real, parameter :: Dxx = 10.0
 
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
    subroutine get_shock_arrival(LagrCoord, ShockTime)
        real, intent(in) :: LagrCoord
        real, intent(out) :: ShockTime
  
        ShockTime = LagrCoord + dShock
  
    end subroutine get_shock_arrival
    !=====================================================================!
    subroutine get_sde_coeffs(X_I, Time, Timestep, DriftCoeff, DiffCoeff)
         
        real, intent(in) :: X_I(nDim), Time, Timestep
        real, intent(out) :: DriftCoeff(nDim), DiffCoeff(nDim)
  
        real :: rho, rhoFuture, dRhodTau, dS, dSdown, dSup
        ! real :: newLagrCoord
  
        call get_ds(Time, X_I(LagrCoord_), dS)
        rho = 1.0 / dS
  
        ! found that this derivative converges around 0.01
        call get_rho(Time + 0.001, X_I(LagrCoord_), rhoFuture)
        dRhodTau = (rhoFuture - rho) / 0.001
  
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
  
    end subroutine get_sde_coeffs
    !=====================================================================!
    subroutine get_random_shock_location(Time, LagrCoord, S)
        real, intent(out) :: Time, LagrCoord, S
        real :: RandomUniform
        
        call random_number(RandomUniform)
  
        LagrCoord = 0.0 + (10000.0 - 0.0) * RandomUniform
        ! call get_shock_arrival(LagrCoord, Time)
        Time = LagrCoord + 30.0
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

        Timestep = 0.01
    end subroutine compute_timestep
    !=====================================================================! 
    subroutine compute_conversion_factor(Time, LagrCoord1, LagrCoord2, ConversionFactor)
        real, intent(in) :: Time, LagrCoord1, LagrCoord2
        real, intent(out) :: ConversionFactor
        real :: dS1, dS2
        call get_ds(Time, LagrCoord1, dS1)
        call get_ds(Time, LagrCoord2, dS2)

        ConversionFactor = 0.5 * (1/dS1 + 1/dS2)
    end subroutine compute_conversion_factor
    !=====================================================================!
end module PT_ModTestFieldline
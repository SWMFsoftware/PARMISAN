module PT_ModSolver

    use PT_ModRandom, ONLY: get_random_normal

    use PT_ModFieldline, ONLY: get_sde_coeffs, StratoFactor, nDim
    ! use PT_ModTestFieldline, ONLY: get_sde_coeffs, StratoFactor, nDim
    implicit none
    SAVE
    real :: Wk, Sk

contains
    !==============================================================================
    ! subroutine milstein(X_I, Time, Timestep, Weight, Xnew_I)
        !     real, intent(in) :: X_I(nDim), Time
        !     real, intent(out) :: Timestep, Weight, Xnew_I(nDim)

        !     real :: DriftCoeff(nDim), DiffCoeff(nDim), dDiffdX(nDim)
        !     real :: RandomNormal

        !     call get_sde_coeffs_milstein(X_I, Timestep, Weight, DriftCoeff, DiffCoeff, dDiffdX)
            
        !     call get_random_normal(RandomNormal)
        !     Wk = sqrt(Timestep)*RandomNormal

        !     Xnew_I = X_I + DriftCoeff * Timestep + DiffCoeff * Wk + &
        !              dDiffdX * (Wk**2 - Timestep)

    ! end subroutine milstein
    !==============================================================================
    subroutine rk2_sde(X_I, Time, Timestep, Xnew_I)
        ! Runge Kutta 2 Scheme for SDEs
        ! Roberts 2012  "Modify the improved Euler scheme..."
        real, intent(in) :: X_I(nDim), Time, Timestep
        real, intent(out) :: Xnew_I(nDim)

        real :: RandomNormal, RandomUniform
        real :: K1_I(nDim), K2_I(nDim) ! runge-kutta coefficients

        ! the same Wk and Sk needs to be used for both K1, K2
        call get_random_normal(RandomNormal)
        call random_number(RandomUniform)
        
        Sk = StratoFactor * sign(sqrt(Timestep), RandomUniform-0.5)
        Wk = sqrt(Timestep)*RandomNormal - Sk
        
        call calculate_k(X_I, Time, TimeStep, K1_I)

        Wk = sqrt(Timestep)*RandomNormal + Sk

        call calculate_k(X_I + K1_I, Time + Timestep, Timestep, K2_I)

        Xnew_I = X_I + 0.5 * (K1_I + K2_I)

    end subroutine rk2_sde
    !==============================================================================
    subroutine calculate_k(X_I, Time, Timestep, K_I)
        ! Calculate K1 and K2 for second order RK scheme for SDE
        real, intent(in) :: X_I(nDim), Time, Timestep
        real, intent(out) :: K_I(nDim)

        integer :: i
        real :: DriftCoeff(nDim), DiffCoeff(nDim) !, pDiff

        call get_sde_coeffs(X_I, Time, Timestep, DriftCoeff, DiffCoeff)
        do i = 1, nDim
            K_I(i) = Timestep * DriftCoeff(i) + Wk * DiffCoeff(i)
        end do

    end subroutine calculate_k
    !==============================================================================
end module PT_ModSolver
!==============================================================================

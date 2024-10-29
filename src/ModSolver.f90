module PT_ModSolver
    use PT_ModConst, ONLY: fourpi
    
    use PT_ModKappa, ONLY: getK
    use PT_ModShockPara, ONLY: getshock, getU, tmin_data, &
        rmax_data, Mach, s, v_shock, V_sw_mod, r_shock, rMin, p0

    implicit none

contains
!==============================================================================
subroutine predictor_corrector(r, p, t, dt, rNew, pNew)
    ! Original predictor-corrector method
    real, intent(in) :: r, p, t, dt
    real, intent(out) :: rNew, pNew

    real :: U, dUdr, K, dKdr, divU
    real :: rn, xi
    real :: rH, pH

    real, parameter :: OneThird = 1.0/3.0
    
    ! advance to the next position and momentum with
    ! stochastic integration method
    ! first step is the prediction step
    call getU(r, t, U, dUdr)
    call getK(r, p, K, dKdr)

    divU=2.d0*U/r + dUdr

    call random_number(rn)
    xi=-1.d0+2.d0*rn

    rH = r + xi*sqrt(6.d0*K*dt) + U*dt + (dKdr+2.d0*K/r)*dt
    pH = p*(1.d0 - OneThird*divU*dt)

    ! second step is the corrector step
    call getU(rH,t,U,dUdr)
    call getK(rH,pH,K,dKdr)
    divU=2.d0*U/rH+dUdr

    rNew = r + xi*sqrt(6.d0*K*dt) + U*dt + (dKdr+2.d0*K/rH)*dt
    pNew = p*(1.d0 - OneThird*divU*dt)

end subroutine predictor_corrector
!==============================================================================
subroutine rk2_sde(r, p, t, dt, rNew, pNew)
    ! Runge Kutta Scheme for SDEs
    ! Stratonovich Version (S_k = 0)
    ! Roberts 2012  "Modify the improved Euler scheme..."
    real, intent(in) :: r, p, t, dt
    real, intent(out) :: rNew, pNew
    
    real :: k1r, k1p, k2r, k2p ! runge-kutta coefficients

    call calculate_k(r, p, t, dt, k1r, k1p)
    call calculate_k(r + k1r, p+k1p, t+dt, dt, k2r, k2p)

    rNew = r + 0.5 * (k1r + k2r)
    pNew = p + 0.5 * (k1p + k2p)

end subroutine rk2_sde
!==============================================================================
subroutine calculate_k(r, p, t, dt, kr, kp)
    ! Calculate K1 and K2 for second order RK scheme for SDE
    real, intent(in) :: r, p, t, dt
    real, intent(out) :: kr, kp

    real :: rn1, rn2, wk1, wk2
    real :: rDrift, rDiff, pDrift, pDiff

    call get_sde_coeffs(r, p, t, rDrift, rDiff, pDrift, pDiff)

    call get_random_normal(rn1)
    call get_random_normal(rn2)
    wk1 = sqrt(dt)*rn1
    wk2 = sqrt(dt)*rn2
    
    kr = dt * rDrift + wk1 * rDiff   ! different wks?
    kp = dt * pDrift + wk2 * pDiff

end subroutine calculate_k
!==============================================================================
subroutine get_sde_coeffs(r, p, t, rDrift, rDiff, pDrift, pDiff)
    ! Get drift and diffusion coefficients for the Stratonovich
    ! representation - see Pei 2010 for Ito coefficients
    real, intent(in) :: r, p, t
    real, intent(out) :: rDrift, rDiff, pDrift, pDiff

    real :: U, dUdr, K, dKdr

    call getU(r, t, U, dUdr)
    call getK(r, p, K, dKdr)

    rDrift = U + (2.0*K/r) + 0.5*dKdr
    rDiff = sqrt(2.0*K)

    pDrift = -p/3.0 * (2*U/r + dUdr)
    pDiff = 0.0

end subroutine get_sde_coeffs
!==============================================================================
subroutine get_random_normal(rn)
    ! returns random number sampled from normal distribution
    ! with mean = 0 and std = 1

    real, intent(out) :: rn
    real :: ru1, ru2

    ! uniform random numbers over [0,1)
    call random_number(ru1)
    call random_number(ru2)

    ! redistribute to (0, 1] to avoid 0
    ru1 = 1 - ru1
    ru2 = 1 - ru2
    
    ! Box-Muller transformation
    rn = sqrt(-2*log(ru1))*cos(0.5*fourpi*ru2)

end subroutine get_random_normal
!==============================================================================
end module PT_ModSolver
!==============================================================================

  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModKappa
  use PT_ModShockPara, ONLY: Rmin, p0
  implicit none
  SAVE
  real, parameter :: dlt = 0.710
  real :: K0
contains
  !============================================================================
  subroutine set_kappa
    use ModConst, ONLY: cAU
    !--------------------------------------------------------------------------

    K0 = 6.d+19 * ((Rmin/cAU)**1.17)   ! K0 is the value at r=Rmin at E0
  end subroutine set_kappa
  !============================================================================
  subroutine getK(r, Momentum, Kappa, dKappaDr)
    real, intent(in)  :: r, Momentum
    real, intent(out) :: Kappa, dKappaDr
    real :: MomentumFactor
    !--------------------------------------------------------------------------
    Kappa = K0*((r/Rmin)**1.170)
    dKappaDr = 1.170*Kappa/r
    MomentumFactor = (Momentum/p0)**(2.0*dlt)
    Kappa = Kappa*MomentumFactor
    dKappaDr=dKappaDr*MomentumFactor
  end subroutine getK
  !============================================================================
end module PT_ModKappa
!==============================================================================

  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModKappa
  use PT_ModShockPara, ONLY: Rmin, K0, p0
  implicit none
  SAVE
  real, parameter :: dlt = 0.710
contains
  subroutine getK(r, t, p, Kappa, dKappaDr)
    real, intent(in)  :: r, t, p
    real, intent(out) :: Kappa, dKappaDr
    !--------------------------------------------------------------------------
  Kappa = K0*((r/Rmin)**1.170)
  dKappaDr = 1.170*Kappa/r
  Kappa = Kappa*((p/p0)**(2.0*dlt))
  dKappaDr=dKappaDr*((p/p0)**(2.0*dlt))
end subroutine getK
end module PT_ModKappa
!==============================================================================

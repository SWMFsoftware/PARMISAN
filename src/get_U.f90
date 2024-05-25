
subroutine getU(r, t, U, dUdr)
  ! based on the appendix in Giacalone, 2015 (ApJ, vol. 799,
  !     article id. 80, Equations C1-C4)
  use PT_ModShockPara, ONLY: r_shock, s, drSHOCK, Ls, v_shock, getshock
  use PT_ModSWPara,    ONLY: getVsw
  implicit none
  real, intent(in)  :: r,t
  real, intent(out) :: U, dUdr
  real :: U2p,U1p,r_shock_p,coshR2, Vsw
  ! Why do you need this 
  ! integer :: ic
  !----------------------------------------------------------------------------
  ! Why do you need this?
  ! data ic/0/
  ! save
  !
  ! if(ic == 0)then
  !   ic=1
  ! end if

  call getShock( t )
  call getVsw( r, Vsw )

  U2p = ((s-1.d0)*v_shock + Vsw)/s
  U1p = Vsw

  r_shock_p = r_shock - 3.d0*drSHOCK
  if(r > r_shock_p)then
    U = 0.5*(U1p+U2p)+0.5*(U1p-U2p)*dtanh((r-r_shock)/drSHOCK)
    coshR2 = (1.d0/dcosh((r-r_shock)/drSHOCK))**2
    dUdr = 0.5*(U1p-U2p)*coshR2/drSHOCK
  elseif(r > r_shock_p-Ls .and. r <= r_shock_p )then
    U = U2p*((r_shock_p/r)**2)
    dUdr = -2.d0*U2p*(r_shock_p**2)/(r**3)
  else
    U = U2p*((r_shock_p/(r_shock_p-Ls))**2)
    dUdr = 0.d0
  end if
end subroutine getU
!==============================================================================

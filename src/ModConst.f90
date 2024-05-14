  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModConst
  use ModConst, ONLY: cTwoPi, cAU, ceV, ckeV, cMeV, cProtonMass, &
       RsunSi => Rsun ! SI
  implicit none
  SAVE
  ! Constants in CGS
  real, parameter:: Rsun = RsunSi*100.0    ! cm
  real, parameter:: AU = cAU*100.0         ! cm
  real, parameter:: eV = ceV*1.0d7         ! erg
  real, parameter:: keV= ckeV*1.0d7        ! erg
  real, parameter:: MeV= cMeV*1.0d7        ! erg
  real, parameter:: mp = cProtonMass*1000.0! g
  real, parameter:: fourpi = 2.0*cTwoPi
end module PT_ModConst
!==============================================================================

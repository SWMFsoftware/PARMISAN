  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModConst
  use ModNumConst, ONLY: TwoPi => cTwoPi   
  use ModConst, ONLY: AU => cAU,    &
                      eV => ceV,   &
                      keV => ckeV, &
                      MeV => cMeV, &
                      LightSpeed => cLightSpeed, & 
                      ProtonMass => cProtonMass, &
                      Rsun ! SI
  implicit none
  SAVE
  ! Constants in CGS
  real, parameter :: cRsun = Rsun*100.0    ! cm
  real, parameter :: cAU = AU*100.0         ! cm
  real, parameter :: ceV = eV*1.0d7         ! erg
  real, parameter :: ckeV = keV*1.0d7        ! erg
  real, parameter :: cMeV = MeV*1.0d7        ! erg
  real, parameter :: cProtonMass = ProtonMass*1000.0! g
  real, parameter :: cFourPi = 2.0*TwoPi
  real, parameter :: cLightSpeed = LightSpeed*100.0 ! cm/s
  real, parameter :: cProtonRestEnergy = cProtonMass * cLightSpeed**2
end module PT_ModConst
!==============================================================================

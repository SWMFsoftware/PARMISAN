  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModConst
  use ModNumConst, ONLY: TwoPi => cTwoPi   
  use ModConst, ONLY: cAU, ceV, ckeV, cMeV, &
                      cLightSpeed, cProtonMass, &
                      cElectronCharge, cRsun => rSun 

  implicit none
  SAVE
  ! Constants in SI
  real, parameter :: cFourPi = 2.0*TwoPi
  real, parameter :: cPi = TwoPi / 2.0
  real, parameter :: cProtonRestEnergy = cProtonMass * cLightSpeed**2

end module PT_ModConst
!==============================================================================

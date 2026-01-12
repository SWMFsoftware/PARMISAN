  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModSize

  implicit none

  private ! except
  public :: nDim, nParticle, nVertexMax

  ! Dimensionality
  integer, parameter:: nDim = 2

  ! Max possible index of a particle on a line set by Config.pl
  integer, parameter:: nVertexMax  = 4000

  ! number of particles per field line
  integer, parameter:: nParticle   = 4000

end module PT_ModSize
!==============================================================================

  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModPlot
  
  use PT_ModConst
  use PT_ModProc, ONLY : iProc, iComm, iError
  
  implicit none
  SAVE
  
  character(len=*), parameter :: OutputDir = 'PT/IO2/'
  character(len=*), parameter :: SplitFile ='split_levels.dat'
  character(len=*), parameter :: TimeBinFile ='time_bin.dat'
  character(len=*), parameter :: EnergyBinFile ='energy_bin.dat'
  character(len=*), parameter :: rCrossFile ='crossings.dat'
  
  integer :: nTimeBin, nEnergybin, nCrossings
  real :: tBinMin, tBinMax, eBinMin, eBinMax, crossMin, crossMax
  
  ! Array for storing the fluxes at several spherical surfaces:
  real, allocatable :: W_III(:,:,:)
  
  character(LEN=:), allocatable, dimension(:) :: NameOutputFile_I
  
  ! Bins
  real, allocatable :: Ebin_I(:), TimeBin_I(:), rCross_I(:)
  
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: CON_stop

    character(len=*), intent(in):: NameCommand ! From PARAM.in
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#PLOT')

       call read_var('nCrossings', nCrossings)
       call read_var('crossMin', crossMin)
       call read_var('crossMax', crossMax)
       call read_var('nTimeBin', nTimeBin )
       call read_var('tBinMin', tBinMin)
       call read_var('tBinMax', tBinMax)
       call read_var('nEnergybin', nEnergybin)
       call read_var('eBinMin', eBinMin)
       call read_var('eBinMax', eBinMax)

      ! Convert to seconds
      tBinMin = tBinMin * 3600.
      tBinMax = tBinMax * 3600.

      ! Convert to keV
      eBinMin = eBinMin * ckeV
      eBinMax = eBinMax * ckeV
      ! Convert to cm
      crossMin = crossMin * cRsun
      crossMax = crossMax * cRsun

    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine set_bins()

    integer :: iBin
    real :: dLogE, dT, dR

    allocate(W_III(nTimeBin, nEnergybin, nCrossings)); W_III = 0.0
    allocate(TimeBin_I(nTimeBin))
    allocate(Ebin_I(nEnergybin))
    allocate(rCross_I(nCrossings))

    dT = (tBinMax - tBinMin) / (nTimeBin-1)
    do iBin = 1, nTimeBin
      TimeBin_I(iBin) = tBinMin + dT * (iBin-1) 
    end do

    dLogE = (log10(eBinMax) - log10(eBinMin)) / (nEnergybin-1)
    do iBin = 1, nEnergybin
       Ebin_I(iBin) = 10.0**(log10(eBinMin)+ dLogE*(iBin-1))
    end do

    dR = (crossMax - crossMin) / (nCrossings-1)
    do iBin = 1, nCrossings
      rCross_I(iBin) = crossMin + dR * (iBin - 1)
    end do

    if(iProc==0) then
       open(801,file=OutputDir//EnergyBinFile,status='unknown')
       write(801,'(1000e15.6)')Ebin_I/ckeV
       close(801)
       open(801,file=OutputDir//TimeBinFile,status='unknown')
       write(801,'(1000e15.6)')TimeBin_I
       close(801)
       open(801,file=OutputDir//rCrossFile,status='unknown')
       write(801,'(1000e15.6)')rCross_I/cRsun
       close(801)
    end if

  end subroutine set_bins
  !============================================================================
  subroutine put_flux_contribution(rStart, rEnd, Time, Energy, Contribution)
    ! Distance along field line at the beginning and at the end of the time step
    real, intent(in) :: rStart, rEnd
    ! Time and energy, to determine the bin to put the contribution:
    real, intent(in) :: Time, Energy
    real, intent(in) :: Contribution
    ! Misc:
    logical :: IsCrossed_I(nCrossings)
    ! Pixel numbers:
    integer :: iDt, iDe
    !--------------------------------------------------------------------------
    IsCrossed_I = (rEnd - rCross_I)*(rStart - rCross_I) < 0.0
    if(.not.any(IsCrossed_I))RETURN
    iDt = minloc(abs(TimeBin_I - Time), DIM=1)
    iDe = minloc(abs(Ebin_I - Energy),  DIM=1)
    where(IsCrossed_I)&
         W_III(iDt,iDe,:) = W_III(iDt,iDe,:) + Contribution
  end subroutine put_flux_contribution
  !============================================================================
  subroutine save_fluxes

    use ModMpi, ONLY: MPI_reduce_real_array, MPI_SUM
    integer :: iTimeBin, iCross
    character(len = 50) :: CrossingFile

    ! Save fluxes stored at different radial distances
    !--------------------------------------------------------------------------
    call MPI_reduce_real_array(W_III, nTimeBin*nEnergybin*nCrossings, MPI_SUM, 0,&
         iComm, iError)
    
    if(iProc/=0) RETURN
    
    do iCross = 1, nCrossings
      write(CrossingFile, '(A13, I0)') 'fluxes_cross_', iCross
      CrossingFile = adjustl(CrossingFile)
       open(901,file=OutputDir//trim(CrossingFile),&
            status='unknown',action="READWRITE")
       do iTimeBin = 1, nTimeBin
          write(901,'(1000e15.6)')W_III(iTimeBin,:,iCross)
       end do
       close(901)
    end do
  end subroutine save_fluxes
end module PT_ModPlot
!==============================================================================

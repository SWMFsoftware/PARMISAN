  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModFluxes
  use PT_ModConst
  use PT_ModProc, ONLY : iProc
  implicit none
  SAVE
  character(len=*), parameter ::  NamePath = &
       '/Users/xhchen/My_work/Matlab/SEP_LaborDay/mpi/data/'
  character(len=*), parameter :: file1 = 'Cross8_'
  character(len=*), parameter :: file2 = 'Cross10_'
  character(len=*), parameter :: file3 = 'Cross12_'
  character(len=*), parameter :: file4 = 'Cross14_'
  character(len=*), parameter :: name_sl ='split_levels.dat'
  character(len=*), parameter :: name_pt ='prms_vs_t.dat'
  character(len=*), parameter :: Time_bin ='t_bin.dat'
  character(len=*), parameter :: Energy_bin ='E_bin.dat'
  character(len=*), parameter :: name_dist = 'dist'
  character(len=100)          :: NameFile
  integer, parameter :: nEbin = 1000, nTimeBin = 1000
  ! Arrays for storing the fluxes at the face:
  real :: w1(nTimeBin, nEbin),w2(nTimeBin, nEbin), &
       w3(nTimeBin, nEbin),w4(nTimeBin, nEbin)
  ! Bins
  real :: Ebin_I(nEbin),TimeBin_I(nTimeBin)
contains
  subroutine set_bins(TimeMin, TimeMax, Emin, Emax)

    real, intent(in) :: TimeMin, TimeMax, Emin, Emax
    integer :: iBin
    !--------------------------------------------------------------------------
    
    do iBin = 1, nTimeBin
       TimeBin_I(iBin) = TimeMin + (TimeMax - TimeMin)*(iBin - 1)/nTimeBin
    end do
    
    do iBin = 1, nEbin
       Ebin_I(iBin) = exp(log(Emin)+(log(Emax)-log(Emin)) &
            /nEbin*(iBin-1) )
    end do
    if(iProc==0) then
       open(801,file=NamePath//Energy_bin,status='unknown')
       write(801,'(1000e15.6)')Ebin_I/keV
       close(801)
       open(802,file=NamePath//Time_bin,status='unknown')
       write(802,'(1000e15.6)')TimeBin_I
       close(802)
    end if
  end subroutine set_bins
  !============================================================================
  subroutine save_fluxes

    integer :: iBin
    !--------------------------------------------------------------------------

    ! Save fluxes stored at different surfaces
    write(NameFile,'(a,i4.4,a)')NamePath//file1, iProc, '.dat'
    open(901,file=trim(NameFile),status='unknown',action="READWRITE")
    write(NameFile,'(a,i4.4,a)')NamePath//file2, iProc, '.dat'
    open(902,file=trim(NameFile),status='unknown',action="READWRITE")
    write(NameFile,'(a,i4.4,a)')NamePath//file3, iProc, '.dat'
    open(903,file=trim(NameFile),status='unknown',action="READWRITE")
    write(NameFile,'(a,i4.4,a)')NamePath//file4, iProc, '.dat'
    open(904,file=trim(NameFile),status='unknown',action="READWRITE")
    
    do iBin = 1, nTimeBin
       write(901,'(1000e15.6)')w1(iBin,:)
       write(902,'(1000e15.6)')w2(iBin,:)
       write(903,'(1000e15.6)')w3(iBin,:)
       write(904,'(1000e15.6)')w4(iBin,:)
    end do
    close(901)
    close(902)
    close(903)
    close(904)
  end subroutine save_fluxes
  !============================================================================
end module PT_ModFluxes
!==============================================================================

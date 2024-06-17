  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModPlot
  use PT_ModConst
  use PT_ModProc, ONLY : iProc, iComm, iError
  implicit none
  SAVE
  character(len=*), parameter ::  NamePath = 'PT/IO2/'
  character(len=*), parameter :: name_sl ='split_levels.dat'
  character(len=*), parameter :: name_pt ='prms_vs_t.dat'
  character(len=*), parameter :: Time_bin ='t_bin.dat'
  character(len=*), parameter :: Energy_bin ='E_bin.dat'
  character(len=*), parameter :: name_dist = 'dist'
  character(len=100)          :: NameFile
  integer :: nTimeBin = 1000, nEbin = 1000, nPlotFile = 4
  ! Array for storing the fluxes at several spherical surfaces:
  real, allocatable :: W_III(:,:,:)
  ! Radii of the said surfaces:
  real, allocatable :: rPlot_I(:)
  character(LEN=:), allocatable, dimension(:) :: NamePlotFile_I
  ! Bins
  real, allocatable :: Ebin_I(:),TimeBin_I(:)
  ! Diagnostics:
  integer, parameter :: ipmx=1000, itdmx=500
  integer :: np(itdmx)
  real :: dN(itdmx,ipmx)
  ! Constants set by set_diagostics
  real :: dpop_diag, dpl_diag, p_diag_min, dp_diag, dtd, den_ep
  integer :: ntd
contains
  !============================================================================
  subroutine read_param(NameCommand)

    use ModReadParam, ONLY: read_var
    use ModUtilities, ONLY: CON_stop

    character(len=*), intent(in):: NameCommand ! From PARAM.in
    ! Loop variable
    integer :: iPlot
    character(len=*), parameter:: NameSub = 'read_param'
    !--------------------------------------------------------------------------
    select case(NameCommand)
    case('#PLOT')
       call read_var('nTimeBin' ,nTimeBin )
       call read_var('nEbin'    ,nEbin    )
       call read_var('nPlotFile',nPlotFile)
       allocate(rPlot_I(nPlotFile))
       allocate(character(LEN=13) :: NamePlotFile_I(nPlotFile))
       do iPlot = 1, nPlotFile
          call read_var('NamePlotFile_I(:)', NamePlotFile_I(iPlot))
          call read_var('rPlot_I(:)'       , rPlot_I(iPlot))
       end do
       ! Convert to CGS
       rPlot_I = rPlot_I*Rsun
    case default
       call CON_stop(NameSub//' Unknown command '//NameCommand)
    end select
  end subroutine read_param
  !============================================================================
  subroutine write_shock_file(TimeMin, TimeMax, dTime)

    use PT_ModShockPara, ONLY: r_shock, v_shock, s, Mach, V_sw_mod, p0, &
         getShock
    use PT_ModKappa,     ONLY: getK
    real, intent(in) :: TimeMin, TimeMax, dTime
    real :: Time, r, Kappa, dKappadR
    !--------------------------------------------------------------------------
    Time = TimeMin - dTime
    open(44, file = NamePath//name_pt, status='unknown')
    do
       Time = Time + dTime
       if(Time > TimeMax) EXIT
       call getShock(Time)
       r=r_shock
       call getK(r, p0, Kappa, dKappaDr)
       write(44,'(6f18.7,2e15.5)')Time/3600.0, Mach, s, &
            v_shock/1.0d+5, V_sw_mod/1.0d+5, r_shock/Rsun, Kappa, dKappaDr
    end do
    close(44)
  end subroutine write_shock_file
  !============================================================================
  subroutine set_bins(TimeMin, TimeMax, Emin, Emax)

    real, intent(in) :: TimeMin, TimeMax, Emin, Emax
    integer :: iBin

    !--------------------------------------------------------------------------
    if(.not.allocated(rPlot_I))then
       nPlotFile = 4
       allocate(rPlot_I(nPlotFile))
       rPlot_I = [8.0, 10.0, 12.0, 14.0]*Rsun
       allocate(character(LEN=13) :: NamePlotFile_I(nPlotFile))
       NamePlotFile_I = ['Cross8.dat   ', 'Cross10.dat  ',&
            'Cross12.dat  ', 'Cross14.dat  ']
    end if
    allocate(W_III(nTimeBin, nEbin, nPlotFile)); W_III = 0.0
    allocate(TimeBin_I(nTimeBin))
    do iBin = 1, nTimeBin
       TimeBin_I(iBin) = TimeMin + (TimeMax - TimeMin)*(iBin - 1)/nTimeBin
    end do
    allocate(Ebin_I(nEbin))
    do iBin = 1, nEbin
       Ebin_I(iBin) = exp(log(Emin)+(log(Emax)-log(Emin)) &
            /nEbin*(iBin-1) )
    end do
    if(iProc==0) then
       open(801,file=NamePath//Energy_bin,status='unknown')
       write(801,'(1000e15.6)')Ebin_I/keV
       close(801)
       open(801,file=NamePath//Time_bin,status='unknown')
       write(801,'(1000e15.6)')TimeBin_I
       close(801)
    end if
  end subroutine set_bins
  !============================================================================
  subroutine put_flux_contribution(rStart, rEnd, Time, Energy, Contribution)
    ! Helisopheric distance at the beginning and at the end of the time step
    real, intent(in) :: rStart, rEnd
    ! Time and energy, to determine the bin to put the contribution:
    real, intent(in) :: Time, Energy
    real, intent(in) :: Contribution
    ! Misc:
    logical :: IsCrossed_I(nPlotFile)
    ! Pixel numbers:
    integer :: iDt, iDe
    !--------------------------------------------------------------------------
    IsCrossed_I = (rEnd - rPlot_I)*(rStart - rPlot_I) < 0.0
    if(.not.any(IsCrossed_I))RETURN
    idt = minloc(abs(TimeBin_I - Time), DIM=1)
    idE = minloc(abs(Ebin_I - Energy),  DIM=1)
    where(IsCrossed_I)&
         W_III(iDt,iDe,:) = W_III(iDt,iDe,:) + Contribution
  end subroutine put_flux_contribution
  !============================================================================
  subroutine save_fluxes

    use ModMpi, ONLY: MPI_reduce_real_array, MPI_SUM
    integer :: iBin, iPlot

    ! Save fluxes stored at different surfaces
    !--------------------------------------------------------------------------
    call MPI_reduce_real_array(W_III, nTimeBin*nEbin*nPlotFile, MPI_SUM, 0,&
         iComm, iError)
    if(iProc/=0)RETURN
    do iPlot = 1, nPlotFile
       open(901,file=NamePath//trim(NamePlotFile_I(iPlot)),&
            status='unknown',action="READWRITE")
       do iBin = 1, nTimeBin
          write(901,'(1000e15.6)')W_III(iBin,:,iPlot)
       end do
       close(901)
    end do
  end subroutine save_fluxes
  !============================================================================
  subroutine set_diagnostics(E0)
    real, intent(in) :: E0
    ! Loop variables
    integer :: itd, ip
    !------diagnostic initializations -------
    !--------------------------------------------------------------------------
    dpop_diag = 0.04
    dpl_diag = log(1.d0+dpop_diag)
    p_diag_min = sqrt(2.d0*E0*mp)
    dtd=60.d0
    ntd=0
    do itd = 1,itdmx
       np(itd)=0
       do ip = 1,ipmx
          dN(itd,ip)=0.d0
       end do
    end do
    den_ep = 2.4d-3
  end subroutine set_diagnostics
  !============================================================================
  subroutine put_diagnostic_contribution(Time, Momentum, Weight)
    real, intent(in) :: Time, Momentum, Weight

    real    :: ap, fp, fpc,rtd,ftd,ftdc
    integer :: ip,ip1,pmx,itd,itd1
    !--------------------------------------------------------------------------

    ap=log(Momentum/p_diag_min)/dpl_diag+0.50000000001
    ip=ap
    ip1=ip+1
    fp=ap-ip
    fpc=1.d0-fp
    ip=min(ipmx,max(1,ip))
    ip1=min(ipmx,max(1,ip1))
    rtd=Time/dtd + 0.5000000000001
    itd=rtd
    itd1=itd+1
    ftd=rtd-itd
    ftdc=1.d0-ftd
    itd=min(itdmx,max(1,itd))
    itd1=min(itdmx,max(1,itd1))
    ntd=max(ntd,itd)
    np(itd)=max(np(itd),ip)
    dN(itd,ip)  = dN(itd,ip)  + Weight*fpc*ftdc
    dN(itd,ip1) = dN(itd,ip1) + Weight*fp*ftdc
    dN(itd1,ip) = dN(itd1,ip) + Weight**fpc*ftd
    dN(itd1,ip1)= dN(itd1,ip1)+ Weight**fp*ftd
  end subroutine put_diagnostic_contribution
  !============================================================================
  subroutine write_diagnostics
    ! Total flux and normalipation factor
    real :: fac, dNs, f, j
    ! Loop variables
    integer:: itd, ip
    ! Energy, momentum
    real :: E, pp
    !--------------------------------------------------------------------------
    write( NameFile,'(a,i4.4,a)')NamePath//name_dist, iProc, '.dat'
    open(14,file=trim(NameFile),status='unknown')
    dNs=0.d0
    do itd = 1,ntd
       do ip = 1,np(itd)
          dNs=dNs+dN(itd,ip)
       end do
    end do
    if(dNs > 0.d0)then
       fac=1.d0/(fourpi*dNs)
    else
       fac=0.d0
    end if
    if(dNs > 0.d0)write(14,*)fac,ntd
    do itd = 1,ntd
       write(14,*)itd,np(itd)
       do ip = 1,np(itd)
          !         p=p_diag_min*exp( dpl_diag*real(ip-1))
          pp=p_diag_min*exp( dpl_diag*(real(ip)-0.5) )
          E=(pp**2)/(2.d0*mp)
          dp_diag=dpop_diag*pp
          ! distribution function:
          f=dN(itd,ip)/((pp**2)*dp_diag)*den_ep*fac
          j=f*(pp**2)                               ! diff. intensity
          write(14,*)E/MeV, dmax1(1.d-20,j*MeV)
          ! do this again, shifted by a bit to make it histogram style
          !         p=p_diag_min*exp( dpl_diag*real(ip))
          !         E=(p**2)/(2.d0*mp)
          !         dp_diag=dpop_diag*p
          !         f=dN(itd,ip)/((p**2)*dp_diag)*den_ep*fac
          !         j=f*(p**2)
          !         write(14,*)E/MeV,dmax1(1.d-20,j*MeV)
       end do
    end do
    close(14)
  end subroutine write_diagnostics
  !============================================================================
end module PT_ModPlot
!==============================================================================

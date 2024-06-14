  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModPlot
  use PT_ModConst
  use PT_ModProc, ONLY : iProc
  implicit none
  SAVE
  character(len=*), parameter ::  NamePath = 'PT/IO2'
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
  real :: w1(nTimeBin, nEbin), w2(nTimeBin, nEbin), &
       w3(nTimeBin, nEbin), w4(nTimeBin, nEbin)
  ! Bins
  real :: Ebin_I(nEbin),TimeBin_I(nTimeBin)
  ! Diagnostics:
  integer, parameter :: ipmx=1000, itdmx=500
  integer :: np(itdmx)
  real :: dN(itdmx,ipmx)
  ! Constants set by set_diagostics
  real :: dpop_diag, dpl_diag, p_diag_min, dp_diag, dtd, den_ep
  integer :: ntd
contains
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
    w1 = 0.0; w2 = 0.0; w3 = 0.0; w4 = 0.0
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
  subroutine put_flux_contribution(rStart, rEnd, Time, Energy, Contribution)
    ! Helisopheric distance at the beginning and at the end of the time step
    real, intent(in) :: rStart, rEnd
    ! Time and energy, to determine the bin to put the contribution:
    real, intent(in) :: Time, Energy
    real, intent(in) :: Contribution
    ! Pixel numbers:
    integer :: idt, ide
    !--------------------------------------------------------------------------

    if (((rEnd - 8.d0*Rsun)*(rStart - 8.d0*Rsun))<0.0) then
       idt = minloc(abs(TimeBin_I - Time),DIM=1)
       idE = minloc(abs(Ebin_I - Energy),DIM=1)
       w1(idt,idE) = w1(idt,idE) + Contribution
    end if
    if (((rEnd - 10.d0*Rsun)*(rStart - 10.d0*Rsun))<0.0) then
       idt = minloc(abs(TimeBin_I - Time),DIM=1)
       idE = minloc(abs(Ebin_I - Energy),DIM=1)
       w2(idt,idE) = w2(idt,idE) + Contribution
    end if
    if (((rEnd - 12.d0*Rsun)*(rStart - 12.d0*Rsun))<0.0) then
       idt = minloc(abs(TimeBin_I - Time),DIM=1)
       idE = minloc(abs(Ebin_I - Energy),DIM=1)
       w3(idt,idE) = w3(idt,idE) + Contribution
    end if
    if (((rEnd - 14.d0*Rsun)*(rStart - 14.d0*Rsun))<0.0) then
       idt = minloc(abs(TimeBin_I - Time),DIM=1)
       idE = minloc(abs(Ebin_I - Energy),DIM=1)
       w4(idt,idE) = w4(idt,idE) + Contribution
    end if
  end subroutine put_flux_contribution
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

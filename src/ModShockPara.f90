  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModShockPara
  use PT_ModConst, ONLY: Rsun
  implicit none
  SAVE
  real :: Rmin
  ! Momentum at the injection energy
  real :: p0
  real :: r_shock, v_shock, s, Mach, V_sw_mod
  ! Maximum and minimum time in the data set
  real :: tmax_data, tmin_data
  ! Maximum and minimum shock wave radius
  real :: Rmax_data,Rmin_data
  logical :: DoReadShockFile = .true.
  real, PRIVATE :: time_A(2000), Mach_A(2000), v_shock_A(2000), s_A(2000), &
       v_sw_mod_A(2000), r_shock_A(2000)
  integer, PRIVATE :: n
contains
  !============================================================================
  subroutine read_shock
    character(len=*), parameter :: name = 'PSP_para.dat', &
         path = '/Users/xhchen/My_work/Matlab/SEP_LaborDay/mpi/'
    ! Variables to read from the shock parameter file, before conversion
    real :: th, rshRsun, M, vshkms, Vswkms, shock_pos_rel_sun, ss
    ! Loop and io variables
    integer :: i, io
    !--------------------------------------------------------------------------
    if(.not.DoReadShockFile)RETURN
    DoReadShockFile = .false.
    open(10,file=path//name,status='old')
    i=0
    do
       read(10,*,iostat=io) th,rshRsun,M,vshkms,Vswkms,ss
       if (io /= 0) then
          ! Why not close(10)???
          close(10)
          EXIT
       end if
      shock_pos_rel_sun = rshRsun*Rsun
      i = i+1
      ! time in seconds
      time_A(i) = th*3600.d0
      ! shock radius (from center of Sun) in cm
      r_shock_A(i) = shock_pos_rel_sun
      ! shock Mach number
      Mach_A(i) = M
      ! shock speed in cm/s
      v_shock_A(i) = vshkms*1.d5
      ! solar wind speed at the shock front
      v_sw_mod_A(i) = Vswkms*1.d+5
      s_A(i) = max(ss, 1.0)
      ! Store the number of lines in the file
      n=i
    end do
    tmax_data = time_A(n)
    tmin_data = time_A(1)
    Rmax_data = r_shock_A(n)
    Rmin_data = r_shock_A(1)
    do i = 2, n
      v_shock_A(i) = (r_shock_A(i) - r_shock_A(i-1))/ &
      (time_A(i) - time_A(i-1))
    end do
    v_shock_A(1) = v_shock_A(2)
  end subroutine read_shock
  !============================================================================
  subroutine getShock(t)
    real, intent(in) :: t
    ! Why FF? FormFactor?
    real :: time1, time2, FF
    integer :: i
    !--------------------------------------------------------------------------
    ! interpolate
    do i = 1,n
       if(time_A(i) <= t) CYCLE
       if(i == 1)then
          Mach = Mach_A(1)
          s = s_A(1)
          v_shock = v_shock_A(1)
          r_shock = r_shock_A(1)
          V_sw_mod = v_sw_mod_A(1)
       else
          time1 = time_A(i-1)
          time2 = time_A(i)
          FF = (t-time1)/(time2-time1)
          Mach = Mach_A(i-1) + FF*(Mach_A(i)-Mach_A(i-1))
          s = s_A(i-1) + FF*(s_A(i)-s_A(i-1))
          v_shock = v_shock_A(i-1) + FF*(v_shock_A(i) - v_shock_A(i-1))
          r_shock = r_shock_A(i-1) + FF*(r_shock_A(i) - r_shock_A(i-1))
          V_sw_mod = V_sw_mod_A(i-1) + FF*(V_sw_mod_A(i) - V_sw_mod_A(i-1))
       end if
       RETURN
    end do
    ! Extrapolate? Time exceeds time_A(n)
    i = n
    time1=time_A(i-1)
    time2=time_A(i)
    FF = (t-time1)/(time2-time1)
    Mach = Mach_A(i-1) + FF*(Mach_A(i)-Mach_A(i-1))
    s = s_A(i-1) + FF*(s_A(i)-s_A(i-1))
    v_shock = v_shock_A(i-1) + FF*(v_shock_A(i)-v_shock_A(i-1))
    r_shock = r_shock_A(i-1) + FF*(r_shock_A(i)-r_shock_A(i-1))
    V_sw_mod = V_sw_mod_A(i-1) + FF*(V_sw_mod_A(i)-V_sw_mod_A(i-1))
  end subroutine getShock
  !============================================================================
  subroutine getVsw(r,VwP)
    real, intent(in)  :: r
    real, intent(out) :: VwP
    real :: rs1,rs2,FF
    integer :: i
    !--------------------------------------------------------------------------

    do i = 1,n
       if(r >= r_shock_A(i))CYCLE
       if(i == 1)then
          VwP=v_sw_mod_A(1)
       else
          rs1=r_shock_A(i-1)
          rs2=r_shock_A(i)
          FF = (r-rs1)/(rs2-rs1)
          VwP=V_sw_mod_A(i-1) + FF*(V_sw_mod_A(i)-V_sw_mod_A(i-1))
       end if
       RETURN
    end do
    VwP=v_sw_mod_A(n)
  end subroutine getVsw
  !============================================================================
  subroutine getU(r, t, U, dUdr)
    ! based on the appendix in Giacalone, 2015 (ApJ, vol. 799,
    !     article id. 80, Equations C1-C4)
    real, intent(in)  :: r,t
    real, intent(out) :: U, dUdr
    real :: U2p, U1p ,r_shock_p, coshR2, Vsw
    ! divergence free region in the downstream (optional)
    real, parameter ::  Ls = Rsun
    ! Shock wave width
    real, parameter :: drSHOCK = 4.d-5*Rsun

    !--------------------------------------------------------------------------
    call getShock( t )
    call getVsw( r, Vsw )

    U2p = ((s - 1.0)*v_shock + Vsw)/s
    U1p = Vsw

    r_shock_p = r_shock - 3.0*drSHOCK
    if(r > r_shock_p)then
       U = 0.50*(U1p + U2p) + 0.50*(U1p - U2p)*tanh((r - r_shock)/drSHOCK)
       coshR2 = (1.0/cosh((r - r_shock)/drSHOCK))**2
       dUdr = 0.50*(U1p - U2p)*coshR2/drSHOCK
    elseif(r > r_shock_p - Ls)then
       U = U2p*((r_shock_p/r)**2)
       dUdr = -2.0*U2p*(r_shock_p**2)/(r**3)
    else
       U = U2p*((r_shock_p/(r_shock_p - Ls))**2)
       dUdr = 0
    end if
  end subroutine getU
  !============================================================================
end module PT_ModShockPara
!==============================================================================

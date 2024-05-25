  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModShockPara
  use PT_ModSWPara, ONLY: n, r_shock_A, v_sw_mod_A
  implicit none
  SAVE
  ! Minimal heliocentric distance?
  real :: Rmin
  ! Shock wave width?
  real :: drSHOCK
  ! What is this?
  real :: Vw
  ! What is this?
  real :: Ls
  ! What is this?
  real :: K0
  ! What is this?
  real :: p0
  ! What is this?
  real :: dlt
  real :: r_shock, v_shock, s, Mach, V_sw_mod
  ! Maximum and minimum time in the data set
  real :: tmax_data, tmin_data
  ! Maximum shock wave radius?
  real :: Rmax_data
  logical :: DoReadShockFile = .true.
  real, PRIVATE :: time_A(2000),Mach_A(2000),v_shock_A(2000),s_A(2000)
contains
  subroutine read_shock
    use PT_ModConst, ONLY: Rsun
    character(len=*), parameter :: name = 'PSP_para.dat', &
         path = '/Users/xhchen/My_work/Matlab/SEP_LaborDay/mpi/'
    ! Variables to read from the shock parameter file, before conversion
    real :: th, rshRsun, M, vshkms, Vswkms, shock_pos_rel_sun, ss
    ! What is this? Not used
    ! real :: x_sh_Rsun,y_sh_Rsun,z_sh_Rsun
    ! Loop and io variables
    integer :: i, io 
    !--------------------------------------------------------------------------
    DoReadShockFile = .false.
    open(10,file=trim(path)//trim(name),status='old')
    i=0
    do
       read(10,*,iostat=io) th,rshRsun,M,vshkms,Vswkms,ss
       if (io /= 0) then
          ! Why not close(10)???
          EXIT
       end if
      shock_pos_rel_sun = rshRsun*Rsun
      i=i+1
      ! time in seconds
      time_A(i)=th*3600.d0
      ! shock radius (from center of Sun) in cm
      r_shock_A(i)=shock_pos_rel_sun
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
    do i = 2, n
      v_shock_A(i)=(r_shock_A(i)-r_shock_A(i-1))/ &
      (time_A(i)-time_A(i-1))
    end do
    v_shock_A(1)=v_shock_A(2)
  end subroutine read_shock
  !============================================================================
  subroutine getShock(t)
    real, intent(in) :: t
    real :: time1,time2,FF
    integer :: i
    !--------------------------------------------------------------------------
    if(DoReadShockFile)call read_shock
    ! interpolate
    do i = 1,n
       if(time_A(i) <= t) CYCLE
       if(i == 1)then
          Mach=Mach_A(1)
          s = s_A(1)
          v_shock=v_shock_A(1)
          r_shock=r_shock_A(1)
          V_sw_mod=v_sw_mod_A(1)
       else
          time1=time_A(i-1)
          time2=time_A(i)
          FF = (t-time1)/(time2-time1)
          Mach = Mach_A(i-1) + FF*(Mach_A(i)-Mach_A(i-1))
          s = s_A(i-1) + FF*(s_A(i)-s_A(i-1))
          v_shock = v_shock_A(i-1) + FF*(v_shock_A(i)-v_shock_A(i-1))
          r_shock = r_shock_A(i-1) + FF*(r_shock_A(i)-r_shock_A(i-1))
          V_sw_mod = V_sw_mod_A(i-1) + FF*(V_sw_mod_A(i)-V_sw_mod_A(i-1))
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
end module PT_ModShockPara
!==============================================================================

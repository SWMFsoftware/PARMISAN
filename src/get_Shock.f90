
subroutine getShock(t)
  use PT_ModShockPara
  use PT_ModSWPara, ONLY: n, r_shock_A, v_sw_mod_A
  use PT_ModConst, ONLY:Rsun
  implicit none
  real :: time_A(2000),Mach_A(2000),v_shock_A(2000), &
  th,rshRsun,M,vshkms,time1,time2,FF,t, &
  Vswkms,x_sh_Rsun,y_sh_Rsun,z_sh_Rsun,shock_pos_rel_sun, &
  ss,s_A(2000)
  integer :: ic,i,io
  character(len=:), allocatable ::  name,path
  !----------------------------------------------------------------------------
  data ic/0/
  save
  name = 'PSP_para.dat'
  path = '/Users/xhchen/My_work/Matlab/SEP_LaborDay/mpi/'
  if(ic == 0)then
    ic=1
    open(10,file=trim(path)//trim(name),status='old')
    i=0
    do
      read(10,*,iostat=io) th,rshRsun,M,vshkms,Vswkms,ss
      if (io /= 0) then
        EXIT
      end if
      shock_pos_rel_sun = rshRsun*Rsun
      i=i+1
      time_A(i)=th*3600.d0               ! time in seconds
      r_shock_A(i)=shock_pos_rel_sun     ! shock radius (from center of Sun) in cm
      Mach_A(i) = M                     ! shock Mach number
      v_shock_A(i) = vshkms*1.d5         ! shock speed in cm/s
      v_sw_mod_A(i) = Vswkms*1.d+5       ! solar wind speed at the shock front
      if (ss < 1.d0) then
        ss = 1.d0
      end if
      s_A(i) = ss
    end do
    n=i
    tmax_data = time_A(n)
    tmin_data = time_A(1)
    Rmax_data = r_shock_A(n)
    do i = 2,n
      v_shock_A(i)=(r_shock_A(i)-r_shock_A(i-1))/ &
      (time_A(i)-time_A(i-1))
    end do
    v_shock_A(1)=v_shock_A(2)
  end if

  ! interpolate
  do i = 1,n
   if(time_A(i) > t) EXIT
  end do
  if(i == 1)then
    Mach=Mach_A(1)
    s = s_A(1)
    v_shock=v_shock_A(1)
    r_shock=r_shock_A(1)
    V_sw_mod=v_sw_mod_A(1)
    RETURN
  end if
  time1=time_A(i-1)
  time2=time_A(i)
  FF = (t-time1)/(time2-time1)
  Mach = Mach_A(i-1) + FF*(Mach_A(i)-Mach_A(i-1))
  s = s_A(i-1) + FF*(s_A(i)-s_A(i-1))
  v_shock = v_shock_A(i-1) + FF*(v_shock_A(i)-v_shock_A(i-1))
  r_shock = r_shock_A(i-1) + FF*(r_shock_A(i)-r_shock_A(i-1))
  V_sw_mod = V_sw_mod_A(i-1) + FF*(V_sw_mod_A(i)-V_sw_mod_A(i-1))
end subroutine getShock
!==============================================================================


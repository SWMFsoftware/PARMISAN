
subroutine getVsw(r,VwP)
  use PT_ModShockPara
  use PT_ModSWPara
  use PT_ModConst, ONLY:Rsun
  implicit none
  real :: r,VwP
  real :: rs1,rs2,FF
  integer :: i

  do i = 1,n
     if(r.ge.r_shock_A(i))CYCLE
     if(i.eq.1)then
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

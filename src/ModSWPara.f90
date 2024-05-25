  !  Copyright (C) 2002 Regents of the University of Michigan,
  !  portions used with permission
  !  For more information, see http://csem.engin.umich.edu/tools/swmf
module PT_ModSWPara
  implicit none
  SAVE
  real :: v_sw_mod_A(2000),r_shock_A(2000)
  integer :: n
contains
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
end module PT_ModSWPara
!==============================================================================

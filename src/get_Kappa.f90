subroutine getK(r,t,p,K,dKdr)
  use PT_ModShockPara
  implicit none
  real :: r,t,p,K,dKdr
  save

  K = K0*((r/Rmin)**1.17)
  dKdr = 1.17*K/r
  K=K*((p/p0)**(2.d0*dlt))
  dKdr=dKdr*((p/p0)**(2.d0*dlt))

return
end subroutine getK
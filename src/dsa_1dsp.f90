program dsa_1D_spherical

  ! solves the Parker transport equation using stochastic integration method
  ! 1D radial geometry is assumed
  !
  ! plasma velocity (radial) and diffusion coefficient are pre-defined
  ! (kinematic)
  !
  ! particle splitting is used
  !
  ! this version adapts the particle advance using the "predictor corrector"
  ! approach
  !
  ! this version assumes the source to have a 1/r^2 dependence
  !
  ! this is an mpi version

  use ModMpi
  use PT_ModConst
  use PT_ModPlot, ONLY : set_bins, write_shock_file, NamePath, name_sl, &
       put_flux_contribution, set_diagnostics, put_diagnostic_contribution,&
       write_diagnostics, save_fluxes
  use PT_ModShockPara, ONLY: getshock, getU, read_shock, tmin_data, &
       rmax_data, Mach, s, v_shock, V_sw_mod, r_shock, rMin, p0
  use PT_ModKappa, ONLY: getK, set_kappa
  use PT_ModProc
  implicit none
  logical ::  UseSplit = .true.

  ! real :: rinj What is this?
  real :: r, rs, rs0, t, p, U, dUdr, K, dKdr
  real :: E0, dt, E, w
  real :: rp,pp,divU,xi,rH,pH,rn
  real, parameter :: OneThird = 1.0/3.0
  real :: E_split_lev(100),E_split_lev_min,dEoE_split,dEl_split, &
       E_split_levL,r_save_split(100),t_save_split(100), &
       p_save_split(100),weight,ri,E_split_lev_max
  real :: tmax,tmin,Rmax,Rshmax
  integer :: time_current,time0 ! ,idt,idE
  integer :: n,npart,n_split_levels,lev,iSplit,iSplitCounter,isp_max, &
       lev_save_split(100)
  integer :: iError,seed(630)
  ! mpi initialization routines
  !----------------------------------------------------------------------------
  call mpi_init(iError)
  call mpi_comm_size( mpi_comm_world, nProc, iError)
  call mpi_comm_rank( mpi_comm_world, iProc, iError)
  seed(630) = 123 + iProc
  call random_seed( put=seed )
  call random_number( xi )

  ! definititions, simulation paramters (cgs units)
  Rmin = 1.1d0*Rsun
  dt = 0.01
  E0 = 50.d0*keV                ! source energy
  p0 = sqrt(2.d0*mp*E0)
  ! particle splitting
  n_split_levels = 40            ! total number of energy levels
  E_split_lev_min = 1.d0*MeV    ! energy of first split level
  E_split_lev_max = 20000.d0*MeV  ! energy of last split level
  isp_max=80             ! max number of split particles for a "mother"
  E_split_lev(1)=E_split_lev_min
  E_split_lev(n_split_levels+1)=E_split_lev_max
  dEL_split=log(E_split_lev(n_split_levels+1)/E_split_lev(1)) &
       /real(n_split_levels)

  if(iProc==0)then
     open(45,file=NamePath//name_sl,status='unknown')
     write(45,*)1, E_split_lev(1)/MeV
  end if
  do lev = 2,n_split_levels
     E_split_lev(lev) = E_split_lev(1)*exp(deL_split*(real(lev)-1))
     if(iProc==0)write(45,*)lev,E_split_lev(lev)/MeV
  end do
  if(iProc==0)close(45)
  call read_shock
  call set_kappa
  tmax=1.3*3600.d0
  tmin=tmin_data
  if(iProc==0)call write_shock_file(tmin, tmax, dt)
  Rmax  = Rmax_data
  Rshmax= 16.d0*Rsun
  ! this stuff, between the dashes, is for the diagnostics used
  !------------------------------------------
  call set_diagnostics(E0)

  ! set up bin matrix
  call set_bins(tmin, tmax, Emin = keV, Emax = 1000.0*MeV)

  time0 = mpi_wtime()

  ! number of particles (unsplit, original particles = "mothers")
  npart = 1000000

  ! particle loop
  do n = 1,npart
     if(iProc==0.and.mod(n,100)==0)write(*,*)n
     iSplit=0
     SPLITLOOP:do
        if( iSplit==0 )then
           ! to get a source that falls as 1/r^2, we pick a time randomly
           ! between 0 and tmax and assume the particle is released
           ! at the shock.
           ! Because the source flux depends on the ambient density (1/r^2)
           ! and the bulk flow speed (constant), and the shock surface
           ! expands as r^2, the approach gives the desired result.
           call random_number(rn)
           t = tmin + rn*(tmax - tmin)             ! initial time
           call getShock(t)
           r = r_shock                  ! initial position (at shock)
           rs = r
           ! What is this? rinj = r
           !        r = Rmin + v_shock*t         ! initial position (at shock)
           p = p0                       ! initial momentum
           weight = 1.d0                ! initial weight
           lev = 1
           iSplitCounter = 0
        else
           iSplitCounter = iSplitCounter + 1
           ! Bug (?) fix. Previous version: if(iSplitCounter > isp) continue
           if(iSplitCounter > iSplit) EXIT SPLITLOOP
           t=t_save_split(iSplitCounter)     ! initial time of split particle
           r=r_save_split(iSplitCounter)     ! initial position of   "    "
           rs = r
           ! What is this? rinj = r
           p=p_save_split(iSplitCounter)     ! initial momentum of   "    "
           lev=lev_save_split(iSplitCounter)
        end if
        ! time loop
        TIMELOOP: do
           t = t + dt
           ! stop particle if max time reached
           if(t > tmax) EXIT

           ! current shock location
           !         r_shock=Rmin+v_shock*t
           call getShock(t)

           ! save particle position
           rp=r
           pp=p

           ! advance to the next position and momentum with
           ! stochastic integration method
           ! first step is the prediction step
           call getU(rp,t,U,dUdr)
           call getK(rp,pp,K,dKdr)
           divU=2.d0*U/rp+dUdr
           call random_number(rn)
           xi=-1.d0+2.d0*rn
           rH = rp + xi*sqrt(6.d0*K*dt) + U*dt + (dKdr+2.d0*K/rp)*dt
           pH = pp*(1.d0 - OneThird*divU*dt)

           ! second step is the corrector step
           call getU(rH,t,U,dUdr)
           call getK(rH,pH,K,dKdr)
           divU=2.d0*U/rH+dUdr
           r = rp + xi*sqrt(6.d0*K*dt) + U*dt + (dKdr+2.d0*K/rH)*dt
           p = pp*(1.d0 - OneThird*divU*dt)

           ! particle energy (in ergs)
           E=(p**2)/(2.0*mp)
           w=weight*(2.d0**(-real(lev-1)))
           rs0 = rs
           rs = r
           if (abs(rs-rs0) > 0.5*Rsun) then
              rs0 = rs
           end if
           call put_flux_contribution(rs0, rs, t, E, w)

           time_current =  mpi_wtime() - time0
           if (time_current > (60))then
              call save_fluxes
              call mpi_finalize(iError)
              stop
           end if

           ! outside the box (this is equivalent to "absorbing"
           ! boundary conditions)
           if(r<Rmin .or. r > Rmax .or. r_shock > Rshmax) EXIT TIMELOOP
           if(t > tmax) EXIT TIMELOOP

           ! diagnostics go here, between the dashed lines.
           !----------------BEGIN DIAGNOSTICS ---------------------
           ! this bit bins in momentum those "particles" that are just
           ! behind the shock
           if(r > r_shock-Rsun .and. r<r_shock)&
                call put_diagnostic_contribution(t, p, w)
           !-------------------END DIAGNOSTICS ----------------------------

           ! particle splitting
           if(UseSplit)then
              if(E > E_split_lev(lev).and.lev<n_split_levels &
                   .and.iSplit <= isp_max)then
                 lev=lev+1
                 iSplit=iSplit+1
                 t_save_split(iSplit)=t
                 r_save_split(iSplit)=r
                 p_save_split(iSplit)=p
                 lev_save_split(iSplit)=lev
              end if
           end if
        end do TIMELOOP
        if( lev==1 ) EXIT SPLITLOOP
        ! need to go back and do all split particles
     end do SPLITLOOP
     ! every 100 "mother" particles, dump the diagnostics
     ! here the spectrum at the shock when it crossed 1AU is determined
     ! this is the differential intensity.  The units are 1/(cm^2 s sr MeV).
     ! the normalization is such that the density of energetic particles
     ! computed from this spectrum is n_ep (an input).
     if(mod(n,1)==0)call write_diagnostics
  end do ! particle loop
  ! finialize mpi routine
  call save_fluxes
  call mpi_finalize(iError)
end program dsa_1D_spherical
!==============================================================================


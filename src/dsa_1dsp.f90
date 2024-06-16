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
  use ModKind
  use ModMpi
  use ModReadParam, ONLY: read_file, read_init
  use ModUtilities, ONLY: remove_file, touch_file
  use PT_ModTime,   ONLY: iIter, init_time  => init, PTTime
  use PT_ModMain,   ONLY: &
       IsLastRead, IsStandAlone,          &
       TimeMax, nIterMax, CpuTimeMax,     &
       PT_read_param => read_param,       &
       PT_check      => check,            &
       PT_initialize => initialize,       &
       PT_run        => run,              &
       PT_finalize   => finalize
  use PT_ModGrid,   ONLY: init_stand_alone, init_grid=>init
  use PT_ModConst
  use PT_ModPlot, ONLY : set_bins, write_shock_file, NamePath, name_sl, &
       put_flux_contribution, set_diagnostics, put_diagnostic_contribution,&
       write_diagnostics, save_fluxes
  use PT_ModShockPara, ONLY: getshock, getU, read_shock, tmin_data, &
       rmax_data, Mach, s, v_shock, V_sw_mod, r_shock, rMin, p0
  use PT_ModKappa, ONLY: getK, set_kappa
  use PT_ModParticle, ONLY: nParticleMax, init_particles,&
       nSplitLev, eSplitLev_I
  use PT_ModProc, ONLY: iProc, nProc, iComm, iError

  implicit none

  integer      :: iSession = 1
  real(Real8_) :: CpuTimeStart, CpuTime
  logical      :: IsFirstSession = .true.
  logical ::  UseSplit = .true.

  ! real :: rinj What is this?
  real :: r, rs, rs0, t, p, U, dUdr, K, dKdr
  real :: E0, dt, E, w
  real :: rp,pp,divU,xi,rH,pH,rn
  real, parameter :: OneThird = 1.0/3.0
  real :: r_save_split(100),t_save_split(100), &
       p_save_split(100),weight,ri
  real :: tmax,tmin,Rmax,Rshmax
  integer :: n,lev,iSplit,iSplitCounter,isp_max, &
       lev_save_split(100)
  ! mpi initialization routines
  !----------------------------------------------------------------------------
  call MPI_INIT( iError )
  ! Assign communicator
  iComm = MPI_COMM_WORLD
  call MPI_COMM_RANK( iComm, iProc, iError )
  call MPI_COMM_SIZE( iComm, nProc, iError )

  ! Initialize time which is used to check CPU time
  CpuTimeStart = MPI_WTIME()
  ! Mark the run as a stand alone
  IsStandAlone = .true.
  ! Read PARAM.in file. Provide default restart file for #RESTART
  call read_file('PARAM.in',iComm)

  SESSIONLOOP: do
     call read_init('  ', iSessionIn=iSession)

     if(iProc==0)&
          write(*,*)'----- Starting Session ',iSession,' ------'
     ! Set and check input parameters for this session
     call PT_read_param  ! Identical to SP_set_param('READ')
     call PT_check       ! Similar to SP_set_param('CHECK'), but see init_time

     if(IsFirstSession)then
        ! Time execution (timing parameters set by SP_read_param)
        call init_grid        ! Similar to SP_set_param('GRID')
        call init_stand_alone ! Distinctions from SWMF (CON_bline) version
        call PT_initialize    ! Similar to SP_init_session, NO origin points
        call init_time        ! StartTime, StartTimeJulian from StartTime_I
        ! DataInputTime=0, PTTime=0 unless set in PARAM.in or restart.H
        ! In SWMF, PTTime is set in PT_set_param('CHECK'), DataInput time is
        ! set either in coupling or in reading the MHD data from file.
        ! definititions, simulation paramters (cgs units)
        Rmin = 1.1d0*Rsun
        dt = 0.01
        E0 = 50.d0*keV                ! source energy
        p0 = sqrt(2.d0*mp*E0)

        isp_max=80             ! max number of split particles for a "mother"
        call init_particles
        call read_shock
        call set_kappa
        tmax=720.0
        tmin=tmin_data
        if(iProc==0)call write_shock_file(tmin, tmax, dt)
        Rmax  = Rmax_data
        Rshmax= 16.d0*Rsun
        ! this stuff, between the dashes, is for the diagnostics used
        !------------------------------------------
        call set_diagnostics(E0)

        ! set up bin matrix
        call set_bins(tmin, tmax, Emin = keV, Emax = 1000.0*MeV)
     end if

     ! particle loop
     do n = 1, nParticleMax
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
              !        r = Rmin + v_shock*t   ! initial position (at shock)
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

              CpuTime =  mpi_wtime() - CpuTimeStart
              if (CpuTime > CpuTimeMax)then
                 call save_fluxes
                 call MPI_finalize(iError)
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
                 if(E > eSplitLev_I(lev).and.lev<nSplitLev &
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
     if(IsLastRead)EXIT SESSIONLOOP
     if(iProc==0) &
          write(*,*)'----- End of Session   ',iSession,' ------'
     iSession       = iSession + 1
     IsFirstSession = .false.
  end do SESSIONLOOP
  ! finialize mpi routine
  call save_fluxes
  call MPI_finalize(iError)
end program dsa_1D_spherical
!==============================================================================


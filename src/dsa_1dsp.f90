program dsa_1D_spherical

  ! solves the Parker transport equation using stochastic integration method
  ! 1D radial geometry is assumed
  !
  ! plasma velocity (radial) and diffusion coefficient are pre-defined (kinematic)
  !
  ! particle splitting is used
  !
  ! this version adapts the particle advance using the "predictor corrector" approach
  !
  ! this version assumes the source to have a 1/r^2 dependence
  !
  ! this is an mpi version

  use ModMpi
  use PT_ModConst
  use PT_ModShockPara
  use PT_ModProc
  implicit none
  character(len=4) :: tmp
  character(len=*), parameter ::  SPLIT = 'y'
  character(len=*), parameter ::  NamePath = &
       '/Users/xhchen/My_work/Matlab/SEP_LaborDay/mpi/data/'
  character(len=*), parameter :: file1 = 'Cross8_'
  character(len=*), parameter :: file2 = 'Cross10_'
  character(len=*), parameter :: file3 = 'Cross12_'
  character(len=*), parameter :: file4 = 'Cross14_'
  character(len=*), parameter :: name_sl ='split_levels.dat'
  character(len=*), parameter :: name_pt ='prms_vs_t.dat'
  character(len=*), parameter :: Time_bin ='t_bin.dat'
  character(len=*), parameter :: Energy_bin ='E_bin.dat'
  character(len=*), parameter :: name_dist = 'dist'
  character(len=100)          ::  NameFile
  real :: r,rs,rinj,rs0,t,p,U,dUdr,K,dKdr
  real :: E0,dt,E
  real :: E_bin(1000),t_bin(1000),w1(1000,1000),w2(1000,1000), &
       w3(1000,1000),w4(1000,1000),t_wall,Emax,Emin
  real :: rp,pp,divU,OneThird,xi,rH,pH,rn
  real :: E_split_lev(100),E_split_lev_min,dEoE_split,dEl_split, &
       E_split_levL,r_save_split(100),t_save_split(100), &
       p_save_split(100),weight,ri,E_split_lev_max
  real :: dN(500,1000),dpop_diag,dpl_diag,p_diag_min,dp_diag, &
       ap,fp,fpc,dNs,den_ep,fac,f,j,w,dtd,rtd,ftd,ftdc
  real :: tmax,tmin,Rmax,Rshmax
  integer :: time_current,time0,n_tbin,n_Ebin,i_tbin,i_Ebin,idt,idE
  integer :: n,npart,n_split_levels,lev,isp,ispcnt,isp_max, &
       lev_save_split(100)
  integer :: ipmx,ip,ip1,pmx,np(500),itd,itd1,itdmx,ntd
  integer :: mpierr,seed(630)
  ! mpi initialization routines
  !----------------------------------------------------------------------------
  call mpi_init( mpierr )
  call mpi_comm_size( mpi_comm_world, nProc, mpierr )
  call mpi_comm_rank( mpi_comm_world, iProc, mpierr )
  seed(630) = 123 + iProc
  call random_seed( put=seed )
  call random_number( xi )

  ! constants (cgs units)
  OneThird = 1.d0/3.d0

  ! definititions, simulation paramters (cgs units)
  Rmin = 1.1d0*Rsun
  Vw = 4.d+7
  Ls = Rsun
  dt = 0.01
  drSHOCK = 4.d-5*Rsun
  E0 = 50.d0*keV                ! source energy
  p0 = dsqrt(2.d0*mp*E0)
  K0 = 6.d+19 * ((Rmin/AU)**1.17)   ! K0 is the value at r=Rmin at E0
  dlt = 0.71d0
  t_wall = 5*60

  write(NameFile,'(a,i4.4,a)')NamePath//file1, iProc, '.dat'
  open(901,file=trim(NameFile),status='unknown',action="READWRITE")
  write(NameFile,'(a,i4.4,a)')NamePath//file2, iProc, '.dat'
  open(902,file=trim(NameFile),status='unknown',action="READWRITE")
  write(NameFile,'(a,i4.4,a)')NamePath//file3, iProc, '.dat'
  open(903,file=trim(NameFile),status='unknown',action="READWRITE")
  write(NameFile,'(a,i4.4,a)')NamePath//file4, iProc, '.dat'
  open(904,file=trim(NameFile),status='unknown',action="READWRITE")

  ! particle splitting
  n_split_levels = 40            ! total number of energy levels
  E_split_lev_min = 1.d0*MeV    ! energy of first split level
  E_split_lev_max = 20000.d0*MeV  ! energy of last split level
  isp_max=80             ! max number of split particles for a "mother"
  E_split_lev(1)=E_split_lev_min
  E_split_lev(n_split_levels+1)=E_split_lev_max
  dEL_split=dlog(E_split_lev(n_split_levels+1)/E_split_lev(1)) &
       /real(n_split_levels)

  if(iProc==0)then
     open(45,file=trim(NamePath)//trim(name_sl),status='unknown')
     write(45,*)1,E_split_lev(1)/MeV
  end if
  do lev = 2,n_split_levels
     E_split_lev(lev) = E_split_lev(1)*exp(deL_split*(real(lev)-1))
     if(iProc==0)write(45,*)lev,E_split_lev(lev)/MeV
  end do
  if(iProc==0)close(45)
  call getShock(t)
  tmax=1.3*3600.d0
  tmin=tmin_data
  t=tmin-dt
  if(iProc==0) then
     open(44,file=trim(NamePath)//trim(name_pt),status='unknown')
  end if
  do
     t=t+dt
     if(t > tmax) EXIT
     call getShock(t)
     r=r_shock
     p=p0
     call getK(r,t,p,K,dKdr)
     if(iProc==0)write(44,'(6f18.7,2e15.5)')t/3600.d0,Mach,s, &
          v_shock/1.d+5,V_sw_mod/1.d+5,r_shock/Rsun,K,dKdr
  end do
  if(iProc==0)close(44)
  Rmax=Rmax_data
  Rshmax=16.d0*Rsun

  ! this stuff, between the dashes, is for the diagnostics used
  !------diagnostic initializations -------
  dpop_diag = 0.04
  dpl_diag = dlog(1.d0+dpop_diag)
  p_diag_min = dsqrt(2.d0*100.d0*keV*mp)
  ipmx=1000
  itdmx=500
  dtd=60.d0
  ntd=0
  do itd = 1,itdmx
     np(itd)=0
     do ip = 1,ipmx
        dN(itd,ip)=0.d0
     end do
  end do
  den_ep = 2.4d-3
  !------------------------------------------

  ! set up bin matrix
  n_tbin = 1000
  n_Ebin = 1000
  do i_tbin = 1,n_tbin
     t_bin(i_tbin) = tmin+(tmax-tmin)*(i_tbin-1)/n_tbin
  end do
  Emax = 1000.d0*MeV
  Emin = 1.d0*keV
  do i_Ebin = 1,n_Ebin
     E_bin(i_Ebin) = dexp(dlog(Emin)+(dlog(Emax)-dlog(Emin)) &
          /n_Ebin*(i_Ebin-1))
  end do
  if(iProc==0) then
     open(801,file=trim(NamePath)//trim(Energy_bin),status='unknown')
     write(801,'(1000e15.6)')E_bin/keV
     close(801)
     open(802,file=trim(NamePath)//trim(Time_bin),status='unknown')
     write(802,'(1000e15.6)')t_bin
     close(802)
  end if

  time0 = mpi_wtime()

  ! number of particles (unsplit, original particles = "mothers")
  npart=1000000

  ! particle loop
  do n = 1,npart
     if(iProc==0.and.mod(n,100)==0)write(*,*)n
     isp=0
     do
        if( isp==0 )then
           ! to get a source that falls as 1/r^2, we pick a time randomly
           ! between 0 and tmax and assume the particle is released
           ! at the shock.
           ! Because the source flux depends on the ambient density (1/r^2)
           ! and the bulk flow speed (constant), and the shock surface
           ! expands as r^2, the approach gives the desired result.
           call random_number(rn)
           t = tmin+rn*(tmax-tmin)             ! initial time
           call getShock(t)
           r = r_shock                  ! initial position (at shock)
           rs = r
           rinj = r
           !        r = Rmin + v_shock*t         ! initial position (at shock)
           p = p0                       ! initial momentum
           weight = 1.d0                ! initial weight
           lev = 1
           ispcnt = 0
        else
           ispcnt=ispcnt+1
           if(ispcnt > isp) continue
           t=t_save_split(ispcnt)     ! initial time of split particle
           r=r_save_split(ispcnt)     ! initial position of   "    "
           rs = r
           rinj = r
           p=p_save_split(ispcnt)     ! initial momentum of   "    "
           lev=lev_save_split(ispcnt)
        end if

        ! time loop
        do
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
           call getK(rp,t,pp,K,dKdr)
           divU=2.d0*U/rp+dUdr
           call random_number(rn)
           xi=-1.d0+2.d0*rn
           rH = rp + xi*dsqrt(6.d0*K*dt) + U*dt + (dKdr+2.d0*K/rp)*dt
           pH = pp*(1.d0-OneThird*divU*dt)

           ! second step is the corrector step
           call getU(rH,t,U,dUdr)
           call getK(rH,t,pH,K,dKdr)
           divU=2.d0*U/rH+dUdr
           r = rp + xi*dsqrt(6.d0*K*dt) + U*dt + (dKdr+2.d0*K/rH)*dt
           p = pp*(1.d0-OneThird*divU*dt)

           ! particle energy (in ergs)
           E=(p**2)/(2.d0*mp)

           rs0 = rs
           rs = r
           if (abs(rs-rs0) > 0.5*Rsun) then
              rs0 = rs
           end if

           if (((rs-8.d0*Rsun)*(rs0-8.d0*Rsun))<0.0) then
              idt = minloc(abs(t_bin-t),DIM=1)
              idE = minloc(abs(E_bin-E),DIM=1)
              w1(idt,idE) = w1(idt,idE)+w
           end if
           if (((rs-10.d0*Rsun)*(rs0-10.d0*Rsun))<0.0) then
              idt = minloc(abs(t_bin-t),DIM=1)
              idE = minloc(abs(E_bin-E),DIM=1)
              w2(idt,idE) = w2(idt,idE)+w
           end if
           if (((rs-12.d0*Rsun)*(rs0-12.d0*Rsun))<0.0) then
              idt = minloc(abs(t_bin-t),DIM=1)
              idE = minloc(abs(E_bin-E),DIM=1)
              w3(idt,idE) = w3(idt,idE)+w
           end if
           if (((rs-14.d0*Rsun)*(rs0-14.d0*Rsun))<0.0) then
              idt = minloc(abs(t_bin-t),DIM=1)
              idE = minloc(abs(E_bin-E),DIM=1)
              w4(idt,idE) = w4(idt,idE)+w
           end if

           time_current =  mpi_wtime()-time0
           if (time_current > (60))then
              do i_tbin = 1,n_tbin
                 write(901,'(1000e15.6)')w1(i_tbin,:)
                 write(902,'(1000e15.6)')w2(i_tbin,:)
                 write(903,'(1000e15.6)')w3(i_tbin,:)
                 write(904,'(1000e15.6)')w4(i_tbin,:)
              end do
              close(901)
              close(902)
              close(903)
              close(904)
              call mpi_finalize( mpierr )
              stop
           end if

           ! outside the box (this is equivalent to "absorbing"
           ! boundary conditions)
           if(r<Rmin .or. r > Rmax .or. r_shock > Rshmax) EXIT
           if(t > tmax) EXIT

           ! diagnostics go here, between the dashed lines.
           !----------------BEGIN DIAGNOSTICS ---------------------
           ! this bit bins in momentum those "particles" that are just
           ! behind the shock
           if(r > r_shock-Rsun .and. r<r_shock)then
              ap=dlog(p/p_diag_min)/dpl_diag+0.50000000001
              ip=ap
              ip1=ip+1
              fp=ap-ip
              fpc=1.d0-fp
              ip=min0(ipmx,max0(1,ip))
              ip1=min0(ipmx,max0(1,ip1))
              rtd=t/dtd + 0.5000000000001
              itd=rtd
              itd1=itd+1
              ftd=rtd-itd
              ftdc=1.d0-ftd
              itd=min0(itdmx,max0(1,itd))
              itd1=min0(itdmx,max0(1,itd1))
              ntd=max0(ntd,itd)
              np(itd)=max0(np(itd),ip)
              w=weight*(2.d0**(-real(lev-1)))
              dN(itd,ip)=dN(itd,ip)+w*fpc*ftdc
              dN(itd,ip1)=dN(itd,ip1)+w*fp*ftdc
              dN(itd1,ip)=dN(itd1,ip)+w*fpc*ftd
              dN(itd1,ip1)=dN(itd1,ip1)+w*fp*ftd
           end if
           !-------------------END DIAGNOSTICS ----------------------------

           ! particle splitting
           if( SPLIT=='Y' .or. SPLIT=='y')then
              if(E > E_split_lev(lev).and.lev<n_split_levels &
                   .and.isp <= isp_max)then
                 lev=lev+1
                 isp=isp+1
                 t_save_split(isp)=t
                 r_save_split(isp)=r
                 p_save_split(isp)=p
                 lev_save_split(isp)=lev
              end if
           end if

        end do
        ! end of time loop

        if( lev==1 ) EXIT   ! need to go back and do all split particles
     end do
     ! every 100 "mother" particles, dump the diagnostics
     ! here the spectrum at the shock when it crossed 1AU is determined
     ! this is the differential intensity.  The units are 1/(cm^2 s sr MeV).
     ! the normalization is such that the density of energetic particles
     ! computed from this spectrum is n_ep (an input).
     if(mod(n,1)==0)then
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
              !         p=p_diag_min*dexp( dpl_diag*real(ip-1))
              p=p_diag_min*dexp( dpl_diag*(real(ip)-0.5) )
              E=(p**2)/(2.d0*mp)
              dp_diag=dpop_diag*p
              ! distribution function:
              f=dN(itd,ip)/((p**2)*dp_diag)*den_ep*fac
              j=f*(p**2)                               ! diff. intensity
              write(14,*)E/MeV,dmax1(1.d-20,j*MeV)
              ! do this again, shifted by a bit to make it histogram style
              !         p=p_diag_min*dexp( dpl_diag*real(ip))
              !         E=(p**2)/(2.d0*mp)
              !         dp_diag=dpop_diag*p
              !         f=dN(itd,ip)/((p**2)*dp_diag)*den_ep*fac
              !         j=f*(p**2)
              !         write(14,*)E/MeV,dmax1(1.d-20,j*MeV)
           end do
        end do
        close(14)
     end if
  end do ! particle loop

  ! finialize mpi routine
  close(901)
  close(902)
  close(903)
  close(904)
  call mpi_finalize( mpierr )
  stop
end program dsa_1D_spherical
!==============================================================================


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
    use PT_ModShockPara, ONLY: getshock, getU, tmax_data, tmin_data, &
    rmax_data, rmin_data, Mach, s, v_shock, V_sw_mod, r_shock, rMin, &
    p0, K0
    use PT_ModKappa, ONLY: getK
    use PT_ModProc
    implicit none
    logical ::  UseSplit = .true.
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
    real :: r(1000000),t(1000000),p(1000000)
    real :: E0,dt,E,U,dUdr,K,dKdr
    real :: E_bin(1000),t_bin(1000),w1(1000,1000),w2(1000,1000), &
    w3(1000,1000),w4(1000,1000)
    real :: tp,rp,pp,divU,xi,rH,pH,rn
    real, parameter :: OneThird = 1.0/3.0
    real :: dN(500,1000),dpop_diag,dpl_diag,p_diag_min,dp_diag, &
    ap,fp,fpc,dNs,den_ep,fac,f,j,w,dtd,rtd,ftd,ftdc
    real :: tmax,tmin,Rmin,Rmax,Emax,Emin,t_wall
    integer :: nParticle,dnParticle,nParticleMax,iParticle,itime, &
    dt_couple,dt_diff,n_diff
    integer :: time_current,time0,n_tbin,n_Ebin,i_tbin,i_Ebin,idt,idE
    integer :: n,npart,ipmx,ip,ip1,pmx,np(500),itd,itd1,itdmx,ntd
    integer :: mpierr,seed(630)
    ! mpi initialization routines
    !----------------------------------------------------------------------------
    call mpi_init( mpierr )
    call mpi_comm_size( mpi_comm_world, nProc, mpierr )
    call mpi_comm_rank( mpi_comm_world, iProc, mpierr )
    seed(630) = 123 + iProc
    call random_seed( put=seed )
    call random_number( xi )

    ! definititions, simulation paramters (cgs units)
    dt_diff = 0.05    ! diffusion time step (Second)
    dt_couple = 30.0 ! MHD coupling time step (Second)
    n_diff = int(dt_couple/dt_diff) ! Diffusions at each coupling time step
    E0 = 50.d0*keV                ! source energy
    p0 = sqrt(2.d0*mp*E0)
    K0 = 6.d+19 * ((Rmin/AU)**1.17)   ! K0 is the value at r=Rmin at E0
    t_wall = 5*60

    write(NameFile,'(a,i4.4,a)')NamePath//file1, iProc, '.dat'
    open(901,file=trim(NameFile),status='unknown',action="READWRITE")
    write(NameFile,'(a,i4.4,a)')NamePath//file2, iProc, '.dat'
    open(902,file=trim(NameFile),status='unknown',action="READWRITE")
    write(NameFile,'(a,i4.4,a)')NamePath//file3, iProc, '.dat'
    open(903,file=trim(NameFile),status='unknown',action="READWRITE")
    write(NameFile,'(a,i4.4,a)')NamePath//file4, iProc, '.dat'
    open(904,file=trim(NameFile),status='unknown',action="READWRITE")

    call getShock(0.0)
    tmax=tmax_data
    tmin=tmin_data
    Rmax=Rmax_data
    Rmin=Rmin_data

    ! this stuff, between the dashes, is for the diagnostics used
    !------diagnostic initializations -------
    dpop_diag = 0.04
    dpl_diag = log(1.d0+dpop_diag)
    p_diag_min = sqrt(2.d0*E0*mp)
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
    Emax = 2000.d0*MeV
    Emin = 1.d0*keV
    do i_Ebin = 1,n_Ebin
      E_bin(i_Ebin) = exp(log(Emin)+(log(Emax)-log(Emin)) &
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


    nParticleMax = 1000000 ! Total particles
    nParticle = 0 
    dnParticle = int(nParticleMax*dt_couple/(tmax-tmin)) ! Injected particles at each time step
    TIMELOOP: do t_couple = tmin,tmax,dt_couple
      if(iProc==0.and.mod(t,60)==0)write(*,*)int(t_couple/60) ! For monitoring the progress 
      ! Inject new particles
      do iParticle = nParticle+1,nParticle+dnParticle
        call random_number(rn)   
        t(iParticle) = t_couple+dt_couple*rn
        call getShock(t(iParticle))
        r(iParticle) = r_shock   ! initial position (at shock)
        p(iParticle) = p0        ! initial momentum
      end do
      nParticle = nParticle+dnParticle
      ParticleLOOP: do iParticle = 1,nParticle
        ! outside the box (this is equivalent to "absorbing" boundary conditions)
        if (r(iParticle)<Rmin .or. r(iParticle)>Rmax) continue
        itime = 1
        DIFFUSIONLOOP: do
          ! save particle status
          rp = r(iParticle)
          pp = p(iParticle)
          tp = t(iParticle)
          ! advance to the next position and momentum with
          ! stochastic integration method
          ! first step is the prediction step
          call getU(rp,tp,U,dUdr)
          call getK(rp,tp,pp,K,dKdr)
          divU=2.d0*U/rp+dUdr
          call random_number(rn)
          xi=-1.d0+2.d0*rn
          rH = rp+xi*sqrt(6.d0*K*dt_diff)+U*dt_diff+(dKdr+2.d0*K/rp)*dt_diff
          pH = pp*(1.d0-OneThird*divU*dt_diff)
          ! second step is the corrector step
          call getU(rH,tp,U,dUdr)
          call getK(rH,tp,pH,K,dKdr)
          divU=2.d0*U/rH+dUdr
          r(iParticle) = rp+xi*sqrt(6.d0*K*dt_diff)+U*dt_diff+(dKdr+2.d0*K/rH)*dt_diff
          p(iParticle) = pp*(1.d0-OneThird*divU*dt_diff)
          t(iParticle) = t(iParticle)+itime*dt_diff
          E=(p(iParticle)**2)/(2.0*mp) ! particle energy (in ergs)
          if (t(iParticle).ge.(t_couple+dt_couple)) EXIT DIFFUSIONLOOP
          if (r(iParticle)<Rmin .or. r(iParticle)>Rmax) EXIT DIFFUSIONLOOP
          ! Bin Particles at some chosen radius, e.g., SC position
          if (((r(iParticle)-14.d0*Rsun)*(rp-14.d0*Rsun))<0.0) then
            idt = minloc(abs(t_bin-t_real),DIM=1)
            idE = minloc(abs(E_bin-E),DIM=1)
            w1(idt,idE) = w1(idt,idE)+1
          end if
          itime = itime+1
        end do DIFFUSIONLOOP
      end do PARTICLELOOP
      ! diagnostics go here, between the dashed lines.
      !----------------BEGIN DIAGNOSTICS ---------------------
      ! this bit bins in momentum those "particles" that are just
      ! behind the shock
      if(r > r_shock-Rsun .and. r<r_shock)then
        ap=log(p/p_diag_min)/dpl_diag+0.50000000001
        ip=ap
        ip1=ip+1
        fp=ap-ip
        fpc=1.d0-fp
        ip=min0(ipmx,max0(1,ip))
        ip1=min0(ipmx,max0(1,ip1))
        rtd=t_real/dtd + 0.5000000000001
        itd=rtd
        itd1=itd+1
        ftd=rtd-itd
        ftdc=1.d0-ftd
        itd=min0(itdmx,max0(1,itd))
        itd1=min0(itdmx,max0(1,itd1))
        ntd=max0(ntd,itd)
        np(itd)=max0(np(itd),ip)
        w=1
        dN(itd,ip)=dN(itd,ip)+w*fpc*ftdc
        dN(itd,ip1)=dN(itd,ip1)+w*fp*ftdc
        dN(itd1,ip)=dN(itd1,ip)+w*fpc*ftd
        dN(itd1,ip1)=dN(itd1,ip1)+w*fp*ftd
      end if
      !-------------------END DIAGNOSTICS ----------------------------

      ! Save data with a wall-time
      time_current =  mpi_wtime()-time0
      if (time_current > t_wall)then
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
            !         p=p_diag_min*exp( dpl_diag*real(ip-1))
            p=p_diag_min*exp( dpl_diag*(real(ip)-0.5) )
            E=(p**2)/(2.d0*mp)
            dp_diag=dpop_diag*p
            ! distribution function:
            f=dN(itd,ip)/((p**2)*dp_diag)*den_ep*fac
            j=f*(p**2)                               ! diff. intensity
            write(14,*)E/MeV,dmax1(1.d-20,j*MeV)
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
      end if

    end do TIMELOOP

    ! Save data
    do i_tbin = 1,n_tbin
      write(901,'(1000e15.6)')w1(i_tbin,:)
      write(902,'(1000e15.6)')w2(i_tbin,:)
      write(903,'(1000e15.6)')w3(i_tbin,:)
      write(904,'(1000e15.6)')w4(i_tbin,:)
    end do
    ! finialize mpi routine
    close(901)
    close(902)
    close(903)
    close(904)
    call mpi_finalize( mpierr )
    stop
end program dsa_1D_spherical
!==============================================================================


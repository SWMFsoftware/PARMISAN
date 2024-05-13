      program dsa_1D_spherical

c solves the Parker transport equation using stochastic integration method
c 1D radial geometry is assumed
c
c plasma velocity (radial) and diffusion coefficient are pre-defined (kinematic)
c
c particle splitting is used
c
c this version adapts the particle advance using the "predictor corrector" approach
c
c this version assumes the source to have a 1/r^2 dependence
c
c this is an mpi version
c

      use mpi

      implicit none
      character*1 SPLIT
      character*4 tmp
      character*120 fname, path
      character*100 name1,name2,name3,name4,name5,name6,name7,name8,
     & name9,name10,name11,name12,name13
      real*8 r,rs,rinj,rs0,t,p,U,dUdr,K,dKdr
      real*8 Rsun,AU,eV,keV,MeV,mp,fourpi
      real*8 E0,Rmin,K0,p0,dt,dlt,E
      real*8 E_bin(1000),t_bin(1000),w1(1000,1000),w2(1000,1000),
     & w3(1000,1000),w4(1000,1000),t_wall,Emax,Emin
      real*8 rp,pp,divU,OneThird,xi,rH,pH,xrand
      real*8 r_shock,v_shock,Vw,drSHOCK,Ls,s,Mach,V_sw_mod
      real*8 E_split_lev(100),E_split_lev_min,dEoE_split,dEl_split,
     &     E_split_levL,r_save_split(100),t_save_split(100),
     &     p_save_split(100),weight,ri,E_split_lev_max
      real*8 dN(500,1000),dpop_diag,dpl_diag,p_diag_min,dp_diag,
     &     ap,fp,fpc,dNs,den_ep,fac,f,j,w,dtd,rtd,ftd,ftdc
      real*8 tmax,tmin,Rmax,tmax_data,tmin_data,Rmax_data,Rshmax
      integer*4 time_current,time0,n_tbin,n_Ebin,i_tbin,i_Ebin,idt,idE
      integer*4 n,npart,n_split_levels,lev,isp,ispcnt,isp_max,
     &     lev_save_split(100)
      integer*4 ipmx,ip,ip1,pmx,np(500),itd,itd1,itdmx,ntd
      integer*4 mype,npe,mpierr,seed(33)

      common /params/ r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt,
     &     Rsun,Mach,V_sw_mod
      common /max_sim_time/ tmax_data,tmin_data,Rmax_data

c mpi initialization routines
      call mpi_init( mpierr )
      call mpi_comm_size( mpi_comm_world, npe, mpierr )
      call mpi_comm_rank( mpi_comm_world, mype, mpierr )
      seed(33) = 123 + mype
      call random_seed( put=seed )
      call random_number( xi )

c constants (cgs units)
      Rsun = 6.96d+10
      AU = 1.496d+13
      eV = 1.6022d-12
      keV = 1000.d0*eV
      MeV = 1.d+6*eV
      mp = 1.6726d-24
      OneThird = 1.d0/3.d0
      fourpi = 8.d0*dasin(1.d0)

c definititions, simulation paramters (cgs units)
      Rmin = 1.1d0*Rsun
      Vw = 4.d+7
c      Ls = 0.025*AU
      Ls = Rsun
      dt = 0.01
      drSHOCK = 4.d-5*Rsun
      E0 = 50.d0*keV                ! source energy
      p0 = dsqrt(2.d0*mp*E0)
      K0 = 6.d+19 * ((Rmin/AU)**1.17)   ! K0 is the value at r=Rmin at E0
      dlt = 0.71d0                   
      t_wall = 5*60
      path =  '/pscratch/sd/x/xhchen/PT1D/Laborday/PSP/data/'
      name1 = 'Cross8_'
      name2 = 'Cross10_'
      name3 = 'Cross12_'
      name4 = 'Cross14_'



      if(mype.le.9)then
         write(tmp,'(i1)')mype
         fname=trim(name1)//tmp(1:1)//'.dat'
         open(901,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name2)//tmp(1:1)//'.dat'
         open(902,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name3)//tmp(1:1)//'.dat'
         open(903,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name4)//tmp(1:1)//'.dat'
         open(904,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
      elseif( mype.ge.10 .and. mype.lt.100)then
         write(tmp,'(i2)')mype
         fname=trim(name1)//tmp(1:2)//'.dat'
         open(901,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name2)//tmp(1:2)//'.dat'
         open(902,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name3)//tmp(1:2)//'.dat'
         open(903,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name4)//tmp(1:2)//'.dat'
         open(904,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
      elseif( mype.ge.100 .and. mype.lt.1000)then
        write(tmp,'(i3)')mype
         fname=trim(name1)//tmp(1:3)//'.dat'
         open(901,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name2)//tmp(1:3)//'.dat'
         open(902,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name3)//tmp(1:3)//'.dat'
         open(903,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name4)//tmp(1:3)//'.dat'
         open(904,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
      else
         write(tmp,'(i4)')mype
         fname=trim(name1)//tmp(1:4)//'.dat'
         open(901,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name2)//tmp(1:4)//'.dat'
         open(902,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name3)//tmp(1:4)//'.dat'
         open(903,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
         fname=trim(name4)//tmp(1:4)//'.dat'
         open(904,file=trim(path)//trim(fname),status='unknown'
     &     ,action="READWRITE")
      end if

c particle splitting
      SPLIT = 'y'
      n_split_levels = 40            ! total number of energy levels
      E_split_lev_min = 1.d0*MeV    ! energy of first split level
      E_split_lev_max = 20000.d0*MeV  ! energy of last split level
      isp_max=80             ! max number of split particles for a "mother"
      E_split_lev(1)=E_split_lev_min
      E_split_lev(n_split_levels+1)=E_split_lev_max
      dEL_split=dlog(E_split_lev(n_split_levels+1)/E_split_lev(1))
     &        /  dfloat(n_split_levels)
     
      name11 ='split_levels.dat'
      if(mype.eq.0)then
       open(45,file=trim(path)//trim(name11),status='unknown')
       write(45,*)1,E_split_lev(1)/MeV
      end if
      do 10 lev=2,n_split_levels
       E_split_lev(lev) = E_split_lev(1)*dexp(deL_split*(dfloat(lev)-1))
       if(mype.eq.0)write(45,*)lev,E_split_lev(lev)/MeV
 10   continue
      if(mype.eq.0)close(45)

c maximum simulation time and where the shock is at this time 
c which defines the max. source location
c      Rmax = 1.4*AU
c      tmax = (Rmax-Rmin)/v_shock_0
c      tmax = 12.d0*3600.d0
c      Rmax = 5.d0*AU

       call getShock(t)
       tmax=1.3*3600.d0
       tmin=tmin_data
       t=tmin-dt
       name12 ='prms_vs_t.dat'
       if(mype.eq.0) then
        open(44,file=trim(path)//trim(name12),status='unknown')
       end if
 777   t=t+dt
        if(t.gt.tmax)goto 778
        call getShock(t)
        r=r_shock
        p=p0
        call getK(r,t,p,K,dKdr)
        if(mype.eq.0)write(44,'(6f18.7,2e15.5)')t/3600.d0,Mach,s,
     &     v_shock/1.d+5,V_sw_mod/1.d+5,r_shock/Rsun,K,dKdr
        goto 777
 778   continue
       if(mype.eq.0)close(44)
       Rmax=Rmax_data
       Rshmax=16.d0*Rsun

       

c this stuff, between the dashes, is for the diagnostics used
c------diagnostic initializations -------
      dpop_diag = 0.04
      dpl_diag = dlog(1.d0+dpop_diag)
      p_diag_min = dsqrt(2.d0*100.d0*keV*mp)
      ipmx=1000
      itdmx=500
      dtd=60.d0
      ntd=0
      do 120 itd=1,itdmx
      np(itd)=0
      do 12 ip=1,ipmx
      dN(itd,ip)=0.d0
 12    continue
 120    continue
      den_ep = 2.4d-3
c------------------------------------------

c set up bin matrix
      n_tbin = 1000
      n_Ebin = 1000
      do i_tbin = 1,n_tbin
        t_bin(i_tbin) = tmin+(tmax-tmin)*(i_tbin-1)/n_tbin
      end do
      Emax = 1000.d0*MeV
      Emin = 1.d0*keV
      do i_Ebin = 1,n_Ebin
        E_bin(i_Ebin) = dexp(dlog(Emin)+(dlog(Emax)-dlog(Emin))
     &       /n_Ebin*(i_Ebin-1))
      end do
      if(mype.eq.0) then
        name5 ='E_bin.dat'
        open(801,file=trim(path)//trim(name5),status='unknown')
        write(801,'(1000e15.6)')E_bin/keV
        close(801)
        name6 ='t_bin.dat'
        open(802,file=trim(path)//trim(name6),status='unknown')
        write(802,'(1000e15.6)')t_bin
        close(802)
      end if


      time0 = mpi_wtime()


c number of particles (unsplit, original particles = "mothers")
      npart=1000000

c particle loop
      do 100 n=1,npart
       if(mype.eq.0.and.mod(n,100).eq.0)print*,n
       isp=0
 15    continue
       if( isp.eq.0 )then
c to get a source that falls as 1/r^2, we pick a time randomly 
c between 0 and tmax and assume the particle is released at the shock.  
c Because the source flux depends on the ambient density (1/r^2) and 
c the bulk flow speed (constant), and the shock surface expands as r^2, 
c the approach gives the desired result.
        t = tmin+xrand()*(tmax-tmin)             ! initial time
        call getShock(t)
        r = r_shock                  ! initial position (at shock)
        rs = r
        rinj = r
c        r = Rmin + v_shock*t         ! initial position (at shock)
        p = p0                       ! initial momentum
        weight = 1.d0                ! initial weight
        lev = 1
        ispcnt = 0
       else
        ispcnt=ispcnt+1
        if(ispcnt.gt.isp)goto 35
        t=t_save_split(ispcnt)     ! initial time of split particle
        r=r_save_split(ispcnt)     ! initial position of   "    "
        rs = r
        rinj = r
        p=p_save_split(ispcnt)     ! initial momentum of   "    "
        lev=lev_save_split(ispcnt)
       end if


c time loop
 20    t = t + dt

c stop particle if max time reached
         if(t.gt.tmax)goto 30

c current shock location
c         r_shock=Rmin+v_shock*t
         call getShock(t)

c save particle position
         rp=r
         pp=p

c advance to the next position and momentum with stochastic integration method
c first step is the prediction step
         call getU(rp,t,U,dUdr) 
         call getK(rp,t,pp,K,dKdr)
         divU=2.d0*U/rp+dUdr
         xi=-1.d0+2.d0*xrand()
         rH = rp + xi*dsqrt(6.d0*K*dt) + U*dt + (dKdr+2.d0*K/rp)*dt
         pH = pp*(1.d0-OneThird*divU*dt)

c second step is the corrector step
         call getU(rH,t,U,dUdr) 
         call getK(rH,t,pH,K,dKdr)
         divU=2.d0*U/rH+dUdr
         r = rp + xi*dsqrt(6.d0*K*dt) + U*dt + (dKdr+2.d0*K/rH)*dt
         p = pp*(1.d0-OneThird*divU*dt)

c particle energy (in ergs)
         E=(p**2)/(2.d0*mp)

        rs0 = rs
        rs = r
        if (abs(rs-rs0).gt.0.5*Rsun) then
            rs0 = rs
        end if





        if (((rs-8.d0*Rsun)*(rs0-8.d0*Rsun)).lt.0.0) then
            idt = minloc(abs(t_bin-t),DIM=1)
            idE = minloc(abs(E_bin-E),DIM=1)
            w1(idt,idE) = w1(idt,idE)+w
        end if
        if (((rs-10.d0*Rsun)*(rs0-10.d0*Rsun)).lt.0.0) then
            idt = minloc(abs(t_bin-t),DIM=1)
            idE = minloc(abs(E_bin-E),DIM=1)
            w2(idt,idE) = w2(idt,idE)+w
        end if
        if (((rs-12.d0*Rsun)*(rs0-12.d0*Rsun)).lt.0.0) then
            idt = minloc(abs(t_bin-t),DIM=1)
            idE = minloc(abs(E_bin-E),DIM=1)
            w3(idt,idE) = w3(idt,idE)+w
        end if
        if (((rs-14.d0*Rsun)*(rs0-14.d0*Rsun)).lt.0.0) then
            idt = minloc(abs(t_bin-t),DIM=1)
            idE = minloc(abs(E_bin-E),DIM=1)
            w4(idt,idE) = w4(idt,idE)+w
        end if



        time_current =  mpi_wtime()-time0
        if (time_current.gt.(10*3600))then
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
        end if





c outside the box (this is equivalent to "absorbing" boundary conditions)
         if(r.lt.Rmin .or. r.gt.Rmax .or. r_shock.gt.Rshmax)goto 30
         if(t.gt.tmax)goto 30

c diagnostics go here, between the dashed lines. 
c----------------BEGIN DIAGNOSTICS ---------------------
c this bit bins in momentum those "particles" that are just
c behind the shock
           if(r.gt.r_shock-Rsun .and. r.lt.r_shock)then
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
            w=weight*(2.d0**(-dfloat(lev-1)))
            dN(itd,ip)=dN(itd,ip)+w*fpc*ftdc
            dN(itd,ip1)=dN(itd,ip1)+w*fp*ftdc
            dN(itd1,ip)=dN(itd1,ip)+w*fpc*ftd
            dN(itd1,ip1)=dN(itd1,ip1)+w*fp*ftd
          end if
c-------------------END DIAGNOSTICS ----------------------------

c particle splitting
         if( SPLIT.eq.'Y' .or. SPLIT.eq.'y')then
          if(E.gt.E_split_lev(lev).and.lev.lt.n_split_levels
     &         .and.isp.le.isp_max)then
             lev=lev+1
             isp=isp+1
             t_save_split(isp)=t
             r_save_split(isp)=r
             p_save_split(isp)=p
             lev_save_split(isp)=lev
          end if
         end if 

         goto 20

 30    continue
c end of time loop

       if( lev.gt.1 )goto 15   ! need to go back and do all split particles

c every 100 "mother" particles, dump the diagnostics
c here the spectrum at the shock when it crossed 1AU is determined
c this is the differential intensity.  The units are 1/(cm^2 s sr MeV).
c the normalization is such that the density of energetic particles computed
c from this spectrum is n_ep (an input). 
        name13 = 'dist'
       if(mod(n,1).eq.0)then
        if(mype.le.9)then
         write(tmp,'(i1)')mype
         fname=trim(name13)//tmp(1:1)//'.dat'
         open(14,file=trim(path)//trim(fname),status='unknown')
        elseif( mype.ge.10 .and. mype.lt.100)then
         write(tmp,'(i2)')mype
         fname=trim(name13)//tmp(1:2)//'.dat'
         open(14,file=trim(path)//trim(fname),status='unknown')
        elseif( mype.ge.100 .and. mype.lt.1000)then
         write(tmp,'(i3)')mype
         fname=trim(name13)//tmp(1:3)//'.dat'
         open(14,file=trim(path)//trim(fname),status='unknown')
        else
         write(tmp,'(i4)')mype
         fname=trim(name13)//tmp(1:4)//'.dat'
         open(14,file=trim(path)//trim(fname),status='unknown')
        end if

        dNs=0.d0
        do 360 itd=1,ntd
        do 36 ip=1,np(itd)
        dNs=dNs+dN(itd,ip)
 36     continue
 360     continue
        if(dNs.gt.0.d0)then
         fac=1.d0/(fourpi*dNs)
        else
         fac=0.d0
        end if
        if(dNs.gt.0.d0)write(14,*)fac,ntd
        do 38 itd=1,ntd
        write(14,*)itd,np(itd)
        do 37 ip=1,np(itd)
c         p=p_diag_min*dexp( dpl_diag*dfloat(ip-1))
         p=p_diag_min*dexp( dpl_diag*(dfloat(ip)-0.5) )
         E=(p**2)/(2.d0*mp)
         dp_diag=dpop_diag*p
         f=dN(itd,ip)/((p**2)*dp_diag)*den_ep*fac     ! distribution function
         j=f*(p**2)                               ! diff. intensity
         write(14,*)E/MeV,dmax1(1.d-20,j*MeV)
c do this again, shifted by a bit to make it histogram style
c         p=p_diag_min*dexp( dpl_diag*dfloat(ip))
c         E=(p**2)/(2.d0*mp)
c         dp_diag=dpop_diag*p
c         f=dN(itd,ip)/((p**2)*dp_diag)*den_ep*fac
c         j=f*(p**2)
c         write(14,*)E/MeV,dmax1(1.d-20,j*MeV)
 37     continue
 38     continue
        close(14)
       end if
       


 35    continue

 100  continue

c finialize mpi routine
      close(901)
      close(902)
      close(903)
      close(904)
      call mpi_finalize( mpierr )

      stop
      end

      subroutine getU(r,t,U,dUdr)

c based on the appendix in Giacalone, 2015 (ApJ, vol. 799, 
c     article id. 80, Equations C1-C4)
  
      implicit none
      real*8 r,t,U,dUdr
      real*8 U2p,U1p,r_shock_p,coshR2
      real*8 r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt,Mach,
     &     V_sw_mod,Vsw
      real*8 Rsun
      integer*4 ic
      common /params/ r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt,
     &     Rsun,Mach,V_sw_mod
      data ic/0/
      save

      if(ic.eq.0)then
       ic=1
      end if

      call getShock( t )
      call getVsw( r, Vsw )

      U2p = ((s-1.d0)*v_shock + Vsw)/s
      U1p = Vsw

      r_shock_p = r_shock - 3.d0*drSHOCK
      if(r.gt.r_shock_p)then
        U = 0.5*(U1p+U2p)+0.5*(U1p-U2p)*dtanh((r-r_shock)/drSHOCK)
        coshR2 = (1.d0/dcosh((r-r_shock)/drSHOCK))**2
        dUdr = 0.5*(U1p-U2p)*coshR2/drSHOCK
      elseif(r.gt.r_shock_p-Ls .and. r.le.r_shock_p )then
        U = U2p*((r_shock_p/r)**2)
        dUdr = -2.d0*U2p*(r_shock_p**2)/(r**3)
      else
        U = U2p*((r_shock_p/(r_shock_p-Ls))**2)
        dUdr = 0.d0
      end if


      return
      end

      subroutine getK(r,t,p,K,dKdr)

c assumed K ir proportional to r**2
  
      implicit none
      real*8 r,t,p,K,dKdr
      real*8 r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt,Mach,
     &     V_sw_mod
      real*8 Rsun
      common /params/ r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt,
     &     Rsun,Mach,V_sw_mod
      save

      K = K0*((r/Rmin)**1.17)
      dKdr = 1.17*K/r

c momentum dependence
      K=K*((p/p0)**(2.d0*dlt))
      dKdr=dKdr*((p/p0)**(2.d0*dlt))


      return
      end

      function xrand()

      real*8 xrand,rn
      save

      call random_number( rn )
      xrand=rn

      return
      end

      subroutine getShock(t)

      implicit none

      real*8 t,r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt
      real*8 time_A(2000),r_shock_A(2000),Mach_A(2000),v_shock_A(2000),
     &    th,rshRsun,M,vshkms,tmax_data,time1,time2,FF,Mach,Rsun,
     &    Vswkms,x_sh_Rsun,y_sh_Rsun,z_sh_Rsun,shock_pos_rel_sun,
     &    V_sw_mod,V_sw_mod_A(2000),tmin_data,ss,s_A(2000),Rmax_data
      integer*4 ic,i,n
      character*30 name
      character*120 path
      common /params/ r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt,
     &     Rsun,Mach,V_sw_mod
      common /max_sim_time/ tmax_data,tmin_data,Rmax_data
      common /v_sw_array/ v_sw_mod_A,r_shock_A,n
      data ic/0/
      save
      name = 'PSP_para.dat'
      path = '/pscratch/sd/x/xhchen/PT1D/Laborday/PSP/'
      if(ic.eq.0)then
       ic=1
       open(10,file=trim(path)//trim(name),
     &     status='old')
       i=0
 10    read(10,*,end=11)th,rshRsun,M,vshkms,Vswkms,ss
        shock_pos_rel_sun = rshRsun*Rsun
        i=i+1
        time_A(i)=th*3600.d0               ! time in seconds
        r_shock_A(i)=shock_pos_rel_sun     ! shock radius (from center of Sun) in cm
        Mach_A(i) = M                     ! shock Mach number
        v_shock_A(i) = vshkms*1.d5         ! shock speed in cm/s
        v_sw_mod_A(i) = Vswkms*1.d+5       ! solar wind speed at the shock front
        if (ss.lt.1.d0) then
            ss = 1.d0
        end if
        s_A(i) = ss
        goto 10
 11    n=i
       tmax_data = time_A(n)
       tmin_data = time_A(1)
       Rmax_data = r_shock_A(n)
       do 12 i=2,n
        v_shock_A(i)=(r_shock_A(i)-r_shock_A(i-1))
     &     /(time_A(i)-time_A(i-1))
 12    continue
       v_shock_A(1)=v_shock_A(2)
      end if


c interpolate
      do 20 i=1,n
       if(time_A(i).gt.t)goto 21
 20   continue
 21   if(i.eq.1)then
       Mach=Mach_A(1)
       s = s_A(1)
       v_shock=v_shock_A(1)
       r_shock=r_shock_A(1)
       V_sw_mod=v_sw_mod_A(1)
       return
      end if
      time1=time_A(i-1)
      time2=time_A(i)

      FF = (t-time1)/(time2-time1)

      Mach = Mach_A(i-1) + FF*(Mach_A(i)-Mach_A(i-1))
      s = s_A(i-1) + FF*(s_A(i)-s_A(i-1))
      v_shock = v_shock_A(i-1) + FF*(v_shock_A(i)-v_shock_A(i-1))
      r_shock = r_shock_A(i-1) + FF*(r_shock_A(i)-r_shock_A(i-1))
      V_sw_mod = V_sw_mod_A(i-1) + FF*(V_sw_mod_A(i)-V_sw_mod_A(i-1))

      return
      end


      subroutine getVsw(r,VwP)

      implicit none
      real*8 r,VwP
      real*8 r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt,
     &     Rsun,Mach,V_sw_mod
      real*8 rs1,rs2,FF
      real*8 r_shock_A(2000),v_sw_mod_A(2000)
      integer*4 i,n
      common /params/ r_shock,v_shock,Vw,drSHOCK,Ls,s,K0,Rmin,p0,dlt,
     &     Rsun,Mach,V_sw_mod
      common /v_sw_array/ v_sw_mod_A,r_shock_A,n

      do 10 i=1,n
       if(r.lt.r_shock_A(i))goto 11
 10   continue
 11   if(i.eq.1)then
       VwP=v_sw_mod_A(1)
      elseif(i.eq.n)then
       VwP=v_sw_mod_A(n)
      else
       rs1=r_shock_A(i-1)
       rs2=r_shock_A(i)
       FF = (r-rs1)/(rs2-rs1) 
       VwP=V_sw_mod_A(i-1) + FF*(V_sw_mod_A(i)-V_sw_mod_A(i-1))     
      end if

      return
      end 

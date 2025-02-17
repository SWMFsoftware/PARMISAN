module PT_ModRandom
    use ModMpi
    use PT_ModConst, ONLY: cFourPi

    implicit none
    
    logical, public :: UseInputSeedFile
contains

    !============================================================================
    subroutine read_param(NameCommand)

        use ModReadParam, ONLY: read_var
        use ModUtilities, ONLY: CON_stop

        character(len=*), intent(in):: NameCommand ! From PARAM.in
        character(len=*), parameter:: NameSub = 'read_param'
        !--------------------------------------------------------------------------
        select case(NameCommand)
        case('#INPUTSEED')
            call read_var("UseInputSeedFile", UseInputSeedFile)

        case default
            call CON_stop(NameSub//' Unknown command '//NameCommand)
        end select
    end subroutine read_param
    !============================================================================
    subroutine init_random_seed()
        use PT_ModProc, ONLY: iProc, nProc
        use ModUtilities, ONLY: CON_stop

        integer, allocatable :: Seed_I(:)
        integer :: nSeed, fileUnit, istat
        character(len=*), parameter:: NameSub = 'init_random_seed'

        ! get size of seed (compiler-dependent)
        call random_seed(size = nSeed)
        allocate(Seed_I(nSeed))

        ! use dev/urandom to set seed for each processor
        open(newunit=fileUnit, file="/dev/urandom", access="stream", &
             form="unformatted", action="read", status="old", iostat=istat)
        if (istat == 0) then
            read(fileUnit) Seed_I
            close(fileUnit)
            call random_seed(put = Seed_I)
        else
            call CON_stop(NameSub//' /dev/urandom not found')
        end if

    end subroutine init_random_seed
    !============================================================================
    subroutine read_seed_file()
        ! reads in seed file to set the PRNG seed for each processor
        use PT_ModProc, ONLY: iProc, nProc

        integer, allocatable :: Seed_I(:)
        integer :: nSeed, io

        ! get size of seed array - compiler-dependent- 
        call random_seed(size = nSeed)
        ! first value in data file is processor id, then seed array of size n
        allocate(Seed_I(nSeed+1))

        ! open input file
        open(8, file = 'PT/Param/seed.in', action = 'read')

        do
            ! read line of seed file
            read(8, *, iostat = io) Seed_I(:)
            ! if processor IDs match, set the PRNG seed
            if (iProc == Seed_I(1)) then
                close(8)
                ! initialize seed
                call random_seed(put = Seed_I(2:nSeed+1))
                exit
            end if
            ! if no line matches, throw an error to screen 
            ! seed is set randomly
            ! need better way of making this error clear?
            if(io.ne.0) then
                close(8)
                write(*,*) 'Error in reading seed file for iProc: ', iProc, ' Setting random seed.'
                call init_random_seed
                exit           
            end if
        end do

    end subroutine read_seed_file
    !============================================================================
    subroutine save_seed()

        use PT_ModProc, ONLY: iProc, nProc, iComm, iError

        integer :: nSeed, i
        integer, allocatable :: Seed_I(:)
        integer, allocatable :: AllSeed_IP(:, :)

        ! get size of seed array - compiler-dependent
        call random_seed(size=nSeed)

        ! allocate seed array
        allocate(Seed_I(1:nSeed+1))

        ! first index is processor id, the rest is the seed
        Seed_I(1) = iProc
        call random_seed(get=Seed_I(2:nSeed+1))

        ! create output array
        allocate(AllSeed_IP(nSeed+1, nProc))
        
        ! combine all processor's seeds
        call MPI_Gather(Seed_I, nSeed+1, MPI_INTEGER, AllSeed_IP, &
                        nSeed+1, MPI_INTEGER, 0, iComm, iError)
        
        if(iProc.ne.0) return

        ! write to file
        open(8, file = 'PT/IO2/seed.out', status='new', action='write')
        do i = 1, nProc
            write(8, *) AllSeed_IP(:, i)
        end do
        close(8)

    end subroutine save_seed
    !============================================================================
    subroutine get_random_normal(RandNormal1)
        ! returns random number sampled from normal distribution
        ! with mean = 0 and std = 1

        real, intent(out) :: RandNormal1 !, RandNormal2
        real :: RandUniform1, RandUniform2

        ! uniform random numbers over [0,1)
        call random_number(RandUniform1)
        call random_number(RandUniform2)

        ! redistribute to (0, 1] to avoid 0
        RandUniform1 = 1 - RandUniform1
        RandUniform2 = 1 - RandUniform2
        
        ! Box-Muller transformation

        ! two independent random variable with standard normal distribution
        ! only need one for 1-D version

        ! this limits the random variable to rn < ~6 
        ! physically limiting the size of the random diffusive process
        
        RandNormal1 = sqrt(-2*log(RandUniform1))*cos(0.5*cFourPi*RandUniform2)
        ! RandNormal2 = sqrt(-2*log(RandUniform1))*sin(0.5*fourpi*RandUniform2)

    end subroutine get_random_normal
!================================================================================
end module PT_ModRandom
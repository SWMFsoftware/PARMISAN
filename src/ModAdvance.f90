!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
!==============================================================================
module PT_ModAdvance

    implicit none
    SAVE

contains
    !--------------------------------------------------------------------------
    subroutine advance(TimeLimit)  
        
        use PT_ModGrid,      only: nLine, Used_B, iShock_IB, Shock_, &
                                   ShockOld_
        use PT_ModShock,     only: DoOutputShock
        use PT_ModParticle,  only: advance_particles, &
                                   inject_particles
        use PT_ModFieldline, only: set_fieldline, advect_fieldline, save_fieldline_data
        use PT_ModPlot,      only: save_distribution_function
        use PT_ModTime,      only: iIter, PTTime, DataInputTime
        use ModMpi
        use PT_ModProc,      only: iComm, iError, iProc

        real, intent(in) :: TimeLimit
        integer :: iLine, iProgress, iShock, iShockOld, nProgress, iShockNew
        real :: TotalDt, DtProgress, NextTimeStep, Alpha, LagrInject

        ! Total timestep between MHD states
        TotalDt = DataInputTime - PTTime

        ! Loop over magnetic field lines
        LINE: do iLine = 1, nLine

            if(.not.Used_B(iLine)) CYCLE LINE

            ! set fieldline number in ModFieldline
            call set_fieldline(iLine)

            ! indices of shock during previous and current MHD state
            iShock = iShock_IB(Shock_, iLine)
            iShockOld = iShock_IB(ShockOld_, iLine)

            ! how many vertices (lagrangian coordinates) the shock moved
            nProgress = max(1, iShock - iShockOld)
            iShockOld = min(iShockOld, iShock-1)
            if(iProc.eq.0) write(*,*) 'ShockOld, ShockNew, dShock: ', &
                                      iShockOld, iShock, iShock - iShockOld
            
            ! timestep for each subinterval where the shock moves one 
            ! lagrangian coordinate
            DtProgress = TotalDt / nProgress
        
            ! loop over subintervals where the shock moves one 
            ! lagrangian coordinate
            PROGRESS: do iProgress = 1, nProgress
                
                ! new shock vertex
                iShockNew = iShockOld + iProgress
                ! interpolation weight
                Alpha = iProgress / real(nProgress)
                ! final time for this subinterval
                NextTimeStep = PTTime + (iProgress*DtProgress)
                ! inject particles 5 lagrangian coordinates upstream of shock
                LagrInject = iShockNew + 5
                ! advect shock
                call advect_fieldline(Alpha, iShockNew, NextTimeStep)

                if(DoOutputShock) call save_fieldline_data(iProgress, NextTimeStep)

                ! inject particles at shock location at start of subinterval
                call inject_particles(iLine, NextTimeStep - DtProgress, LagrInject)

                ! advance particles in time until NextTimeStep
                call advance_particles(iLine, min(TimeLimit, NextTimeStep), TimeLimit)

                if(NextTimeStep.ge.TimeLimit) EXIT PROGRESS

            end do PROGRESS
            
            call MPI_BARRIER(iComm, iError)
            call save_distribution_function(iIter + 1, TimeLimit)

        end do LINE

    end subroutine advance
    !--------------------------------------------------------------------------
end module PT_ModAdvance
!==============================================================================
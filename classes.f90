
module classes

    implicit none
    
    type solid_parameters
    
    end type
    
    type fluid_parameters
    
    end type
    
    type, public :: Phase
        real(8) :: beta
        character :: phase_type(6)
        real(8), allocatable :: K, composition
    contains
        procedure :: update_variables_x => null()
        procedure :: update_variables_T => null()
        procedure :: calc_fugcoef => null()
    end type
    
    type, public, extends(Phase) :: FluidPhase
        
    end type
    
    type, public, extends(FluidPhase) :: SolidPhase
        real(8), allocatable :: tau
    end type
    
    contains
    
    subroutine load_phase(phase_object, phase_type_out)
        character :: phase_type(:)
        
        if(phase_type_out .eq. "vapor" .or. phase_type_out .eq. "liquid") then
            phase_object = FluidPhase()
            phase_object%update_variables_x => FFFFFFFFFFFFF
            phase_object%update_variables_T => FFFFFFFFFFFFF
            phase_object%calc_fugcoef => FFFFFFFFFFFFF
        else if(phase_type_out .eq. "solid") then
            phase_object = SolidPhase()
            phase_object%update_variables_x => FFFFFFFFFFFFF
            phase_object%update_variables_T => FFFFFFFFFFFFF
            phase_object%calc_fugcoef => FFFFFFFFFFFFF
        else
            print*,
            print*, "Error while loading phase: No match for phase_type."
            print*,
            stop
        end if
        
        phase_object%phase_type = phase_type_out
        allocate(phase_object%K(Ncomp),phase_object%composition(Ncomp))
    end subroutine
    
end module


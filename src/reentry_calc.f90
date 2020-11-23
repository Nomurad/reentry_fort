module mod_reentry_params
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    real(real64), parameter :: rho0         = 1.225d0  ! kg/m^3
    real(real64), parameter :: Scale_height = 7.163d0  ! kk
    real(real64), parameter :: r_earth      = 6378d0
    real(real64), parameter :: g            = 9.8d-3
    real(real64), parameter :: mu           = 3.986004418d5
    real(real64), parameter :: omg_e        = 7.292d5

    type :: reentry_parameters
        real(real64) :: rho0         = rho0
        real(real64) :: Scale_height = Scale_height
        real(real64) :: r_earth      = r_earth
        real(real64) :: g            = g
        real(real64) :: mu           = mu
        real(real64) :: omg_e        = omg_e
    end type reentry_parameters
    type(reentry_parameters), protected, public :: reentry_params

end module mod_reentry_params


module mod_reentry_calc
    use, intrinsic :: iso_fortran_env
    use mod_calcutil
    use mod_vector
    use mod_reentry_params
    implicit none

    type t_SpaceCraft
        real(real64) weight
    end type t_SpaceCraft

    type t_ReentryCalc
        real(real64) :: height = 0d0
        real(real64) :: r = 0d0
        real(real64) :: V = 0d0
        real(real64) :: dV = 0d0

        real(real64) :: theta = 0d0
        real(real64) :: phi = 0d0
        real(real64) :: gamma = 0d0
        real(real64) :: d_gamma = 0d0
        real(real64) :: psi = 0d0
        real(real64) :: d_kai = 0d0
        real(real64) :: alpha = 0d0
        real(real64) :: d_alpha = 0d0
        real(real64) :: delta = 0d0
        real(real64) :: d_delta = 0d0

        type(t_Vector) r_vec
        type(t_Vector) V_vec
        type(t_Vector) dV_vec

        contains
            procedure, public :: atmos_density
            procedure :: set_r_vec_arr
            procedure :: set_r_vec_vector
            generic :: set_r_vec => set_r_vec_arr, set_r_vec_vector 
    end type

    interface t_ReentryCalc
        module procedure init
    end interface

    contains
    type(t_ReentryCalc) function init(craft) result(res)
        class(t_SpaceCraft), intent(in) :: craft

       
    end function init

    real(real64) function atmos_density(this, h) result(res)
        class(t_ReentryCalc), intent(in) :: this
        real(real64), intent(in) :: h 
        res = reentry_params%rho0 * exp(-h/reentry_params%Scale_height)

    end function atmos_density

    real(real64) function gravity(this, r) result(res)
        class(t_ReentryCalc), intent(in) :: this
        real(real64), intent(in) :: r
        res = reentry_params%mu / (r**2)

    end function gravity

    subroutine set_r_vec_arr(this, r_vec)
        class(t_ReentryCalc), intent(inout) :: this
        real(real64), intent(in) :: r_vec(3)
        real(real64) x, y, z 
        
        x = r_vec(1)
        y = r_vec(2)
        z = r_vec(3)
        this%r = norm(r_vec) 
        this%theta = atan2(y, x)
        if(this%theta < 0d0) then
            this%theta = 2*pi + this%theta
        end if

        this%phi = acos(norm([x,y,0])/this%r) 

    end subroutine set_r_vec_arr

    subroutine set_r_vec_vector(this, r_vec)
        class(t_ReentryCalc), intent(inout) :: this
        class(t_Vector), intent(in) :: r_vec
        real(real64) x, y, z 
        
        x = r_vec%vec(1)
        y = r_vec%vec(2)
        z = r_vec%vec(3)
        this%r = norm(r_vec%vec) 
        this%theta = atan2(y, x)
        if(this%theta < 0d0) then
            this%theta = 2*pi + this%theta
        end if

        this%phi = acos(norm([x,y])/this%r) 

    end subroutine set_r_vec_vector

    subroutine set_v_vec_arr(

end module mod_reentry_calc



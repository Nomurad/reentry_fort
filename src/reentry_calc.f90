module mod_reentry_calc
    use, intrinsic :: iso_fortran_env
    use mod_calcutil
    use mod_vector
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
    type(reentry_parameters), public :: reentry_params
    

    contains

end module mod_reentry_calc


#ifdef _debug
program debug
    use, intrinsic :: iso_fortran_env 
    use mod_reentry_calc 
    implicit none 
    real(real64) :: v(3)
    real(real64) :: a(3,3), b(4,3)
    integer :: i, j
    
    reentry_params%rho0 = 0d0
    print *, reentry_params
end program
#endif 

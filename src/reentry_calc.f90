module mod_reentry_params
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    real(real64), parameter :: rho0         = 1.225d0  ! kg/m^3
    real(real64), parameter :: Scale_height = 7.163d0  ! kk
    real(real64), parameter :: r_earth      = 6378d0
    real(real64), parameter :: g            = 9.8d-3
    real(real64), parameter :: mu           = 3.986004418d5
    real(real64), parameter :: omg_e        = 7.292d-5

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


module mod_spacecraft
    use, intrinsic :: iso_fortran_env 
    use mod_vector
    implicit none 
    private

    type, public :: t_SpaceCraft
        real(real64) weight
        real(real64) ref_area
        type(t_Vector) :: posi
        real(real64) :: C_L = 0d0
        real(real64) :: C_D = 1d0
        real(real64) :: ballistic
        contains
            procedure, public :: set_ballistic 
    end type t_SpaceCraft

    interface t_SpaceCraft 
        module procedure init_craft
    end interface t_SpaceCraft 

    contains

    function init_craft(posi, weight, ref_area) result(craft)
        type(t_SpaceCraft) :: craft
        type(t_Vector), intent(in) :: posi
        real(real64), intent(in) :: weight
        real(real64), intent(in) :: ref_area
        
        craft%posi = posi 
        craft%weight = weight 
        craft%ref_area = ref_area 
        call craft%set_ballistic()

    end function init_craft

    subroutine set_ballistic(this, b, res)
        class(t_SpaceCraft), intent(inout) :: this 
        real(real64), optional, intent(in) :: b 
        real(real64), optional, intent(out) :: res

        if(present(b)) then 
            this%ballistic = b 
        else 
            this%ballistic = this%weight / (1d0*this%ref_area)
        end if

        if(present(res)) then
            res = this%ballistic
        end if

    end subroutine set_ballistic

end module mod_spacecraft


module mod_reentry_calc
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: iso_c_binding
    use mod_calcutil
    use mod_vector
    use mod_reentry_params
    use mod_spacecraft
    implicit none
    private

    public :: reentry_params

    type, public :: trj_results
        real(real64), allocatable :: r_hist(:)
        real(real64), allocatable :: theta_hist(:)
        real(real64), allocatable :: phi_hist(:)
        real(real64), allocatable :: rho_hist(:)
        real(real64), allocatable :: V_hist(:)
        real(real64), allocatable :: gamma_hist(:)
        real(real64), allocatable :: psi_hist(:)
        real(real64), allocatable :: theta_dot_hist(:)
        real(real64), allocatable :: phi_dot_hist(:)
    end type

    type, public :: t_ReentryCalc
        type(t_SpaceCraft) :: craft
        type(t_Vector)    r_vec
        type(t_Vector)    V_vec
        type(t_Vector)    dV_vec
        type(t_Vector)    posi
        type(t_Vector)     polor
        type(trj_results) results

        real(real64) :: height = 0d0
        ! real(real64), pointer :: r 
        ! real(real64), pointer :: V
        ! real(real64), pointer :: dV
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

        real(real64) Lift, Drag
        real(real64) La, Da

        contains
            procedure, public :: atmos_density, gravity
            procedure :: set_r_vec_arr, set_r_vec_vector
            procedure :: set_v_vec_arr, set_v_vec_vector
            procedure :: add_v_scalar, add_v_vec
            procedure :: calc_Lift, calc_Drag
            procedure :: diff_r, diff_phi, diff_theta
            procedure :: diff2_r, diff2_theta, diff2_phi
            procedure :: trj_calc_3d
            generic, public :: add_v => add_v_scalar, add_v_vec
            generic, public :: set_r_vec => set_r_vec_arr, set_r_vec_vector 
            generic, public :: set_v_vec => set_v_vec_arr, set_v_vec_vector 
    end type t_ReentryCalc

    
    interface t_ReentryCalc
        module procedure init
    end interface

    interface trj_results 
        module procedure hist_alloc 
    end interface

    contains
        function hist_alloc(len)
            type(trj_results) :: hist_alloc
            integer, intent(in) :: len

            allocate(hist_alloc%r_hist(len), source=-1d0)
            allocate(hist_alloc%theta_hist(len), source=1000d0)
            allocate(hist_alloc%phi_hist(len))
            allocate(hist_alloc%rho_hist(len))
            allocate(hist_alloc%V_hist(len))
            allocate(hist_alloc%gamma_hist(len))
            allocate(hist_alloc%psi_hist(len))
            allocate(hist_alloc%theta_dot_hist(len))
            allocate(hist_alloc%phi_dot_hist(len))

        end function hist_alloc

        function init(craft, h) result(this)
            type(t_ReentryCalc) :: this 
            class(t_SpaceCraft), intent(in) :: craft
            real(real64), optional, intent(in) :: h 
            real(real64) rho

            this%craft = craft
            if(present(h)) then 
                this%height = h
            else 
                this%height = 120d0
            end if 

            this%r_vec = t_Vector(3) 
            this%V_vec = t_Vector(3)
            this%dV_vec = t_Vector(3)
            
            this%r_vec = craft%posi 
            ! this%r = reentry_params%r_earth + this%height 
            this%r = norm(this%r_vec) 
            this%V = norm(this%V_vec)
            this%dV = norm(this%dV_vec)
            this%theta = 0d0
            this%phi = 0d0

            rho = this%atmos_density(this%height)
            this%Lift = this%calc_Lift(rho)
            this%Drag = this%calc_Drag(rho)
            this%La = this%Lift/this%craft%weight 
            this%Da = this%Drag/this%craft%weight 

            this%posi = this%craft%posi  ! vecのpointerは同じ場所を指す
        end function init

        real(real64) function atmos_density(this, h) result(rho)
            class(t_ReentryCalc), intent(in) :: this
            real(real64), intent(in) :: h 

            if(h > 120d0) return
            
            rho = reentry_params%rho0 * exp(-h/reentry_params%Scale_height)

        end function atmos_density

        real(real64) function gravity(this, r) result(res)
            class(t_ReentryCalc), intent(in) :: this
            real(real64), intent(in) :: r

            ! if(r-reentry_params%r_earth > 120d0) return
            res = reentry_params%mu / (r**2)
            ! print *, "g = ", res

        end function gravity

        subroutine set_r_vec_arr(this, r_vec)
            class(t_ReentryCalc), intent(inout) :: this
            real(real64), target, intent(in) :: r_vec(3)
            real(real64), pointer :: x, y, z 
            
            x => r_vec(1)
            y => r_vec(2)
            z => r_vec(3)
            this%r = norm(r_vec) 
            this%theta = atan2(y, x)
            if(this%theta < 0d0) then
                this%theta = 2*pi + this%theta
            end if

            this%phi = acos(norm([x,y])/this%r) 

        end subroutine set_r_vec_arr

        subroutine set_r_vec_vector(this, r_vec)
            class(t_ReentryCalc), intent(inout) :: this
            class(t_Vector), intent(in) :: r_vec
            ! real(real64) x, y, z 
            
            call this%set_r_vec_arr(r_vec%vec)

        end subroutine set_r_vec_vector

        subroutine set_v_vec_arr(this, v_vec)
            class(t_ReentryCalc), intent(inout) :: this
            real(real64), target, intent(in) :: v_vec(3)
            real(real64), pointer :: x, y, z 
            real(real64) :: r, v
            type(t_Vector) :: rv, vv
            integer i

            x => v_vec(1)
            y => v_vec(2)
            z => v_vec(3)
            forall(i=1:3) this%V_vec%vec(i) = v_vec(i)
            this%V = norm(v_vec)
            rv = this%r_vec
            vv = this%V_vec
            r = this%r 
            v = this%v

            ! gamma = 経路角 
            this%gamma = pi/2d0 - acos(rv%dot_product(vv)/(r*v))

        end subroutine set_v_vec_arr 

        subroutine set_v_vec_vector(this, v_vec)
            class(t_ReentryCalc), intent(inout) :: this
            class(t_Vector), intent(in) :: v_vec 

            call this%set_v_vec_arr(v_vec%vec)

        end subroutine set_v_vec_vector 

        subroutine add_v_scalar(this, v_scalar)
            class(t_ReentryCalc), intent(inout) :: this 
            real(real64), target, intent(in) :: v_scalar
            ! real(real64), pointer :: x, y, z
            ! integer i
            real(real64) dummy

            this%V = this%V + v_scalar
            dummy = diff_r(this)
            dummy = diff_theta(this)
            dummy = diff_phi(this)

        end subroutine

        subroutine add_v_vec(this, v_vec)
            class(t_ReentryCalc), intent(inout) :: this 
            real(real64), target, intent(in) :: v_vec(3) 
            real(real64), pointer :: x, y, z 
            integer i

        end subroutine

        real(real64) function calc_Lift(this, rho) result(Lift)
            class(t_ReentryCalc), intent(inout) :: this
            real(real64), intent(in) :: rho 
            real(real64) V, S, C_L 

            V = this%V 
            S = this%craft%ref_area 
            C_L = this%craft%C_L
            Lift = rho*(V**2)*S*C_L/2.0d0
            this%Lift = Lift 

        end function calc_Lift

        real(real64) function calc_Drag(this, rho) result(Drag)
            class(t_ReentryCalc), intent(inout) :: this
            real(real64), intent(in) :: rho 
            real(real64) V, S, C_D 

            V = this%V 
            S = this%craft%ref_area 
            C_D = this%craft%C_D
            Drag = rho*(V**2)*S*C_D/2.0d0
            this%Drag = Drag 

        end function calc_Drag

        real(real64) function diff_r(this)
            class(t_ReentryCalc), intent(inout) :: this
            diff_r = this%V * sin(this%gamma)
        end function diff_r

        real(real64) function diff_theta(this) 
            class(t_ReentryCalc), intent(inout) :: this
            diff_theta = this%V * cos(this%gamma) * cos(this%psi) / this%r
        end function diff_theta

        real(real64) function diff_phi(this) 
            class(t_ReentryCalc), intent(inout) :: this
            diff_phi = (this%V/this%r) * cos(this%gamma) * sin(this%psi)
        end function diff_phi

        real(real64) function diff2_r(this, dots, F_thrust) 
            class(t_ReentryCalc), intent(inout) :: this
            real(real64), intent(in) :: dots(3)
            real(real64), optional, value :: F_thrust
            real(real64) r_dot, theta_dot, phi_dot
            real(real64) r, phi, L, D, m, gamma_, delta, omg_e, g 
            real(real64) F_r, F_t
            real(real64) cosphi2
            logical is_exist_thrust

            is_exist_thrust = (present(F_thrust))
            if(.not. is_exist_thrust) F_thrust = 0d0

            r_dot     = dots(1)
            theta_dot = dots(2)
            phi_dot   = dots(3)

            r      = this%r
            phi    = this%phi
            L      = this%Lift
            D      = this%Drag
            m      = this%craft%weight
            gamma_ = this%gamma
            delta  = this%delta
            omg_e  = reentry_params%omg_e
            g      = this%gravity(r)
            
            F_r     = (L/m)*cos(delta)*cos(gamma_) - (D/m)*sin(gamma_) + F_thrust/m
            cosphi2 = cos(phi)**2
            diff2_r = r*(phi_dot**2 + cosphi2*(theta_dot + omg_e)**2) + F_r - g

        end function diff2_r

        real(real64) function diff2_theta(this, dots, F_thrust) 
            class(t_ReentryCalc), intent(inout) :: this
            real(real64), intent(in) :: dots(3)
            real(real64) r_dot, theta_dot, phi_dot
            real(real64), optional, value :: F_thrust
            real(real64) r, phi, L, D, m, gamma_, psi, delta, omg_e, g 
            real(real64) F_theta
            real(real64) phi_dot2
            logical is_exist_thrust

            is_exist_thrust = (present(F_thrust))
            if(is_exist_thrust .eqv. .false.) F_thrust = 0d0

            r_dot     = dots(1)
            theta_dot = dots(2)
            phi_dot   = dots(3)

            r      = this%r
            phi    = this%phi
            L      = this%Lift
            D      = this%Drag
            m      = this%craft%weight
            gamma_ = this%gamma
            psi    = this%psi
            delta  = this%delta
            omg_e  = reentry_params%omg_e
            !g     = this%gravity(r)
            phi_dot2 = phi_dot**2
            
            F_theta     = (L/m)*sin(delta)*sin(psi) - ((L/m)*cos(delta)*sin(gamma_) + (D/m)*cos(gamma_))*cos(psi)
            diff2_theta = -2*r_dot*theta_dot*cos(phi) + 2*r*theta_dot*phi_dot*sin(phi) &
                        & - 2*omg_e*(r_dot*cos(phi)-r_dot*phi_dot*sin(phi)) + F_theta + F_thrust/m
            diff2_theta = diff2_theta/(r*cos(phi))

        end function diff2_theta

        real(real64) function diff2_phi(this, dots, F_thrust) 
            class(t_ReentryCalc), intent(inout) :: this
            real(real64), intent(in) :: dots(3)
            real(real64) r_dot, theta_dot, phi_dot
            real(real64), optional, value :: F_thrust
            real(real64) r, phi, L, D, m, gamma_, psi, delta, omg_e, g 
            real(real64) F_phi
            logical is_exist_thrust

            is_exist_thrust = (present(F_thrust))
            if(is_exist_thrust .eqv. .false.) F_thrust = 0d0

            r_dot     = dots(1)
            theta_dot = dots(2)
            phi_dot   = dots(3)

            r      = this%r
            phi    = this%phi
            L      = this%Lift
            D      = this%Drag
            m      = this%craft%weight
            gamma_ = this%gamma
            psi    = this%psi
            delta  = this%delta
            omg_e  = reentry_params%omg_e
            !g     = this%gravity(r)
            
            F_phi     = -(L/m)*sin(delta)*cos(psi) - ((L/m)*cos(delta)*sin(gamma_) + (D/m)*cos(gamma_))*sin(psi)
            diff2_phi = -2*r_dot*phi_dot - r*sin(phi)*cos(phi)*((theta_dot+omg_e)**2) + F_phi + F_thrust/m
            diff2_phi = diff2_phi/r

        end function diff2_phi

        subroutine trj_calc_3d(this, x0, t, F_t, res)
            class(t_ReentryCalc), target, intent(inout) :: this
            real(real64), intent(in) :: x0(3), t(:)
            ! real(real64), intent(in) :: dv_s(:)
            real(real64), intent(in) :: F_t(:, :)
            type(trj_results), intent(out) :: res
            real(real64) r0, theta0, phi0
            real(real64) r, theta, phi
            real(real64) r_dot, theta_dot, phi_dot
            real(real64) dots(3)
            real(real64) rho, r_earth
            integer iter_max, i
            real(real64) dt
            real(real64) dummy
            real(real64) F(3)
            real(real64), allocatable :: F_t2(:, :)
            ! type(t_Vector) k(4), l(4)
            real(real64) r_k(4), theta_k(4), phi_k(4)
            real(real64) r_l(4), theta_l(4), phi_l(4)

            res        = trj_results(size(t))
            r0         = x0(1)
            theta0     = x0(2)
            phi0       = x0(3)

            r          = r0
            theta      = theta
            phi        = phi0
            
            this%r     = r
            this%theta = theta
            this%phi   = phi

            ! forall(i = 1:4) k(i) = t_Vector(3)
            
            r_dot      = this%diff_r()
            theta_dot  = this%diff_theta()
            phi_dot    = this%diff_phi()

            r_earth = reentry_params%r_earth
            rho = this%atmos_density(r0 - r_earth)

            iter_max = size(t)
            ! allocate(F_t2(3, iter_max), source=0d0)
            ! F_t2 = F_t
            ! F = [-1.0d-3, 0.0d0, 0.0d0]
            ! F_t2(1,:) = F(1)
            ! F_t2(2,:) = F(2)/this%r
            ! F_t2(3,:) = F(3)/this%r
            
            this%results = trj_results(iter_max)
            dt = t(2) - t(1)

            print *, i, r
            print *, r_dot, r*theta_dot
            do i = 1, iter_max
                ! call this%add_v(dv_s(i))
                rho       = this%atmos_density(r - r_earth)
                dummy     = this%calc_Drag(rho)
                dummy     = this%calc_Lift(rho)
                r_dot     = this%diff_r()
                theta_dot = this%diff_theta()
                phi_dot   = this%diff_phi()
                
                dots       = [r_dot, theta_dot, phi_dot]
                r_k(1)     = r_dot 
                theta_k(1) = theta_dot 
                phi_k(1)   = phi_dot 
                r_l(1)     = this%diff2_r(dots, F_t(1,i))
                theta_l(1) = this%diff2_theta(dots, F_t(2,i))
                phi_l(1)   = this%diff2_phi(dots, F_t(3,i))

                dots       = [r_dot+(r_k(1)/2d0), theta_dot+(theta_k(1)/2d0), phi_dot+(phi_k(1)/2d0)]
                r_k(2)     = r_dot     + r_l(1)/2d0
                theta_k(2) = theta_dot + theta_l(1)/2d0
                phi_k(2)   = phi_dot   + phi_l(1)/2d0
                r_l(2)     = this%diff2_r(dots, F_t(1,i))
                theta_l(2) = this%diff2_theta(dots, F_t(2,i))
                phi_l(2)   = this%diff2_phi(dots, F_t(3,i))

                dots       = [r_dot+(r_k(2)/2d0), theta_dot+(theta_k(2)/2d0), phi_dot+(phi_k(2)/2d0)]
                r_k(3)     = r_dot     + r_l(2)/2d0
                theta_k(3) = theta_dot + theta_l(2)/2d0
                phi_k(3)   = phi_dot   + phi_l(2)/2d0
                r_l(3)     = this%diff2_r(dots, F_t(1,i))
                theta_l(3) = this%diff2_theta(dots, F_t(2,i))
                phi_l(3)   = this%diff2_phi(dots, F_t(3,i))

                dots       = [r_dot+r_k(3), theta_dot+theta_k(3), phi_dot+phi_k(3)]
                r_k(4)     = r_dot     + r_l(3)
                theta_k(4) = theta_dot + theta_l(3)
                phi_k(4)   = phi_dot   + phi_l(3)
                r_l(4)     = this%diff2_r(dots, F_t(1,i))
                theta_l(4) = this%diff2_theta(dots, F_t(2,i))
                phi_l(4)   = this%diff2_phi(dots, F_t(3,i))

                r         = r     + dt*(r_k(1) + 2d0*r_k(2) + 2d0*r_k(3) + r_k(4))/6d0
                theta     = theta + dt*(theta_k(1) + 2d0*theta_k(2) + 2d0*theta_k(3) + theta_k(4))/6d0
                phi       = phi   + dt*(phi_k(1) + 2d0*phi_k(2) + 2d0*phi_k(3) + phi_k(4))/6d0
                r_dot     = r_dot     + dt*(r_l(1) + 2d0*r_l(2) + 2d0*r_l(3) + r_l(4))/6d0
                theta_dot = theta_dot + dt*(theta_l(1) + 2d0*theta_l(2) + 2d0*theta_l(3) + theta_l(4))/6d0
                phi_dot   = phi_dot   + dt*(phi_l(1) + 2d0*phi_l(2) + 2d0*phi_l(3) + phi_l(4))/6d0

                this%r     = r
                this%theta = theta
                this%phi   = phi
                call this%set_v_vec([r_dot, r*theta_dot, r*phi_dot])
                ! this%V = norm([r_dot, r*theta_dot, r*phi_dot])
                ! this%V = norm(this%V_vec)
                this%psi   = atan((r*phi_dot)/(r*theta_dot))
                this%gamma = atan(r_dot / norm([(r*phi_dot), (r*theta_dot)]))

                if ((mod(i, 10) == 0d0) .or. (i == 1)) then
                    write(*,"(a, 10(g15.8:))", advance="no") achar(z"0d"), dt*i, r-reentry_params%r_earth, r_dot
                    ! write(*,"(a, 10(g15.8:))", advance="no") achar(z"0a"), dt*i, r-reentry_params%r_earth, r_dot
                    ! write(*,"(i7,3f15.8)") i, r-reentry_params%r_earth, this%V, F_t(:,i)
                end if

                if(r <= reentry_params%r_earth) then 
                    print *, ""
                    print "('r=', f15.3, ' / iter: ', i10)", r, i 
                    print *, "landing."
                    this%results%r_hist         =  this%results%r_hist(:i-1)
                    this%results%theta_hist     =  this%results%theta_hist(:i-1)
                    this%results%phi_hist       =  this%results%phi_hist(:i-1)
                    this%results%V_hist         =  this%results%V_hist(:i-1)
                    this%results%gamma_hist     =  this%results%gamma_hist(:i-1)
                    this%results%psi_hist       =  this%results%psi_hist(:i-1)
                    this%results%rho_hist       =  this%results%rho_hist(:i-1)
                    this%results%theta_dot_hist =  this%results%theta_dot_hist(:i-1)
                    this%results%phi_dot_hist   =  this%results%phi_dot_hist(:i-1)
                    exit
                end if

                this%results%r_hist(i)         = r 
                this%results%theta_hist(i)     = theta 
                this%results%phi_hist(i)       = phi 
                this%results%V_hist(i)         = this%V 
                this%results%gamma_hist(i)     = this%gamma 
                this%results%psi_hist(i)       = this%psi 
                this%results%rho_hist(i)       = rho 
                this%results%theta_dot_hist(i) = theta_dot
                this%results%phi_dot_hist(i)   = phi_dot

            end do
            if(i >= iter_max) print *, "landing failed..."
            res = this%results

        end subroutine trj_calc_3d
                           
end module mod_reentry_calc



subroutine reentry_from_py(posi, weight, ref_area, initV, gamma, psi, tsize, ts, F_t)
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: iso_c_binding
    use mod_vector
    use mod_reentry_params
    use mod_spacecraft
    use mod_reentry_calc
    implicit none 
    real(real64), intent(in) :: posi(3), weight, ref_area
    real(real64), intent(in) :: initV
    real(real64), intent(in) :: gamma, psi
    integer(c_int), intent(in) :: tsize
    real(real64), intent(in) :: ts(tsize)
    ! real(real64), intent(in) :: F_t(3,tsize)
    real(real64), intent(in) :: F_t(tsize)
    ! real(real64), intent(inout) :: ress
    ! type(trj_results), intent(inout) :: res
    type(trj_results) res
    real(real64), allocatable :: outputs(:,:)

    type(t_Vector) craft_vec
    type(t_SpaceCraft) craft 
    type(t_ReentryCalc) reentry 
    real(real64) dt
    integer i
    real(real64) gravity

    print *, "start calc!"
    ! ress = 1d0
    craft_vec = t_Vector(3)
    craft_vec%vec = posi
    craft = t_SpaceCraft(craft_vec, weight, ref_area)
    reentry = t_ReentryCalc(craft)
    
    gravity = reentry%gravity(reentry%r)

    call reentry%set_r_vec(posi)
    ! call reentry%set_v_vec(posi)
    reentry%V = initV
    reentry%gamma = gamma 
    reentry%psi = psi
    dt = ts(2) - ts(1)
    ! debug
#ifdef _debug
    print *, "----- fortran ------"
    print *, "dt =", dt
    print *, "g = ", gravity, reentry%gravity(6378d0)
    ! print *, "posi = ", craft%posi
    print *, "posi = ", craft_vec
    print *, "R =", reentry%r
    print *, "V = ", reentry%v
#endif

    ! call reentry%trj_calc_3d([reentry_params%r_earth+400d0, 0d0, 0d0], ts, F_t, res)
    call reentry%trj_calc_3d([reentry%r, reentry%theta, reentry%phi], ts, F_t, res)
    ! print *, "time:", size(res%r_hist)*dt, "[s]"

    allocate(outputs(11, ubound(res%r_hist, dim=1)))
    outputs(1,:) = ts(:ubound(res%r_hist, dim=1)) 
    outputs(2,:) = res%r_hist 
    outputs(3,:) = res%theta_hist
    outputs(4,:) = res%phi_hist
    outputs(5,:) = res%rho_hist
    outputs(6,:) = res%V_hist
    outputs(7,:) = res%gamma_hist 
    outputs(8,:) = res%psi_hist 
    outputs(9,:) = res%r_dot_hist 
    outputs(10,:) = res%theta_dot_hist 
    outputs(11,:) = res%phi_dot_hist 
    open(20, file="trj_calc_rslt.dat", form="unformatted")
    write(20) outputs
    ! do i = 1, ubound(res%r_hist, dim=1)
    !     write(20) ts(i), res%r_hist(i)
    ! end do
    close(20)
    
end subroutine
program debug 
    use mod_vector 
    use mod_calcutil 
    use mod_spacecraft
    use mod_reentry_calc 
    implicit none 
    type(t_Vector) :: vector
    type(t_SpaceCraft) :: mycraft
    type(t_ReentryCalc) :: reentry 
    type(t_Vector) :: tmp, v(4)
    integer i
    integer, parameter :: tsize = 1000000
    real(real64) ts(tsize)
    real(real64) :: t = 1d-3
    type(trj_results) res

    print *, "reentry parameters:"
    print *, reentry_params

    vector = t_Vector(3)

    vector%vec = [reentry_params%r_earth + 120d0, 0d0, 0d0]
    mycraft = t_SpaceCraft(vector, 100d0, pi)
    reentry = t_ReentryCalc(mycraft)
    print *, "position vec"
    print *, reentry%posi%vec
    
    reentry%V = 7.74d0
    reentry%gamma = deg2rad(-1.6d0)
    reentry%psi = deg2rad(0d0)
    print *, "V = ", reentry%V
    forall(i=1:tsize) ts(i) = t*i
    call reentry%trj_calc_3d([reentry_params%r_earth + 120d0, 0d0, 0d0], ts, res)
    print *, "data size: ", size(res%r_hist), res%r_hist(ubound(res%r_hist))

end program

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
    integer nunit, ios
    character(64) :: img
    integer, parameter :: tsize = 10000000
    real(real64) ts(tsize)
    real(real64) dv_s(tsize)
    real(real64) F_t(3, tsize)
    ! real(real64) F_t(tsize)
    real(real64) :: t = 1d-4
    type(trj_results) res

    print *, "reentry parameters:"
    print *, reentry_params

    vector = t_Vector(3)

    vector%vec = [reentry_params%r_earth + 400d0, 0d0, 0d0]
    mycraft = t_SpaceCraft(vector, 100d0, pi)
    reentry = t_ReentryCalc(mycraft)
    print *, "position vec"
    print *, reentry%posi%vec
    
    reentry%V = 7.74d0
    ! reentry%gamma = deg2rad(-1.6d0)
    reentry%gamma = deg2rad(0d0)
    reentry%psi = deg2rad(0d0)
    print *, "V = ", reentry%V
    forall(i=1:tsize) ts(i) = t*i
    dv_s = 0d0
    ! F_t(:) = -1d-3   ! F

    ! call reentry%trj_calc_3d([reentry_params%r_earth + 120d0, 0d0, 0d0], ts, F_t, res)
    call reentry%trj_calc_3d([reentry_params%r_earth + 400d0, 0d0, 0d0], ts, F_t, res)
    print *, "data size: ", size(res%r_hist)
    print *, "last height: ", res%r_hist(ubound(res%r_hist))
    print *, "time: ", size(res%r_hist)*t, "[s]"

    print *, "移動距離(水平) => ", reentry_params%r_earth*(res%theta_hist(ubound(res%r_hist))-res%theta_hist(1))

    open(newunit=nunit, file="trj.csv", iostat=ios, iomsg=img)
    write(nunit, "(4(a,','), a)") "r_hist", "t", "fr", "ftheta", "fphi"  ! header
    do i = 1, ubound(res%r_hist, dim=1)
        write(nunit, "(4(f12.7,','), f12.7)") res%r_hist(i), ts(i), F_t(:,i)
    end do

end program

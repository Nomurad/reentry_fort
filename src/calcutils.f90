module mod_calcutil
    use, intrinsic :: iso_fortran_env 
    implicit none 
    real(real64), parameter :: pi = 4.0d0*atan(1.0d0)

    contains

    real(real64) function get_pi() result(res)
        res = pi 
    end function get_pi

    real(real64) function deg2rad(deg) result(rad)
        real(real64), intent(in) :: deg 

        rad = deg*(pi/180.0d0)

    end function deg2rad 

    real(real64) function rad2deg(rad) result(deg)
        real(real64), intent(in) :: rad 

       deg = rad*(pi/180.0d0)

    end function rad2deg 

    function rot_matrix(theta, axis)
        real(real64) :: rot_matrix(3,3)
        real(real64), intent(in) :: theta 
        integer, intent(in) :: axis
        integer :: i 
        real(real64) :: q 

        rot_matrix = 0
        q = theta 
        if(axis == 1) then
            rot_matrix(1, :) = [ 1d0,    0d0 ,   0d0 ]  ! 1行目
            rot_matrix(2, :) = [ 0d0,  cos(q), sin(q)]  ! 2行目
            rot_matrix(3, :) = [ 0d0, -sin(q), cos(q)]  ! 3行目
        else if(axis == 2) then 
            rot_matrix(1, :) = [ cos(q), 0d0, -sin(q)]
            rot_matrix(2, :) = [   0d0 , 1d0,    0d0 ]
            rot_matrix(3, :) = [ sin(q), 0d0,  cos(q)]
        else if(axis == 3) then 
            rot_matrix(1, :) = [ cos(q), sin(q), 0d0 ]
            rot_matrix(2, :) = [-sin(q), cos(q), 0d0 ]
            rot_matrix(3, :) = [   0d0 ,   0d0 , 1d0 ]
        else
            do i = 1, 3
                rot_matrix(i,i) = 1d0
            end do
        end if

    end function rot_matrix

    real(real64) function norm(vec) result(res)
        real(real64), intent(in) :: vec(:)
        integer i 
        real(real64) tmp(:)

        tmp = [(vec(i)**2, i=1,size(vec)]
        print *, tmp 
        res = sqrt(sum(tmp))

    end function norm
        

end module mod_calcutil


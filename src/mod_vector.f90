module mod_vector
    use, intrinsic :: iso_fortran_env
    implicit none 
    private

    type, public :: t_Vector
        integer :: ndim = 3
        real(real64), pointer, public:: vec(:) => null()
        
        contains
            procedure, public :: getdim
            procedure :: formatted_write
            procedure :: add => vecadd
            procedure, public :: dot_product => vec_dot_product
            procedure, public :: cross_product => vec_cross_product
            generic :: operator(+) => add
            generic :: write(formatted) => formatted_write 
    end type t_Vector

    interface t_Vector 
        module procedure init_vec1, init_vec2
    end interface

    contains
    type(t_Vector) function init_vec1(x, y, z) result(res)
        real(real64), intent(in) :: x, y
        real(real64), intent(in), optional :: z
        real(real64), allocatable :: tmpvec(:)

        if(present(z)) then 
            res%ndim = 3 
        else 
            res%ndim = 2
        end if

        allocate(tmpvec(res%getdim()))
        allocate(res%vec(res%getdim()))
        res%vec(1) = x
        res%vec(2) = y
        if(present(z)) res%vec(3) = z

    end function init_vec1

    type(t_Vector) function init_vec2(ndim) result(res)
        integer, intent(in) :: ndim

        res%ndim = ndim
        allocate(res%vec(ndim))

    end function init_vec2

    integer function getdim(this) result(res)
        implicit none
        class(t_Vector), intent(in) :: this
        
        res = (this%ndim)

    end function getdim

    logical function is_eq_dimention(vec1, vec2) result(res)
        implicit none 
        class(t_Vector), intent(in) :: vec1
        class(t_Vector), intent(in) :: vec2
        
        res = (vec1%getdim() == vec2%getdim())

    end function is_eq_dimention

    type(t_Vector) function vecadd(vec1, vec2) result(res)
        implicit none 
        class(t_Vector), intent(in) :: vec1
        class(t_Vector), intent(in) :: vec2

        if(.not. is_eq_dimention(vec1, vec2)) return 
        res = t_Vector(vec1%getdim())

        ! do i = 1, vec1%ndim 
        !     res%vec(i) = vec1%vec(i) + vec2%vec(i)
        ! end do
        res%vec = vec1%vec + vec2%vec
        
    end function vecadd

    real(real64) function vec_dot_product(vec1, vec2) result(res)
        implicit none 
        class(t_Vector), intent(in) :: vec1
        class(t_Vector), intent(in) :: vec2

        if(.not. is_eq_dimention(vec1, vec2)) return
        res = dot_product(vec1%vec, vec2%vec)
    end function vec_dot_product

    type(t_Vector) function vec_cross_product(vec1, vec2) result(res)
        implicit none
        class(t_Vector), intent(in) :: vec1
        class(t_Vector), intent(in) :: vec2
        real(real64), pointer :: v1(:), v2(:), vr(:)

        if(.not. is_eq_dimention(vec1, vec2)) return 
        if(vec1%ndim == 2) return

        res = t_Vector(vec1%getdim())
        v1 => vec1%vec 
        v2 => vec2%vec 
        vr => res%vec 
        
        vr(1) = v1(2)*v2(3) - v1(3)*v2(2)
        vr(2) = v1(3)*v2(1) - v1(1)*v2(3)
        vr(3) = v1(1)*v2(2) - v1(2)*v2(1)

    end function vec_cross_product
    
    subroutine formatted_write(this,unit,iotype,vlist,iostat,iomsg)
        ! write(*,*) t_Vector でv_Vector%vecを出力するための関数
        class(t_Vector),INTENT(IN) :: this
        INTEGER,INTENT(IN) :: unit
        CHARACTER(*),INTENT(IN) :: iotype
        INTEGER,INTENT(IN) :: vlist(:)
        INTEGER,INTENT(OUT) :: iostat
        CHARACTER(*),INTENT(INOUT) :: iomsg

        if(.not. associated(this%vec)) then 
            iostat = 3
            iomsg = "this vector didn't allocated."
        end if
        ! write(unit,*) unit, iotype, vlist, iomsg
        write(unit, "(3f12.7)", iostat=iostat) this%vec
    end subroutine

end module mod_vector


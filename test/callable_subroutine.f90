subroutine test(n,A) bind(c)
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: n 
    real(8), intent(inout) :: A(n,n)
    real(8), allocatable :: AT(:,:)
    real(8), allocatable :: B(:)
    real(8) :: dat = 100d0
    integer(c_int) i, j

    AT = (A)

    write(*, *) "** write from fortran"
    write(*, *) n 
    do i = 1, n
        write(*, *) (AT(i, j), j=1,n)
    end do

    AT = AT * 10
    A = AT

    allocate(B(n*n))
    open(20, file="input.dat", form="unformatted")
    read(20) B
    close(20)
    print *, "read dat"
    write(*, "(3f10.5)") reshape(B, [3,3])
    print *, "** end fortran"

    open(20, file="output.dat", form="unformatted")
    write(20) AT
    close(20)

end subroutine


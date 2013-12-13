
! To compile:
! gfortran -o test_insol.x ../coord/interp1D.f90 insolation.f90 test_insol.f90

program test_insol 

    use interp1D 
    use insolation 

    integer, parameter :: nlat = 90 
    double precision   :: lats(nlat)
    integer :: i 

    do i = 1, nlat
        lats(i) = -89.d0 + (i-1)*180.d0/dble(nlat)
    end do 

    ! Initialize insolation module data 
    call INI_SINSOL("input") 

end program 
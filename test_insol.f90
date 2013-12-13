
! To compile:
! gfortran -o test_insol.x ../coord/interp1D.f90 insolation.f90 test_insol.f90

program test_insol 

    use interp1D 
    use insolation 

    integer, parameter :: nlat = 90
    integer, parameter :: nday = 360  
    double precision   :: lats(nlat), days(nday), insol(nlat,nday)

    integer :: i, day

    do i = 1, nlat
        lats(i) = -89.d0 + (i-1)*180.d0/dble(nlat)
    end do 

    ! Calculate insolation at these latitudes
    do day = 1, nday 
        days(day) = day 
        insol(:,day) = calc_insol_daily(day,lats,0.d0)
    end do 

end program 
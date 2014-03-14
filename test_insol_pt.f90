
! To compile:
! gfortran -o test_insol.x insolation.f90 test_insol_pt.f90

program test_insol 

    use insolation 

    integer          :: nday, day
    double precision :: lat, time_bp, insol

    nday    = 360 
    lat     = 65.d0 
    time_bp = 0.d0 

    ! Test daily calculations for a given latitude and time
    do day = 1, nday 
        insol = calc_insol_day(day,lat,time_bp)
        write(*,"(i10,f10.2)") day, insol
    end do 

end program 
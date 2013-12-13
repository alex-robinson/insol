
! To compile:
! gfortran -o test_insol.x ../coord/interp1D.f90 insolation.f90 test_insol.f90

program test_insol 

    use interp1D 
    use insolation 
    use ncio 

    integer, parameter :: nlat = 91
    integer, parameter :: nday = 360  
    double precision   :: lats(nlat), days(nday), insol(nlat,nday)
    integer :: i, day
    character(len=256) :: filename 

    do i = 1, nlat
        lats(i) = -90.d0 + (i-1)*180.d0/dble(nlat)
    end do 

    ! Calculate insolation at these latitudes
    do day = 1, nday 
        days(day) = day 
        insol(:,day) = calc_insol_daily(day,lats,0.d0)
    end do 

    filename = "insol_0BP.nc"

    call nc_create(filename)
    call nc_write_dim(filename,"lat",x=lats)
    call nc_write_dim(filename,"day",x=days)
    call nc_write(filename,insol,"insol",dim1="lat",dim2="day")

end program 
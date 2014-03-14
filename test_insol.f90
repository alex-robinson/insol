
! To compile:
! gfortran -o test_insol.x ../coord/interp1D.f90 insolation.f90 test_insol.f90
! gfortran -o test_insol.x ../ncio/ncio3.f90 sinsol_orbit.f test_insol.f90

program test_insol 

    use insolation 
    use ncio 

    implicit none 

    integer, parameter :: nlat = 181
    integer, parameter :: nday = 360  
    double precision   :: lats(nlat), days(nday), insol(nlat,nday)
    integer :: i, day
    character(len=256) :: filename 
    double precision :: lats2D(37,61), insol2D(37,61)
    double precision :: insol_pt 

    ! Load 2D latitudes to calculate 2D insolation
    filename = "../gridding/output/GRL-50KM_TOPO.nc"
    call nc_read(filename,"lat2D",lats2D)

    filename = "insol_0BP.nc"
    call nc_create(filename)
    call nc_write_dim(filename,"lat",x=lats)
    call nc_write_dim(filename,"day",x=days)
    call nc_write_dim(filename,"xc",x=0,nx=37)
    call nc_write_dim(filename,"yc",x=0,nx=61)
    call nc_write(filename,"lat2D",lats2D,dim1="xc",dim2="yc")
    
    do i = 1, nlat
        lats(i) = -90.d0 + (i-1)*180.d0/dble(nlat)
    end do 

    ! Test 1D lat calculations
    do day = 1, nday 
        days(day) = day 
        insol(:,day) = calc_insol_day(day,lats,0.d0)
    end do 
    call nc_write(filename,"insol",insol,dim1="lat",dim2="day")

    ! Test 2D lat calculations 
    !do i = 1, 10
    do day = 1, nday 
        days(day) = day 
        insol2D   = calc_insol_day(day,lats2D,0.d0)
        call nc_write(filename,"insol2D",insol2D,dim1="xc",dim2="yc",dim3="day", &
                      start=[1,1,day],count=[37,61,1])
    end do
    !end do

    ! Test point lat calculations
    do day = 1, nday 
        insol_pt = calc_insol_day(day,65.d0,0.d0)
!         write(*,"(i10,f10.2)") day, insol_pt 
    end do 

end program 
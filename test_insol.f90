
program test_insol 

    use insolation 
    use ncio 

    implicit none 

    integer, parameter :: nlat = 181
    integer, parameter :: nday = 360 
    integer, parameter :: nx = 37, ny = 61 
    double precision   :: lats(nlat), days(nday), insol(nlat,nday)
    integer :: i, day
    character(len=256) :: filename 
    double precision :: lats2D(nx,ny), insol2D(nx,ny)
    double precision :: xc(nx), yc(ny)
    double precision :: insol_pt 

    ! Define vectors of days, lats 
    do day = 1, nday 
        days(day) = day 
    end do 

    do i = 1, nlat
        lats(i) = -90.d0 + (i-1)*181.d0/dble(nlat)
    end do 

    ! Load a 2D array of lats from a stereographic projection 
    filename = "test_data/GRL-50KM_TOPO.nc"
    call nc_read(filename,"xc",xc)
    call nc_read(filename,"yc",yc)
    call nc_read(filename,"lat2D",lats2D)

    ! Initially load module with a simple call
    ! (this reads the orbital values from the input files)
    insol_pt = calc_insol_day(day,65.d0,0.d0)

    ! Test point lat calculations (65 N, year 0)
    write(*,"(a12,a12)") "day", "S (W/m2)" 
    do day = 1, nday 
        insol_pt = calc_insol_day(day,65.d0,0.d0)
        write(*,"(i12,f12.2)") day, insol_pt 
    end do 
    
    ! Prepare NetCDF test output file with dimension variables
    filename = "test_data/insol_0BP.nc"
    call nc_create(filename)
    call nc_write_dim(filename,"lat",x=lats)
    call nc_write_dim(filename,"day",x=days)
    call nc_write_dim(filename,"xc",x=xc)
    call nc_write_dim(filename,"yc",x=yc)
    
    ! Write the 2D latitude field to file
    call nc_write(filename,"lat2D",lats2D,dim1="xc",dim2="yc")

    ! Test 1D vector lat calculations (90 S to 90 N, year 0) 
    do day = 1, nday 
        insol(:,day) = calc_insol_day(day,lats,0.d0)
    end do 
    call nc_write(filename,"insol",insol,dim1="lat",dim2="day")

    ! Test 2D array lat calculations (stereographic projection, year 0) 
    do day = 1, nday 
        insol2D   = calc_insol_day(day,lats2D,0.d0)
        call nc_write(filename,"insol2D",insol2D,dim1="xc",dim2="yc",dim3="day", &
                      start=[1,1,day],count=[37,61,1])
    end do

    write(*,*)
    stop 

end program 


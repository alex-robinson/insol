
! To compile:
! gfortran -o test_insol.x insolation.f90 test_insol_pt.f90

program test_dni 

    use irradiance  
    use ncio 

    integer  :: i, j, k, nx, ny, nt 
    real(wp), allocatable :: lon(:) 
    real(wp), allocatable :: lat(:) 
    real(wp), allocatable :: time(:) 
    real(wp), allocatable :: ghi(:,:,:) 
    real(wp), allocatable :: zenith_angle(:,:,:) 
    real(wp), allocatable :: dni(:,:,:) 

    character(len=512) :: filename_in 
    character(len=512) :: filename_out  

    filename_in  = "data/era5/southspain/era5_southspain_2019.nc" 
    filename_out = "data/era5/southspain/era5_southspain_2019_dni.nc"
    
    nx = nc_size(filename_in,"longitude")
    ny = nc_size(filename_in,"latitude")
    nt = nc_size(filename_in,"time") 

    allocate(lon(nx))
    allocate(lat(ny))
    allocate(time(nt))
    allocate(ghi(nx,ny,nt)) 
    allocate(zenith_angle(nx,ny,nt)) 
    allocate(dni(nx,ny,nt))

    call nc_read(filename_in,"longitude",lon)
    call nc_read(filename_in,"latitude",lat)
    call nc_read(filename_in,"time",time)
    call nc_read(filename_in,"ssrd",ghi)

    time = time - minval(time) 
    ghi  = ghi / 3600.0_wp 

    
    write(*,*) "range(lon):  ", minval(lon),  maxval(lon)
    write(*,*) "range(lat):  ", minval(lat),  maxval(lat)
    write(*,*) "range(time): ", minval(time), maxval(time)
    write(*,*) "range(ghi):  ", minval(ghi),  maxval(ghi)

    zenith_angle = 0.0_wp 
    dni = 0.0_wp 

    ! i = 2 
    ! j = 2 

    !     do k = 1, nt 
    !         call calc_zenith_angle(zenith_angle(i,j,k),time(k),lat(j),lon(i))
    !         call calc_dni_disc(dni(i,j,k),ghi(i,j,k),time(k),lat(j),lon(i))
    !     end do


    do i = 1, nx 
    do j = 1, ny 

        do k = 1, nt 
            call calc_dni_disc(dni(i,j,k),ghi(i,j,k),time(k),lat(j),lon(i))
        end do

    end do 
    end do 


    ! Write to output file 

    call nc_create(filename_out)
    call nc_write_dim(filename_out,"longitude",x=lon)
    call nc_write_dim(filename_out,"latitude",x=lat)
    call nc_write_dim(filename_out,"time",x=time)
    
    
    call nc_write(filename_out,"zenith_angle",  zenith_angle,   dim1="longitude",dim2="latitude",dim3="time")
    call nc_write(filename_out,"ghi",           ghi,            dim1="longitude",dim2="latitude",dim3="time")
    call nc_write(filename_out,"dni",           dni,            dim1="longitude",dim2="latitude",dim3="time")


end program test_dni
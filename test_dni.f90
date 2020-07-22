
! To compile:
! gfortran -o test_insol.x insolation.f90 test_insol_pt.f90

program test_dni 

    use irradiance  
    use ncio 

    integer  :: i, j, k, nx, ny, nt 
    real(wp), allocatable :: lon(:) 
    real(wp), allocatable :: lat(:) 
    real(wp), allocatable :: time(:) 
    real(wp), allocatable :: tisr(:,:,:) 
    real(wp), allocatable :: ghi(:,:,:) 
    real(wp), allocatable :: spres(:,:,:) 
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
    allocate(tisr(nx,ny,nt))
    allocate(ghi(nx,ny,nt)) 
    allocate(spres(nx,ny,nt))
    allocate(zenith_angle(nx,ny,nt)) 
    allocate(dni(nx,ny,nt))

    call nc_read(filename_in,"longitude",lon)
    call nc_read(filename_in,"latitude",lat)
    call nc_read(filename_in,"time",time)
    call nc_read(filename_in,"tisr",tisr)       ! TOA Incident solar radiation
    call nc_read(filename_in,"ssrd",ghi)        ! Surface global horizontal irradiance
    call nc_read(filename_in,"sp",spres)        ! Surface pressure

    time = time - minval(time) 
    tisr = tisr / 3600.0_wp                     ! [J/m2] => [W/m2]
    ghi  = ghi  / 3600.0_wp                     ! [J/m2] => [W/m2]
    
    write(*,*) "range(lon):  ", minval(lon),  maxval(lon)
    write(*,*) "range(lat):  ", minval(lat),  maxval(lat)
    write(*,*) "range(time): ", minval(time), maxval(time)
    write(*,*) "range(tisr): ", minval(tisr),  maxval(tisr)
    write(*,*) "range(ghi):  ", minval(ghi),  maxval(ghi)
    write(*,*) "range(sp):   ", minval(spres),  maxval(spres)

    zenith_angle = 0.0_wp 
    dni          = 0.0_wp 

    do i = 1, nx 
    do j = 1, ny 

        do k = 1, nt
            call calc_zenith_angle(zenith_angle(i,j,k),time(k),lat(j),lon(i)) 
            !call calc_dni_disc(dni(i,j,k),time(k),lat(j),lon(i),ghi(i,j,k),spres(i,j,k),tisr(i,j,k))
            call calc_dni_disc(dni(i,j,k),time(k),lat(j),lon(i),ghi(i,j,k),-1.0_wp,-1.0_wp)
            
        end do

    end do 
    end do 


    ! Write to output file 

    call nc_create(filename_out)
    call nc_write_dim(filename_out,"longitude",x=lon)
    call nc_write_dim(filename_out,"latitude",x=lat)
    call nc_write_dim(filename_out,"time",x=time)
    
    
    call nc_write(filename_out,"zenith_angle",  zenith_angle,   dim1="longitude",dim2="latitude",dim3="time")
    call nc_write(filename_out,"tisr",          tisr,           dim1="longitude",dim2="latitude",dim3="time")
    call nc_write(filename_out,"ghi",           ghi,            dim1="longitude",dim2="latitude",dim3="time")
    call nc_write(filename_out,"spres",         spres,          dim1="longitude",dim2="latitude",dim3="time")
    call nc_write(filename_out,"dni",           dni,            dim1="longitude",dim2="latitude",dim3="time")


end program test_dni
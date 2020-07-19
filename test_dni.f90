
! To compile:
! gfortran -o test_insol.x insolation.f90 test_insol_pt.f90

program test_dni 

    use irradiance  
    use ncio 

    integer          :: nday, day
    double precision :: lat, time_bp, insol

    integer  :: nt 
    real(wp), allocatable :: time(:) 
    real(wp), allocatable :: ghi(:) 

    character(len=512) :: filename_in 
    character(len=512) :: filename_out  

    filename_in = "data/era5/download_2019.nc" 

    nt = nc_size(filename_in,"time") 

    allocate(time(nt))
    allocate(ghi(nt)) 

    call nc_read(filename_in,"time",time)
    call nc_read(filename_in,"ssrd",ghi,start=[2,2,1],count=[1,1,nt])

    time = time - minval(time) 
    ghi  = ghi / 86400.0_wp 
    
    write(*,*) "range(time): ", minval(time), maxval(time)
    write(*,*) "range(ghi):  ", minval(ghi), maxval(ghi)

end program test_dni
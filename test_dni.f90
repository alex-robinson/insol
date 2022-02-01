
! To compile:
! gfortran -o test_insol.x insolation.f90 test_insol_pt.f90

program test_dni 

    use irradiance  
    use ncio 

    
    ! Run a diagnostic test with specific parameter values for testing
    ! (hard-coded inside routine below) 
    !call calc_point_dni_test()


    ! Run a calculation of NREL station data for a given year
    call calc_nrel_station(fldr="data/nrel/NRELNSRDBViewerData",stn="848945",year=1998)

contains 

    subroutine calc_point_dni_test()
        ! Diagnostic testing, e.g., to compare with values in NREL Excel file smalldisc.xls
        ! https://www.nrel.gov/grid/solar-resource/disc.html

        implicit none

        ! Local variables 
        real(wp) :: lat, lon, time_zone
        real(wp) :: hour
        real(wp) :: zenith_angle 
        real(wp) :: pres 
        real(wp) :: ghi
        real(wp) :: dni

        ! === Specify input variable values ===
        
        lon       = -75.0
        lat       =  40.0
        time_zone = -5.0 
        hour      = 48.0
        ghi       = -5.0
        pres      = 840.0 * 100.0  ! [mb] => [Pa]

        ! === Calculate Zenith Angle and DNI ===
        ! (with print_diagnostics=.TRUE., zenith_angle and dni calculation info will be printed to screen)

        !call calc_zenith_angle(zenith_angle,hour,lat,lon,time_zone,print_diagnostics=.TRUE.)
        call calc_dni_disc(dni,hour,lat,lon,time_zone,ghi,pres,-1.0_wp,print_diagnostics=.TRUE.)
        
        stop " === Done. === "

        return

    end subroutine calc_point_dni_test

    subroutine calc_nrel_station(fldr,stn,year)

        implicit none

        character(len=*), intent(IN) :: fldr
        character(len=*), intent(IN) :: stn
        integer, intent(IN) :: year
        
        ! Local variables 
        character(len=512) :: filepath_info
        character(len=512) :: filepath_in
        character(len=512) :: filepath_out
        character(len=4)   :: year_str
        integer  :: id 
        real(wp) :: lat, lon, time_zone, elev
        real(wp) :: hour
        real(wp) :: zenith_angle 
        real(wp) :: pres 
        real(wp) :: ghi
        real(wp) :: dni

        ! Define filenames 
        write(year_str,"(i4)") year 
        filepath_info = trim(fldr)//"/"//trim(stn)//"-info.txt"
        filepath_in   = trim(fldr)//"/"//trim(stn)//"-"//year_str//".txt"
        filepath_out  = trim(fldr)//"/"//trim(stn)//"-dni-"//year_str//".txt"
        
        write(*,*) "filepath_info: ", trim(filepath_info)
        write(*,*) "filepath_in:   ", trim(filepath_in)
        write(*,*) "filepath_out:  ", trim(filepath_out)
        
        stop 

        call calc_zenith_angle(zenith_angle,hour,lat,lon,time_zone,print_diagnostics=.FALSE.) 
        call calc_dni_disc(dni,hour,lat,lon,time_zone,ghi,pres,-1.0_wp,print_diagnostics=.FALSE.)
            
        return

    end subroutine calc_nrel_station

    subroutine calc_era5_dni()

        implicit none 

        ! Local variables 
        character(len=56) :: loc_now 
        character(len=56) :: year_now 
        integer :: narg 

        character(len=512) :: filename_in 
        character(len=512) :: filename_out  

        integer  :: i, j, k, nx, ny, nt 
        real(wp) :: time_zone 
        real(wp), allocatable :: lon(:) 
        real(wp), allocatable :: lat(:) 
        real(wp), allocatable :: time(:)
        real(wp), allocatable :: hour(:) 
        real(wp), allocatable :: tisr(:,:,:) 
        real(wp), allocatable :: ghi(:,:,:) 
        real(wp), allocatable :: fdir(:,:,:) 
        real(wp), allocatable :: pres(:,:,:) 
        real(wp), allocatable :: zenith_angle(:,:,:) 
        real(wp), allocatable :: dni(:,:,:) 

        ! ===== Load location and year from command line arguments =====
        narg = command_argument_count()

        if (narg .ne. 2) then 
            write(*,*) "test_dni:: Error: location and year must be provided as an argument."
            stop 
        end if 

        call get_command_argument(1,loc_now)
        call get_command_argument(2,year_now)

        ! === Define filenames === 

        !filename_in  = "data/era5/southspain/era5_southspain_2019.nc" 
        !filename_out = "data/era5/southspain/era5_southspain_2019_dni.nc"
        
        filename_in  = "data/era5/"//trim(loc_now)//"/era5_"//trim(loc_now)//"_"//trim(year_now)//".nc" 
        filename_out = "data/era5/"//trim(loc_now)//"/era5-dni_"//trim(loc_now)//"_"//trim(year_now)//".nc" 
        
        ! === Get data dimensions === 

        nx = nc_size(filename_in,"longitude")
        ny = nc_size(filename_in,"latitude")
        nt = nc_size(filename_in,"time") 

        allocate(lon(nx))
        allocate(lat(ny))
        allocate(time(nt))
        allocate(hour(nt))
        allocate(tisr(nx,ny,nt))
        allocate(ghi(nx,ny,nt)) 
        allocate(fdir(nx,ny,nt)) 
        allocate(pres(nx,ny,nt))
        allocate(zenith_angle(nx,ny,nt)) 
        allocate(dni(nx,ny,nt))

        ! === Read in dataset === 

        call nc_read(filename_in,"longitude",lon)
        call nc_read(filename_in,"latitude",lat)
        call nc_read(filename_in,"time",time)
        call nc_read(filename_in,"tisr",tisr)       ! TOA Incident solar radiation
        call nc_read(filename_in,"ssrd",ghi)        ! Surface global horizontal irradiance == [GHI]
        call nc_read(filename_in,"fdir",fdir)       ! Surface direct horizontal irradiance == [DNI*cos(zenith_angle)]
        call nc_read(filename_in,"sp",pres)         ! Surface pressure [Pa]

        hour = time - minval(time) + 1.0            ! First hour == 1.0
        tisr = tisr / 3600.0_wp                     ! [J/m2] => [W/m2]
        ghi  = ghi  / 3600.0_wp                     ! [J/m2] => [W/m2]
        fdir = fdir / 3600.0_wp                     ! [J/m2] => [W/m2] 

        ! Assume input data is already in UTC
        time_zone = 0.0_wp 

        write(*,*) "range(lon):  ", minval(lon),  maxval(lon)
        write(*,*) "range(lat):  ", minval(lat),  maxval(lat)
        write(*,*) "range(time): ", minval(time), maxval(time)
        write(*,*) "range(tisr): ", minval(tisr), maxval(tisr)
        write(*,*) "range(ghi):  ", minval(ghi),  maxval(ghi)
        write(*,*) "range(fdir): ", minval(fdir), maxval(fdir)
        write(*,*) "range(sp):   ", minval(pres), maxval(pres)


        ! === Calculate Zenith Angle and DNI ===

        zenith_angle = 0.0_wp 
        dni          = 0.0_wp 

        do i = 1, nx 
        do j = 1, ny 

            do k = 1, nt
                call calc_zenith_angle(zenith_angle(i,j,k),hour(k),lat(j),lon(i),time_zone,print_diagnostics=.FALSE.) 
                call calc_dni_disc(dni(i,j,k),hour(k),lat(j),lon(i),time_zone,ghi(i,j,k),pres(i,j,k),-1.0_wp,print_diagnostics=.FALSE.)
                !call calc_dni_disc(dni(i,j,k),hour(k),lat(j),lon(i),ghi(i,j,k),-1.0_wp,-1.0_wp,print_diagnostics=.FALSE.)
                
            end do

        end do 
        end do 


        ! === Write to output file ===

        call nc_create(filename_out)
        call nc_write_dim(filename_out,"longitude",x=lon)
        call nc_write_dim(filename_out,"latitude",x=lat)
        call nc_write_dim(filename_out,"time",x=time)
        
        call nc_write(filename_out,"zenith_angle",  zenith_angle,   dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"tisr",          tisr,           dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"ghi",           ghi,            dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"fdir",          fdir,           dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"sp",            pres,           dim1="longitude",dim2="latitude",dim3="time")
        call nc_write(filename_out,"dni",           dni,            dim1="longitude",dim2="latitude",dim3="time")

        return

    end subroutine calc_era5_dni

end program test_dni
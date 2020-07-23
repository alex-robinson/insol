module irradiance 

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    ! Missing value and aliases
    real(wp), parameter :: MISSING_VALUE_DEFAULT = real(-9999.0,wp)
    real(wp), parameter :: MISSING_VALUE = MISSING_VALUE_DEFAULT
    real(wp), parameter :: MV = MISSING_VALUE_DEFAULT
    
    ! Error distance (very large), error index, and smallest number epsilon 
    real(wp), parameter :: ERR_DIST = real(1E8,wp) 
    integer,    parameter :: ERR_IND  = -1 
    real(wp), parameter :: tol_underflow = real(1E-15,wp)

    ! Mathematical constants
    real(wp), parameter :: pi  = real(2._dp*acos(0.0_dp),wp)
    real(wp), parameter :: degrees_to_radians = real(pi / 180._dp,wp)  ! Conversion factor between radians and degrees
    real(wp), parameter :: radians_to_degrees = real(180._dp / pi,wp)  ! Conversion factor between degrees and radians
    
    private 
    public :: sp, dp, wp, MISSING_VALUE, MV, pi, degrees_to_radians, radians_to_degrees
    public :: calc_dni_disc
    public :: calc_zenith_angle

contains 

    subroutine calc_dni_disc(dni,hour,lat,lon,ghi,pres,tisr)
        ! Direct Insolation Simulation Code (DISC) model
        ! Maxwell (1987), SERI Technical Report 
        ! https://www.nrel.gov/grid/solar-resource/disc.html
        !def getDNI(ghi=0, hour=0, lat=0, lon=0, press=1013.25, zone=8, debug=False):
        
        implicit none 

        real(wp), intent(OUT) :: dni        ! [W/m2]
        real(wp), intent(IN)  :: hour       ! [Hour of the year]
        real(wp), intent(IN)  :: lat        ! [Degrees North]
        real(wp), intent(IN)  :: lon        ! [Degrees East]
        real(wp), intent(IN)  :: ghi        ! [W/m2]
        real(wp), intent(IN)  :: pres       ! [Pa]
        real(wp), intent(IN)  :: tisr       ! [W/m2] TOA incident solar radiation (horizontal)

        ! Local variables 
        real(wp) :: hour_of_year
        integer  :: day_of_year
        real(wp) :: latitude
        real(wp) :: longitude
        real(wp) :: day_angle
        real(wp) :: time_zone 
        real(wp) :: pressure 

        real(wp) :: I0                  ! Insolation/irradiance top-of-atmosphere (extraterrestial radiation)
        real(wp) :: I0_horizontal       ! Horizontal TOA irrandiance
        real(wp) :: dec 
        real(wp) :: eqt  
        real(wp) :: zenith_angle
        real(wp) :: cos_zenith_angle
        real(wp) :: m_air, kt, kn, knc, del_kn
        real(wp) :: a, b, c 
        real(wp) :: kt2, kt3  

        real(wp), parameter :: S0 = 1370.0_wp 
        real(wp), parameter :: standard_pressure = 101325_wp  

        ! Assign input arguments to local variables 
        hour_of_year = hour
        day_of_year  = int((hour_of_year - 1) / 24) + 1
        latitude     = lat
        longitude    = lon
        day_angle    = (2.0*pi) * (day_of_year - 1) / 365.0_wp

        ! Assume that all input data is in UTC time zone 
        time_zone = 0.0_wp 

        if (pres .gt. 0.0_wp) then 
            ! Use pressure provided as an argument 

            pressure = pres 

        else  
            ! Assume standard pressure for now [Pa]
            pressure = standard_pressure

        end if 
        
        ! Calculate zenith angle and its cosine
        call calc_zenith_angle(zenith_angle,hour_of_year,latitude,longitude)
        cos_zenith_angle = cos(zenith_angle*degrees_to_radians) 

        ! Calculate air mass 
        if (zenith_angle < 80) then 
            m_air = 1.0 / (cos_zenith_angle + 0.15 / pow(93.885 - zenith_angle, 1.253)) * (pressure / standard_pressure)
        else
            m_air = 0.0
        end if

        ! Normal insolation at the top of the atmosphere
        ! following Spencer 1971, ref Maxwell 1987)
        I0 = S0 * (1.00011 + 0.034221 * cos(day_angle) + 0.00128 * sin(day_angle) + &
              0.000719 * cos(2.0 * day_angle) + 0.000077 * sin(2.0 * day_angle))

        if (tisr .ge. 0.0_wp) then 
            ! Use horizontal incident irrandiance from input 

            I0_horizontal = tisr 

        else 
            ! Calculate from parameterized irradiance 

            I0_horizontal = (cos_zenith_angle * I0)

        end if
        
        ! Calculate clearness index kt:
        ! Global horizontal irradiance over horizontal insolation TOA 
        if (m_air > 0) then 
            kt = ghi / I0_horizontal
        else
            kt = 0.0
        end if
        
        if (kt .gt. 1.0) then 
            write(*,*) "kt > 1 !" 
            write(*,*) "kt = ", kt, ghi, I0, I0_horizontal, pressure, cos_zenith_angle, zenith_angle
            stop 
        end if

        ! Calculate powers of kt 
        kt2 = kt*kt 
        kt3 = kt2*kt 

        if (kt .gt. 0.6) then 

            a = -5.743 +  21.77*kt -  27.49*kt2 + 11.56*kt3
            b =  41.40 - 118.50*kt +  66.05*kt2 + 31.90*kt3
            c = -47.01 + 184.20*kt - 222.00*kt2 + 73.81*kt3

        else 

            a =  0.512 - 1.560*kt + 2.286*kt2 - 2.222*kt3
            b =  0.370 + 0.962*kt
            c = -0.280 + 0.932*kt - 2.048*kt2

        end if 

        del_kn = a + b * exp(c * m_air)
        knc    = 0.886 - 0.122*m_air + 0.0121*m_air**2.0 &
                        - 0.000653*m_air**3.0 + 0.000014*m_air**4.0

        kn     = knc - del_kn 

        if (kt .gt. 0.0 .and. I0*kn .ge. 0.0) then 
            dni = I0*kn
        else 
            dni = 0.0_wp 
        end if

        !write(*,*) "testing: ", hour_of_year, zenith_angle, m_air, kt, a, b, c, kn, knc, dni
        !stop 

        return 

    end subroutine calc_dni_disc


    subroutine calc_zenith_angle(zenith_angle,hour,lat,lon)

        implicit none 

        real(wp), intent(OUT) :: zenith_angle 
        real(wp), intent(IN)  :: hour
        real(wp), intent(IN)  :: lat
        real(wp), intent(IN)  :: lon
        
        ! Local variables 
        real(wp) :: hour_of_year
        integer  :: day_of_year
        real(wp) :: latitude
        real(wp) :: longitude
        real(wp) :: day_angle
        real(wp) :: time_zone 
         
        real(wp) :: dec 
        real(wp) :: eqt 
        real(wp) :: hour_angle 
        real(wp) :: cos_zenith_angle 

        real(wp), parameter :: S0 = 1370.0_wp 

        ! Assign input arguments to local variables 
        hour_of_year = hour
        day_of_year  = int((hour_of_year - 1) / 24) + 1
        latitude     = lat
        longitude    = lon
        !day_angle    = (2.0*pi) * min( (day_of_year - 1) / 365.0_wp, 1.0)
        day_angle    = (2.0*pi) * (day_of_year - 1) / 365.0_wp

        ! Assume that all input data is in UTC time zone 
        time_zone = 0.0_wp 

        dec = (0.006918 - 0.399912 * cos(day_angle) + 0.070257 * sin(day_angle) -   &
              0.006758 * cos(2 * day_angle) + 0.000907 * sin(2 * day_angle) -       &
              0.002697 * cos(3 * day_angle) + 0.00148 * sin(3 * day_angle)) * (180.0 / 3.14159)
        eqt = (0.000075 + 0.001868 * cos(day_angle) - 0.032077 * sin(day_angle) -   &
              0.014615 * cos(2 * day_angle) - 0.040849 * sin(2 * day_angle)) * (229.18)
        
        hour_angle = 15 * (hour_of_year - 12 - 0.5 + eqt / 60 + ((longitude - time_zone * 15) * 4) / 60)
        
        ! Calculate cosine of zenith angle 
        cos_zenith_angle =  cos(dec*degrees_to_radians)  * &
                            cos(latitude*degrees_to_radians) * &
                            cos(hour_angle*degrees_to_radians) &
                          + sin(dec*degrees_to_radians) * &
                            sin(latitude*degrees_to_radians)

        ! Finally, calculate the zenith angle [degrees]
        zenith_angle = acos(cos_zenith_angle) * radians_to_degrees

        return 

    end subroutine calc_zenith_angle

    function pow(var,power) result(var_pow)

        implicit none 

        real(wp), intent(IN) :: var
        real(wp), intent(IN) :: power 
        real(wp) :: var_pow 

        var_pow = var**power 

        return 

    end function pow
    
end module irradiance

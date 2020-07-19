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
    

contains 

    subroutine calc_dni_point(dni,ghi,hour,lat,lon,pres,zone)
        !def getDNI(ghi=0, hour=0, lat=0, lon=0, press=1013.25, zone=8, debug=False):
        
        implicit none 

        real(wp), intent(OUT) :: dni 
        real(wp), intent(IN)  :: ghi
        real(wp), intent(IN)  :: hour
        real(wp), intent(IN)  :: lat
        real(wp), intent(IN)  :: lon
        real(wp), intent(IN)  :: pres
        integer,  intent(IN)  :: zone
        
        ! Local variables 
        real(wp) :: hour_of_year
        integer  :: day_of_year
        real(wp) :: latitude
        real(wp) :: longitude
        real(wp) :: day_angle
        real(wp) :: time_zone 
        real(wp) :: pressure 

        real(wp) :: etr 
        real(wp) :: dec 
        real(wp) :: eqt 
        real(wp) :: hour_angle 
        real(wp) :: zenith_angle
        real(wp) :: am, kt, kn, knc
        real(wp) :: a, b, c 

         

        real(wp), parameter :: S0 = 1370.0_wp 

        ! Assign input arguments to local variables 
        hour_of_year = hour
        day_of_year  = int((hour_of_year - 1) / 24) + 1
        latitude     = lat
        longitude    = lon
        day_angle    = (2.0*pi) * min( (day_of_year - 1) / 365.0_wp, 1.0)
        time_zone    = zone
        pressure     = pres

        etr = S0 * (1.00011 + 0.034221 * cos(day_angle) + 0.00128 * sin(day_angle) + &
              0.000719 * cos(2.0 * day_angle) + 0.000077 * sin(2.0 * day_angle))
        dec = (0.006918 - 0.399912 * cos(day_angle) + 0.070257 * sin(day_angle) -  &
              0.006758 * cos(2.0 * day_angle) + 0.000907 * sin(2.0 * day_angle) -  &
              0.002697 * cos(3.0 * day_angle) + 0.00148 * sin(3.0 * day_angle)) * radians_to_degrees
        eqt = (0.000075 + 0.001868 * cos(day_angle) - 0.032077 * sin(day_angle) -  &
              0.014615 * cos(2 * day_angle) - 0.040849 * sin(2.0 * day_angle)) * (229.18)
        hour_angle = 15 * (hour_of_year - 12 - 0.5 + eqt / 60 + ((longitude - time_zone * 15) * 4) / 60)
        zenith_angle = acos(cos(dec*degrees_to_radians) * cos(latitude*degrees_to_radians) * cos(hour_angle*degrees_to_radians) +  &
              sin(dec*degrees_to_radians) * sin(latitude*degrees_to_radians)) * radians_to_degrees
        if (zenith_angle < 80) then 
            am = 1.0 / (cos(zenith_angle*degrees_to_radians) + 0.15 / pow(93.885 - zenith_angle, 1.253)) * (pressure / 1013.25)
        else
            am = 0.0
        end if

        if (am > 0) then 
            kt = ghi / (cos(zenith_angle*degrees_to_radians) * etr)
        else
            kt = 0.0
        end if
        

        a = 0.0
        if (kt > 0) then 
            if (kt > 0.6) then 
                a = -5.743 + 21.77 * kt - 27.49 * pow(kt, 2.0_wp) + 11.56 * pow(kt, 3.0_wp)
            else if (kt < 0.6) then
                a = 0.512 - 1.56 * kt + 2.286 * pow(kt, 2.0_wp) - 2.222 * pow(kt, 3.0_wp)
            end if
        end if

        b = 0.0
        if (kt > 0) then
            if (kt > 0.6) then
                b = 41.4 - 118.5 * kt + 66.05 * pow(kt, 2.0_wp) + 31.9 * pow(kt, 3.0_wp)
            else if (kt < 0.6) then
                b = 0.37 + 0.962 * kt
            end if 
        end if 

        c = 0.0
        if (kt > 0) then
            if (kt > 0.6) then
                c = -47.01 + 184.2 * kt - 222 * pow(kt, 2.0_wp) + 73.81 * pow(kt, 3.0_wp)
            else if (kt < 0.6) then
                c = -0.28 + 0.932 * kt - 2.048 * pow(kt, 2.0_wp)
            end if
        end if 

        kn = 0.0
        if (kt > 0) then
            kn = a + b * exp(c * am)
        end if 

        knc = 0.0
        if (kt > 0) then 
            knc = 0.886 - 0.122 * am + 0.0121 * pow(am, 2.0_wp) - 0.000653 * pow(am, 3.0_wp) + 0.000014 * pow(am, 4.0_wp)
        end if

        if (kt > 0 .and. etr * (knc - kn) >= 0) then 
            dni = etr * (knc - kn)
        else 
            dni = 0.0_wp 
        end if

        return 

    end subroutine calc_dni_point


    function pow(var,power) result(var_pow)

        implicit none 

        real(wp), intent(IN) :: var
        real(wp), intent(IN) :: power 
        real(wp) :: var_pow 

        var_pow = var**power 

        return 

    end function pow
    
end module irradiance

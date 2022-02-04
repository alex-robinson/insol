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
    integer,  parameter :: ERR_IND  = -1 
    real(wp), parameter :: tol_underflow = real(1E-15,wp)

    ! Mathematical constants
    real(wp), parameter :: pi  = real(2._dp*acos(0.0_dp),wp)
    real(wp), parameter :: degrees_to_radians = real(pi / 180._dp,wp)  ! Conversion factor between radians and degrees
    real(wp), parameter :: radians_to_degrees = real(180._dp / pi,wp)  ! Conversion factor between degrees and radians
    
    ! Earth constants
    real(wp), parameter :: S0 = 1370.0_wp   ! Solar constant 

    private 
    public :: sp, dp, wp, MISSING_VALUE, MV, pi, degrees_to_radians, radians_to_degrees
    public :: calc_dni_60min
    public :: calc_dni_disc
    public :: calc_zenith_angle

contains 
    
    subroutine calc_dni_60min(dni,hour,lat,lon,time_zone,ghi,pres,tisr,print_diagnostics)

        implicit none

        real(wp), intent(OUT) :: dni        ! [W/m2]
        real(wp), intent(IN)  :: hour       ! [Hour of the year]
        real(wp), intent(IN)  :: lat        ! [Degrees North]
        real(wp), intent(IN)  :: lon        ! [Degrees East]
        real(wp), intent(IN)  :: time_zone  ! [Hours difference to UTC]
        real(wp), intent(IN)  :: ghi        ! [W/m2]
        real(wp), intent(IN)  :: pres       ! [Pa]
        real(wp), intent(IN)  :: tisr       ! [W/m2] TOA incident solar radiation (horizontal)
        logical,  intent(IN)  :: print_diagnostics

        ! Local variables 
        integer :: k, nmin
        real(wp), allocatable :: hour_vals(:)
        real(wp), allocatable :: dni_vals(:)


        ! Do 60 minutes of calculations by default
        nmin = 60

        allocate(hour_vals(nmin))
        allocate(dni_vals(nmin))

        dni_vals = 0.0_wp 

        do k = 1, nmin 
            hour_vals(k) = hour + (-(nmin/2.0_wp) + real(k,wp)) /60_wp
            call calc_dni_disc(dni_vals(k),hour_vals(k),lat,lon,time_zone,ghi,pres,tisr,print_diagnostics)

            !write(*,*) "calc_dni_60min: ", hour_vals(k), dni_vals(k)

        end do

        ! Return the mean DNI value for this hour
        dni = sum(dni_vals)/real(nmin,wp)

        return

    end subroutine calc_dni_60min

    subroutine calc_dni_disc(dni,hour,lat,lon,time_zone,ghi,pres,tisr,print_diagnostics)
        ! Direct Insolation Simulation Code (DISC) model
        ! Maxwell (1987), SERI Technical Report 
        ! https://www.nrel.gov/grid/solar-resource/disc.html
        !def getDNI(ghi=0, hour=0, lat=0, lon=0, press=1013.25, zone=8, debug=False):
        
        implicit none 

        real(wp), intent(OUT) :: dni        ! [W/m2]
        real(wp), intent(IN)  :: hour       ! [Hour of the year]
        real(wp), intent(IN)  :: lat        ! [Degrees North]
        real(wp), intent(IN)  :: lon        ! [Degrees East]
        real(wp), intent(IN)  :: time_zone  ! [Hours difference to UTC]
        real(wp), intent(IN)  :: ghi        ! [W/m2]
        real(wp), intent(IN)  :: pres       ! [Pa]
        real(wp), intent(IN)  :: tisr       ! [W/m2] TOA incident solar radiation (horizontal)
        logical,  intent(IN)  :: print_diagnostics

        ! Local variables 
        real(wp) :: hour_of_year
        real(wp) :: hour_of_day
        integer  :: day_of_year
        real(wp) :: latitude
        real(wp) :: longitude
        real(wp) :: day_angle 
        real(wp) :: pressure 

        real(wp) :: I0                      ! Insolation/irradiance top-of-atmosphere (extraterrestial radiation)
        real(wp) :: I0_horizontal           ! Horizontal TOA irrandiance
        real(wp) :: dec 
        real(wp) :: eqt  
        real(wp) :: zenith_angle
        real(wp) :: cos_zenith_angle
        real(wp) :: m_air_sl, m_air, kt, kn, knc, del_kn
        real(wp) :: a, b, c 
        real(wp) :: kt2, kt3  

        real(wp) :: dhi 

        real(wp), parameter :: standard_pressure = 101325_wp  

        ! Limit kt to a reasonable value, usually between 0 and 1 (see comments below).
        real(wp), parameter :: kt_max = 1.0_wp 


        ! Assign input arguments to local variables 
        hour_of_year = hour                         ! Note: first hour of year at 00:00:00 UTC == 0
        hour_of_day  = mod(hour_of_year,24.0_wp)    ! Just for diagnostic output
        day_of_year  = floor(hour_of_year/24) + 1
        latitude     = lat
        longitude    = lon
        day_angle    = (2.0*pi) * (day_of_year-1) / 365.0_wp

        if (pres .gt. 0.0_wp) then 
            ! Use pressure provided as an argument 

            pressure = pres 

        else  
            ! Assume standard pressure for now [Pa]
            pressure = standard_pressure

        end if 
        
        ! Calculate zenith angle and its cosine
        call calc_zenith_angle(zenith_angle,hour_of_year,latitude,longitude, &
                                                        time_zone,print_diagnostics)
        cos_zenith_angle = cos(zenith_angle*degrees_to_radians) 

        ! Calculate air mass 
        if (zenith_angle < 87) then 
            
            ! Use Bird & Hulstrom coefficient values (DEFAULT)
            ! Bird, R., Hulstrom, R., A Simplified Clear Sky Model for Direct and Diffuse Insolation on Horizontal Surfaces, 1981
            ! NREL SERI/TP-642-761
            m_air_sl = 1.0 / (cos_zenith_angle + 0.15 * (93.885-zenith_angle)**(-1.253))
        
            ! Use Kasten / Sandia Version
            ! Kasten, F. Discussions on the Relative Air Mass, Light Research Technology 25, 129.
            !m_air_sl = 1.0 / (cos_zenith_angle + 0.5057 * (96.080-zenith_angle)**(-1.634))
            
        else
            m_air_sl = 0.0
        end if

        ! Scale by pressure too to get actual estimated air mass
        m_air = m_air_sl * (pressure / standard_pressure)

        ! Insolation at the top of the atmosphere
        ! following Spencer (1971), ref Maxwell (1987)
        I0 = S0 * (1.00011 + 0.034221 * cos(day_angle) + 0.00128 * sin(day_angle) + &
              0.000719 * cos(2.0 * day_angle) + 0.000077 * sin(2.0 * day_angle))

if (.FALSE.) then
    ! Calculate clearness index inline 

        if (tisr .ge. 0.0_wp) then 
            ! Use TOA horizontal incident irrandiance from input 

            I0_horizontal = tisr 

        else 
            ! Calculate from parameterized irradiance 

            I0_horizontal = (cos_zenith_angle * I0)

        end if
        
        ! Calculate clearness index kt:
        kt = ghi / I0_horizontal

        ! Limit kt to a reasonable value. As the zenith_angle
        ! approaches 90 degrees (ie, at sunrise), then 
        ! I0_horizontal gets reduced and the result can
        ! be  kt >> 1. In this case, the factor del_kn below
        ! can give an infinity value. Thus, it is better to
        ! limit kt to a reasonable range like 0 to 1.
        kt = min(kt,kt_max)
else

    ! Use subroutine to calculate clearness index

        call calc_clearness_index(kt,hour_of_year,zenith_angle,ghi,kt_max)

        I0_horizontal = missing_value 

end if
        


select case("disc")


case("disc")

        ! Calculate powers of kt 
        kt2 = kt*kt 
        kt3 = kt2*kt 

        ! Calculate coefficients A, B and C 
        if (kt .gt. 0.6) then 

            a = -5.743 +  21.77*kt -  27.49*kt2 + 11.56*kt3
            b =  41.40 - 118.50*kt +  66.05*kt2 + 31.90*kt3
            c = -47.01 + 184.20*kt - 222.00*kt2 + 73.81*kt3

        else if (kt .gt. 0.0) then

            a =  0.512 - 1.560*kt + 2.286*kt2 - 2.222*kt3
            b =  0.370 + 0.962*kt
            c = -0.280 + 0.932*kt - 2.048*kt2

        else 

            a = 0.0
            b = 0.0 
            c = 0.0 

        end if 

        ! Calculate the clear-sky atmospheric transmissivity 

        if (kt .gt. 0.0) then 
            knc    = 0.886 - 0.122*m_air + 0.0121*m_air**2.0 &
                            - 0.000653*m_air**3.0 + 0.000014*m_air**4.0
            
            del_kn = a + b * exp(c * m_air)
            
            kn     = max(knc - del_kn,0.0_wp)

            ! Finally calculate DNI 
            dni = I0*kn 

        else 

            knc    = 0.0 
            del_kn = 0.0 
            kn     = 0.0 

            ! Set DNI to zero
            dni = 0.0 

        end if 

case("brl")

    ! AJR: This doesn't seem to work yet!!!!
    ! It seems that dividing by cos(zenith_angle) may 
    ! give problems. Needs further investigation. 

        write(*,*) "calc_dni:: BRL method not ready yet!" 
        stop 



        ! Calculate diffuse horizontal irradiance (dhi)
        call calc_dhi_brl(dhi,kt,ghi)

        ! Calculate dni from other components 
        call calc_dni_from_diffuse(dni,ghi,dhi,zenith_angle)

end select

    
        ! Further corrections
        if (ghi .lt. 0.0 .or. dni .lt. 0.0) then 

            dni = 0.0

        end if 

        if (print_diagnostics) then

            write(*,*) "=== calc_dni_disc ==="
            write(*,*)
            
            write(*,*) "hour_of_year     = ", hour_of_year 
            write(*,*) "hour_of_day      = ", hour_of_day 
            write(*,*) "lat              = ", lat 
            write(*,*) "lon              = ", lon 
            write(*,*) "time_zone        = ", time_zone 
            write(*,*) "zenith_angle     = ", zenith_angle 
            write(*,*) "cos_zenith_angle = ", cos_zenith_angle 
            write(*,*) "pressure         = ", pressure
            write(*,*) "pressure ratio   = ", (pressure / standard_pressure)

            write(*,*) "m_air            = ", m_air 
            write(*,*) "ghi              = ", ghi
            write(*,*) "I0               = ", I0
            write(*,*) "I0_horiz         = ", I0_horizontal
            write(*,*) "kt               = ", kt 
            write(*,*) "a                = ", a 
            write(*,*) "b                = ", b 
            write(*,*) "c                = ", c 
            write(*,*) "del_kn           = ", del_kn
            write(*,*) "kn               = ", kn 
            write(*,*) "knc              = ", knc 
            write(*,*) "dni              = ", dni
            write(*,*) 

        end if 

        return 

    end subroutine calc_dni_disc

    subroutine calc_dni_boland2013(dni,ghi)
        ! Logistic model for DNI calculation
        ! Boland et al. (2013, Renewable and Sustainable Energy Reviews)
        ! https://www.sciencedirect.com/science/article/pii/S1364032113005637

        implicit none

        real(wp), intent(OUT) :: dni 
        real(wp), intent(IN)  :: ghi 

        ! Local variables 
        real(wp) :: p 
        real(wp) :: day_angle 
        real(wp) :: cos_zenith_angle
        real(wp) :: I0 
        real(wp) :: I0_horizontal
        real(wp) :: kt 

        real(wp), parameter :: N = 0.006
        real(wp), parameter :: M = 4.38 

        real(wp), parameter :: b1 = 7.75
        real(wp), parameter :: b2 = 1.185
        real(wp), parameter :: b3 = 1.05
        real(wp), parameter :: b4 = 0.004
        real(wp), parameter :: b5 = -0.003



        
        !p = -b1*kt - b2*ast - b3*alpha - b4*kt_day + b5*psi 

        dni = (N*M) / (N+(M-N)*exp(p))

        return

    end subroutine calc_dni_boland2013

    subroutine calc_dni_from_diffuse(dni,ghi,dhi,zenith_angle)

        implicit none

        real(wp), intent(OUT) :: dni 
        real(wp), intent(IN)  :: ghi
        real(wp), intent(IN)  :: dhi 
        real(wp), intent(IN)  :: zenith_angle

        ! Boland et al (2013), Eq. 7 
        dni = (ghi - dhi) / cos(zenith_angle*pi/180.0_wp)

        return

    end subroutine calc_dni_from_diffuse

    subroutine calc_dhi_brl(dhi,kt,ghi)
        ! Logistic model for DHI calculation
        ! Ridley et al. (2010, Renewable Energy)
        ! https://www.sciencedirect.com/science/article/pii/S0960148109003012

        implicit none

        real(wp), intent(OUT) :: dhi 
        real(wp), intent(IN)  :: kt 
        real(wp), intent(IN)  :: ghi

        ! Local variables 
        real(wp) :: p 
        real(wp) :: diffuse_frac 

        ! Simple logistic function that depends only on kt 
        ! Ridley et al (2010), Eq. 3

        p   = -5.0033 + 8.6025*kt 
        diffuse_frac = 1.0_wp / (1.0_wp + exp(p))

        dhi = diffuse_frac*ghi

        return

    end subroutine calc_dhi_brl

    subroutine calc_clearness_index(kt,hour_of_year,zenith_angle,ghi,kt_max)

        implicit none

        real(wp), intent(OUT) :: kt
        real(wp), intent(IN)  :: hour_of_year
        real(wp), intent(IN)  :: zenith_angle
        real(wp), intent(IN)  :: ghi
        real(wp), intent(IN)  :: kt_max 

        ! Local variables 
        real(wp) :: day_of_year
        real(wp) :: day_angle 
        real(wp) :: cos_zenith_angle
        real(wp) :: I0 
        real(wp) :: I0_horizontal

        ! Calculate some information
        day_of_year  = floor(hour_of_year/ 24) + 1
        day_angle    = (2.0*pi) * (day_of_year-1) / 365.0_wp

        ! Calculate insolation at the top of the atmosphere
        ! following Spencer (1971), ref Maxwell (1987)
        I0 = S0 * (1.00011 + 0.034221 * cos(day_angle) + 0.00128 * sin(day_angle) + &
              0.000719 * cos(2.0 * day_angle) + 0.000077 * sin(2.0 * day_angle))

        ! Calculate global horizontal irradiance at TOA
        cos_zenith_angle = cos(zenith_angle*degrees_to_radians) 
        I0_horizontal    = (cos_zenith_angle * I0)

        ! Calculate clearness index kt:
        ! Global horizontal irradiance over horizontal insolation TOA 
        kt = ghi / I0_horizontal

        ! Limit kt to reasonable values (kt_max ~ 1.0)
        kt = min(kt,kt_max)

        return

    end subroutine calc_clearness_index

    subroutine calc_clearness_index_day()

        implicit none


        return

    end subroutine calc_clearness_index_day
    
    subroutine calc_zenith_angle(zenith_angle,hour,lat,lon,time_zone,print_diagnostics)

        implicit none 

        real(wp), intent(OUT) :: zenith_angle 
        real(wp), intent(IN)  :: hour
        real(wp), intent(IN)  :: lat
        real(wp), intent(IN)  :: lon
        real(wp), intent(IN)  :: time_zone      ! Relative to UTC (0 == UTC)
        logical,  intent(IN)  :: print_diagnostics

        ! Local variables 
        real(wp) :: hour_of_year
        real(wp) :: hour_of_day
        integer  :: day_of_year
        real(wp) :: latitude
        real(wp) :: longitude
        real(wp) :: day_angle

        real(wp) :: dec 
        real(wp) :: eqt 
        real(wp) :: hour_angle 
        real(wp) :: cos_zenith_angle 
        
        ! Assign input arguments to local variables 
        hour_of_year = hour                         ! Note: first hour of year at 00:00:00 UTC == 0
        hour_of_day  = mod(hour_of_year,24.0_wp)    ! Just for diagnostic output
        day_of_year  = floor(hour_of_year/24) + 1
        latitude     = lat
        longitude    = lon
        day_angle    = (2.0*pi) * (day_of_year-1) / 365.0_wp

        call calc_hour_angle(hour_angle,hour_of_year,longitude,time_zone)

        dec = (0.006918 - 0.399912 * cos(day_angle) + 0.070257 * sin(day_angle) -   &
              0.006758 * cos(2 * day_angle) + 0.000907 * sin(2 * day_angle) -       &
              0.002697 * cos(3 * day_angle) + 0.00148 * sin(3 * day_angle)) * (180.0 / 3.14159)
        
        ! Calculate cosine of zenith angle 
        cos_zenith_angle =  cos(dec*degrees_to_radians)  * &
                            cos(latitude*degrees_to_radians) * &
                            cos(hour_angle*degrees_to_radians) &
                          + sin(dec*degrees_to_radians) * &
                            sin(latitude*degrees_to_radians)

        ! Finally, calculate the zenith angle [degrees]
        zenith_angle = acos(cos_zenith_angle) * radians_to_degrees

        if (print_diagnostics) then

            write(*,*) "== calc_zenith_angle =="
            write(*,*)
            write(*,*) "hour_of_year     = ", hour_of_year 
            write(*,*) "hour_of_day      = ", hour_of_day 
            write(*,*) "time_zone        = ", time_zone 
            write(*,*) "lat              = ", lat 
            write(*,*) "lon              = ", lon 
            write(*,*) "day_of_year      = ", day_of_year 
            write(*,*) "day_angle        = ", day_angle 
            write(*,*) "hour_angle       = ", hour_angle 
            write(*,*) "dec              = ", dec 
            write(*,*) "eqt              = ", eqt 
            write(*,*) "cos_zenith_angle = ", cos_zenith_angle 
            write(*,*) "zenith_angle     = ", zenith_angle
            write(*,*)  
        end if 

        return 

    end subroutine calc_zenith_angle

    subroutine calc_hour_angle(hour_angle,hour_of_year,longitude,time_zone)

        implicit none 

        real(wp), intent(OUT) :: hour_angle 
        real(wp), intent(IN)  :: hour_of_year
        real(wp), intent(IN)  :: longitude
        real(wp), intent(IN)  :: time_zone      ! Relative to UTC (0 == UTC)

        ! Local variables 
        integer  :: day_of_year
        real(wp) :: day_angle

        real(wp) :: eqt 

        ! Assign input arguments to local variables 
        day_of_year  = floor(hour_of_year/24) + 1
        day_angle    = (2.0*pi) * (day_of_year-1) / 365.0_wp

        eqt = (0.000075 + 0.001868 * cos(day_angle) - 0.032077 * sin(day_angle) -   &
              0.014615 * cos(2 * day_angle) - 0.040849 * sin(2 * day_angle)) * (229.18)
        
        hour_angle = 15.0 * (hour_of_year - 12.0 - 0.5 + eqt / 60.0 + ((longitude - time_zone * 15.0) * 4.0) / 60.0)
        
        return 

    end subroutine calc_hour_angle

end module irradiance

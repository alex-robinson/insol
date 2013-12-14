
module insolation

    use interp1D 

    implicit none 

    integer,  parameter :: dp  = kind(1.0d0)

    real (dp), parameter :: pi = 2.d0*acos(0.d0)
    real (dp), parameter :: deg_to_rad = pi/180.0_dp
        
    integer, parameter                   :: n_orbit = 6001 
    real (dp), dimension(n_orbit) :: ECCM, PERM, XOBM 
    real (dp)                     :: ECCP, XOBP, PERP 

    ! Vector for generic insol calculations at predefined latitudes (-90:90)
    integer, parameter          :: NLAT0 = 91
    real (dp), dimension(NLAT0) :: LATS0, INSOL0

    integer :: init_sinsol = 0 

    interface calc_insol_day 
        module procedure calc_insol_day_1D, calc_insol_day_2D 
    end interface 

    private
    public :: calc_insol_day 

contains 

    function calc_insol_day_1D(day,lats,time_bp,S0,day_year) result(insol)
        ! Given day of year, latitudes and time before present (1950),
        ! Calculate the daily insolation values at each latitude 

        implicit none 

        integer   :: day
        real (dp) :: lats(:) 
        real (dp) :: time_bp
        real (dp) :: insol(size(lats))

        integer, parameter :: nh = 24 

        real (dp) :: ECC, XOBCH, TPERI, ZAVEXPE
        real (dp) :: PZEN3, PRAE, PCLOCK, PYTIME 
        real (dp) :: PDISSE(nh), PZEN1(nh), PZEN2(nh)

        integer :: j, h 

        real (dp), optional :: S0
        real (dp)           :: S0_value 

        integer, optional   :: day_year
        integer             :: day_max 

        ! If necessary, initialize input data
        if ( init_sinsol .eq. 0 ) then
            call INI_SINSOL("input")
            init_sinsol = 1
        end if
        
        ! Define the solar constant
        S0_value = 1365.0_dp
        if (present(S0)) S0_value = S0  

        ! Define year length in days 
        day_max = 360
        if (present(day_year)) day_max = day_year


        ! Call Berger subroutine to calculate 
        ! orbital parameters for current year 
        call BERGOR(time_bp,ECC,XOBCH,TPERI,ZAVEXPE)

        ! Get hourly zenith angle values for input into daily insol function
        PYTIME = dble(day)*2.0_dp*pi/dble(day_max)
        do h = 1, nh
            PCLOCK=dble(h)*2.0_dp*pi/dble(nh) 
            call ORBIT(ECC,XOBCH,TPERI,ZAVEXPE, &
                       PCLOCK,PYTIME,PDISSE(h),PZEN1(h),PZEN2(h),PZEN3,PRAE)

!             write(*,*) h, PDISSE(h), PZEN1(h), PZEN2(h)
        end do 
        
        ! First calculate daily insolation at predefined latitude values
        do j = 1, NLAT0
            INSOL0(j) = calc_insol_day_internal(LATS0(j), &
                                         PDISSE,PZEN1,PZEN2,S0_value)
        end do 

        ! Use spline interpolation to get insolation at desired latitudes
        insol = interp_spline(LATS0,INSOL0,lats)
        !insol = INSOL0 

        return 

    end function calc_insol_day_1D

    function calc_insol_day_2D(day,lats,time_bp,S0,day_year) result(insol2D)

        integer   :: day
        real (dp) :: lats(:,:) 
        real (dp) :: time_bp
        real (dp) :: insol2D(size(lats,1),size(lats,2))
        real (dp), allocatable :: lats1D(:), insol1D(:)
        integer   :: n 

        real (dp), optional :: S0
        integer, optional   :: day_year

        n = size(lats,1)*size(lats,2)
        allocate(lats1D(n),insol1D(n))

        lats1D = reshape(lats,[n])
        insol1D = calc_insol_day_1D(day,lats1D,time_bp,S0,day_year)
        insol2D = reshape(insol1D,[size(lats,1),size(lats,2)])

        return 

    end function calc_insol_day_2D

    function calc_insol_day_internal(lat,PDISSE,PZEN1,PZEN2,S0) result(solarm)
        ! Given latitude and sun hourly angle,
        ! Calculate the daily insolation value 

        implicit none 

        integer   :: day, h, nh   
        real (dp), intent(IN) :: lat, S0  
        real (dp), intent(IN) :: PDISSE(:), PZEN1(:), PZEN2(:)
        real (dp) :: coszm, cosp, cosn, S 
        real (dp) :: solarm 

        ! How many hours in the day?
        nh = size(PDISSE)

        ! Daily insolation is calculated by hourly integration for each day
        solarm = 0.0_dp
        coszm  = 0.0_dp 

        do h = 1, nh 

            cosp = PZEN1(h)*dsin(lat*deg_to_rad) + PZEN2(h)*dcos(lat*deg_to_rad)
            cosn = max(cosp,0.0_dp)

            S = S0*cosn*PDISSE(h)
            solarm = solarm + S
            coszm  = coszm + S*cosn

        end do

        ! Get daily average insolation and zenith angle 
        solarm = solarm/24.0_dp
        if (solarm .gt. 0.0_dp) then
            coszm = coszm/(solarm*24.0_dp)
        else
            coszm = 0.0_dp 
        endif

        return 

    end function calc_insol_day_internal

    subroutine INI_SINSOL(input_dir)  
    !********************************************************************  
    !
    !  Purpose:  Calculation of solar insolation and averaged zenite
    !            angle for different times
    ! 
    !  By: A.Ganopolski 
    !  Last modification: 28.04.2011 
    !                                                CLIMBER-2.4
    ! Global variables: ECCM(6001),PERM(6001),XOBM(6001)
        real :: nnf
        character (len=*)   :: input_dir
        character (len=512) :: filen, filep

        integer :: n 

        ! Load Berger orbital parameters
!        open (1,file=INP/orbit_par_5000.txt')

!        do n=1,5001
!            read (1,'(I5,3F10.5)') nn ,ECCM(n), PERM(n), XOBM(n)
!        end do 

        ! Load Laskar2004 input parameters instead of Berger   
        ! Note: year 0: n=5001 ; year -5000 ka: n=1; year 1000 ka: n=6001

        filen = trim(input_dir)//'/LA2004.INSOLN.5Ma.txt'
        filep = trim(input_dir)//'/LA2004.INSOLP.1Ma.txt'

        ! ==== PAST (NEGATIVE) TIMES ====
        write(*,*) "Loading orbital parameters from file: "//trim(filen)
        open (1,file=trim(filen))

        do n=1,5001
            read (1,*) nnf, ECCM(5001-n+1), XOBM(5001-n+1), PERM(5001-n+1)
        enddo

        close (1)

        ! ==== FUTURE (POSITIVE) TIMES ====
        write(*,*) "Loading orbital parameters from file: "//trim(filep)
        open (1,file=trim(filep))

        do n=1,1001
            read (1,*) nnf, ECCM(5000+n), XOBM(5000+n), PERM(5000+n)
        enddo

        close (1)

        ! Convert from radians to degrees      
        XOBM = XOBM * 180.0/pi
        PERM = PERM * 180.0/pi

        ! Also initialize the default latitude vector for insol calculations
        ! (actual lat values are interpolated from these generic calcs)
        do n = 1, NLAT0
            LATS0(n) = -90.0_dp + (n-1)*180.0_dp/dble(NLAT0-1)
        end do 

        return

    end subroutine INI_SINSOL

    subroutine ORBIT(ECC,XOBCH,TPERI,ZAVEXPE,  &
                     PCLOCK,PYTIME,PDISSE,PZEN1,PZEN2,PZEN3,PRAE)
    !   *ORBIT* - FOR SOLAR ORBITAL PARAMETERS.
!
!       J.F.GELEYN     E.C.M.W.F.     03/06/82.
!
!       CHANGED TO 18000YBP-CONDITIONS BY LAUTENSCHLAGER
!                       MPI F. METEOR. HAMBURG  4.87
!       CHANGED TO OPTIONAL ORBITAL PARAMETERS BY GRAF, MPI 1/91
!       CHANGED AGAIN TO A MORE EXACT SOLUTION BY HOFFMANN, MPI 7/93
!       REVISED BY LORENZ, UNIVERSITY OF BREMEN, 10/93
!       -- last change: 94-08-03
!       -- Converted to FORTRAN90 by A. Robinson, 12/2013
!
!       PURPOSE.
!       --------
!
!            THIS ROUTINE COMPUTES THE SOLAR CONSTANT NORMALISED BY ITS
!       ANNUAL MEAN AND THREE ORBITAL PARAMETERS DEPENDING ON THE TIME
!       OF THE DAY AS WELL AS OF THE YEAR (BOTH IN RADIANS). TIME SCALE
!       ORIGIN IS 01/01 00 GMT OF THE DESIGNED TIME DISC. THE DIFFERENT
!       TIME DISCS ARE SCALED ON 03/21 12 GMT AS THE TIME OF THE VERNAL
!       EQUINOX. ALSO RETURNED IS A CONSTANT FOR THE EFFECT OF THE EARTH'S
!       CURVATURE ON THE COSINE OF THE SOLAR ZENITH ANGLE.
!
!       THE ORBITAL PARAMETERS CAN BE TAKEN FROM PROGRAM BERGOR, WRITTEN
!       BY LORENZ, UNIVERSITY OF BREMEN, 8/94
!
!       INTERFACE.
!       ----------
!
!            *ORBIT* IS CALLED FROM *PHYSC* AT THE FIRST LATITUDE ROW.
!            THERE ARE SEVEN DUMMY ARGUMENTS: *PCLOCK* IS THE TIME OF THE
!       DAY
!                                             *PYTIME* IS THE TIME OF THE
!       YEAR (BOTH IN RADIANS).
!                                             *PDISSE* IS THE RATIO OF THE
!       SOLAR CONSTANT TO ITS ANNUAL MEAN
!                                             *PZEN1*, *PZEN2* AND *PZEN3*
!       ARE ZENITHAL PARAMETERS.
!                                             *PRAE* IS THE RATIO OF THE
!       HEIGHT OF THE EQUIVALENT ATMOSPHERE TO THE RADIUS OF THE EARTH.
!
!       METHOD.
!       -------
!
!            STAIGHTFORWARD. INTERMEDIATE VARIABLES ARE THE SOLAR
!       DECLINATION AND THE EQUATION OF TIME.
!
!       EXTERNALS.
!       ----------
!
!            NONE.
!
!       REFERENCE.
!       ----------
!
!            NONE.
!
!
!       DATA STATEMENTS.
!       ----------------
!
!            *ZCDIS*, *ZCDEC* AND *ZCEQT* ARE ARRAYS FOR A SECOND ORDER
!       *FOURIER DEVELOPMENT FOR RESPECTIVELY: SOLAR CONSTANT, SOLAR
!       DECLINATION AND EQUATION OF TIME.
!       THEY ARE VALID FOR TODAY'S ORBITAL PARAMETERS ONLY (LORENZ, 11/93).
!
!       DIMENSION ZCDIS(5),ZCDEC(5),ZCEQT(5)
!       DATA ZCDIS/+1.000110,+0.034221,+0.001280,+0.000719,+0.000077/
!       DATA ZCDEC/+0.006918,-0.399912,+0.070257,-0.006758,+0.000907/
!       DATA ZCEQT/+0.000075,+0.001868,-0.032077,-0.014615,-0.040849/
!
!       *ZRAE* IS THE VALUE FOR *PRAE*.
!
        real (dp), parameter :: ZRAE = 0.1277e-02_dp   

        real(dp) :: TTROP 

        real (dp) :: BTIME, ECC, XOBCH, TPERI, ZAVEXPE
        real (dp) :: PCLOCK, PYTIME, PDISSE, PZEN1, PZEN2, PZEN3, PRAE

        real (dp) :: ZCLOCK, ZYTIME

        real (dp) :: time, eold, enew, eps, epc, zeps, cose, E
        real (dp) :: ZDISSE, zsqecc, ztgean, znu, zlambda, xobche
        real (dp) :: zsinde, zdecli, ZZEN1, ZZEN2, ZZEN3 
        integer   :: niter 

!       1. PRELIMINARY SETTINGS
!          --------------------

        ZCLOCK=PCLOCK
        ZYTIME=PYTIME

!       TROPICAL YEAR

        TTROP=360.0_dp 
!
!       DAY OF VERNAL EQUINOX: DEFINED ON 21. OF MARCH, 12.00 GMT
!       (BEGIN OF YEAR: 01. OF JANUARY, 0.00 GMT: ZYTIME=0.0, THEREAFTER
!       E. G. 2. OF JANUARY 12.00 GMT: TIME IN DAYS FROM BEGIN OF YEAR = 1.5)
!       - NOT USED IN THE EXACT FORMULATION, *ZAVEXPE* USED INSTEAD
!       TVEREX=80.5              ! TROPICAL YEAR: 360 DAYS
!       TVEREX=79.5              ! TROPICAL YEAR: 365 DAYS

!        ***  ORBITAL PARAMETERS : VALUES FOR DIFFERENT TIME DISCS
!             THESE VALUES CAN BE TAKEN FROM PROGRAM BERGOR (S. LORENZ)
!
!       2. COMPUTATIONS
!          ------------

!       ZC1YT=COS(ZYTIME)
!       ZS1YT=SIN(ZYTIME)
!       ZC2YT=ZC1YT**2-ZS1YT**2
!       ZS2YT=2.*ZS1YT*ZC1YT
!
!       CALCULATE ECCENTRI!   ANOMALY:
!
!       LESS EXACT METHOD OF H. GRAF: SIN(E)=E
!       E=(ZYTIME - 2*PI*TPERI/TTROP) / (1.-ECC)
!
!       EXACT CALCULATION WITH KEPLER'S LAW (ELLIPSE)
!       (MONIN, 1986 : AN INTRODUCTION TO THE THEORY OF CLIMATE)
!       USE FIRST DERIVATIVE (NEWTON) TO CALCULATE SOLUTION OF
!       EQUATION OF ECCENTRI!   ANOMALY *E*
!
        time  = ZYTIME-2*pi*TPERI/TTROP 
        eold  = time/(1.-ecc)
        enew  = time
        eps   = 1.e-6_dp 
        niter = 0

        do while ( dabs(zeps) .gt. eps ) 
            zeps = eold - enew
            if (niter .ge. 30) then 
                write(*,*) ' SUBROUTINE *ORBIT* - eccentri!   anomaly not found!'
                write(*,*) ' ERROR IN   *ORBIT* -- STOP'
                stop
            end if 

            niter = niter+1
            eold  = enew
            cose  = cos(enew)
            enew  = (time+ecc*(sin(enew)-enew*cose))/(1.0_dp-ecc*cose)
        end do 

        E=enew
        ZDISSE=(1.0_dp/(1.0_dp-ECC*cos(E)))**2

!       Change in the calculation of the declination.
!       Used are not the formulas from Paltridge/Platt as in the ECHAM model
!       with fixed constants for contemporanious conditions
!       but the exact equations from Monin (s.a. - Keplers law)
!       are solved. Day of vernal equinox is fixed for a 360 day year on the
!       21. Maerz noon, with start of year at 01/01/00 GMT (s.a.).

        zsqecc=sqrt((1+ecc)/(1-ecc))
        ztgean=tan(E/2)

!       znu: true anomaly (actual angle of Earth's position from perihelion)
!       zlambda: true longitude of the Earth (actual angle from vernal equinox)

        znu=2.*atan(zsqecc*ztgean)
        zlambda=znu+zavexpe
        xobche=xobch/180.*pi
        zsinde=sin(xobche)*sin(zlambda)
        zdecli=asin(zsinde)

        ZZEN1=sin(zdecli)
        ZZEN2=cos(zdecli)*cos(zdecli)
        ZZEN3=cos(zdecli)*sin(zdecli)

        ! 3. RETURN
        !    ------

        PDISSE=ZDISSE
        PZEN1=ZZEN1
        PZEN2=ZZEN2
        PZEN3=ZZEN3
        PRAE=ZRAE
 
        return

    end subroutine ORBIT 

    subroutine BERGOR(BTIME,ECC,XOB,TPER,ZAN)
    !**** *BERGOR* - calculates earth's orbital parameters for ECHAM use
    !
    !   S. J. Lorenz   University of Bremen    08-03-94.
    !     -- last change:  31-07-96
    !
    !   PURPOSE.
    !   --------
    !
    !      This programm uses a simple algorithm from A. Berger to compute
    !   earth's orbital parameters of any interesting time disc. The
    !   original parameters are modified to the values needed by new
    !   subroutine "ORBIT" (by S. Lorenz) of AGCM "ECHAM" of MPI/DKRZ
    !   in Hamburg.
    !
    !
    !   INTERFACE.
    !   ----------
    !      none.
    !
    !   METHOD.
    !   -------
    !      The program BERGOR is selfexplanatory.
    !
    !   REFERENCE.
    !   ----------
    !      Berger, A.L., 1978: Long-term variations of daily insolation and
    !         Quaternary climati!   changes, J. Atm. Sci., 35, 2362-2367.

        real (dp) :: BTIME, ECC, XOB, TPER, ZANE
        real (dp) :: ang, PER, tgnu, zan, sqecc, tian
        real (dp) :: days, e, epc

        ! Get pre-calculated orbital parameters with subroutine ORBIT_PAR
        call ORBIT_PAR(BTIME,PER,ECC,XOB)     

        !  compute new orbital parameters for *ORBIT.NEW*:
        !  -----------------------------------------------

        !  calculate time of perihelion in year:
        !  MONIN, A.S.: An Introduction to the Theory of Climate (Reidel Publ.)
        !   - PER: perihelion from aut. equinox in degrees (BERGER)
        !   - ECC: eccentricity of Earth's orbit (BERGER)

        ang = PER - 180.0_dp
        if (ang .lt. 0.0_dp) ang = ang+360.0_dp
        zan   = ang / 180.0_dp*pi                !   =zavexpe
        tgnu  = tan(zan/2.0_dp)
        sqecc = dsqrt((1.0_dp-ECC)/(1.0_dp+ECC))
        e     = 2.0_dp*atan(sqecc*tgnu)          !   Eccentri!   Anomaly

        !   time angle in radians of perihelion from vernal equinox

        tian = e - ECC*sin(e)
        days = tian*180.0_dp/pi                  !   360 day year only
        if (days .lt. 0.0_dp) days = days+360.0_dp !   days from ver.eq. to perh.

        !   time in days from begin of year: vernal equinox fixed at 3/21, 12 GMT
        !    = 80.5 days in 360 day year

        TPER = days+80.5_dp
        if (TPER .gt. 360.0_dp) TPER = TPER-360.0_dp

        epc = ECC*100.0_dp
      
        return

    end subroutine BERGOR 

    subroutine ORBIT_PAR(T,PER,ECC,XOB)
    ! INPUT : T - Time in years (negative for the past, reference year 1950 a. D.)                
    ! OUTPUT: Earth's orbital parameters
    !           PER  - Longitude of Perihelion (measured from vernal equinox,
    !                  relative to observation from the earth, i.e. angle from
    !                  autumnal equinox to perihelion)
    !           EC!    - Eccentricity
    !           XOB  - Obliquity
    !           T    - Attention! T (input) is multiplied by 1000.
    !
    ! Global variables: ECCM(6001),PERM(6001),XOBM(6001)
    !                   ECCP,XOBP,PERP 

        real (dp) :: T, PER, ECC, XOB
        integer :: N1, N2  
        real (dp) :: TIME, ef, PERM1, PERM2 

        N1=5001+T/1000.0_dp 
        N2=MIN0(N1+1,6001)

        TIME=T/1000
        ef=int(TIME)-TIME
        if (ef.eq.0.) ef=1.0_dp 

        ECC=ef*ECCM(N1)+(1.0_dp-ef)*ECCM(N2)
        XOB=ef*XOBM(N1)+(1.0_dp-ef)*XOBM(N2)

        PERM1=PERM(N1)
        PERM2=PERM(N2)
        if (PERM2.lt.PERM1) PERM2=PERM2+360.0_dp
        PER=ef*PERM1+(1.0_dp-ef)*PERM2
        if (PER.gt.360.0_dp) PER=PER-360.0_dp

        ECCP=ECC
        XOBP=XOB
        PERP=PER

        return
      
      end subroutine ORBIT_PAR 

end module 
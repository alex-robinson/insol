*********************************************************************  
                      SUBROUTINE INI_SINSOL(input_dir)  
*********************************************************************  
*
*  Purpose:  Calculation of solar insolation and averaged zenite
*            angle for different times
* 
*  By: A.Ganopolski 
*  Last modification: 28.04.2011 
*                                                        CLIMBER-2.4
******************************************************************** 
      COMMON  /orbitpar/ ECCM(6001),PERM(6001),XOBM(6001)
********************************************************************
      REAL nnf
      character (len=*) input_dir
      character (len=512) filen, filep

C     ! Load Berger orbital parameters
c      open (1,file=INP/orbit_par_5000.txt')

C      do n=1,5001
C       read (1,'(I5,3F10.5)') nn ,ECCM(n), PERM(n), XOBM(n)
C      enddo 

C     ! Load Laskar2004 input parameters instead of Berger   
C     ! Note: year 0: n=5001 ; year -5000 ka: n=1; year 1000 ka: n=6001

      filen = trim(input_dir)//'/LA2004.INSOLN.5Ma.txt'
      filep = trim(input_dir)//'/LA2004.INSOLP.1Ma.txt'

      write(*,*) "Loading orbital parameters from file: "//trim(filen)
      open (1,file=trim(filen))

      do n=1,5001
        read (1,*) nnf,
     >       ECCM(5001-n+1), XOBM(5001-n+1), PERM(5001-n+1)
      enddo
      
      close (1)
      
      write(*,*) "Loading orbital parameters from file: "//trim(filep)
      open (1,file=trim(filep))

      do n=1,1001
        read (1,*) nnf,
     >       ECCM(5000+n), XOBM(5000+n), PERM(5000+n)
      enddo
      
      close (1)


C     ! Convert from radians to degrees      
      XOBM = XOBM * 180.0/3.141592653589793
      PERM = PERM * 180.0/3.141592653589793
      
      return
      end

*********************************************************************
       SUBROUTINE ORBIT(ECC,XOBCH,TPERI,ZAVEXPE,
     >            PCLOCK,PYTIME,PDISSE,PZEN1,PZEN2,PZEN3,PRAE)
*********************************************************************
c
C**** *ORBIT* - FOR SOLAR ORBITAL PARAMETERS.
C
C     J.F.GELEYN     E.C.M.W.F.     03/06/82.
C
C     CHANGED TO 18000YBP-CONDITIONS BY LAUTENSCHLAGER
C                     MPI F. METEOR. HAMBURG  4.87
C     CHANGED TO OPTIONAL ORBITAL PARAMETERS BY GRAF, MPI 1/91
C     CHANGED AGAIN TO A MORE EXACT SOLUTION BY HOFFMANN, MPI 7/93
C     REVISED BY LORENZ, UNIVERSITY OF BREMEN, 10/93
c -- last change: 94-08-03
C
C     PURPOSE.
C     --------
C
C          THIS ROUTINE COMPUTES THE SOLAR CONSTANT NORMALISED BY ITS
C     ANNUAL MEAN AND THREE ORBITAL PARAMETERS DEPENDING ON THE TIME
C     OF THE DAY AS WELL AS OF THE YEAR (BOTH IN RADIANS). TIME SCALE
C     ORIGIN IS 01/01 00 GMT OF THE DESIGNED TIME DISC. THE DIFFERENT
C     TIME DISCS ARE SCALED ON 03/21 12 GMT AS THE TIME OF THE VERNAL
C     EQUINOX. ALSO RETURNED IS A CONSTANT FOR THE EFFECT OF THE EARTH'S
C     CURVATURE ON THE COSINE OF THE SOLAR ZENITH ANGLE.
C
C     THE ORBITAL PARAMETERS CAN BE TAKEN FROM PROGRAM BERGOR, WRITTEN
C     BY LORENZ, UNIVERSITY OF BREMEN, 8/94
C
C**   INTERFACE.
C     ----------
C
C          *ORBIT* IS CALLED FROM *PHYSC* AT THE FIRST LATITUDE ROW.
C          THERE ARE SEVEN DUMMY ARGUMENTS: *PCLOCK* IS THE TIME OF THE
C     DAY
C                                           *PYTIME* IS THE TIME OF THE
C     YEAR (BOTH IN RADIANS).
C                                           *PDISSE* IS THE RATIO OF THE
C     SOLAR CONSTANT TO ITS ANNUAL MEAN
C                                           *PZEN1*, *PZEN2* AND *PZEN3*
C     ARE ZENITHAL PARAMETERS.
C                                           *PRAE* IS THE RATIO OF THE
C     HEIGHT OF THE EQUIVALENT ATMOSPHERE TO THE RADIUS OF THE EARTH.
C
C     METHOD.
C     -------
C
C          STAIGHTFORWARD. INTERMEDIATE VARIABLES ARE THE SOLAR
C     DECLINATION AND THE EQUATION OF TIME.
C
C     EXTERNALS.
C     ----------
C
C          NONE.
C
C     REFERENCE.
C     ----------
C
C          NONE.
C
C
C*    DATA STATEMENTS.
C     ---- -----------
C
C          *ZCDIS*, *ZCDEC* AND *ZCEQT* ARE ARRAYS FOR A SECOND ORDER
C     *FOURIER DEVELOPMENT FOR RESPECTIVELY: SOLAR CONSTANT, SOLAR
C     DECLINATION AND EQUATION OF TIME.
C     THEY ARE VALID FOR TODAY'S ORBITAL PARAMETERS ONLY (LORENZ, 11/93).
C
C     DIMENSION ZCDIS(5),ZCDEC(5),ZCEQT(5)
C     DATA ZCDIS/+1.000110,+0.034221,+0.001280,+0.000719,+0.000077/
C     DATA ZCDEC/+0.006918,-0.399912,+0.070257,-0.006758,+0.000907/
C     DATA ZCEQT/+0.000075,+0.001868,-0.032077,-0.014615,-0.040849/
C
C     *ZRAE* IS THE VALUE FOR *PRAE*.
C
      DATA ZRAE/+.1277E-02/
C
C     ------------------------------------------------------------------
C
C*         1.     PRELIMINARY SETTING.
C                 ----------- --------
C
  100 CONTINUE

      ZCLOCK=PCLOCK
      ZYTIME=PYTIME

      PI=3.141592654

C     TROPICAL YEAR

      TTROP=360.
C
C     DAY OF VERNAL EQUINOX: DEFINED ON 21. OF MARCH, 12.00 GMT
C     (BEGIN OF YEAR: 01. OF JANUARY, 0.00 GMT: ZYTIME=0.0, THEREAFTER
C     E. G. 2. OF JANUARY 12.00 GMT: TIME IN DAYS FROM BEGIN OF YEAR = 1.5)
C     - NOT USED IN THE EXACT FORMULATION, *ZAVEXPE* USED INSTEAD
C     TVEREX=80.5              ! TROPICAL YEAR: 360 DAYS
C     TVEREX=79.5              ! TROPICAL YEAR: 365 DAYS

C      ***  ORBITAL PARAMETERS : VALUES FOR DIFFERENT TIME DISCS
C           THESE VALUES CAN BE TAKEN FROM PROGRAM BERGOR (S. LORENZ)
C
       
C
C     ------------------------------------------------------------------
C
C*         2.     COMPUTATIONS.
C                 -------------
C
  200 CONTINUE
C
C     ZC1YT=COS(ZYTIME)
C     ZS1YT=SIN(ZYTIME)
C     ZC2YT=ZC1YT**2-ZS1YT**2
C     ZS2YT=2.*ZS1YT*ZC1YT
C
C
C     CALCULATE ECCENTRIC ANOMALY:
C
C     LESS EXACT METHOD OF H. GRAF: SIN(E)=E
C     E=(ZYTIME - 2*PI*TPERI/TTROP) / (1.-ECC)
C
C     EXACT CALCULATION WITH KEPLER'S LAW (ELLIPSE)
C     (MONIN, 1986 : AN INTRODUCTION TO THE THEORY OF CLIMATE)
C     USE FIRST DERIVATIVE (NEWTON) TO CALCULATE SOLUTION OF
C     EQUATION OF ECCENTRIC ANOMALY *E*
C
      time=zytime-2*pi*tperi/ttrop
      eold=time/(1.-ecc)
      enew=time
      eps=1.E-6
      niter=0
  250 continue
      zeps=eold-enew
      if (niter.ge.30) go to 270
      if (abs(zeps).lt.eps) go to 280
      niter=niter+1
      eold=enew
      cose=cos(enew)
      enew=(time+ecc*(sin(enew)-enew*cose))/(1.-ecc*cose)
      go to 250
  270 print*,' SUBROUTINE *ORBIT* - eccentric anomaly not found!'
      print*,' ERROR IN   *ORBIT* -- STOP'
      stop
  280 continue
      E=enew
      ZDISSE=(1./(1.-ECC*COS(E)))**2
C
C     Change in the calculation of the declination.
C     Used are not the formulas from Paltridge/Platt as in the ECHAM model
C     with fixed constants for contemporanious conditions
C     but the exact equations from Monin (s.a. - Keplers law)
C     are solved. Day of vernal equinox is fixed for a 360 day year on the
C     21. Maerz noon, with start of year at 01/01/00 GMT (s.a.).
C
      zsqecc=sqrt((1+ecc)/(1-ecc))
      ztgean=tan(E/2)
c
c     znu: true anomaly (actual angle of Earth's position from perihelion)
c     zlambda: true longitude of the Earth (actual angle from vernal equinox)
c
      znu=2.*atan(zsqecc*ztgean)
      zlambda=znu+zavexpe
      xobche=xobch/180.*pi
      zsinde=sin(xobche)*sin(zlambda)
      zdecli=asin(zsinde)
 
C
      ZZEN1=SIN(ZDECLI)
      ZZEN2=COS(ZDECLI)*COS(ZCLOCK)
      ZZEN3=COS(ZDECLI)*SIN(ZCLOCK)
C
C     ------------------------------------------------------------------
C
C*         3.     RETURN.
C                 -------
C
  300 CONTINUE
 
      PDISSE=ZDISSE
      PZEN1=ZZEN1
      PZEN2=ZZEN2
      PZEN3=ZZEN3
      PRAE=ZRAE
 
      RETURN
      END

*************************************************************************
      SUBROUTINE BERGOR(BTIME,ECC,XOB,TPER,ZAN)
*************************************************************************
c
c**** *BERGOR* - calculates earth's orbital parameters for ECHAM use
c
c     S. J. Lorenz   University of Bremen    08-03-94.
c       -- last change:  31-07-96
c
c     PURPOSE.
c     --------
c
c        This programm uses a simple algorithm from A. Berger to compute
c     earth's orbital parameters of any interesting time disc. The
c     original parameters are modified to the values needed by new
c     subroutine "ORBIT" (by S. Lorenz) of AGCM "ECHAM" of MPI/DKRZ
c     in Hamburg.
c
c
c     INTERFACE.
c     ----------
c        none.
c
c     METHOD.
c     -------
c        The program BERGOR is selfexplanatory.
c
c     REFERENCE.
c     ----------
c        Berger, A.L., 1978: Long-term variations of daily insolation and
c           Quaternary climatic changes, J. Atm. Sci., 35, 2362-2367.
c
c

c  get orbital parameters with subroutine ORBIT_PAR

c    ===============================
      call ORBIT_PAR(BTIME,PER,ECC,XOB)
c    ===============================      

c  compute new orbital parameters for *ORBIT.NEW*:
c  -----------------------------------------------

      pi=3.1415926

c  calculate time of perihelion in year:
c  MONIN, A.S.: An Introdiuction to the Theory of Climate (Reidel Publ.)
c   - per: perihelion from aut. equinox in degrees (BERGER)
c   - ecc: eccentricity of Earth's orbit (BERGER)

      ang=per-180.
      if (ang.lt.0.) ang=ang+360.
      zan=ang/180.*pi                !   =zavexpe
      tgnu=tan(zan/2.)
      sqecc=sqrt((1-ecc)/(1+ecc))
      e=2.*atan(sqecc*tgnu)          !   Eccentric Anomaly

c   time angle in radians of perihelion from vernal equinox

      tian=e-ecc*sin(e)
      days=tian*180./pi              !   360 day year only
      if (days.lt.0.) days=days+360. !   days from ver.eq. to perh.

c   time in days from begin of year: vernal equinox fixed at 3/21, 12 GMT
c    = 80.5 days in 360 day year

      tper=days+80.5
      if (tper.gt.360.) tper=tper-360.

      epc=ecc*100.
      
      
      return
      end

***********************************************************************
      SUBROUTINE ORBIT_PAR(T,PER,ECC,XOB)
***********************************************************************
c     INPUT : T - Time in years (negative for the past, reference year
c                 1950 a. D.)
c     OUTPUT: Earth's orbital parameters
c             PER  - Longitude of Perihelion (measured from vernal equinox,
c                    relative to observation from the earth, i.e. angle from
c                    autumnal equinox to perihelion)
c             ECC  - Eccentricity
c             XOB  - Obliquity
c             T    - Attention! T (input) is multiplied by 1000.
*****************************************************************************
      COMMON /orbitpar/ ECCM(6001),PERM(6001),XOBM(6001)
      common /or/ ECCP,XOBP,PERP 
*****************************************************************************

      N1=5001+T/1000. 
      N2=MIN0(N1+1,6001)
      
      TIME=T/1000
      iTIME=TIME
      ef=iTIME-TIME
      if (ef.eq.0.) ef=1.
      
      ECC=ef*ECCM(N1)+(1.-ef)*ECCM(N2)
      XOB=ef*XOBM(N1)+(1.-ef)*XOBM(N2)
      
      PERM1=PERM(N1)
      PERM2=PERM(N2)
      if (PERM2.lt.PERM1) PERM2=PERM2+360.
      PER=ef*PERM1+(1.-ef)*PERM2
      if (PER.gt.360.) PER=PER-360. 
      
      ECCP=ECC
      XOBP=XOB
      PERP=PER

      RETURN
      END


! To compile:
! gfortran -o test_insol.x ../coord/interp1D.f90 insolation.f90 test_insol.f90
! gfortran -o test_insol.x ../ncio/ncio3.f90 sinsol_orbit.f test_insol.f90

program test_insol 

    use interp1D 
    use insolation 
    use ncio 

    integer, parameter :: nlat = 91
    integer, parameter :: nday = 360  
    double precision   :: lats(nlat), days(nday), insol(nlat,nday)
    integer :: i, day
    character(len=256) :: filename 
    double precision, parameter :: s0 = 1365.d0 
    double precision, parameter :: pi = 2.d0*acos(0.d0)

    integer :: init_sinsol = 0 

    do i = 1, nlat
        lats(i) = -90.d0 + (i-1)*180.d0/dble(nlat)
    end do 

    ! Calculate insolation at these latitudes
    do day = 1, nday 
        days(day) = day 
        insol(:,day) = calc_insol_daily(day,lats,0.d0)
    end do 

!     call sinsol2d(insol,lats,0.0)

    filename = "insol_0BP.nc"

    call nc_create(filename)
    call nc_write_dim(filename,"lat",x=lats)
    call nc_write_dim(filename,"day",x=days)
    call nc_write(filename,insol,"insol",dim1="lat",dim2="day")

contains 

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Subroutine :  s i n s o l 2 d
  ! Purpose    :  Same as climber module sinsol, except it is
  !               based on a 2d grid 
  ! Author     :  Alex Robinson (24. June 2008)
  ! =Input======  
  !  S0        - solar constant, normally 1365.d0 W / m2
  !  BTIME     - astronomical time (=NYRA), present day 1950 = 0.0
  !  LATS      - 2d array of latitude values, in radians!!
  ! =Output======
  !  SOLARM2D  - daily solar radiation at the top of the atmosphere
  !  COSZM2D   - daily averaged solar zenith angle
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   subroutine sinsol2d(solarm,lats,BTIME)

!     implicit none
    
!     integer, parameter :: nl = 91, nk = 360
!     real (8) :: lats(nl)
!     real (8) :: lats0(nl), solarm(nl,nk), coszm(nl,nk)
    
!     integer :: i, j, k, h, jj1, jj0
!     real (8) :: s, cosn, cosp, fi
    
!     real (4) :: BTIME, ECC, XOBCH, TPERI, ZAVEXPE
!     real (4) :: PCLOCK, PYTIME, PDISSE, PZEN1, PZEN2, PZEN3, PRAE

!     if ( init_sinsol .eq. 0 ) then
!       call INI_SINSOL("input")
!       init_sinsol = 1
!     end if
    
!     ! Populate the lats lookup table    
!     lats0 = lats 
    
! !...1) Berger program calculates orbital parameters
! !   ===========================================   
!     call BERGOR(BTIME,ECC,XOBCH,TPERI,ZAVEXPE)
! !   =========================================== 
    
! !...2) Daily insolation is calculated by hourly integration for each day
     
!     ! Loop over interpolation matrix (lats0, solarm, coszm)
!     do j = 1, nl

!         fi=lats0(j)*pi/180.d0
!         do k = 1, nk
!           PYTIME=dble(k)*2.d0*pi/360.d0
!           solarm(j,k)=0.d0
!           coszm(j,k) =0.d0        
!           do h = 1, 24
!             PCLOCK=dble(h)*2.d0*pi/24.d0   
! !     =================================================================          
!             call ORBIT(ECC,XOBCH,TPERI,ZAVEXPE, &
!                  PCLOCK,PYTIME,PDISSE,PZEN1,PZEN2,PZEN3,PRAE)
! !     =================================================================                    
!             cosp=PZEN1*dsin(fi)+PZEN2*dcos(fi)
!             cosn=max(cosp,0.0)
      
!             s=s0*cosn*PDISSE
!             solarm(j,k) = solarm(j,k)+s
!             coszm(j,k)  = coszm(j,k)+s*cosn
!           end do
!         end do

      
! !...  Daily insolation and zenite angle
                                 

!         do k = 1, nk
!           solarm(j,k)=solarm(j,k)/24.d0
!           if (solarm(j,k) .gt. 0.) then
!             coszm(j,k)=coszm(j,k)/(solarm(j,k)*24.d0)
!           else
!             coszm(j,k)=0.d0
!           end if
!         end do
      
!     end do ! end j-loop
    
!     return
  
!   end subroutine sinsol2d 
  
end program 
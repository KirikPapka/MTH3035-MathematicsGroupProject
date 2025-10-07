      ! **************************************
      !   GRID
      ! **************************************

      ! Vertical resolution (m) 
!      ! redundant for log linear levels  
      LOGICAL, PARAMETER :: stretch_grid=.true.  
      !REAL, PARAMETER :: dz_grid=10.0

      ! Surface layer starting level for bl height diagnosis
      INTEGER, PARAMETER :: k_sl=10
      ! Model top (m)
      REAL, PARAMETER :: zztop=4000.0
      ! number of fine grid points.
      INTEGER, PARAMETER :: n_fine=200001
      ! log grid constants (m)
      REAL, PARAMETER :: a_grid=20.0, b_grid=2.0

      ! initial inversion top
      REAL, PARAMETER :: z_i=30.0
      ! Overlying inversion strength K/m
      REAL, PARAMETER :: inv=0.0384


      ! Number of levels
      INTEGER, PARAMETER :: kkp=401

      ! Time step (s)
      REAL :: dt

      ! Frequency of output (s)
      REAL :: freq_out       
      
      ! Length of run (s)
      REAL :: run_length

      NAMELIST/run/ run_length, freq_out, dt


      ! FORCING AND OTHER CONSTANTS
      ! Surface boundary condition option surf_opt
      ! =0 prescribed sst
      ! =1 prescribed flux
      INTEGER, PARAMETER :: surf_opt=1
      ! surface heating rate K/h  for surf_opt=0
      REAL, PARAMETER :: heating=1.5

      ! Maximum surface heat flux
      REAL, PARAMETER :: h_max=173.0
      ! SBL heat flux
      REAL, PARAMETER :: h_sbl=-15.0
      ! Timescale for heat flux (s)
      REAL, PARAMETER :: tau_h=20900.0
      ! Time of morning transition(s)
      REAL, PARAMETER :: t_trans=-3600.0 
        
      ! Geostrophic wind
      REAL :: U_g, V_g
      NAMELIST/geowind/ U_g, V_g
      ! Reference potential temperature (K)
      REAL, PARAMETER :: th_ref = 280.0 
      
      ! Amplitude of exponential initial pot. temp profile
      REAL, PARAMETER :: th_i=4.0
      ! scale height for exponential initial pot. temp profile
!      REAL, PARAMETER :: h_3=120.0
      REAL, PARAMETER :: h_3=300.0
      ! scale height for exponential initial momentum profile
      REAL, PARAMETER :: h_m=120.0       
 
      ! Parameters for a jet profile.
      REAL, PARAMETER :: CC=0.4
      REAL, PARAMETER :: DD=0.57
      ! pi *h_m=hh (for Ekman spiral)
      REAL, PARAMETER :: hh=376.99
      REAL, PARAMETER :: h_2=450.0
            
      REAL,  PARAMETER ::h_beare=300.0
 
      ! Initial surface potential temperature (K)
      REAL, PARAMETER :: th0_init=300.0

      ! Coriolis parameter (s^-1)
      REAL, PARAMETER :: f_0=1.0e-4

      ! Gravitational acceleration (m/s^2)
      REAL, PARAMETER :: g=9.81

      ! Constant neutral mixing length in Derbyshire 99 model (m)
      REAL, PARAMETER :: lambda_const=50.0

      ! Soil equivalent depth (m)
      REAL, PARAMETER :: SED=40.0

      ! Infrared forcing; cooling is positive (Wm^-2)
      REAL, PARAMETER :: IR=20.0

      ! Surface density (kgm^-3)
      REAL, PARAMETER :: rho_0=1.0

      ! Specific heat capacity of air at constant pressure (JK^-1kg^-1)
      REAL, PARAMETER :: c_p= 1005.0

      ! Karman constant
      REAL, PARAMETER :: kar=0.4

      ! Momentum (and thermal) rhoughness length
      REAL, PARAMETER :: z_0=0.1
      
!      ! 
!      REAL, PARAMETER ::

      ! Implicit weight (0.5 implies Crank-Nicolson)
      REAL :: gamma=0.5

      ! Critical Richardson number for bl height
      ! diagnosis
      REAL, PARAMETER :: ri_cr=0.5
      ! Zilitinkevich '72 constant for bl depth
      REAL, PARAMETER :: c_zil=0.7
      ! factor b in Holtslag and Boville '93
      REAL, PARAMETER :: b_hb=8.5
      ! factor a in  Holtslag and Boville '93
      ! factor in non-gradient term
      REAL, PARAMETER :: a_hb=7.2

      ! Non-gradient momentum term
      ! following Brown and Grant '97
      REAL, PARAMETER :: a_brgr=5.4


      ! The factor in defining the momentum CBL velocity
      ! scale in Holtslag and Boville '93
      REAL, PARAMETER :: c_1_hb=0.6  

      ! Option for wind initialisation
      ! 0 = geostrophic initialisation
      ! 1 = exponential in u, 0 in v
      ! 2 = Beare '06 jet
      ! 3 = Ekman spiral
      
      ! Moist Dynamics
      logical, parameter :: moist_dyn=.true.
      real, parameter :: psfc = 101000. !Pa
      real, parameter :: p0 = 1.0e5 ! Pa
      real, parameter :: Rd = 287. ! J/kg
      real, parameter :: Rv = 461. ! J/kg
      real, parameter :: Lv = 2.5e6 ! Latent Heat of vaporisation
      
 
      INTEGER :: wind_opt
      NAMELIST/wind/wind_opt

      

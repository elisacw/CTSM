module SoilBiogeochemDecompCascadeMIMICSplusMod

  !-----------------------------------------------------------------------
  ! DESCRIPTION:

  ! MIMICS+ soil decomposition model based on Aas et al (2023)
  
  ! This module is based on: SoilBiogeochemDecompCascadeMIMICSMod 
  ! Changes
  ! - Creating the option to use mimicsplus as a soil decomposition method
  ! - created folders for century, mimicsplus and mimicsplusFatesCold in: /ctsm/cime_config/testdefs/testmods_dirs/clm/

  ! Additions inside SoilBiogeochemDecompCascadeMIMICSplusMod
  ! - EcM and AM pools
  ! - modified parameter equations for soil microbes
  ! - added parameters / parameter equations for mycorrhiza
  ! - relpacing w_d_o_scalar with moisture function 
  ! - TO DO add mycorrhizal modifier

  ! Parameters
  ! - new parameter file including updated MIMICS parameters & new parameters, saved under /cluster/home/elisacw/parameter_mimicsplus
  ! - added read in method for new created parameters in SoilBiogeochemDecompCascadeMIMICSplusMod


  ! 
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_const_mod                      , only : SHR_CONST_TKFRZ
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : nlevdecomp, ndecomp_pools_max
  use clm_varpar                         , only : i_str_lit, i_met_lit, i_cop_mic, i_oli_mic, i_cwd
  use clm_varpar                         , only : i_phys_som, i_chem_som, i_avl_som, i_ecm_myc, i_am_myc
  use clm_varpar                         , only : i_litr_min, i_litr_max, i_cwdl2
  use clm_varctl                         , only : iulog, spinup_state, anoxia, use_lch4, use_fates
  use clm_varcon                         , only : zsoi
  use decompMod                          , only : bounds_type, subgrid_level_patch
  use spmdMod                            , only : masterproc
  use abortutils                         , only : endrun
  use CNSharedParamsMod                  , only : CNParamsShareInst, nlev_soildecomp_standard 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilStateType                      , only : soilstate_type
  use TemperatureType                    , only : temperature_type 
  use CNVegCarbonFluxType                , only : cnveg_carbonflux_type
  use ch4Mod                             , only : ch4_type
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  use CLMFatesInterfaceMod               , only : hlm_fates_interface_type
  use WaterStateBulkType                 , only : waterstatebulk_type
   
  implicit none
  private
  !

  !                           | MIMICS+  | MIMICS 
  ! 
  ! Metabolic litter          | LITm     | l1        | i_met_lit |
  ! Structural litter         | LITs     | l2        | i_str_lit |
  ! Microbes                  | SAPb     | m1 (MICr) | i_cop_mic |
  ! Microbes                  | SAPf     | m2 (MICk) | i_oli_mic |
  ! Physically protected SOM  | SOMp     | s3        | i_phys_som |
  ! Available SOM             | SOMa     | s1        | i_avl_som |
  ! Chemically protected SOM  | SOMc     | s2        | i_chem_som |
  ! Ectomycorrhiza            | EcM      | - 
  ! Arbuscular Mycorrhiza     | AM       | - 


  !PUBLIC MEMBER FUNCTIONS:
  public :: readParams                          ! Read in parameters from params file
  public :: init_decompcascade_mimicsplus       ! Initialization
  public :: decomp_rates_mimicsplus             ! Figure out decomposition rates
  !public :: decomp_rates_after_FUN
  public :: calc_myc_mining_rates
  public :: calc_myc_mortality
  public :: calc_myc_roi
  public :: cost_FUN
  public :: fun_fluxes_myc_update1
  private :: r_moist                            ! calculates moisture modifier according to CORPSE
  !
  ! !PUBLIC DATA MEMBERS 
  !
  ! !PRIVATE DATA MEMBERS 
  ! next four 2d vars are dimensioned (columns,nlevdecomp)
      
  real(r8), private, allocatable :: desorp(:,:)
  real(r8), private, allocatable :: fphys_m1(:,:)
  real(r8), private, allocatable :: fphys_m2(:,:)
  real(r8), private, allocatable :: p_scalar(:,:) 
  integer, private :: i_l1m1        ! indices of transitions, eg l1m1: metabolic litter -> microbial bacteria
  integer, private :: i_l1m2                                         ! metabolic litter -> microbial fungi
  integer, private :: i_l2m1                                         ! structural litter -> microbial bacteria     
  integer, private :: i_l2m2                                         ! structural litter -> microbial fungi     
  integer, private :: i_s1m1                                         ! available SOM -> microbial bacteria       
  integer, private :: i_s1m2                                         ! available SOM -> microbial fungi      
  integer, private :: i_m1s1                                         ! microbial bacteria -> available SOM     
  integer, private :: i_m1s2                                         ! microbial bacteria -> chemically protected SOM 
  integer, private :: i_m1s3                                         ! microbial bacteria -> physically protected SOM       
  integer, private :: i_m2s1                                         ! microbial fungi -> available SOM       
  integer, private :: i_m2s2                                         ! microbial fungi -> chemically protected SOM 
  integer, private :: i_m2s3                                         ! microbial fungi -> physically protected SOM              
  integer, private :: i_s2s1                                         ! chemically protected SOM -> available SOM    
  integer, private :: i_s3s1                                         ! physically protected SOM -> available SOM     
  integer, private :: i_myc1s1                                       ! ectomycorrhizal fungi -> avaliable SOM  
  integer, private :: i_myc1s2                                       ! ectomycorrhizal fungi -> chemically protected SOM 
  integer, private :: i_myc1s3                                       ! ectomycorrhizal fungi -> physically protected SOM     
  integer, private :: i_myc2s1                                       ! abruscular fungi -> avaliable SOM  
  integer, private :: i_myc2s2                                       ! abruscular fungi -> chemically protected SOM
  integer, private :: i_myc2s3                                       ! abruscular fungi -> physically protected SOM
  integer, private :: i_myc1s1_e ! enzyme production by ectomycorrhizal fungi | ectomycorrhizal fungi -> avaliable SOM  
  integer, private :: i_s2myc1_m ! nitrogen mining | chemically protected SOM -> ectomycorrhizal fungi
  integer, private :: i_s3myc1_m ! nitrogen mining | physically protected SOM -> ectomycorrhizal fungi

        
  ! MICHAELIS MENTEN
  real(r8), private :: rf_l1m1  ! respiration fractions by transition
  real(r8), private :: rf_l1m2
  real(r8), private :: rf_l2m1
  real(r8), private :: rf_l2m2
  real(r8), private :: rf_s1m1
  real(r8), private :: rf_s1m2
  real(r8), private :: vint_l1_m1  ! regression intercepts by transition
  real(r8), private :: vint_l2_m1
  real(r8), private :: vint_s1_m1
  real(r8), private :: vint_l1_m2
  real(r8), private :: vint_l2_m2
  real(r8), private :: vint_s1_m2
  real(r8), private :: kint_l1_m1  ! regression intercepts by transition
  real(r8), private :: kint_l2_m1
  real(r8), private :: kint_s1_m1
  real(r8), private :: kint_l1_m2
  real(r8), private :: kint_l2_m2
  real(r8), private :: kint_s1_m2
  real(r8), private :: vmod_l1_m1  ! vmod = vmod * av from Wieder et al 2015
  real(r8), private :: vmod_l2_m1
  real(r8), private :: vmod_s1_m1
  real(r8), private :: vmod_l1_m2
  real(r8), private :: vmod_l2_m2
  real(r8), private :: vmod_s1_m2
  real(r8), private :: kmod_l1_m1  ! kmod = ak / kmod from Wieder et al 2015
  real(r8), private :: kmod_l2_m1
  real(r8), private :: kmod_s1_m1
  real(r8), private :: kmod_l1_m2
  real(r8), private :: kmod_l2_m2
  real(r8), private :: kmod_s1_m2
  real(r8), private :: vslope_l1_m1  ! regression coefficients by transition
  real(r8), private :: vslope_l2_m1
  real(r8), private :: vslope_s1_m1
  real(r8), private :: vslope_l1_m2
  real(r8), private :: vslope_l2_m2
  real(r8), private :: vslope_s1_m2
  real(r8), private :: kslope_l1_m1  ! regression coefficients by transition
  real(r8), private :: kslope_l2_m1
  real(r8), private :: kslope_s1_m1
  real(r8), private :: kslope_l1_m2
  real(r8), private :: kslope_l2_m2
  real(r8), private :: kslope_s1_m2

  ! Parameters from parameter file
  type, private :: params_type
      real(r8) :: mimicsplus_nue_into_mic  ! microbial N use efficiency for N fluxes
      real(r8) :: mimicsplus_desorpQ10
      real(r8) :: mimicsplus_densdep       ! exponent controling the density dependence of microbial turnover
      real(r8) :: mimicsplus_tau_mod_factor  ! (1 / tauModDenom) from testbed code
      real(r8) :: mimicsplus_tau_mod_min
      real(r8) :: mimicsplus_tau_mod_max
      real(r8) :: mimicsplus_ko_r  ! increase in half-saturation constant for oxidation (SOMc -> SOMa) microbial bacteria
      real(r8) :: mimicsplus_ko_k  ! increase in half-saturation constant for oxidation (SOMc -> SOMa) microbial fungi
      real(r8) :: mimicsplus_cn_r  ! C:N of microbial bacteria
      real(r8) :: mimicsplus_cn_k  ! C:N of microbial fungi
      real(r8) :: mimicsplus_cn_mod_num  ! adjusts microbial CN based on fmet
      real(r8) :: mimicsplus_t_soi_ref  ! reference soil temperature (degC)
      real(r8) :: mimicsplus_initial_Cstocks_depth  ! Soil depth for initial C stocks for a cold-start (m)
      real(r8), allocatable :: mimicsplus_initial_Cstocks(:)  ! Initial C stocks for a cold-start (gC/m3)
      ! The next few vectors are dimensioned by the number of decomposition
      ! transitions that make use of the corresponding parameters, currently
      ! six. The transitions are represented in this order:
      ! l1m1 l2m1 s1m1 l1m2 l2m2 s1m2
      real(r8), allocatable :: mimicsplus_mge(:)  ! Microbial growth efficiency (mg/mg)
      real(r8), allocatable :: mimicsplus_vmod(:)  ! vmod = vmod * av from Wieder et al 2015
      real(r8), allocatable :: mimicsplus_vint(:)  ! regression intercepts (5.47 ln(mg Cs (mg MIC)-1 h-1) )
      real(r8), allocatable :: mimicsplus_vslope(:)  ! regression coeffs (ln(mg Cs (mg MIC)-1 h-1) ¡C-1)
      real(r8), allocatable :: mimicsplus_kmod(:)  ! kmod = ak / kmod from Wieder et al 2015
      real(r8), allocatable :: mimicsplus_kint(:)  ! regression intercepts
      real(r8), allocatable :: mimicsplus_kslope(:)  ! regression coeffs
      ! The next few vectors are dimensioned by the number of parameters with
      ! the same name (eg 2 for mimicsplus_tau_r_p1, mimicsplus_tau_r_p2) used in the
      ! respective formula. In the formulas, we use scalar copies of each
      ! parameter with suffixes _p1, p2, p3, ... to distinguish among them.
      ! See allocate statements below for the size of each of the following
      ! vectors.
      real(r8), allocatable :: mimicsplus_fmet(:)
      real(r8), allocatable :: mimicsplus_p_scalar(:)
      real(r8), allocatable :: mimicsplus_fphys_r(:)
      real(r8), allocatable :: mimicsplus_fphys_k(:)
      real(r8), allocatable :: mimicsplus_fchem_r(:)
      real(r8), allocatable :: mimicsplus_fchem_k(:)
      real(r8), allocatable :: mimicsplus_desorp(:)
      real(r8), allocatable :: mimicsplus_tau_r(:)
      real(r8), allocatable :: mimicsplus_tau_k(:)
      
      !Mycorrhoiza parameter:
      real(r8) :: mimicsplus_k_myc_som   ! 'Turnover rate for Mycorrhiza', 'units': 'h⁻1'
      real(r8) :: mimicsplus_k_mo        ! 'Mycorrhizal decay rate', 'units': 'm²gC⁻¹hr⁻¹'
      real(r8) :: mimicsplus_vmax_myc    ! 'Max. mycorrhizal uptake up inorg. N', 'units': 'g x g⁻¹h⁻¹'
      real(r8) :: mimicsplus_k_m_emyc    ! 'Half saturation constant of ectomycorrhizal uptake of inorg N', 'units': 'gNm⁻²'
      real(r8) :: mimicsplus_mge_ecm     ! 'Microbial growth efficiency (mg/mg) for ectomycorrhiza ', 'units': 'unitless'
      real(r8) :: mimicsplus_mge_am      ! 'Growth efficiency of arbuscular mycorrhiza ', 'units': 'unitless'
      !real(r8) :: mimicsplus_r_myc       ! 'Mycorrhizal modifier', 'units': 'unitless'
      real(r8) :: mimicsplus_fphys_ecm   ! 'Fraction ectomycorrhizal necromass into physically protected soil organic matter pool', 'units': 'unitless'
      real(r8) :: mimicsplus_fchem_ecm   ! 'Fraction ectomycorrhizal necromass into chemically protected soil organic matter pool', 'units': 'unitless'
      real(r8) :: mimicsplus_tau_ecm     ! 'Fraction ectomycorrhizal necromass into avaliable soil organic matter pool', 'units': 'unitless'
      real(r8) :: mimicsplus_fphys_am    ! 'Fraction arbuscular mycorrhizal necromass into physically protected soil organic matter pool', 'units': 'unitless'
      real(r8) :: mimicsplus_fchem_am    ! 'Fraction arbuscular mycorrhizal necromass into chemically protected soil organic matter pool', 'units': 'unitless'
      real(r8) :: mimicsplus_tau_am      ! 'Fraction arbuscular mycorrhizal necromass into avaliable soil organic matter pool', 'units': 'unitless'
      real(r8) :: mimicsplus_cn_myc      ! 'Optimal CN ratio for mycorrhizal fungi', 'units': 'unitless'
  end type params_type
  !
  type(params_type), private :: params_inst

  character(len=*), parameter, private :: sourcefile = &
        __FILE__

  !-----------------------------------------------------------------------

  contains

  ! Parameters are read in from parameter file:
  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
      !
      ! !DESCRIPTION:
      !
      ! !USES:
      use ncdio_pio    , only: file_desc_t,ncd_io
      !
      ! !ARGUMENTS:
      type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
      !
      ! !LOCAL VARIABLES:
      character(len=32)  :: subname = 'readMimicsplusParams'
      character(len=100) :: errCode = 'Error reading MIMICSplus params '
      logical            :: readv   ! has variable been read in or not
      real(r8)           :: tempr   ! temporary to read in constant
      character(len=100) :: tString ! temp. var for reading
      !-----------------------------------------------------------------------

      ! Read off of netcdf file
      tString='mimicsplus_initial_Cstocks_depth'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_initial_Cstocks_depth=tempr

      allocate(params_inst%mimicsplus_initial_Cstocks(ndecomp_pools_max))
      tString='mimicsplus_initial_Cstocks'
      call ncd_io(trim(tString), params_inst%mimicsplus_initial_Cstocks(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_mge(ndecomp_pools_max))
      tString='mimicsplus_mge'
      call ncd_io(trim(tString), params_inst%mimicsplus_mge(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_vmod(ndecomp_pools_max))
      tString='mimicsplus_vmod'
      call ncd_io(trim(tString), params_inst%mimicsplus_vmod(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_vslope(ndecomp_pools_max))
      tString='mimicsplus_vslope'
      call ncd_io(trim(tString), params_inst%mimicsplus_vslope(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_vint(ndecomp_pools_max))
      tString='mimicsplus_vint'
      call ncd_io(trim(tString), params_inst%mimicsplus_vint(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_kmod(ndecomp_pools_max))
      tString='mimicsplus_kmod'
      call ncd_io(trim(tString), params_inst%mimicsplus_kmod(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_kslope(ndecomp_pools_max))
      tString='mimicsplus_kslope'
      call ncd_io(trim(tString), params_inst%mimicsplus_kslope(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_kint(ndecomp_pools_max))
      tString='mimicsplus_kint'
      call ncd_io(trim(tString), params_inst%mimicsplus_kint(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_p_scalar(2))
      tString='mimicsplus_p_scalar'
      call ncd_io(trim(tString), params_inst%mimicsplus_p_scalar(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_desorp(2))
      tString='mimicsplus_desorp'
      call ncd_io(trim(tString), params_inst%mimicsplus_desorp(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_fphys_r(2))
      tString='mimicsplus_fphys_r'
      call ncd_io(trim(tString), params_inst%mimicsplus_fphys_r(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_fphys_k(2))
      tString='mimicsplus_fphys_k'
      call ncd_io(trim(tString), params_inst%mimicsplus_fphys_k(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_fmet(4))
      tString='mimicsplus_fmet'
      call ncd_io(trim(tString), params_inst%mimicsplus_fmet(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_fchem_r(4))
      tString='mimicsplus_fchem_r'
      call ncd_io(trim(tString), params_inst%mimicsplus_fchem_r(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_fchem_k(4))
      tString='mimicsplus_fchem_k'
      call ncd_io(trim(tString), params_inst%mimicsplus_fchem_k(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_tau_r(2))
      tString='mimicsplus_tau_r'
      call ncd_io(trim(tString), params_inst%mimicsplus_tau_r(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      allocate(params_inst%mimicsplus_tau_k(2))
      tString='mimicsplus_tau_k'
      call ncd_io(trim(tString), params_inst%mimicsplus_tau_k(:), 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))

      tString='mimicsplus_nue_into_mic'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_nue_into_mic = tempr

      tString='mimicsplus_tau_mod_factor'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_tau_mod_factor = tempr

      tString='mimicsplus_tau_mod_min'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_tau_mod_min = tempr

      tString='mimicsplus_tau_mod_max'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_tau_mod_max = tempr

      tString='mimicsplus_ko_r'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_ko_r = tempr

      tString='mimicsplus_ko_k'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_ko_k = tempr

      tString='mimicsplus_densdep'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_densdep = tempr

      tString='mimicsplus_desorpQ10'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_desorpQ10 = tempr

      tString='mimicsplus_t_soi_ref'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_t_soi_ref = tempr

      tString='mimicsplus_cn_mod_num'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_cn_mod_num = tempr

      tString='mimicsplus_cn_r'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_cn_r = tempr

      tString='mimicsplus_cn_k'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_cn_k = tempr

      tString='mimicsplus_k_myc_som'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_k_myc_som = tempr

      tString='mimicsplus_k_mo'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_k_mo = tempr

      tString='mimicsplus_vmax_myc'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_vmax_myc = tempr

      tString='mimicsplus_k_m_emyc'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_k_m_emyc = tempr

      tString='mimicsplus_mge_ecm'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_mge_ecm = tempr

      tString='mimicsplus_mge_am'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_mge_am = tempr

      !tString='mimicsplus_r_myc'
      !call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      !if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      !params_inst%mimicsplus_r_myc = tempr

      tString='mimicsplus_fphys_ecm'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_fphys_ecm = tempr

      tString='mimicsplus_fchem_ecm'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_fchem_ecm = tempr

      tString='mimicsplus_tau_ecm'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_tau_ecm = tempr

      tString='mimicsplus_fphys_am'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_fphys_am = tempr

      tString='mimicsplus_fchem_am'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_fchem_am = tempr

      tString='mimicsplus_tau_am'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_tau_am = tempr

      tString='mimicsplus_cn_myc'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_cn_myc = tempr

  end subroutine readParams  

  !-----------------------------------------------------------------------
  subroutine init_decompcascade_mimicsplus(bounds, soilbiogeochem_state_inst, soilstate_inst)
      !
      ! !DESCRIPTION:
      ! initialize rate constants and decomposition pathways following the
      ! decomposition cascade of the MIMICS+ model.
      !
      ! !USES:
      use clm_varcon, only: pct_to_frac
      !
      ! !ARGUMENTS:
      type(bounds_type)               , intent(in)    :: bounds  
      type(soilbiogeochem_state_type) , intent(inout) :: soilbiogeochem_state_inst
      type(soilstate_type)            , intent(in)    :: soilstate_inst
      !

      ! !LOCAL VARIABLES
      !-- properties of each decomposing pool

      !copied from above:
      ! The next few vectors are dimensioned by the number of parameters with
      ! the same name (eg 2 for mimicsplus_tau_r_p1, mimicsplus_tau_r_p2) used in the
      ! respective formula. In the formulas, we use scalar copies of each
      ! parameter with suffixes _p1, p2, p3, ... to distinguish among them.
      ! See allocate statements below for the size of each of the following
      ! vectors.
      real(r8) :: mimicsplus_nue_into_mic
      real(r8) :: mimicsplus_p_scalar_p1     !physical protection scalar
      real(r8) :: mimicsplus_p_scalar_p2
      real(r8) :: mimicsplus_fphys_r_p1      !microbial bacteria -> SOMp
      real(r8) :: mimicsplus_fphys_r_p2
      real(r8) :: mimicsplus_fphys_k_p1      !microbial fungi -> SOMp
      real(r8) :: mimicsplus_fphys_k_p2
      real(r8) :: mimicsplus_desorp_p1       !Desorption rate !ECW parameters and equation needs to be edited
      real(r8) :: mimicsplus_desorp_p2

      real(r8):: speedup_fac                  ! acceleration factor, higher when vertsoilc = .true.

      real(r8) :: clay_frac  ! local copy of cellclay converted from % (fraction)
      integer  :: c, j  ! indices
      !-----------------------------------------------------------------------

      
      associate(                                                                                     &
            nue_decomp_cascade             => soilbiogeochem_state_inst%nue_decomp_cascade_col      , & ! Output: [real(r8)          (:)     ]  N use efficiency for a given transition (gN going into microbe / gN decomposed)

            cellclay                       => soilstate_inst%cellclay_col                           , & ! Input:  [real(r8)          (:,:)   ]  column 3D clay (%)
            
            cascade_donor_pool             => decomp_cascade_con%cascade_donor_pool                 , & ! Output: [integer           (:)     ]  which pool is C taken from for a given decomposition step 
            cascade_receiver_pool          => decomp_cascade_con%cascade_receiver_pool              , & ! Output: [integer           (:)     ]  which pool is C added to for a given decomposition step   
            floating_cn_ratio_decomp_pools => decomp_cascade_con%floating_cn_ratio_decomp_pools     , & ! Output: [logical           (:)     ]  TRUE => pool has fixed C:N ratio                          
            is_microbe                     => decomp_cascade_con%is_microbe                         , & ! Output: [logical           (:)     ]  TRUE => pool is microbial
            is_litter                      => decomp_cascade_con%is_litter                          , & ! Output: [logical           (:)     ]  TRUE => pool is a litter pool                             
            is_soil                        => decomp_cascade_con%is_soil                            , & ! Output: [logical           (:)     ]  TRUE => pool is a soil pool                               
            is_cwd                         => decomp_cascade_con%is_cwd                             , & ! Output: [logical           (:)     ]  TRUE => pool is a cwd pool                                
            initial_cn_ratio               => decomp_cascade_con%initial_cn_ratio                   , & ! Output: [real(r8)          (:)     ]  c:n ratio for initialization of pools                    
            initial_stock                  => decomp_cascade_con%initial_stock                      , & ! Output: [real(r8)          (:)     ]  initial concentration for seeding at spinup              
            initial_stock_soildepth        => decomp_cascade_con%initial_stock_soildepth            , & ! Output: [real(r8)          (:)     ]  soil depth for initial concentration for seeding at spinup              
            is_metabolic                   => decomp_cascade_con%is_metabolic                       , & ! Output: [logical           (:)     ]  TRUE => pool is metabolic material                        
            is_cellulose                   => decomp_cascade_con%is_cellulose                       , & ! Output: [logical           (:)     ]  TRUE => pool is cellulose                                 
            is_lignin                      => decomp_cascade_con%is_lignin                          , & ! Output: [logical           (:)     ]  TRUE => pool is lignin                                    
            spinup_factor                  => decomp_cascade_con%spinup_factor                      , & ! Output: [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
            is_mycorrhiza                  => decomp_cascade_con%is_mycorrhiza                        & ! Output: [logical           (:)     ]  TRUE => pool is mycorrhiza
            )

        allocate(desorp(bounds%begc:bounds%endc,1:nlevdecomp))
        allocate(fphys_m1(bounds%begc:bounds%endc,1:nlevdecomp))
        allocate(fphys_m2(bounds%begc:bounds%endc,1:nlevdecomp))
        allocate(p_scalar(bounds%begc:bounds%endc,1:nlevdecomp))


        !Parameter file has (partly) several values for one parameter, 
        !here is one value assignt to one transition between the pools
        !xxx_p1 = mimicsplus_xxx(1)
        !this xxx_p1 can than later be used in equations

        !------- time-constant coefficients ---------- !
        mimicsplus_nue_into_mic = params_inst%mimicsplus_nue_into_mic
        mimicsplus_p_scalar_p1 = params_inst%mimicsplus_p_scalar(1)
        mimicsplus_p_scalar_p2 = params_inst%mimicsplus_p_scalar(2)
        mimicsplus_fphys_r_p1 = params_inst%mimicsplus_fphys_r(1)
        mimicsplus_fphys_r_p2 = params_inst%mimicsplus_fphys_r(2)
        mimicsplus_fphys_k_p1 = params_inst%mimicsplus_fphys_k(1)
        mimicsplus_fphys_k_p2 = params_inst%mimicsplus_fphys_k(2)
        mimicsplus_desorp_p1 = params_inst%mimicsplus_desorp(1)
        mimicsplus_desorp_p2 = params_inst%mimicsplus_desorp(2)

        ! set respiration fractions for fluxes between microbial compartments
        rf_l1m1 = 1.0_r8 - params_inst%mimicsplus_mge(1)
        rf_l2m1 = 1.0_r8 - params_inst%mimicsplus_mge(2)
        rf_s1m1 = 1.0_r8 - params_inst%mimicsplus_mge(3)
        rf_l1m2 = 1.0_r8 - params_inst%mimicsplus_mge(4)
        rf_l2m2 = 1.0_r8 - params_inst%mimicsplus_mge(5)
        rf_s1m2 = 1.0_r8 - params_inst%mimicsplus_mge(6)

        ! vmod = "old" vmod * av  AND  kmod = ak / "old" kmod
        ! Table B1 Wieder et al. (2015) and MIMICS params file give diff
        ! ak and av values. I used the params file values.
        vmod_l1_m1 = params_inst%mimicsplus_vmod(1)
        vmod_l2_m1 = params_inst%mimicsplus_vmod(2)
        vmod_s1_m1 = params_inst%mimicsplus_vmod(3)
        vmod_l1_m2 = params_inst%mimicsplus_vmod(4)
        vmod_l2_m2 = params_inst%mimicsplus_vmod(5)
        vmod_s1_m2 = params_inst%mimicsplus_vmod(6)
        kmod_l1_m1 = params_inst%mimicsplus_kmod(1)
        kmod_l2_m1 = params_inst%mimicsplus_kmod(2)
        kmod_s1_m1 = params_inst%mimicsplus_kmod(3)
        kmod_l1_m2 = params_inst%mimicsplus_kmod(4)
        kmod_l2_m2 = params_inst%mimicsplus_kmod(5)
        kmod_s1_m2 = params_inst%mimicsplus_kmod(6)
        vslope_l1_m1 = params_inst%mimicsplus_vslope(1)
        vslope_l2_m1 = params_inst%mimicsplus_vslope(2)
        vslope_s1_m1 = params_inst%mimicsplus_vslope(3)
        vslope_l1_m2 = params_inst%mimicsplus_vslope(4)
        vslope_l2_m2 = params_inst%mimicsplus_vslope(5)
        vslope_s1_m2 = params_inst%mimicsplus_vslope(6)
        kslope_l1_m1 = params_inst%mimicsplus_kslope(1)
        kslope_l2_m1 = params_inst%mimicsplus_kslope(2)
        kslope_s1_m1 = params_inst%mimicsplus_kslope(3)
        kslope_l1_m2 = params_inst%mimicsplus_kslope(4)
        kslope_l2_m2 = params_inst%mimicsplus_kslope(5)
        kslope_s1_m2 = params_inst%mimicsplus_kslope(6)
        vint_l1_m1 = params_inst%mimicsplus_vint(1)
        vint_l2_m1 = params_inst%mimicsplus_vint(2)
        vint_s1_m1 = params_inst%mimicsplus_vint(3)
        vint_l1_m2 = params_inst%mimicsplus_vint(4)
        vint_l2_m2 = params_inst%mimicsplus_vint(5)
        vint_s1_m2 = params_inst%mimicsplus_vint(6)
        kint_l1_m1 = params_inst%mimicsplus_kint(1)
        kint_l2_m1 = params_inst%mimicsplus_kint(2)
        kint_s1_m1 = params_inst%mimicsplus_kint(3)
        kint_l1_m2 = params_inst%mimicsplus_kint(4)
        kint_l2_m2 = params_inst%mimicsplus_kint(5)
        kint_s1_m2 = params_inst%mimicsplus_kint(6)

        ! some of these are dependent on the soil texture properties
        ! One-time initializations here.
        ! Time-dep params in subr. decomp_rates_mimicsplus.

        !SOIL INETIALIZATION 
        do c = bounds%begc, bounds%endc
            do j = 1, nlevdecomp
              ! The parameter values currently in the params files always lead to
              ! positive values for the expressions below, so we do not
              ! need to use the max function to limit these expressions.
              ! We apply the min function on cellclay because we are looping over
              ! some non-soil columns here that contain cellclay = 1e36.
              clay_frac = pct_to_frac * &
                          dmin1(100.0_r8, cellclay(c,j))  ! conv. % to fraction
              desorp(c,j) = mimicsplus_desorp_p1 * dexp(mimicsplus_desorp_p2 * clay_frac)  
              fphys_m1(c,j) = dmin1(1.0_r8, mimicsplus_fphys_r_p1 * &
                                            dexp(mimicsplus_fphys_r_p2 * clay_frac))       !parameter equation for microbial necromass SAPb -> SOMp
              fphys_m2(c,j) = dmin1(1.0_r8, mimicsplus_fphys_k_p1 * &
                                            dexp(mimicsplus_fphys_k_p2 * clay_frac))       !parameter equation for microbial necromass SAPf -> SOMp
              p_scalar(c,j) = 1.0_r8 / (mimicsplus_p_scalar_p1 * &
                                      dexp(mimicsplus_p_scalar_p2 * dsqrt(clay_frac)))
            end do
        end do
        initial_stock_soildepth = params_inst%mimicsplus_initial_Cstocks_depth

        !-------------------  list of pools and their attributes  ------------
        i_litr_min = 1                                                                ! metabolic litter pool LITm
        i_met_lit = i_litr_min
        floating_cn_ratio_decomp_pools(i_met_lit) = .true.
        decomp_cascade_con%decomp_pool_name_restart(i_met_lit) = 'litr1'
        decomp_cascade_con%decomp_pool_name_history(i_met_lit) = 'LIT_MET'
        decomp_cascade_con%decomp_pool_name_long(i_met_lit) = 'metabolic litter'
        decomp_cascade_con%decomp_pool_name_short(i_met_lit) = 'L1'
        is_microbe(i_met_lit) = .false.
        is_litter(i_met_lit) = .true.
        is_soil(i_met_lit) = .false.
        is_cwd(i_met_lit) = .false.
        is_mycorrhiza(i_met_lit) = .false. 
        initial_cn_ratio(i_met_lit) = 10._r8  ! 90 in BGC; not used in MIMICS
        initial_stock(i_met_lit) = params_inst%mimicsplus_initial_Cstocks(i_met_lit)
        is_metabolic(i_met_lit) = .true.
        is_cellulose(i_met_lit) = .false.
        is_lignin(i_met_lit) = .false.

        i_str_lit = i_met_lit + 1                                                     ! structural litter pool LITs
        floating_cn_ratio_decomp_pools(i_str_lit) = .true.
        decomp_cascade_con%decomp_pool_name_restart(i_str_lit) = 'litr2'
        decomp_cascade_con%decomp_pool_name_history(i_str_lit) = 'LIT_STR'
        decomp_cascade_con%decomp_pool_name_long(i_str_lit) = 'structural litter'
        decomp_cascade_con%decomp_pool_name_short(i_str_lit) = 'L2'
        is_microbe(i_str_lit) = .false.
        is_litter(i_str_lit) = .true.
        is_soil(i_str_lit) = .false.
        is_cwd(i_str_lit) = .false.
        is_mycorrhiza(i_str_lit) = .false. 
        initial_cn_ratio(i_str_lit) = 10._r8  ! 90 in BGC; not used in MIMICS
        initial_stock(i_str_lit) = params_inst%mimicsplus_initial_Cstocks(i_str_lit)
        is_metabolic(i_str_lit) = .false.
        is_cellulose(i_str_lit) = .true.
        is_lignin(i_str_lit) = .true.

        i_litr_max = i_str_lit
        if (i_litr_min /= 1 .or. i_litr_max < 2 .or. i_litr_max > 3) then
            write(iulog,*) 'Expecting i_litr_min = 1 and i_litr_max = 2 or 3.'
            write(iulog,*) 'See pftconMod, SoilBiogeochemCarbonFluxType, and'
            write(iulog,*) 'clmfates_interfaceMod for ramifications of changing'
            write(iulog,*) 'this assumption.'
            call endrun(msg='ERROR: i_litr_min and/or i_litr_max out of range '// &
              errMsg(sourcefile, __LINE__))
        end if

        i_avl_som = i_str_lit + 1                                                      ! avaliable SOM 
        floating_cn_ratio_decomp_pools(i_avl_som) = .true.
        decomp_cascade_con%decomp_pool_name_restart(i_avl_som) = 'soil1'
        decomp_cascade_con%decomp_pool_name_history(i_avl_som) = 'SOM_AVL'
        decomp_cascade_con%decomp_pool_name_long(i_avl_som) = 'available soil organic matter'
        decomp_cascade_con%decomp_pool_name_short(i_avl_som) = 'S1'
        is_microbe(i_avl_som) = .false.
        is_litter(i_avl_som) = .false.
        is_soil(i_avl_som) = .true.
        is_cwd(i_avl_som) = .false.
        is_mycorrhiza(i_avl_som) = .false. 
        initial_cn_ratio(i_avl_som) = 10._r8  ! cn_s1 in BGC; not used in MIMICS
        initial_stock(i_avl_som) = params_inst%mimicsplus_initial_Cstocks(i_avl_som)
        is_metabolic(i_avl_som) = .false.
        is_cellulose(i_avl_som) = .false.
        is_lignin(i_avl_som) = .false.

        i_chem_som = i_avl_som + 1                                                    ! chemically protected SOM
        floating_cn_ratio_decomp_pools(i_chem_som) = .true.
        decomp_cascade_con%decomp_pool_name_restart(i_chem_som) = 'soil2'
        decomp_cascade_con%decomp_pool_name_history(i_chem_som) = 'SOM_CHEM'
        decomp_cascade_con%decomp_pool_name_long(i_chem_som) = 'chemically protected soil organic matter'
        decomp_cascade_con%decomp_pool_name_short(i_chem_som) = 'S2'
        is_microbe(i_chem_som) = .false.
        is_litter(i_chem_som) = .false.
        is_soil(i_chem_som) = .true.
        is_cwd(i_chem_som) = .false.
        is_mycorrhiza(i_chem_som) = .false. 
        initial_cn_ratio(i_chem_som) = 10._r8  ! cn_s2 in BGC; not used in MIMICS
        initial_stock(i_chem_som) = params_inst%mimicsplus_initial_Cstocks(i_chem_som)
        is_metabolic(i_chem_som) = .false.
        is_cellulose(i_chem_som) = .false.
        is_lignin(i_chem_som) = .false.

        i_phys_som = i_chem_som + 1                                                   ! physically protected SOM
        floating_cn_ratio_decomp_pools(i_phys_som) = .true.
        decomp_cascade_con%decomp_pool_name_restart(i_phys_som) = 'soil3'
        decomp_cascade_con%decomp_pool_name_history(i_phys_som) = 'SOM_PHYS'
        decomp_cascade_con%decomp_pool_name_long(i_phys_som) = 'physically protected soil organic matter'
        decomp_cascade_con%decomp_pool_name_short(i_phys_som) = 'S3'
        is_microbe(i_phys_som) = .false.
        is_litter(i_phys_som) = .false.
        is_soil(i_phys_som) = .true.
        is_cwd(i_phys_som) = .false.
        is_mycorrhiza(i_phys_som) = .false. 
        initial_cn_ratio(i_phys_som) = 10._r8  ! cn_s3 in BGC; not used in MIMICS
        initial_stock(i_phys_som) = params_inst%mimicsplus_initial_Cstocks(i_phys_som)
        is_metabolic(i_phys_som) = .false.
        is_cellulose(i_phys_som) = .false.
        is_lignin(i_phys_som) = .false.

        i_cop_mic = i_phys_som + 1                                                   ! micribial bacteria
        floating_cn_ratio_decomp_pools(i_cop_mic) = .true.
        decomp_cascade_con%decomp_pool_name_restart(i_cop_mic) = 'micr1'
        decomp_cascade_con%decomp_pool_name_history(i_cop_mic) = 'MIC_COP'
        decomp_cascade_con%decomp_pool_name_long(i_cop_mic) = 'copiotrophic microbes'
        decomp_cascade_con%decomp_pool_name_short(i_cop_mic) = 'M1'
        is_microbe(i_cop_mic) = .true.
        is_litter(i_cop_mic) = .false.
        is_soil(i_cop_mic) = .false.
        is_cwd(i_cop_mic) = .false.
        is_mycorrhiza(i_cop_mic) = .false. 
        initial_cn_ratio(i_cop_mic) = 10._r8  ! MIMICS may use this
        initial_stock(i_cop_mic) = params_inst%mimicsplus_initial_Cstocks(i_cop_mic)
        is_metabolic(i_cop_mic) = .false.
        is_cellulose(i_cop_mic) = .false.
        is_lignin(i_cop_mic) = .false.

        i_oli_mic = i_cop_mic + 1                                                     ! microbial fungi
        floating_cn_ratio_decomp_pools(i_oli_mic) = .true.
        decomp_cascade_con%decomp_pool_name_restart(i_oli_mic) = 'micr2'
        decomp_cascade_con%decomp_pool_name_history(i_oli_mic) = 'MIC_OLI'
        decomp_cascade_con%decomp_pool_name_long(i_oli_mic) = 'oligotrophic microbes'
        decomp_cascade_con%decomp_pool_name_short(i_oli_mic) = 'M2'
        is_microbe(i_oli_mic) = .true.
        is_litter(i_oli_mic) = .false.
        is_soil(i_oli_mic) = .false.
        is_cwd(i_oli_mic) = .false.
        is_mycorrhiza(i_oli_mic) = .false. 
        initial_cn_ratio(i_oli_mic) = 10._r8
        initial_stock(i_oli_mic) = params_inst%mimicsplus_initial_Cstocks(i_oli_mic)
        is_metabolic(i_oli_mic) = .false.
        is_cellulose(i_oli_mic) = .false.
        is_lignin(i_oli_mic) = .false.

        i_ecm_myc = i_oli_mic + 1                                                     ! Ectomycorrhiza
        floating_cn_ratio_decomp_pools(i_ecm_myc) = .true.                            
        decomp_cascade_con%decomp_pool_name_restart(i_ecm_myc) = 'myc1'
        decomp_cascade_con%decomp_pool_name_history(i_ecm_myc) = 'MYC_ECM'
        decomp_cascade_con%decomp_pool_name_long(i_ecm_myc) = 'ectomycorrhiza'
        decomp_cascade_con%decomp_pool_name_short(i_ecm_myc) = 'MYC1'
        is_microbe(i_ecm_myc) = .false.
        is_litter(i_ecm_myc) = .false.
        is_soil(i_ecm_myc) = .false.
        is_cwd(i_ecm_myc) = .false.
        is_mycorrhiza(i_ecm_myc) = .true.
        initial_cn_ratio(i_ecm_myc) = 10._r8
        initial_stock(i_ecm_myc) = params_inst%mimicsplus_initial_Cstocks(i_ecm_myc)
        is_metabolic(i_ecm_myc) = .false.
        is_cellulose(i_ecm_myc) = .false.
        is_lignin(i_ecm_myc) = .false.

        i_am_myc = i_ecm_myc + 1                                                       ! Arbuscular Mycorrhiza
        floating_cn_ratio_decomp_pools(i_am_myc) = .true.
        decomp_cascade_con%decomp_pool_name_restart(i_am_myc) = 'myc2'
        decomp_cascade_con%decomp_pool_name_history(i_am_myc) = 'MYC_AM'
        decomp_cascade_con%decomp_pool_name_long(i_am_myc) = 'arbuscular mycorrhiza '
        decomp_cascade_con%decomp_pool_name_short(i_am_myc) = 'MYC2'
        is_microbe(i_am_myc) = .false.
        is_litter(i_am_myc) = .false.
        is_soil(i_am_myc) = .false.
        is_cwd(i_am_myc) = .false. 
        is_mycorrhiza(i_am_myc) = .true. !I invented this new catogory, for mycorrhiza! 
        initial_cn_ratio(i_am_myc) = 10._r8  ! MIMICS may use this
        initial_stock(i_am_myc) = params_inst%mimicsplus_initial_Cstocks(i_am_myc)
        is_metabolic(i_am_myc) = .false.
        is_cellulose(i_am_myc) = .false.
        is_lignin(i_am_myc) = .false.
        

        if (.not. use_fates) then
            ! CWD
            i_cwd = i_am_myc + 1
            floating_cn_ratio_decomp_pools(i_cwd) = .true.
            decomp_cascade_con%decomp_pool_name_restart(i_cwd) = 'cwd'
            decomp_cascade_con%decomp_pool_name_history(i_cwd) = 'CWD'
            decomp_cascade_con%decomp_pool_name_long(i_cwd) = 'coarse woody debris'
            decomp_cascade_con%decomp_pool_name_short(i_cwd) = 'CWD'
            is_microbe(i_cwd) = .false.
            is_litter(i_cwd) = .false.
            is_soil(i_cwd) = .false.
            is_cwd(i_cwd) = .true.
            initial_cn_ratio(i_cwd) = 10._r8  ! 90 in BGC; not used in MIMICS
            initial_stock(i_cwd) = params_inst%mimicsplus_initial_Cstocks(i_cwd)
            is_metabolic(i_cwd) = .false.
            is_cellulose(i_cwd) = .false.
            is_lignin(i_cwd) = .false.
        endif

        speedup_fac = 1._r8

        !lit1,2
        spinup_factor(i_met_lit) = 1._r8
        spinup_factor(i_str_lit) = 1._r8
        !CWD
        if (.not. use_fates) then
            spinup_factor(i_cwd) = max(1._r8, (speedup_fac * CNParamsShareInst%tau_cwd * 0.5_r8 ))
        end if
        !som1,2,3
        spinup_factor(i_avl_som) = 1._r8
        spinup_factor(i_chem_som) = 1._r8  ! BGC used cwd formula above but
        spinup_factor(i_phys_som) = 1._r8  ! ...w the respective tau_s values
        ! micr1,2
        spinup_factor(i_cop_mic) = 1._r8
        spinup_factor(i_oli_mic) = 1._r8

        !ecm, am
        spinup_factor(i_ecm_myc) = 1._r8
        spinup_factor(i_am_myc) = 1._r8

        if ( masterproc ) then
            write(iulog,*) 'Spinup_state ',spinup_state
            write(iulog,*) 'Spinup factors ',spinup_factor
        end if

        !----------------  list of transitions and their time-independent coefficients  ---------------!
        
        i_l1m1 = 1
        decomp_cascade_con%cascade_step_name(i_l1m1) = 'L1M1'
        cascade_donor_pool(i_l1m1) = i_met_lit
        cascade_receiver_pool(i_l1m1) = i_cop_mic
        nue_decomp_cascade(i_l1m1) = mimicsplus_nue_into_mic

        i_l1m2 = 2
        decomp_cascade_con%cascade_step_name(i_l1m2) = 'L1M2'
        cascade_donor_pool(i_l1m2) = i_met_lit
        cascade_receiver_pool(i_l1m2) = i_oli_mic
        nue_decomp_cascade(i_l1m2) = mimicsplus_nue_into_mic

        i_l2m1 = 3
        decomp_cascade_con%cascade_step_name(i_l2m1) = 'L2M1'
        cascade_donor_pool(i_l2m1) = i_str_lit
        cascade_receiver_pool(i_l2m1) = i_cop_mic
        nue_decomp_cascade(i_l2m1) = mimicsplus_nue_into_mic

        i_l2m2 = 4
        decomp_cascade_con%cascade_step_name(i_l2m2) = 'L2M2'
        cascade_donor_pool(i_l2m2) = i_str_lit
        cascade_receiver_pool(i_l2m2) = i_oli_mic
        nue_decomp_cascade(i_l2m2) = mimicsplus_nue_into_mic

        i_s1m1 = 5
        decomp_cascade_con%cascade_step_name(i_s1m1) = 'S1M1'
        cascade_donor_pool(i_s1m1) = i_avl_som
        cascade_receiver_pool(i_s1m1) = i_cop_mic
        nue_decomp_cascade(i_s1m1) = mimicsplus_nue_into_mic

        i_s1m2 = 6
        decomp_cascade_con%cascade_step_name(i_s1m2) = 'S1M2'
        cascade_donor_pool(i_s1m2) = i_avl_som
        cascade_receiver_pool(i_s1m2) = i_oli_mic
        nue_decomp_cascade(i_s1m2) = mimicsplus_nue_into_mic

        i_s2s1 = 7
        decomp_cascade_con%cascade_step_name(i_s2s1) = 'S2S1'
        cascade_donor_pool(i_s2s1) = i_chem_som
        cascade_receiver_pool(i_s2s1) = i_avl_som
        nue_decomp_cascade(i_s2s1) = 1.0_r8

        i_s3s1 = 8
        decomp_cascade_con%cascade_step_name(i_s3s1) = 'S3S1'
        cascade_donor_pool(i_s3s1) = i_phys_som
        cascade_receiver_pool(i_s3s1) = i_avl_som
        nue_decomp_cascade(i_s3s1) = 1.0_r8

        i_m1s1 = 9
        decomp_cascade_con%cascade_step_name(i_m1s1) = 'M1S1'
        cascade_donor_pool(i_m1s1) = i_cop_mic
        cascade_receiver_pool(i_m1s1) = i_avl_som
        nue_decomp_cascade(i_m1s1) = 1.0_r8

        i_m1s2 = 10
        decomp_cascade_con%cascade_step_name(i_m1s2) = 'M1S2'
        cascade_donor_pool(i_m1s2) = i_cop_mic
        cascade_receiver_pool(i_m1s2) = i_chem_som
        nue_decomp_cascade(i_m1s2) = 1.0_r8

        i_m1s3 = 11
        decomp_cascade_con%cascade_step_name(i_m1s3) = 'M1S3'
        cascade_donor_pool(i_m1s3) = i_cop_mic
        cascade_receiver_pool(i_m1s3) = i_phys_som
        nue_decomp_cascade(i_m1s3) = 1.0_r8

        i_m2s1 = 12
        decomp_cascade_con%cascade_step_name(i_m2s1) = 'M2S1'
        cascade_donor_pool(i_m2s1) = i_oli_mic
        cascade_receiver_pool(i_m2s1) = i_avl_som
        nue_decomp_cascade(i_m2s1) = 1.0_r8

        i_m2s2 = 13
        decomp_cascade_con%cascade_step_name(i_m2s2) = 'M2S2'
        cascade_donor_pool(i_m2s2) = i_oli_mic
        cascade_receiver_pool(i_m2s2) = i_chem_som
        nue_decomp_cascade(i_m2s2) = 1.0_r8

        i_m2s3 = 14
        decomp_cascade_con%cascade_step_name(i_m2s3) = 'M2S3'
        cascade_donor_pool(i_m2s3) = i_oli_mic
        cascade_receiver_pool(i_m2s3) = i_phys_som
        nue_decomp_cascade(i_m2s3) = 1.0_r8

        i_myc1s1 = 15                                                       ! EcM to SOMa, necromass
        decomp_cascade_con%cascade_step_name(i_myc1s1) = 'MYC1S1'
        cascade_donor_pool(i_myc1s1) = i_ecm_myc
        cascade_receiver_pool(i_myc1s1) = i_avl_som
        nue_decomp_cascade(i_myc1s1) = 1.0_r8

        i_myc1s2 = 16                                                       ! EcM to SOMc, necromass
        decomp_cascade_con%cascade_step_name(i_myc1s2) = 'MYC1S2'
        cascade_donor_pool(i_myc1s2) = i_ecm_myc
        cascade_receiver_pool(i_myc1s2) = i_chem_som
        nue_decomp_cascade(i_myc1s2) = 1.0_r8

        i_myc1s3 = 17                                                       ! EcM to SOMp, necromass
        decomp_cascade_con%cascade_step_name(i_myc1s3) = 'MYC1S3'
        cascade_donor_pool(i_myc1s3) = i_ecm_myc
        cascade_receiver_pool(i_myc1s3) = i_phys_som
        nue_decomp_cascade(i_myc1s3) = 1.0_r8

        i_myc2s1 = 18                                                       ! AM to SOMa, necromass
        decomp_cascade_con%cascade_step_name(i_myc2s1) = 'MYC2S1'
        cascade_donor_pool(i_myc2s1) = i_am_myc
        cascade_receiver_pool(i_myc2s1) = i_avl_som
        nue_decomp_cascade(i_myc2s1) = 1.0_r8

        i_myc2s2 = 19                                                       ! AM to SOMc, necromass
        decomp_cascade_con%cascade_step_name(i_myc2s2) = 'MYC2S2'
        cascade_donor_pool(i_myc2s2) = i_am_myc
        cascade_receiver_pool(i_myc2s2) = i_chem_som
        nue_decomp_cascade(i_myc2s2) = 1.0_r8

        i_myc2s3 = 20                                                       ! AM to SOMp, necromass
        decomp_cascade_con%cascade_step_name(i_myc1s3) = 'MYC2S3'
        cascade_donor_pool(i_myc2s3) = i_am_myc
        cascade_receiver_pool(i_myc2s3) = i_phys_som
        nue_decomp_cascade(i_myc2s3) = 1.0_r8

      
        if (.not. use_fates) then
            i_cwdl2 = 21
            decomp_cascade_con%cascade_step_name(i_cwdl2) = 'CWDL2'
            cascade_donor_pool(i_cwdl2) = i_cwd
            cascade_receiver_pool(i_cwdl2) = i_str_lit
            nue_decomp_cascade(i_cwdl2) = 1.0_r8
        end if

        deallocate(params_inst%mimicsplus_mge)
        deallocate(params_inst%mimicsplus_vmod)
        deallocate(params_inst%mimicsplus_vint)
        deallocate(params_inst%mimicsplus_vslope)
        deallocate(params_inst%mimicsplus_kmod)
        deallocate(params_inst%mimicsplus_kint)
        deallocate(params_inst%mimicsplus_kslope)
        deallocate(params_inst%mimicsplus_p_scalar)
        deallocate(params_inst%mimicsplus_desorp)
        deallocate(params_inst%mimicsplus_fphys_r)
        deallocate(params_inst%mimicsplus_fphys_k)
        deallocate(params_inst%mimicsplus_initial_Cstocks)

      end associate

  end subroutine init_decompcascade_mimicsplus


  !-----------------------------------------------------------------------
  subroutine decomp_rates_mimicsplus(bounds, num_bgc_soilc, filter_bgc_soilc, &
        num_soilp, filter_soilp, clm_fates, &
        soilstate_inst, temperature_inst, cnveg_carbonflux_inst, &
        ch4_inst, soilbiogeochem_carbonflux_inst, waterstatebulk_inst, &
        soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, idop)
      !
      ! !DESCRIPTION:
      ! Calculate rates and decomposition pathways for the MIMICS+
      ! decomposition cascade model
      !
      ! !USES:
      use clm_time_manager , only : get_average_days_per_year
      use clm_varcon       , only : secspday, secsphr, tfrz, spval
      use clm_varcon       , only : g_to_mg, cm3_to_m3
      use clm_varpar       , only : nlevdecomp_full
      use subgridAveMod    , only : p2c
      use PatchType        , only : patch
      use pftconMod        , only : pftname
      use TillageMod       , only : get_do_tillage
      use TillageMod       , only : get_apply_tillage_multipliers
      !
      ! !ARGUMENTS:
      type(bounds_type)                    , intent(in)    :: bounds          
      integer                              , intent(in)    :: num_soilp       ! number of soil patches in filter
      integer                              , intent(in)    :: filter_soilp(:) ! filter for soil patches
      integer                              , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
      integer                              , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
      type(soilstate_type)                 , intent(in)    :: soilstate_inst
      type(temperature_type)               , intent(in)    :: temperature_inst
      type(cnveg_carbonflux_type)          , intent(in)    :: cnveg_carbonflux_inst
      type(ch4_type)                       , intent(in)    :: ch4_inst
      type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
      type(soilbiogeochem_carbonstate_type), intent(in)    :: soilbiogeochem_carbonstate_inst
      type(soilbiogeochem_state_type)      , intent(in)    :: soilbiogeochem_state_inst
      type(waterstatebulk_type)            , intent(in)    :: waterstatebulk_inst
      type(hlm_fates_interface_type)       , intent(inout) :: clm_fates
      integer, optional                    , intent(in)    :: idop(:) ! patch day of planting
    
      !
      ! !LOCAL VARIABLES:
      real(r8), parameter :: eps = 1.e-6_r8
      real(r8), parameter :: min_modifier = 0.1     ! minimum value in microbial turnover
      real(r8):: frw(bounds%begc:bounds%endc) ! rooting fraction weight
      real(r8), allocatable:: fr(:,:)         ! column-level rooting fraction by soil depth
      real(r8), allocatable:: norm_froot_prof(:,:) ! normalized fine root profile
      real(r8):: min_froot !minimum root fraction in column
      real(r8):: max_froot !maximum root fraction in column
      real(r8):: psi                          ! temporary soilpsi for water scalar
      real(r8):: k_frag                       ! fragmentation rate constant CWD (1/sec)
      real(r8):: moist_mod                    ! moisture modifier to replace w_scalar
      real(r8):: fmet
      real(r8):: favl
      real(r8):: fchem_m1
      real(r8):: fchem_m2
      real(r8):: desorption
      real(r8):: vmax_l1_m1  !
      real(r8):: vmax_l2_m1  !
      real(r8):: vmax_s1_m1  !
      real(r8):: vmax_l1_m2  !
      real(r8):: vmax_l2_m2  !
      real(r8):: vmax_s1_m2  !
      real(r8):: km_l1_m1  !
      real(r8):: km_l2_m1  !
      real(r8):: km_s1_m1  !
      real(r8):: km_l1_m2  !
      real(r8):: km_l2_m2  !
      real(r8):: km_s1_m2  !
      real(r8):: tau_m1  ! Microbial tunrnover rate bacteria
      real(r8):: tau_m2  ! Mircobial turnover rate fungi
      real(r8):: tau_myc ! Mycorrhizal turnover rate for EcM & AM 
      real(r8):: tau_mod
      real(r8):: m1_conc       ! Carbon concentration in microbial barteria pool
      real(r8):: m2_conc       ! Carbon concentration in microbial fungi pool
      real(r8):: myc1_conc     ! Carbon concentration in ectomycorrhiza pool
      real(r8):: myc2_conc     ! Carbon concentration in abruscular mycorrhiza pool
      real(r8):: avl_som_conc  ! Carbon concentration in avaliable SOM pool
      real(r8):: chem_som_conc ! Carbon concentration in chemically protected SOM pool
      real(r8):: phys_som_conc ! Carbon concentration in physically protected SOM pool
      real(r8):: term_1  !
      real(r8):: term_2  !
      real(r8):: t_soi_degC
  !   real(r8):: decomp_depth_efolding        ! (meters) e-folding depth for reduction in decomposition [
      integer :: p, fp, c, fc, j, k, l, s  ! indices
      integer :: pf  ! fates patch index
      integer :: nc  ! clump index
      real(r8):: days_per_year                ! days per year
      real(r8):: depth_scalar(bounds%begc:bounds%endc,1:nlevdecomp) 
      real(r8):: w_d_o_scalars  ! product of w_scalar * depth_scalar * o_scalar
      real(r8):: mino2lim                     !minimum anaerobic decomposition rate
      real(r8):: mimicsplus_fmet_p1
      real(r8):: mimicsplus_fmet_p2
      real(r8):: mimicsplus_fmet_p3
      real(r8):: mimicsplus_fmet_p4
      real(r8):: mimicsplus_fchem_r_p1 !ECW in my parameter file this has only 2 values, in testbed and CLM it has 3
      real(r8):: mimicsplus_fchem_r_p2
      real(r8):: mimicsplus_fchem_r_p3
      real(r8):: mimicsplus_fchem_k_p1
      real(r8):: mimicsplus_fchem_k_p2
      real(r8):: mimicsplus_fchem_k_p3
      real(r8):: mimicsplus_tau_mod_min
      real(r8):: mimicsplus_tau_mod_max
      real(r8):: mimicsplus_tau_mod_factor
      real(r8):: mimicsplus_tau_r_p1
      real(r8):: mimicsplus_tau_r_p2
      real(r8):: mimicsplus_tau_k_p1
      real(r8):: mimicsplus_tau_k_p2
      real(r8):: mimicsplus_ko_r
      real(r8):: mimicsplus_ko_k
      real(r8):: mimicsplus_densdep
      real(r8):: mimicsplus_desorpQ10
      real(r8):: mimicsplus_t_soi_ref
      real(r8):: mimicsplus_cn_mod_num
      real(r8):: mimicsplus_cn_r
      real(r8):: mimicsplus_cn_k

      real(r8):: mimicsplus_fphys_ecm
      real(r8):: mimicsplus_fphys_am
      real(r8):: mimicsplus_fchem_ecm
      real(r8):: mimicsplus_fchem_am
      real(r8):: mimicsplus_tau_ecm
      real(r8):: mimicsplus_tau_am

      real(r8):: mimicsplus_cn_myc
      !real(r8):: mimicsplus_r_myc

      real(r8):: mimicsplus_mge_am
      real(r8):: mimicsplus_mge_ecm
      real(r8):: mimicsplus_k_m_emyc
      real(r8):: mimicsplus_vmax_myc
      real(r8):: mimicsplus_k_mo
      real(r8):: mimicsplus_k_myc_som

      real(r8):: tau_myc1
      real(r8):: tau_myc2
      real(r8):: fchem_myc1
      real(r8):: fchem_myc2
      real(r8):: fphys_myc1
      real(r8):: fphys_myc2

      real(r8):: fphys_avl       ! desorpion
      real(r8):: fchem_avl       ! oxidation

      !--------------------------------------------

      real(r8):: spinup_geogterm_l1(bounds%begc:bounds%endc) ! geographically-varying spinup term for l1
      real(r8):: spinup_geogterm_l2(bounds%begc:bounds%endc)  ! geographically-varying spinup term for l2
      real(r8):: spinup_geogterm_cwd(bounds%begc:bounds%endc) ! geographically-varying spinup term for cwd
      real(r8):: spinup_geogterm_s1(bounds%begc:bounds%endc) ! geographically-varying spinup term for s1
      real(r8):: spinup_geogterm_s2(bounds%begc:bounds%endc) ! geographically-varying spinup term for s2
      real(r8):: spinup_geogterm_s3(bounds%begc:bounds%endc) ! geographically-varying spinup term for s3
      real(r8):: spinup_geogterm_m1(bounds%begc:bounds%endc)  ! geographically-varying spinup term for m1
      real(r8):: spinup_geogterm_m2(bounds%begc:bounds%endc)  ! geographically-varying spinup term for m2
      real(r8):: spinup_geogterm_myc1(bounds%begc:bounds%endc) ! geographically-varying spinup term for myc1
      real(r8):: spinup_geogterm_myc2(bounds%begc:bounds%endc) ! geographically-varying spinup term for myc2 !ECW check with MATVEY for some reason litter and SOM are also mentioned in: module SoilBiogeochemDecompCascadeBGCMod, but not m1 und m2
      real(r8):: annsum_npp_col_local(bounds%begc:bounds%endc)  ! local annual sum of NPP at the column level
      real(r8):: annsum_npp(bounds%begp:bounds%endp)  ! local annual sum of NPP at the patch level
      real(r8):: annsum_npp_col_scalar  ! annual sum of NPP, scalar in column-level loop

      !-----------------------------------------------------------------------

      associate(                                                           &

            rf_cwdl2       => CNParamsShareInst%rf_cwdl2                  , & ! Input:  [real(r8)         ]  respiration fraction in CWD to litter2 transition (frac)
            minpsi         => CNParamsShareInst%minpsi                    , & ! Input:  [real(r8)         ]  minimum soil suction (mm)
            maxpsi         => CNParamsShareInst%maxpsi                    , & ! Input:  [real(r8)         ]  maximum soil suction (mm)
            soilpsi        => soilstate_inst%soilpsi_col                  , & ! Input:  [real(r8) (:,:)   ]  soil water potential in each soil layer (MPa)          
            watsat         => soilstate_inst%watsat_col                   , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity)  
            t_soisno       => temperature_inst%t_soisno_col               , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       
            o2stress_sat   => ch4_inst%o2stress_sat_col                   , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
            o2stress_unsat => ch4_inst%o2stress_unsat_col                 , & ! Input:  [real(r8) (:,:)   ]  Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
            finundated     => ch4_inst%finundated_col                     , & ! Input:  [real(r8) (:)     ]  fractional inundated area                                
            decomp_cpools_vr => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col , &  ! Input: [real(r8) (:,:,:) ] (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) C pools
            pathfrac_decomp_cascade => soilbiogeochem_carbonflux_inst%pathfrac_decomp_cascade_col, &  ! Output: [real(r8) (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)
            froot_prof     => soilbiogeochem_state_inst%froot_prof_patch  , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots                     
            h2osoi_liq       =>    waterstatebulk_inst%h2osoi_liq_col     , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)             
            h2osoi_ice       =>    waterstatebulk_inst%h2osoi_ice_col     , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2) (new)
            rf_decomp_cascade       => soilbiogeochem_carbonflux_inst%rf_decomp_cascade_col    , & ! Input:  [real(r8)          (:,:,:) ]  respired fraction in decomposition step (frac)
            w_scalar       => soilbiogeochem_carbonflux_inst%w_scalar_col , & ! Output: [real(r8) (:,:)   ]  soil water scalar for decomp                           
            o_scalar       => soilbiogeochem_carbonflux_inst%o_scalar_col , & ! Output: [real(r8) (:,:)   ]  fraction by which decomposition is limited by anoxia   
            cn_col         => soilbiogeochem_carbonflux_inst%cn_col       , & ! Output: [real(r8) (:,:)   ]  C:N ratio
            ligninNratioAvg => soilbiogeochem_carbonflux_inst%litr_lig_c_to_n_col, &  ! Input: [real(r8) (:) ] C:N ratio of litter lignin
            decomp_k       => soilbiogeochem_carbonflux_inst%decomp_k_col , & ! Output: [real(r8) (:,:,:) ]  rate for decomposition (1./sec)
            spinup_factor  => decomp_cascade_con%spinup_factor              & ! Input:  [real(r8)          (:)     ]  factor for AD spinup associated with each pool           
            )

            if (get_do_tillage() .and. .not. present(idop)) then
               call endrun("Do not enable tillage without providing idop to decomp_rate_constants_mimics().")
            end if

        mino2lim = CNParamsShareInst%mino2lim

        days_per_year = get_average_days_per_year()

  !     ! Set "decomp_depth_efolding" parameter
  !     decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding

        ! translate to per-second time constant
        k_frag = 1._r8 / (secspday * days_per_year * CNParamsShareInst%tau_cwd)

        allocate(norm_froot_prof(bounds%begc:bounds%endc, 1:nlevdecomp_full))
        norm_froot_prof(bounds%begc:bounds%endc, 1:nlevdecomp_full)=0._r8

      ! calc ref rate
        if ( spinup_state >= 1 ) then
            do fc = 1,num_bgc_soilc
              c = filter_bgc_soilc(fc)
              !
              if ( abs(spinup_factor(i_met_lit) - 1._r8) .gt. eps) then
                  spinup_geogterm_l1(c) = spinup_factor(i_met_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              else
                  spinup_geogterm_l1(c) = 1._r8
              endif
              !
              if ( abs(spinup_factor(i_str_lit) - 1._r8) .gt. eps) then
                  spinup_geogterm_l2(c) = spinup_factor(i_str_lit) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              else
                  spinup_geogterm_l2(c) = 1._r8
              endif
              !
              if ( .not. use_fates ) then
                  if ( abs(spinup_factor(i_cwd) - 1._r8) .gt. eps) then
                    spinup_geogterm_cwd(c) = spinup_factor(i_cwd) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
                  else
                    spinup_geogterm_cwd(c) = 1._r8
                  endif
              endif
              !
              if ( abs(spinup_factor(i_avl_som) - 1._r8) .gt. eps) then
                  spinup_geogterm_s1(c) = spinup_factor(i_avl_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              else
                  spinup_geogterm_s1(c) = 1._r8
              endif
              !
              if ( abs(spinup_factor(i_chem_som) - 1._r8) .gt. eps) then
                  spinup_geogterm_s2(c) = spinup_factor(i_chem_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              else
                  spinup_geogterm_s2(c) = 1._r8
              endif
              !
              if ( abs(spinup_factor(i_phys_som) - 1._r8) .gt. eps) then
                  spinup_geogterm_s3(c) = spinup_factor(i_phys_som) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              else
                  spinup_geogterm_s3(c) = 1._r8
              endif
              !
              if ( abs(spinup_factor(i_cop_mic) - 1._r8) .gt. eps) then
                  spinup_geogterm_m1(c) = spinup_factor(i_cop_mic) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              else
                  spinup_geogterm_m1(c) = 1._r8
              endif
              !
              if ( abs(spinup_factor(i_oli_mic) - 1._r8) .gt. eps) then
                  spinup_geogterm_m2(c) = spinup_factor(i_oli_mic) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              else
                  spinup_geogterm_m2(c) = 1._r8
              endif
              
              if ( abs(spinup_factor(i_ecm_myc) - 1._r8) .gt. eps) then
                 spinup_geogterm_myc1(c) = spinup_factor(i_ecm_myc) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
              else
                 spinup_geogterm_myc1(c) = 1._r8
              endif
            !
            if ( abs(spinup_factor(i_am_myc) - 1._r8) .gt. eps) then
              spinup_geogterm_myc2(c) = spinup_factor(i_am_myc) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
            else
              spinup_geogterm_myc2(c) = 1._r8
            endif
              !
            end do
        else
            do fc = 1,num_bgc_soilc
              c = filter_bgc_soilc(fc)
              spinup_geogterm_l1(c) = 1._r8
              spinup_geogterm_l2(c) = 1._r8
              spinup_geogterm_cwd(c) = 1._r8
              spinup_geogterm_s1(c) = 1._r8
              spinup_geogterm_s2(c) = 1._r8
              spinup_geogterm_s3(c) = 1._r8
              spinup_geogterm_m1(c) = 1._r8
              spinup_geogterm_m2(c) = 1._r8
              spinup_geogterm_myc1(c) = 1._r8 
              spinup_geogterm_myc2(c) = 1._r8 
            end do
        endif

        ! root fraction is used for calculating microbial tunrover rate, it functions as a sort of depth profile
        ! get root fraction from patch to column level
        call p2c(bounds, nlevdecomp_full, &
                 froot_prof(bounds%begp:bounds%endp,1:nlevdecomp_full), &
                 norm_froot_prof(bounds%begc:bounds%endc,1:nlevdecomp_full), 'unity')
        ! normalize root fraction
        do fc = 1,num_bgc_soilc
           c = filter_bgc_soilc(fc)
              min_froot = minval(norm_froot_prof(c,1:nlevdecomp))
              max_froot = maxval(norm_froot_prof(c,1:nlevdecomp))
           do j = 1,nlevdecomp
              ! all turnover will be happening in just 1 soil layer
              if ( nlevdecomp .eq. 1 ) then
                 norm_froot_prof(c,j) = 1._r8
              else
                 if (max_froot - min_froot == 0._r8 .or. min_froot == spval .or. max_froot == spval) then
                    norm_froot_prof(c,j) = min_modifier
                 else
                    norm_froot_prof(c,j) = (norm_froot_prof(c,j) - min_froot) / (max_froot - min_froot)
                 endif
              endif
           enddo
         enddo


        !--- time dependent coefficients-----!
        if ( nlevdecomp .eq. 1 ) then

            ! calculate function to weight the temperature and water potential scalars
            ! for decomposition control.  

            ! the following normalizes values in fr so that they
            ! sum to 1.0 across top nlevdecomp levels on a column
            frw(bounds%begc:bounds%endc) = 0._r8
            allocate(fr(bounds%begc:bounds%endc,nlev_soildecomp_standard))
            do j=1,nlev_soildecomp_standard
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  frw(c) = frw(c) + col%dz(c,j)
               end do
            end do
            do j = 1,nlev_soildecomp_standard
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if (frw(c) /= 0._r8) then
                     fr(c,j) = col%dz(c,j) / frw(c)
                  else
                     fr(c,j) = 0._r8
                  end if
               end do
            end do

            ! calculate the rate constant scalar for soil water content.
            ! Uses the log relationship with water potential given in
            ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
            ! a comparison of models. Ecology, 68(5):1190-1200.
            ! and supported by data in
            ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
            ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

            do j = 1,nlev_soildecomp_standard
              do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if (j==1) w_scalar(c,:) = 0._r8
                  psi = min(soilpsi(c,j),maxpsi)
                  ! decomp only if soilpsi is higher than minpsi
                  if (psi > minpsi) then
                    w_scalar(c,1) = w_scalar(c,1) + (log(minpsi/psi)/log(minpsi/maxpsi))*fr(c,j)
                  end if
              end do
            end do

            ! Calculate ANOXIA
            ! anoxia = .true. when (use_lch4)

            if (anoxia) then

              do j = 1,nlev_soildecomp_standard
                  do fc = 1,num_bgc_soilc
                    c = filter_bgc_soilc(fc)

                    if (j==1) o_scalar(c,:) = 0._r8

                    o_scalar(c,1) = o_scalar(c,1) + fr(c,j) * max(o2stress_unsat(c,j), mino2lim)
                  end do
              end do
            else
              o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
            end if

            deallocate(fr)

        else

            ! calculate the rate constant scalar for soil water content.
            ! Uses the log relationship with water potential given in
            ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
            ! a comparison of models. Ecology, 68(5):1190-1200.
            ! and supported by data in
            ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
            ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

            do j = 1,nlevdecomp
              do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  psi = min(soilpsi(c,j),maxpsi)
                  ! decomp only if soilpsi is higher than minpsi
                  if (psi > minpsi) then
                    w_scalar(c,j) = (log(minpsi/psi)/log(minpsi/maxpsi))
                  else
                    w_scalar(c,j) = 0._r8
                  end if
              end do
            end do

            ! Calculate ANOXIA
            ! anoxia = .true. when (use_lch4)

            if (anoxia) then
              do j = 1,nlevdecomp
                  do fc = 1,num_bgc_soilc
                    c = filter_bgc_soilc(fc)

                    o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim)
                  end do
              end do
            else
              o_scalar(bounds%begc:bounds%endc,1:nlevdecomp) = 1._r8
            end if

        end if

        ! Term that reduces decomposition rate at depth
        ! Placeholder. For now depth_scalar = 1.
        do j = 1, nlevdecomp
            do fc = 1, num_bgc_soilc
              c = filter_bgc_soilc(fc)
              ! Using fixed e-folding depth as in
              ! SoilBiogeochemDecompCascadeBGCMod.F90
              ! depth_scalar(c,j) = exp(-zsoi(j) / decomp_depth_efolding)
              depth_scalar(c,j) = 1.0_r8
            end do
        end do 

        ! TODO @ekluzek suggested possibly making the Left Hand Sides into arrays
        ! and I wonder in that case whether to skip these assignments altogether
        ! and use the Right Hand Sides directly
        mimicsplus_fmet_p1 = params_inst%mimicsplus_fmet(1)
        mimicsplus_fmet_p2 = params_inst%mimicsplus_fmet(2)
        mimicsplus_fmet_p3 = params_inst%mimicsplus_fmet(3)
        mimicsplus_fmet_p4 = params_inst%mimicsplus_fmet(4)
        mimicsplus_fchem_r_p1 = params_inst%mimicsplus_fchem_r(1)
        mimicsplus_fchem_r_p2 = params_inst%mimicsplus_fchem_r(2)
        mimicsplus_fchem_r_p3 = params_inst%mimicsplus_fchem_r(3)
        mimicsplus_fchem_k_p1 = params_inst%mimicsplus_fchem_k(1)
        mimicsplus_fchem_k_p2 = params_inst%mimicsplus_fchem_k(2)
        mimicsplus_fchem_k_p3 = params_inst%mimicsplus_fchem_k(3)
        mimicsplus_tau_mod_min = params_inst%mimicsplus_tau_mod_min
        mimicsplus_tau_mod_max = params_inst%mimicsplus_tau_mod_max
        mimicsplus_tau_mod_factor = params_inst%mimicsplus_tau_mod_factor
        mimicsplus_tau_r_p1 = params_inst%mimicsplus_tau_r(1)
        mimicsplus_tau_r_p2 = params_inst%mimicsplus_tau_r(2)
        mimicsplus_tau_k_p1 = params_inst%mimicsplus_tau_k(1)
        mimicsplus_tau_k_p2 = params_inst%mimicsplus_tau_k(2)

        !assign parameters 
        !some parameters have several values in the parameter file:
        !then it has to be specified (x) which one is taken, otherwise just assigned
        mimicsplus_fphys_ecm = params_inst%mimicsplus_fphys_ecm
        mimicsplus_fphys_am = params_inst%mimicsplus_fphys_am
        mimicsplus_fchem_ecm = params_inst%mimicsplus_fchem_ecm
        mimicsplus_fchem_am = params_inst%mimicsplus_fchem_am
        mimicsplus_tau_ecm = params_inst%mimicsplus_tau_ecm
        mimicsplus_tau_am = params_inst%mimicsplus_tau_am

        mimicsplus_cn_myc = params_inst%mimicsplus_cn_myc
        !mimicsplus_r_myc = params_inst%mimicsplus_r_myc

        mimicsplus_mge_am = params_inst%mimicsplus_mge_am
        mimicsplus_mge_ecm = params_inst%mimicsplus_mge_ecm
        mimicsplus_k_m_emyc = params_inst%mimicsplus_k_m_emyc
        mimicsplus_vmax_myc = params_inst%mimicsplus_vmax_myc
        mimicsplus_k_mo = params_inst%mimicsplus_k_mo
        mimicsplus_k_myc_som = params_inst%mimicsplus_k_myc_som
        !-----------------------------------------------------------------------------

        mimicsplus_ko_r = params_inst%mimicsplus_ko_r
        mimicsplus_ko_k = params_inst%mimicsplus_ko_k
        mimicsplus_densdep = params_inst%mimicsplus_densdep
        mimicsplus_desorpQ10 = params_inst%mimicsplus_desorpQ10
        mimicsplus_t_soi_ref = params_inst%mimicsplus_t_soi_ref
        mimicsplus_cn_mod_num = params_inst%mimicsplus_cn_mod_num
        mimicsplus_cn_r = params_inst%mimicsplus_cn_r
        mimicsplus_cn_k = params_inst%mimicsplus_cn_k

        ! If FATES-MIMICS, then use FATES copy of annsum_npp.
        ! The FATES copy of annsum_npp is available when use_lch4 = .true., so
        ! we limit FATES-MIMICS to if (use_lch4).
        fates_if: if (use_fates) then
            lch4_if: if (use_lch4) then

              ! Loop over p to get FATES copy of annsum_npp
              nc = bounds%clump_index
              do fp = 1, num_soilp

                  p = filter_soilp(fp)
                  c = patch%column(p)

                  pf = p - col%patchi(c)
                  s  = clm_fates%f2hmap(nc)%hsites(c)
                  annsum_npp(p) = clm_fates%fates(nc)%bc_out(s)%annsum_npp_pa(pf)

                  ! Initialize local column-level annsum_npp before averaging
                  annsum_npp_col_local(c) = 0._r8

              end do  ! p loop

              ! Calculate the column-level average !MVD from patch to colum wide use for root frac      
              ! Is this not done automatically? on top for if fates is used and below for without fates         

              call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
                  annsum_npp(bounds%begp:bounds%endp), &
                  annsum_npp_col_local(bounds%begc:bounds%endc))
            else
              call endrun(msg='ERROR: soil_decomp_method = MIMICSWieder2015 or MIMICSplusAas2023 '// &
              'will work with use_fates = .true. only if use_lch4 = .true. '// &
              errMsg(sourcefile, __LINE__))
            end if lch4_if
        end if fates_if

        ! calculate rates for all litter and som pools
        do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)

            if (use_fates) then
              annsum_npp_col_scalar = max(0._r8, annsum_npp_col_local(c))
            else
              annsum_npp_col_scalar = max(0._r8, cnveg_carbonflux_inst%annsum_npp_col(c))
            end if

            ! Time-dependent params from Wieder et al. 2015 & testbed code

            ! Set limits on Lignin:N to keep fmet > 0.2
            ! Necessary for litter quality in boreal forests with high cwd flux
            ! TODO Check for high-freq variations in ligninNratioAvg. To avoid,
            !      replace pool_to_litter terms with ann or other long term mean
            !      in CNVegCarbonFluxType.
            fmet = mimicsplus_fmet_p1 * (mimicsplus_fmet_p2 - mimicsplus_fmet_p3 * &
            min(mimicsplus_fmet_p4, ligninNratioAvg(c)))
            tau_mod = min(mimicsplus_tau_mod_max, max(mimicsplus_tau_mod_min, &
            sqrt(mimicsplus_tau_mod_factor * annsum_npp_col_scalar)))

            ! These two get used in SoilBiogeochemPotentialMod.F90
            ! cn(c,i_cop_mic), cn(c,i_oli_mic) are CN_r, CN_k in the testbed code
            ! cn_r, cn_k are CNr, CNk in the testbed code
            ! https://github.com/wwieder/biogeochem_testbed_1.1/blob/Testbed_CN/SOURCE_CODE/mimicsplus_cycle_CN.f90#L1613
            cn_col(c,i_cop_mic) = mimicsplus_cn_r * sqrt(mimicsplus_cn_mod_num / fmet)
            cn_col(c,i_oli_mic) = mimicsplus_cn_k * sqrt(mimicsplus_cn_mod_num / fmet)

            
            ! Parameter equation microbial necromass to chemically protected SOM
            fchem_m1 = min(1._r8, max(0._r8, mimicsplus_fchem_r_p1 * &
                    exp(mimicsplus_fchem_r_p2 * fmet)))                       !parameter equation for microbial necromass SAPb -> SOMc
            fchem_m2 = min(1._r8, max(0._r8, mimicsplus_fchem_k_p1 * &
                    exp(mimicsplus_fchem_k_p2 * fmet)))                       !parameter equation for microbial necromass SAPb -> SOMc

            do j = 1,nlevdecomp
              ! vmax ends up in units of per hour but is expected
              ! in units of per second, so convert here; alternatively
              ! place the conversion once in w_d_o_scalars
              ! Table B1 Wieder et al. 2015 & MIMICS params file give diff
              ! kslope. I used the params file value(s).
              t_soi_degC = t_soisno(c,j) - tfrz

              vmax_l1_m1 = exp(vslope_l1_m1 * t_soi_degC + vint_l1_m1) * &
                          vmod_l1_m1 / secsphr
              vmax_l1_m2 = exp(vslope_l1_m2 * t_soi_degC + vint_l1_m2) * &
                          vmod_l1_m2 / secsphr
              vmax_l2_m1 = exp(vslope_l2_m1 * t_soi_degC + vint_l2_m1) * &
                          vmod_l2_m1 / secsphr
              vmax_l2_m2 = exp(vslope_l2_m2 * t_soi_degC + vint_l2_m2) * &
                          vmod_l2_m2 / secsphr
              vmax_s1_m1 = exp(vslope_s1_m1 * t_soi_degC + vint_s1_m1) * &
                          vmod_s1_m1 / secsphr
              vmax_s1_m2 = exp(vslope_s1_m2 * t_soi_degC + vint_s1_m2) * &
                          vmod_s1_m2 / secsphr

              km_l1_m1 = exp(kslope_l1_m1 * t_soi_degC + kint_l1_m1) * &
                        kmod_l1_m1
              km_l1_m2 = exp(kslope_l1_m2 * t_soi_degC + kint_l1_m2) * &
                        kmod_l1_m2
              km_l2_m1 = exp(kslope_l2_m1 * t_soi_degC + kint_l2_m1) * &
                        kmod_l2_m1
              km_l2_m2 = exp(kslope_l2_m2 * t_soi_degC + kint_l2_m2) * &
                        kmod_l2_m2
              km_s1_m1 = exp(kslope_s1_m1 * t_soi_degC + kint_s1_m1) * &
                        kmod_s1_m1 * p_scalar(c,j)
              km_s1_m2 = exp(kslope_s1_m2 * t_soi_degC + kint_s1_m2) * &
                        kmod_s1_m2 * p_scalar(c,j)

              ! Desorption a function of soil temperature and
              ! Q10 = 1.1 w/ reference temperature of 25C.
              ! Expected in units of per second, so convert; alternatively
              ! place the conversion once in w_d_o_scalars
              desorption = (desorp(c,j) / secsphr) * mimicsplus_desorpQ10 * &
                          exp((t_soi_degC - mimicsplus_t_soi_ref) / 10.0_r8)

              ! Microbial concentration with necessary unit conversions
              ! mgC/cm3 = gC/m3 * (1e3 mg/g) / (1e6 cm3/m3)
              m1_conc = (decomp_cpools_vr(c,j,i_cop_mic) / col%dz(c,j)) * & 
                        g_to_mg * cm3_to_m3
              m2_conc = (decomp_cpools_vr(c,j,i_oli_mic) / col%dz(c,j)) * &
                        g_to_mg * cm3_to_m3

              ! Mycorrhizal concentration with necerssary unit conversions
                        ! mgC/cm3 = gC/m3 * (1e3 mg/g) / (1e6 cm3/m3)
              myc1_conc = (decomp_cpools_vr(c,j,i_ecm_myc) / col%dz(c,j)) * & !ECW is this what goes into ROI function (see end of module)
                        g_to_mg * cm3_to_m3
              myc2_conc = (decomp_cpools_vr(c,j,i_am_myc) / col%dz(c,j)) * &
                        g_to_mg * cm3_to_m3

             ! Soil organic matter concentration with necerssary unit conversions
                        ! mgC/cm3 = gC/m3 * (1e3 mg/g) / (1e6 cm3/m3)
              avl_som_conc = (decomp_cpools_vr(c,j,i_avl_som) / col%dz(c,j)) * & !ECW This is for Baskaram N mining
                        g_to_mg * cm3_to_m3
              chem_som_conc = (decomp_cpools_vr(c,j,i_chem_som) / col%dz(c,j)) * &
                        g_to_mg * cm3_to_m3
              phys_som_conc = (decomp_cpools_vr(c,j,i_phys_som) / col%dz(c,j)) * &
                        g_to_mg * cm3_to_m3

            
              tau_myc = mimicsplus_k_myc_som / secsphr ! Turnover rate mycorrhiza with unit conversions (hourly -> second)

              ! Product of w_scalar * depth_scalar * o_scalar
              w_d_o_scalars = w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
              moist_mod = r_moist(h2osoi_liq(c,j),watsat(c,j), h2osoi_ice(c,j), col%dz(c,j)) 
              
              ! Turnover rate microbes with unit conversions (hourly -> second)
              ! modifies turnover based on soil temperature and depth (norm_froot_prof)
              if (t_soi_degC < 0) then
                 tau_m1 = mimicsplus_tau_r_p1 * exp(mimicsplus_tau_r_p2 * fmet) * min_modifier / &
                        secsphr
                 tau_m2 = mimicsplus_tau_k_p1 * exp(mimicsplus_tau_k_p2 * fmet) * min_modifier / &
                        secsphr
              else
                 tau_m1 = mimicsplus_tau_r_p1 * exp(mimicsplus_tau_r_p2 * fmet) * max(norm_froot_prof(c,j),min_modifier) / &
                        secsphr
                 tau_m2 = mimicsplus_tau_k_p1 * exp(mimicsplus_tau_k_p2 * fmet) * max(norm_froot_prof(c,j),min_modifier) / &
                        secsphr
              end if

              ! RATE calculatons: depending on how much C is avaliable, the rate changes 
              ! rate = the part of the flux equation which is multiplied by the donor pool
              ! C amount from the donor pool is only added later, because it will still change over time

              ! decomp_k used in SoilBiogeochemPotentialMod.F90
              ! also updating pathfrac terms that vary with time
              term_1 = vmax_l1_m1 * m1_conc / (km_l1_m1 + m1_conc)                     ! rate: metabolic litter -> microbial bacteria
              term_2 = vmax_l1_m2 * m2_conc / (km_l1_m2 + m2_conc)                     ! rate: metabolic litter -> microbial fungi
              decomp_k(c,j,i_met_lit) = (term_1 + term_2) * moist_mod !w_d_o_scalars
              if (term_1 + term_2 /= 0._r8) then
                  pathfrac_decomp_cascade(c,j,i_l1m1) = term_1 / (term_1 + term_2)
                  pathfrac_decomp_cascade(c,j,i_l1m2) = term_2 / (term_1 + term_2)
              else
                  pathfrac_decomp_cascade(c,j,i_l1m1) = 0._r8
                  pathfrac_decomp_cascade(c,j,i_l1m2) = 0._r8
              end if
              
              term_1 = vmax_l2_m1 * m1_conc / (km_l2_m1 + m1_conc)                      !rate: structural litter -> microbial bacteria
              term_2 = vmax_l2_m2 * m2_conc / (km_l2_m2 + m2_conc)                      !rate: structural litter -> microbial fungi
              decomp_k(c,j,i_str_lit) = (term_1 + term_2) * moist_mod !w_d_o_scalars
              if (term_1 + term_2 /= 0._r8) then
                  pathfrac_decomp_cascade(c,j,i_l2m1) = term_1 / (term_1 + term_2)
                  pathfrac_decomp_cascade(c,j,i_l2m2) = term_2 / (term_1 + term_2)
              else
                  pathfrac_decomp_cascade(c,j,i_l2m1) = 0._r8
                  pathfrac_decomp_cascade(c,j,i_l2m2) = 0._r8
              end if

            !Microbes can access SOMa- (not SOMc and SOMp)
              term_1 = vmax_s1_m1 * m1_conc / (km_s1_m1 + m1_conc)                      ! rate: avaliable SOM -> microbial bacteria
              term_2 = vmax_s1_m2 * m2_conc / (km_s1_m2 + m2_conc)                      ! rate: avaliable SOM -> microbial fungi
              decomp_k(c,j,i_avl_som) = (term_1 + term_2) * moist_mod !w_d_o_scalars               !decomp_k calculates the rate of decomposition per second based on: c soil coulum filter | j soil layer | and the carbon in the pool
              if (term_1 + term_2 /= 0._r8) then
                  pathfrac_decomp_cascade(c,j,i_s1m1) = term_1 / (term_1 + term_2)
                  pathfrac_decomp_cascade(c,j,i_s1m2) = term_2 / (term_1 + term_2)
              else
                  pathfrac_decomp_cascade(c,j,i_s1m1) = 0._r8
                  pathfrac_decomp_cascade(c,j,i_s1m2) = 0._r8
              end if

              
              ! Desorption and Oxidation processes microbes & mycorrhiza
              ! Here are Fluxes C11 & C12, further down under mining rates are C25 & C26, they need to be added to SOMa pool still!
              !!!CHECK if depth_scalar(c,j) is the same as dz (in Elin):
              ! dz is used in mycorrhizal equations, because parameter values come form Bashkaran et al. need to be transformed into C/m^3 
              fphys_avl     = desorption * depth_scalar(c,j)                                   !Rate for microbial desorption SOMp -> SOMa
              decomp_k(c,j,i_phys_som) = fphys_avl                                             


              term_1 = vmax_l2_m1 * m1_conc / (mimicsplus_ko_r * km_l2_m1 + m1_conc)     
              term_2 = vmax_l2_m2 * m2_conc / (mimicsplus_ko_k * km_l2_m2 + m2_conc)     
              ! The right hand side is OXIDAT in the testbed (line 1145)
              fchem_avl = (term_1 + term_2) * moist_mod !w_d_o_scalars                          !Rate for microbial oxidation SOMc -> SOMa
              decomp_k(c,j,i_chem_som) = (fchem_avl)                           

                  

              !Parameter equations microbial necromass

              decomp_k(c,j,i_cop_mic) = tau_m1 * &
                    m1_conc**(mimicsplus_densdep - 1.0_r8) * moist_mod !w_d_o_scalars    !decomposition rate for microbial bacteria
              favl = min(1.0_r8, max(0.0_r8, 1.0_r8 - fphys_m1(c,j) - fchem_m1))         !parameter equation for microbial necromass SAPb -> SOMa
              pathfrac_decomp_cascade(c,j,i_m1s1) = favl                                 !stores parameter to be transferred in next module
              pathfrac_decomp_cascade(c,j,i_m1s2) = fchem_m1                             !stores parameter to be transferred in next module

              decomp_k(c,j,i_oli_mic) = tau_m2 * &
                    m2_conc**(mimicsplus_densdep - 1.0_r8) * moist_mod !w_d_o_scalars    !decomposition rate for microbial fungi
              favl = min(1.0_r8, max(0.0_r8, 1.0_r8 - fphys_m2(c,j) - fchem_m2))         !parameter equation for microbial necromass SAPf -> SOMa
              pathfrac_decomp_cascade(c,j,i_m2s1) = favl                                 !stores parameter to be transferred in next module
              pathfrac_decomp_cascade(c,j,i_m2s2) = fchem_m2                             !stores parameter to be transferred in next module

              
              !ECW parameter for fraction of mycorrjizal necromass to SOM pools

              tau_myc1 = min(1._r8, max(0._r8, mimicsplus_tau_ecm))
              fchem_myc1 = min(1._r8, max(0._r8, mimicsplus_fchem_ecm))
              fphys_myc1 = min(1._r8, max(0._r8, mimicsplus_fphys_ecm))
              tau_myc2 = min(1._r8, max(0._r8, mimicsplus_tau_am))
              fchem_myc2 = min(1._r8, max(0._r8, mimicsplus_fchem_am))
              fphys_myc2 = min(1._r8, max(0._r8, mimicsplus_fphys_am))
              
              decomp_k(c,j,i_ecm_myc) = tau_myc * &
                    myc1_conc**((mimicsplus_densdep - 1.0_r8) * moist_mod) !w_d_o_scalars)     !decomposition rate for ecm !not sure
              pathfrac_decomp_cascade(c,j,i_myc1s1) = tau_myc1 !ECW should I add the enzyme flux here
              pathfrac_decomp_cascade(c,j,i_myc1s2) = fchem_myc1
              pathfrac_decomp_cascade(c,j,i_myc1s3) = fphys_myc1

              decomp_k(c,j,i_am_myc) = tau_myc * &
                    myc2_conc**((mimicsplus_densdep - 1.0_r8) * moist_mod) !w_d_o_scalars)     !decomposition rate for am
              pathfrac_decomp_cascade(c,j,i_myc2s1) = tau_myc2
              pathfrac_decomp_cascade(c,j,i_myc2s2) = fchem_myc2
              pathfrac_decomp_cascade(c,j,i_myc2s3) = fphys_myc2




              ! Same for cwd but only if fates not enabled; fates handles cwd on
              ! its own structure
              ! TODO This shows how BGC applies the spinup coefficients
              if (.not. use_fates) then
                  decomp_k(c,j,i_cwd) = k_frag * moist_mod !w_d_o_scalars  ! * spinup_geogterm_cwd(c)
               end if

               ! Tillage
               if (get_do_tillage()) then
                  call get_apply_tillage_multipliers(idop, c, j, decomp_k(c,j,:))
               end if
            end do
        end do

        deallocate(norm_froot_prof)

        ! pathfrac terms not calculated in the previous loop
        pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = 1.0_r8
        pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = 1.0_r8
        pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m1s3) = fphys_m1(bounds%begc:bounds%endc,1:nlevdecomp)
        pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m2s3) = fphys_m2(bounds%begc:bounds%endc,1:nlevdecomp)
        if (.not. use_fates) then
            pathfrac_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = 1.0_r8
            rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
        end if
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m1) = rf_l1m1
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l1m2) = rf_l1m2
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m1) = rf_l2m1
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_l2m2) = rf_l2m2
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1m1) = rf_s1m1
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s1m2) = rf_s1m2
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s2s1) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_s3s1) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m1s1) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m1s2) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m1s3) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m2s1) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m2s2) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_m2s3) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_myc1s1) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_myc1s2) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_myc1s3) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_myc2s1) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_myc2s2) = 0.0_r8
        rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_myc2s3) = 0.0_r8

      ! rf_decomp_cascade(bounds%begc:bounds%endc,1:nlevdecomp,i_myc1s1_e) = 0.0_r8

        !Track down order of things called in moudule that calls decomp_rate, when is what called (eg. FUN)

    end associate

end subroutine decomp_rates_mimicsplus



! subroutine decomp_rates_after_FUN (bounds, num_bgc_soilc, filter_bgc_soilc, &
!    num_soilp, filter_soilp, &
!    soilstate_inst, cnveg_carbonflux_inst, &
!    soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenflux_inst, &
!    soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst,& 
!    soilbiogeochem_state_inst &
!    )
! ! 
! ! DESCRIPTION
! !! add a new subroutine decomp_rates_after_FUN, will be called in FUN 
! ! col variables, do we have corresponding patch variables
! ! patch variables have to be p2c

! ! USES:
!    use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
!    use SoilBiogeochemNitrogenStateType  , only : soilbiogeochem_nitrogenstate_type
!    use subgridAveMod                    , only : p2c

! ! ARGUMENTS:
!    type(bounds_type)                    , intent(in)    :: bounds          
!    integer                              , intent(in)    :: num_soilp       ! number of soil patches in filter
!    integer                              , intent(in)    :: filter_soilp(:) ! filter for soil patches
!    integer                              , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
!    integer                              , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
!    type(soilstate_type)                 , intent(in)    :: soilstate_inst
!    type(cnveg_carbonflux_type)          , intent(in)    :: cnveg_carbonflux_inst
!    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
!    type(soilbiogeochem_carbonstate_type), intent(in)    :: soilbiogeochem_carbonstate_inst
!    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
!    type(soilbiogeochem_carbonstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
!    type(soilbiogeochem_state_type)      , intent(in)    :: soilbiogeochem_state_inst


!    !LOCAL VARIABLES:
!    integer   :: fn                    ! number of values in pft filter
!    integer   :: fp                    ! lake filter pft index 
!    integer   :: fc                    ! lake filter column index
!    integer   :: p, c                  ! pft index
!    integer   :: g, l                  ! indices
!    integer   :: j, i, k               ! soil/snow level index
!    real(r8), parameter :: eps = 1.e-6_r8
!    real(r8):: tau_myc ! Mycorrhizal turnover rate for EcM & AM 
!    real(r8):: myc1_conc     ! Carbon concentration in ectomycorrhiza pool
!    real(r8):: myc2_conc     ! Carbon concentration in abruscular mycorrhiza pool
!    real(r8):: avl_som_conc  ! Carbon concentration in avaliable SOM pool
!    real(r8):: chem_som_conc ! Carbon concentration in chemically protected SOM pool
!    real(r8):: phys_som_conc ! Carbon concentration in physically protected SOM pool
!    real(r8):: term_1  !
!    real(r8):: term_2  !
!    real(r8):: t_soi_degC

!    associate ( &
!       ! Carbon
!       decomp_cpools_vr => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col , &  ! Input: [real(r8) (:,:,:) ] (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) C pools
!       pathfrac_decomp_cascade => soilbiogeochem_carbonflux_inst%pathfrac_decomp_cascade_col, &  ! Output: [real(r8) (:,:,:) ]  what fraction of C leaving a given pool passes through a given transition (frac)
!       rf_decomp_cascade       => soilbiogeochem_carbonflux_inst%rf_decomp_cascade_col    , & ! Input:  [real(r8)          (:,:,:) ]  respired fraction in decomposition step (frac)
        
!       decomp_npools_vr      => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col         

!       )

!       ! prep
!       tau_myc = mimicsplus_k_myc_som / secsphr ! Turnover rate mycorrhiza with unit conversions (hourly -> second)
!       tau_myc1 = min(1._r8, max(0._r8, mimicsplus_tau_ecm))
!       fchem_myc1 = min(1._r8, max(0._r8, mimicsplus_fchem_ecm))
!       fphys_myc1 = min(1._r8, max(0._r8, mimicsplus_fphys_ecm))
!       tau_myc2 = min(1._r8, max(0._r8, mimicsplus_tau_am))
!       fchem_myc2 = min(1._r8, max(0._r8, mimicsplus_fchem_am))
!       fphys_myc2 = min(1._r8, max(0._r8, mimicsplus_fphys_am))
!       ! main column loop
!       do fc = 1,num_bgc_soilc
!          c = filter_bgc_soilc(fc)

!               ! Mycorrhizal concentration with necerssary unit conversions
!                         ! mgC/cm3 = gC/m3 * (1e3 mg/g) / (1e6 cm3/m3)
!          myc1_conc = (decomp_cpools_vr(c,j,i_ecm_myc) / col%dz(c,j)) * & !ECW is this what goes into ROI function (see end of module)
!          g_to_mg * cm3_to_m3
!          myc2_conc = (decomp_cpools_vr(c,j,i_am_myc) / col%dz(c,j)) * &
!          g_to_mg * cm3_to_m3

!          ! Soil organic matter concentration with necerssary unit conversions
!          ! mgC/cm3 = gC/m3 * (1e3 mg/g) / (1e6 cm3/m3)
!          avl_som_conc = (decomp_cpools_vr(c,j,i_avl_som) / col%dz(c,j)) * & !ECW This is for Baskaram N mining
!          g_to_mg * cm3_to_m3
!          chem_som_conc = (decomp_cpools_vr(c,j,i_chem_som) / col%dz(c,j)) * &
!          g_to_mg * cm3_to_m3
!          phys_som_conc = (decomp_cpools_vr(c,j,i_phys_som) / col%dz(c,j)) * &
!          g_to_mg * cm3_to_m3


!       ! FLUXES from SOMC 

!       ! FLUXES from SOMP

!       ! FLUXES from ECM
      
!       ! FLUXES from AM

!       ! CALCULATE MORTALITY FLUXES


!       ! CALCULATE FLUXES CARBON FLUXES

!       !
!       end do
!    end associate
! end subroutine decomp_rates_after_FUN

subroutine calc_myc_mortality(cpool_myc, npool_myc, m_fr,fc_myc2som, fn_myc2som)
   ! DESCRIPTION:
   ! Calculates mortality of mycorrhiza, based on: carbon content in mycorrhizal pool, mycorrhizal turnover rate and 
   ! fraction of necromass going from mycorrhizal pools into soil organic matter pools.
   
   ! USES

   ! ARGUMENTS
   real(r8), intent(in) :: cpool_myc     ! Carbon pool of mycorrhiza [gC/m3]
   real(r8), intent(in) :: npool_myc     ! Nitrogen pool of mycorrhiza [gN/m3]
   real(r8), intent(in) :: m_fr          ! mortality fraction
   real(r8), intent(inout) :: fc_myc2som ! carbon flux to SOM pool [gC/m3/s]
   real(r8), intent(inout) :: fn_myc2som ! nitrogen flux [gN/m3/s]

   ! LOCAL VARIABLES
   real(r8) :: secphr = 60._r8 * 60.0_r8
   real(r8), parameter :: small_flux = 1.e-14_r8

   fc_myc2som = cpool_myc * (params_inst%mimicsplus_k_myc_som / secphr ) * m_fr
   if (cpool_myc > small_flux) then 
      fn_myc2som = fc_myc2som * (npool_myc / cpool_myc)
   else
      fn_myc2som = 0.0_r8
      fc_myc2som = 0.0_r8
   endif

end subroutine calc_myc_mortality
   
   subroutine calc_myc_roi(cpool_myc, npool_myc,cpool_somp,cpool_soma,cpool_somc, &
                           npool_somp, npool_somc, sminn, myc_type, dz, big_roi, roi)
      !
      ! DESCRIPTION:
      ! Calculates ROI (Return Of Investment) of nitrogen per carbon invested from vegetation to mycorrhiza

      ! Contains:
      ! - Mycorrhizal modifer (controls nitrogen flux to vegetaion, based on carbon flux to mycorrhiza)
      ! - Connection to inorganic nitrogen cycle (mycorrhiza receive inorganic nitrogen)
      ! 
      ! !USES:
      use clm_varcon       , only : secspday, secsphr, tfrz, spval
            
      !
      ! !ARGUMENTS:
     real(r8), intent(in) :: cpool_myc    ! Carbon pool of mycorrhiza [gC/m3] ()
     real(r8), intent(in) :: npool_myc    ! Nitrogen pool of mycorrhiza [gN/m3] ()
     real(r8), intent(in) :: cpool_somp   ! physically protected SOM pool [gC/m3]
     real(r8), intent(in) :: cpool_soma   ! available SOM pool [gC/m3]
     real(r8), intent(in) :: cpool_somc   ! chemically protected SOM pool [gC/m3] 
     real(r8), intent(in) :: npool_somp   ! physically protected SOM pool [gN/m3]
     real(r8), intent(in) :: npool_somc   ! chemically protected SOM pool [gN/m3]]
     real(r8), intent(in) :: sminn        ! soil mineral nitrogen (NO3+NH4) [gN/m3/s]
     integer, intent(in)  :: myc_type     ! type of mycorrhiza EcM=1, AM=2
     real(r8), intent(in) :: dz           ! layer thickness [m]
     real(r8), intent(in) :: big_roi      ! large number
     real(r8), intent(inout) :: roi       ! RoI nitrogen per carbon invested from vegetation to mycorrhiza [gN/gC]
     ! Local variables
     !real(r8) :: turnover, r_myc, N_inorg, N_inorg_myc, C_mining_s3, N_mining_s3, C_mining_s2, N_mining_s2, N_acquired
     
     real(r8) :: fc_myc2somp,fc_myc2soma,fc_myc2somc  ! carbon fluxes myc to som (necrormass) [gC/m2/s]
     real(r8) :: fc_somp2soma,fc_somc2soma            ! carbon fluxes som to som due to mining [gC/m2/s]
     real(r8) :: fn_myc2somp,fn_myc2soma,fn_myc2somc  ! nitrogen fluxes myc to som (necrormass) [gN/m2/s]
     real(r8) :: fn_mining_somc,fn_mining_somp        ! nitrogen fluxes to som to myc (mining + scavenging) [gN/m2/s]
     real(r8) :: fn_smin2myc                          ! nitrogen flux from mineral soil to myc [gN/m3/s2]

     !real(r8)            :: r_myc ! mycorrhizal modifier (constrains mining and scavaging) [-]
     real(r8), parameter :: small_flux = 1.e-14_r8
     real(r8), parameter :: eps = 0.5

     fn_smin2myc = (params_inst%mimicsplus_vmax_myc / secsphr) * sminn  * &
                    (cpool_myc / (cpool_myc + params_inst%mimicsplus_k_m_emyc / dz)) !* r_myc
     ! Initialize mining rates
     fc_somp2soma = 0.0_r8
     fc_somc2soma = 0.0_r8
     fn_mining_somp = 0.0_r8
     fn_mining_somc = 0.0_r8
     if (myc_type == 1) then
      ! check units
      call calc_myc_mining_rates(dz, cpool_somp,cpool_myc, npool_myc,fc_somp2soma,fn_mining_somp) !,r_myc)
      call calc_myc_mining_rates(dz, cpool_somc,cpool_myc, npool_myc,fc_somc2soma,fn_mining_somc) !,r_myc)
     endif
     
 
     ! Calculate ROI
     if (params_inst%mimicsplus_k_myc_som < small_flux) then
         call endrun(msg='k_myc_som in the parameter file is too small')
     endif
     if (myc_type == 1) then 
      if (cpool_myc > 0.0_r8) then
      roi = (fn_smin2myc + fn_mining_somc + fn_mining_somp) / (params_inst%mimicsplus_k_myc_som / secsphr ) * &
             params_inst%mimicsplus_mge_ecm / cpool_myc
             if (roi == 0.0_r8) then
               roi = 1.0_r8/big_roi
             end if 
      else
         roi = 1.0_r8/big_roi
      endif
     else
      if (cpool_myc > 0.0_r8) then
         roi = fn_smin2myc/ (params_inst%mimicsplus_k_myc_som / secsphr ) * &
                params_inst%mimicsplus_mge_am / cpool_myc
         if (roi == 0.0_r8) then
            roi = 1.0_r8/big_roi
         end if 
      else
         roi = big_roi
      endif

     endif
 
   end subroutine calc_myc_roi


   subroutine cost_FUN(fc_veg2myc,cpool_myc, npool_myc,cpool_somp,cpool_soma,cpool_somc, &
      npool_somp, npool_somc, sminn, sminfrc_no3, myc_type, dz, big_cost, cost_myc_no3, cost_myc_nh4)
   ! DESCRIPTION:
   ! Calculates the cost in gC (from vegetation) per gN (from mycorrhiza).
   ! Mycorrhizal pools have ac constant C:N ratio, these are fullfilled by uptake fluxes from inorganic nitrogen (to EcM & AM), and mining (to EcM). 
   ! The Nitrogen that is left after C:N ratio is fullfilled can go to vegetation (uptake - demand equation)
   ! Under limited nitrogen soil ocnditions, the carbon use efficiency can be lowered
   ! !ARGUMENTS:
      real(r8), intent(in) :: fc_veg2myc   ! Carbon flux from plant to mycorrhiza
      real(r8), intent(in) :: cpool_myc    ! Carbon pool of mycorrhiza [gC/m3] ()
      real(r8), intent(in) :: npool_myc    ! Nitrogen pool of mycorrhiza [gN/m3] ()
      real(r8), intent(in) :: cpool_somp   ! physically protected SOM pool [gC/m3]
      real(r8), intent(in) :: cpool_soma   ! available SOM pool [gC/m3]
      real(r8), intent(in) :: cpool_somc   ! chemically protected SOM pool [gC/m3] 
      real(r8), intent(in) :: npool_somp   ! physically protected SOM pool [gN/m3]
      real(r8), intent(in) :: npool_somc   ! chemically protected SOM pool [gN/m3]]
      real(r8), intent(in) :: sminn        ! soil mineral nitrogen (NO3+NH4) [gN/m2]
      real(r8), intent(in) :: sminfrc_no3  ! fraction of soil mineral nitrogen NO3
      integer, intent(in)  :: myc_type     ! type of mycorrhiza EcM=1, AM=2
      real(r8), intent(in) :: dz           ! layer thickness [m]
      real(r8), intent(in) :: big_cost     ! large cost
      real(r8), intent(out) :: cost_myc_no3    ! cost function for mycorrhiza in FUN UNITS!!!!
      real(r8), intent(out) :: cost_myc_nh4    ! cost function for mycorrhiza in FUN UNITS!!!!

   !LOCAL VARIABLES
   real(r8), parameter :: secphr = 60.0_r8 * 60.0_r8
   real(r8), parameter :: f_enz = 0.1_r8                ! [-]Fraction of C from vegetation to EcM, that goes into SOMa for mining
   real(r8), parameter :: f_growth = 0.5_r8            ! [-] Fraction of mycorrhizal N uptake that needs to stay within the fungi (not given to plant)
   real(r8), parameter :: small_Value = 1.0e-6_r8

   real(r8) :: fn_myc2veg                           ! nitrogen fluxes mycorrhiza to vegetation
   real(r8) :: fn_smin2myc                          ! nitrogen flux from mineral soil to myc [gN/m2/s]
   real(r8) :: fc_somp2soma,fc_somc2soma            ! carbon fluxes som to som due to mining [gC/m2/s]
   real(r8) :: fn_mining_somc,fn_mining_somp        ! nitrogen fluxes to som to myc (mining + scavenging) [gN/m2/s]

   real(r8) ::     N_demand_myc   ![gN/m3 h]
   real(r8) ::     N_uptake_myc   ![gN/m3 h]
   real(r8) ::     c_use_eff      ! carbon use efficiency [-]

   if (myc_type == 1) then
      c_use_eff = params_inst%mimicsplus_mge_ecm
   else
      c_use_eff = params_inst%mimicsplus_mge_am
   endif

   fn_smin2myc = (params_inst%mimicsplus_vmax_myc / secphr) * sminn / dz * &
                 (cpool_myc / (cpool_myc + params_inst%mimicsplus_k_m_emyc)) !ECW double check thats ecphr is correct here

   ! Initialize mining rates
   fc_somp2soma = 0.0_r8
   fc_somc2soma = 0.0_r8
   fn_mining_somp = 0.0_r8
   fn_mining_somc = 0.0_r8
   if (myc_type == 1) then
    call calc_myc_mining_rates(dz, cpool_somp,cpool_myc, npool_myc,fc_somp2soma,fn_mining_somp)
    call calc_myc_mining_rates(dz, cpool_somc,cpool_myc, npool_myc,fc_somc2soma,fn_mining_somc)
   endif


   ! Nitrogen flux from mycorrhiza to vegetation
   if (myc_type == 1) then 
      !fn_myc_veg = (fn_smin2myc + fn_somc2myc_ecm + fn_somp2myc_ecm) - c_use_eff * fc_veg2myc * (1.0_r8 - f_enz) / params_inst%mimicsplus_cn_myc
      N_uptake_myc = fn_smin2myc + fn_mining_somc + fn_mining_somp
      N_demand_myc = c_use_eff * fc_veg2myc * (1.0_r8 - f_enz) / params_inst%mimicsplus_cn_myc
      if (N_uptake_myc > N_demand_myc) then
         fn_myc2veg = N_uptake_myc - N_demand_myc
      else
         fn_myc2veg = (1-f_growth) * N_uptake_myc
         ! here should CUE go, but we will do this in the update phase
      end if
   else
         N_uptake_myc = fn_smin2myc
         N_demand_myc = c_use_eff * fc_veg2myc / params_inst%mimicsplus_cn_myc
         if (N_uptake_myc > N_demand_myc) then
            fn_myc2veg = N_uptake_myc - N_demand_myc
         else
            fn_myc2veg = (1-f_growth) * N_uptake_myc
            ! here should CUE go, but we will do this in the update phase
         end if
   end if

   if (fn_myc2veg > small_Value) then
      if (sminfrc_no3 > small_Value) then
         cost_myc_no3 = fc_veg2myc / (fn_myc2veg) / sminfrc_no3
         if ( (1.0_r8 - sminfrc_no3)> small_Value) then
            cost_myc_nh4 = fc_veg2myc / (fn_myc2veg) / (1.0_r8 - sminfrc_no3)
         else
            cost_myc_nh4 = big_cost
         end if
      else
         cost_myc_no3 = big_cost
         cost_myc_nh4 = fc_veg2myc / (fn_myc2veg) / (1.0_r8 - sminfrc_no3)
      endif
   else 
      cost_myc_no3 = big_cost
      cost_myc_nh4 = big_cost
   end if
 
   
   end subroutine cost_FUN



  subroutine calc_myc_mining_rates(dz, cpool_som,cpool_myc, npool_myc, fc_som2soma,fn_mining_som)

   ! DESCRIPTION:
   ! Calculates mining of ectomycorrhizal fungi for nitrogen in soil organic matter pools.
   ! (only!) Ectomycorrhizal fungi can a allocate fraction of the incoming carbon from vegetation (not in this subroutine) to the avaliable SOM pool.
   ! They mine for nitrogen in the chemmically and physically protected SOM pools and create nitrogen fluxes from theses to EcM.
   ! During this process carbon is release from the chemmically and physically protected SOM pools which enters the avaliable SOM pool.

   ! USES:

   ! ARGUMENTS:
   real(r8), intent(in) :: cpool_som          ! SOM pool [gC/m3]
   real(r8), intent(in) :: cpool_myc          ! Carbon pool of mycorrhiza [gC/m3]
   real(r8), intent(in) :: npool_myc          ! Nitrogen pool of mycorrhiza [gC/m3]
   real(r8), intent(in) :: dz                 ! layer thickness [m]
   real(r8), intent(inout) :: fc_som2soma     ! carbon flux to available SOM pool [gC/m3/s]
   real(r8), intent(inout) :: fn_mining_som   ! nitrogen mining flux [gN/m3/s]

   ! LOCAL VARIABLES:
   real(r8) :: secphr = 60.0_r8 * 60.0_r8
   real(r8), parameter :: small_flux = 1.e-14_r8
   
   ! SOM carbon flux
   fc_som2soma = (params_inst%mimicsplus_k_mo / secphr) * dz * cpool_myc * cpool_som
   ! Nitrogen mining flux
   if (cpool_myc > small_flux) then
     fn_mining_som = fc_som2soma * (npool_myc / cpool_myc )
   else 
     fn_mining_som = 0.0_r8
     fc_som2soma = 0.0_r8
   endif

  end subroutine calc_myc_mining_rates

  subroutine fun_fluxes_myc_update1 (cpool_myc, npool_myc, cpool_somp, cpool_soma,cpool_somc, &
                                    npool_somp, npool_somc, sminn, sminfrc_no3, myc_type, dz, &
                                    fc_veg2myc_no3, fc_veg2myc_nh4, &
                                    no3_from_nonmyc,nh4_from_nonmyc, &
                                    no3_from_myc, nh4_from_myc, &
                                    fn_smin_no3_2myc, fn_smin_nh4_2myc, c_myc_resp,c_myc_growth,n_myc_growth,c_ecm_enz, &
                                    fc_somp2soma, fc_somc2soma, fn_mining_somp,fn_mining_somc)

   ! ! DESCRIPTION:
   ! Updating Nitrogen and Carbon fluxes into mycorrhizal pools, to let them grow
   !
   ! ! USES:
   use clm_time_manager, only: get_step_size_real

   !
   ! !ARGUMENTS:

   real(r8), intent(in) :: cpool_myc    ! Carbon pool of mycorrhiza            [gC/m3]
   real(r8), intent(in) :: npool_myc    ! Nitrogen pool of mycorrhiza          [gN/m3]
   real(r8), intent(in) :: cpool_somp   ! physically protected SOM pool        [gC/m3]
   real(r8), intent(in) :: cpool_soma   ! available SOM pool                   [gC/m3]
   real(r8), intent(in) :: cpool_somc   ! chemically protected SOM pool        [gC/m3] 
   real(r8), intent(in) :: npool_somp   ! physically protected SOM pool        [gN/m3]
   real(r8), intent(in) :: npool_somc   ! chemically protected SOM pool        [gN/m3]
   real(r8), intent(in) :: sminn        ! soil mineral nitrogen (NO3+NH4)      [gN/m2]
   real(r8), intent(in) :: sminfrc_no3  ! fraction of soil mineral nitrogen NO3 [-]
   integer,  intent(in) :: myc_type     ! type of mycorrhiza EcM=1, AM=2        [-]
   real(r8), intent(in) :: dz           ! layer thickness                       [m]
 
      ! FUN fluxes 
   real(r8), intent(inout) :: fc_veg2myc_no3           ! Carbon flux from plant to mycorrhiza                        [gC/m3] 
   real(r8), intent(inout) :: fc_veg2myc_nh4           ! Carbon flux from plant to mycorrhiza                        [gC/m2]
   real(r8), intent(inout) :: no3_from_nonmyc          ! nitrogen flux from non mycorrhizal pool to vegetation (NO3) [gN/m3]
   real(r8), intent(inout) :: nh4_from_nonmyc          ! nitrogen flux from non mycorrhizal pool to vegetation (NH4) [gN/m3]
   real(r8), intent(inout) :: no3_from_myc             ! nitrogen flux from mycorrhizal pool to vegetation (NO3)     [gN/m3]
   real(r8), intent(inout) :: nh4_from_myc             ! nitrogen flux from mycorrhizal pool to vegetation (NH4)     [gN/m3]

   real(r8), intent(inout) :: fn_smin_no3_2myc         ! nitrogen flux from inorganic NO3 pool to mycorrhiza [gN/m2]
   real(r8), intent(inout) :: fn_smin_nh4_2myc         ! nitrogen flux from inorganic NH4 pool to mycorrhiza [gN/m2]

     ! Mycorrhiza internal fluxes
   real(r8), intent(inout) :: c_myc_resp               ! carbon respiration flux for AM mycorrhiza     [gC/m3/s]
   real(r8), intent(inout) :: c_myc_growth             ! carbon growth flux for AM mycorrhiza          [gC/m3/s]
   real(r8), intent(inout) :: n_myc_growth             ! nitrogen growth flux for AM mycorrhiza        [gN/m3/s]
   real(r8), intent(inout), optional :: c_ecm_enz      ! carbon enzyme production flux for ECM mycorrhiza [gC/m3/s]

   ! Mining fluxes
   real(r8), intent(inout), optional :: fc_somp2soma    ! carbon fluxes somp to soma due to mining [gC/m3/s]
   real(r8), intent(inout), optional :: fc_somc2soma    ! carbon fluxes somc to soma due to mining [gC/m3/s]
   real(r8) ,intent(inout), optional :: fn_mining_somp  ! nitrogen fluxes to somp to myc (mining + scavenging) [gN/m3/s]
   real(r8) ,intent(inout), optional :: fn_mining_somc  ! nitrogen fluxes to somc to myc (mining + scavenging) [gN/m3/s]

   ! LOCAL VARIABLES
   real(r8) :: secphr = 60.0_r8 * 60.0_r8
   real(r8)            :: dt                            ! timestep size (seconds)
   real(r8), parameter :: f_enz = 0.1_r8                ! [-]Fraction of C from vegetation to EcM, that goes into SOMa for mining
   real(r8), parameter :: f_growth = 0.5_r8             ! [-] Fraction of mycorrhizal N uptake that needs to stay within the fungi (not given to plant)
   real(r8) :: fc_veg2myc                               ! carbon flux vegetation to mycorrhiza []
   real(r8) :: fn_myc2veg                               ! nitrogen fluxes mycorrhiza to vegetation  [gN/m3/s]
   real(r8) :: l_c_ecm_enz, l_fc_somp2soma, l_fc_somc2soma, l_fn_mining_somp, l_fn_mining_somc !locals of optional args
   real(r8) :: N_uptake_myc                             ! [gN/m3/s]
   real(r8) :: N_demand_myc                             ! [gN/m3/s]
   real(r8) :: c_use_eff                                ! carbon use efficiency [-]
   real(r8) :: frac_no3_myc, frac_no3_nonmyc            ! [-]
   real(r8) :: fn_smin2myc, fno3_myc2veg,fnh4_myc2veg   ! local fluxes for balancing [gN/m3/s]
   real(r8) :: no3_unpaid, nh4_unpaid                   ! nitrogen that did not actually reach the plant [gN/m2]
   real(r8) :: smin_overflow                            ! How much N mycorrhiza and non mycorrhiza promised to take out of soil and give to plant
   real(r8) :: l_no3_frac
   real(r8), parameter :: small_value =1.0e-7_r8 

   character(len=256) :: message

   dt           = get_step_size_real()

   ! If statements check, if mining is happening. Mining includes: C Enzymes production flux from EcM, 
   ! C fluxes from SOMp and c to SOMa, and N fluxes from SOMp and c to EcM. 
   ! If no mining occurs eg. only AM is present those fluxes are set to 0.
   if (present(c_ecm_enz)) then
      l_c_ecm_enz = c_ecm_enz
   else
      l_c_ecm_enz = 0.0_r8
   end if

   if (present(fc_somp2soma)) then
      l_fc_somp2soma = fc_somp2soma
   else
      l_fc_somp2soma = 0.0_r8
   end if

   if (present(fc_somc2soma)) then
      l_fc_somc2soma = fc_somc2soma
   else
      l_fc_somc2soma = 0.0_r8
   end if
 
   if (present(fn_mining_somp)) then
      l_fn_mining_somp = fn_mining_somp
   else
      l_fn_mining_somp = 0.0_r8
   end if

   if (present(fn_mining_somc)) then
      l_fn_mining_somc = fn_mining_somc
   else
      l_fn_mining_somc = 0.0_r8
   end if

   ! Initialize carbon use efficiencies
   if (myc_type == 1) then
      c_use_eff = params_inst%mimicsplus_mge_ecm
   else
      c_use_eff = params_inst%mimicsplus_mge_am
   endif

   ! Get no3/total N fraction
   if (no3_from_myc + nh4_from_myc > 0.0_r8) then
      frac_no3_myc = no3_from_myc / (no3_from_myc + nh4_from_myc)
   else
      frac_no3_myc = 1.0_r8
   end if
   if (no3_from_nonmyc + nh4_from_nonmyc > 0.0_r8) then
      frac_no3_nonmyc = no3_from_nonmyc / (no3_from_nonmyc + nh4_from_nonmyc)
   else
      frac_no3_nonmyc = 1.0_r8
   end if

   ! Soil mineral nitrogen flux to mycorrhiza, split into NO3 and NH4
   fn_smin2myc = (params_inst%mimicsplus_vmax_myc / secphr) * sminn * &  
   (cpool_myc / (cpool_myc + params_inst%mimicsplus_k_m_emyc)) * dt ! multiplied by dt to get it in same units as no3_from_nonmyc
   fn_smin_no3_2myc = fn_smin2myc * sminfrc_no3
   fn_smin_nh4_2myc = fn_smin2myc * (1.0_r8 - sminfrc_no3)
   if (sminn * sminfrc_no3 <= small_value) then
      fn_smin_no3_2myc = 0.0_r8
      no3_from_nonmyc = 0.0_r8
   else
     ! If there is not enough nitrogen in the soil mineral, NO3 and NH4
      if (fn_smin_no3_2myc + no3_from_nonmyc > sminn * sminfrc_no3) then
         smin_overflow =  (fn_smin_no3_2myc + no3_from_nonmyc) -  sminn * sminfrc_no3! NO3 that is actually in soil - 
         fn_smin_no3_2myc = fn_smin_no3_2myc - smin_overflow * fn_smin_no3_2myc / (sminn * sminfrc_no3 + smin_overflow ) ! mycorrhiza limit uptake from inorganic N pools
         no3_from_nonmyc = no3_from_nonmyc - smin_overflow * no3_from_nonmyc / (sminn * sminfrc_no3 + smin_overflow)    ! non mycorrhiza limit uptake from inorganic N pools
         if (& 
            no3_from_nonmyc < 0.0_r8 .or. fn_smin_no3_2myc < 0.0_r8) then 
            write(iulog,*) 'ERROR: type,myc, nonmyc, sminno3', myc_type, fn_smin_no3_2myc, no3_from_nonmyc, sminn *sminfrc_no3 - smin_overflow, sminfrc_no3 * sminn, fn_smin2myc
            call endrun(                             msg= errMsg(sourcefile,  __LINE__))
         end if
      endif
    end if 
    if (sminn *(1.0_r8 - sminfrc_no3) <= small_value) then
      nh4_from_nonmyc = 0.0_r8
      fn_smin_nh4_2myc = 0.0_r8
   else
      if (fn_smin_nh4_2myc + nh4_from_nonmyc > sminn * (1.0_r8 - sminfrc_no3)) then
         smin_overflow = (fn_smin_nh4_2myc + nh4_from_nonmyc) - sminn * (1.0_r8 - sminfrc_no3)
         fn_smin_nh4_2myc = fn_smin_nh4_2myc - smin_overflow * fn_smin_nh4_2myc / (sminn * (1.0_r8 - sminfrc_no3) + smin_overflow)
         nh4_from_nonmyc = nh4_from_nonmyc - smin_overflow * nh4_from_nonmyc / (sminn * (1.0_r8 - sminfrc_no3) + smin_overflow)
         if ( & 
            nh4_from_nonmyc < 0.0_r8 .or. fn_smin_nh4_2myc < 0.0_r8) then 
            write(iulog,*) 'ERROR: type,myc, nonmyc, sminnh4', myc_type, fn_smin_nh4_2myc, nh4_from_nonmyc, sminn * (1.0_r8 - sminfrc_no3) - smin_overflow, (1.0_r8 - sminfrc_no3) * sminn, fn_smin2myc
            call endrun(                             msg= errMsg(sourcefile,  __LINE__))
         end if
      endif
   endif



   ! use new fraction and sminn if they got limited
   fn_smin2myc = (fn_smin_no3_2myc + fn_smin_nh4_2myc) /dz /dt ! put it back to gN/m3/s to be consistent with mining fluxes
   if(fn_smin2myc > 0.0_r8) then
      l_no3_frac = fn_smin_no3_2myc / fn_smin2myc
   else
      l_no3_frac = 1.0_r8
   endif

   fc_veg2myc=fc_veg2myc_no3 + fc_veg2myc_nh4
   if (myc_type == 1) then                                                ! EcM
      call calc_myc_mining_rates(dz, cpool_somp,cpool_myc, npool_myc,l_fc_somp2soma,l_fn_mining_somp)
      call calc_myc_mining_rates(dz, cpool_somc,cpool_myc, npool_myc,l_fc_somc2soma,l_fn_mining_somc)
      N_demand_myc = c_use_eff * (1 - f_enz) * (fc_veg2myc_no3 + fc_veg2myc_nh4) / params_inst%mimicsplus_cn_myc / dz / dt
      N_uptake_myc = fn_smin2myc + l_fn_mining_somp + l_fn_mining_somc 
      if (N_uptake_myc > N_demand_myc) then
         fn_myc2veg = N_uptake_myc - N_demand_myc                        ! N flux myc -> veg
         n_myc_growth = N_demand_myc                                     ! How much N the need to grow
         l_c_ecm_enz  = fc_veg2myc * f_enz
         ! we add the enzymes here into the EcM pool because we first took them away, but now add them again
         c_myc_growth = n_myc_growth * params_inst%mimicsplus_cn_myc + l_c_ecm_enz
         ! enzyme flux will go to soma pool in the next update routine
         c_myc_resp  = fc_veg2myc - c_myc_growth                          ! C that they don't need to grow
      else ! less N in soil, so we limit N flux to vegetaion and mycorrhiza N demand so their sum is equal to N uptake
         fn_myc2veg = (1-f_growth) * N_uptake_myc
         n_myc_growth = f_growth * N_uptake_myc
         l_c_ecm_enz  = fc_veg2myc * f_enz
         ! we add the enzymes here into the EcM pool because we first took them away, but now add them again
         c_myc_growth = n_myc_growth * params_inst%mimicsplus_cn_myc + l_c_ecm_enz
         c_myc_resp  = fc_veg2myc - c_myc_growth                          ! C that they don't need to grow
      end if

      ! update inout optional arguments
      if (present(c_ecm_enz)) then
         c_ecm_enz = l_c_ecm_enz
      endif
      if (present(fc_somp2soma)) then
         fc_somp2soma = l_fc_somp2soma
      end if
      if (present(fc_somc2soma)) then
         fc_somc2soma = l_fc_somc2soma
      end if
      if (present(fn_mining_somp)) then
         fn_mining_somp = l_fn_mining_somp
      end if
      if (present(fn_mining_somc)) then
         fn_mining_somc = l_fn_mining_somc
      end if

   else                                                                   ! AM 
      N_demand_myc = c_use_eff * (fc_veg2myc_no3 + fc_veg2myc_nh4) / params_inst%mimicsplus_cn_myc / dz / dt
      N_uptake_myc = fn_smin2myc
      if (N_uptake_myc > N_demand_myc) then
         fn_myc2veg = N_uptake_myc - N_demand_myc                        ! N flux myc -> veg
         n_myc_growth = N_demand_myc                                     ! How much N the need to grow
         c_myc_growth = n_myc_growth * params_inst%mimicsplus_cn_myc     ! How much C they need to grow
         c_myc_resp  = fc_veg2myc - c_myc_growth                         ! C that they don't need to grow

      else
         fn_myc2veg = (1-f_growth) * N_uptake_myc
         n_myc_growth = f_growth * N_uptake_myc 
         c_myc_growth = n_myc_growth * params_inst%mimicsplus_cn_myc / fc_veg2myc
         c_myc_resp = fc_veg2myc - c_myc_growth
      end if
   endif

   fno3_myc2veg = fn_myc2veg * l_no3_frac ! Amount of NO3 / NH4 mycorrhiza could give to plant 
   fnh4_myc2veg = fn_myc2veg * (1 - l_no3_frac)
   no3_unpaid = fno3_myc2veg * dz * dt - no3_from_myc ! plants did not actually get all the promised nitrogen (NO3 / NH4)
   nh4_unpaid = fnh4_myc2veg * dz * dt - nh4_from_myc
   !if (no3_unpaid + nh4_unpaid < 0) then
   !   message = "N flux from AM mycorrhiza to plant is greater than mycorrhiza can actually give to the plant"
   !   call endrun(msg=message)
   !else 
   ! Bringing nitrogen fluxes to vegetation back to gN/m2
      no3_from_myc = fno3_myc2veg * dz * dt   ! if we can give more nitrogen than we promised just give it
      nh4_from_myc = fnh4_myc2veg * dz * dt 
   !endif

  end subroutine fun_fluxes_myc_update1
  
  
          
  !Moisture function, based on testbed code: https://github.com/wwieder/biogeochem_testbed/blob/957a5c634b9f2d0b4cdba0faa06b5a91216ace33/SOURCE_CODE/mimics_cycle.f90#L401-L419
  real(r8) function r_moist(h2osoi_liq,watsat, h2osoi_ice, dz) !As in testbed (and CLM) version of MIMICS            
   
   !  !USES
   use clm_varcon, only: denh2o, denice
   
   ! !ARGUMENTS:
   real(r8), intent(in)    :: h2osoi_liq ! liquid water content kg/m2
   real(r8), intent(in)    :: watsat ! porosity m3/m3
   real(r8), intent(in)    :: h2osoi_ice ! ice content kg/m2 
   real(r8), intent(in)    :: dz ! soil layer thickness 

   ! !LOCAL VARIABLES
   real(r8):: wliq !water liquid
   real(r8):: wice !water ice


   !NOTE: This moisture function represent both inhibition bc. very dry conditions, and very wet (anaerobic) conditions. 
   !FROM mimics_cycle.f90 in testbed:
   ! ! Read in soil moisture data as in CORPSE
   !  theta_liq  = min(1.0, casamet%moistavg(npt)/soil%ssat(npt))     ! fraction of liquid water-filled pore space (0.0 - 1.0)
   !  theta_frzn = min(1.0, casamet%frznmoistavg(npt)/soil%ssat(npt)) ! fraction of frozen water-filled pore space (0.0 - 1.0)
   !  air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
   !
   !  if (mimicsbiome%fWFunction .eq. CORPSE) then
   !    ! CORPSE water scalar, adjusted to give maximum values of 1
   !    fW = (theta_liq**3 * air_filled_porosity**2.5)/0.022600567942709
   !    fW = max(0.05, fW)
   wliq = h2osoi_liq / dz * denh2o
   wice = h2osoi_ice / dz * denice
   wliq  = min(1.0_r8, wliq/watsat)     ! fraction of liquid water-filled pore space (0.0 - 1.0)
   wice = min(1.0_r8, wice/watsat)     ! fraction of frozen water-filled pore space (0.0 - 1.0)
   !ECW check equations, find out how theta_l & theta_f are called in the rest of CTSM 

   r_moist = ((wliq**3)*max(0.0_r8, 1.0_r8-wliq-wice)**2.5_r8)/0.022600567942709_r8
   r_moist = max(0.05, r_moist) !ECW This is probably what I will replace w_d_o_scalar with
  end function r_moist

end module SoilBiogeochemDecompCascadeMIMICSplusMod

!Problems for future me:
! - Elin had a function for mycorrhizal mortality, didn't includ ethat yet bc I think it's the same as the parameter



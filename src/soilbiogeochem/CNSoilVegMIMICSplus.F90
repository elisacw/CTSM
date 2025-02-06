module CNSoilVegMIMICSplus

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module connects the Soil decomposition module MIMICS+ (Aas et al. 2023) with the vegetation through the 
  ! symbiosis between mycorrhizal fungi and plants.
  ! Coupling follows Sulman et al. (2019)
  
  ! !USES:
  use clm_varpar                          , only : nlevdecomp
  use clm_varpar                          , only : i_litr_min, i_litr_max, i_cwd
  use clm_varpar                          , only : i_met_lit, i_str_lit, i_phys_som, i_chem_som
  use clm_time_manager                    , only : get_average_days_per_year, get_step_size
  use clm_varcon                          , only : secspday, secsphr, tfrz, spval
  use SoilBiogeochemDecompCascadeConType  , only : mimicsplus_decomp, decomp_method
  use decompMod                           , only : bounds_type
  use PatchType                           , only : patch
  use CNVegstateType                      , only : cnveg_state_type
  use CNVegCarbonStateType                , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType                 , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType              , only : cnveg_nitrogenstate_type
  use CNVegnitrogenfluxType               , only : cnveg_nitrogenflux_type
  use SoilBiogeochemNitrogenFluxType      , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType     , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemCarbonFluxType        , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType       , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemStateType             , only : soilbiogeochem_state_type
  use SoilBiogeochemDecompCascadeMIMICSMod, only : sulman_params => params_inst 
  use WaterStateType                      , only : waterstate_type
  use SoilStateType                       , only : soilstate_type

  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: CN_soil_veg_exchange         !
  private :: active_root_N_uptake         !
  private :: myc_scavenger_N_uptake       !
  private :: myc_miner_N_uptake           !
  private :: miner_decomposition          !

  ! !FUNCTIONS:
  private :: theta_func                   ! Functional dependence on soil moisture / air 
  private :: resp_aerobic                 !
  private :: resp_denitrif                !
  private :: resp_myc                     ! Respiration of mycorrhiza
  private :: max_immobilization_rate      ! 
  private :: Vmax                         !
  private :: Vmax_denitrif                !
  private :: Vmax_myc                     ! Michaelis Menten Kinetics for mycorrhiza

 
  type, public :: symbiont_type

  !real(r8), pointer, private :: C_biomass             (:) ! Carbon biomass
  !real(r8), pointer, private :: N_biomass             (:) ! Nitrogen biomass
  !real(r8), pointer, private :: C_inter_biomass       (:) ! Carbon intermediate pool biomass
  !real(r8), pointer, private :: N_inter_biomass       (:) ! Nitrogen intermediate pool biomass
  
   contains

     ! Public procedures
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: InitCold

     ! Private procedures
     procedure, private :: InitAllocate
     procedure, private :: InitHistory

  end type symbiont_type


  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(symbiont_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
   !
   ! !ARGUMENTS:
   class(symbiont_type) :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:
   integer :: begp, endp
   integer :: begc, endc
   !------------------------------------------------------------------------

   begp = bounds%begp; endp= bounds%endp
   begc = bounds%begc; endc= bounds%endc

   ! allocate(this%ac_phs_patch      (begp:endp,2,1:nlevcan)) ; this%ac_phs_patch      (:,:,:) = nan

   end subroutine InitAllocate

   
   !-----------------------------------------------------------------------
   subroutine InitHistory(this, bounds)
      !
      ! !USES:
      use histFileMod   , only: hist_addfld1d, hist_addfld2d
      !
      ! !ARGUMENTS:
      class(symbiont_type) :: this
      type(bounds_type), intent(in) :: bounds
      real(r8), pointer  :: ptr_1d(:)  ! pointer to 1d patch array
      !
      ! !LOCAL VARIABLES:
      integer :: begp, endp
      !---------------------------------------------------------------------
  
      begp = bounds%begp; endp= bounds%endp
  
      !this%rh_leaf_patch(begp:endp) = spval
      !call hist_addfld1d (fname='RH_LEAF', units='fraction', &
      !   avgflag='A', long_name='fractional humidity at leaf surface', &
      !   ptr_patch=this%rh_leaf_patch, set_spec=spval, default='inactive')
   
   end subroutine InitHistory

   !------------------------------------------------------------------------
   subroutine Restart(this, bounds, ncid, flag)
     !
     ! !USES:
     use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
     use restUtilMod
     !
     ! !ARGUMENTS:
     class(symbiont_type) :: this
     type(bounds_type), intent(in)    :: bounds
     type(file_desc_t), intent(inout) :: ncid   ! netcdf id
     character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
     !
     ! !LOCAL VARIABLES:
     integer :: j,c ! indices
     logical :: readvar      ! determine if variable is on initial file
     !-----------------------------------------------------------------------

   end subroutine Restart

   !-----------------------------------------------------------------------

   subroutine InitCold(this, bounds)
      !
      ! !ARGUMENTS:
      class(symbiont_type)          :: this
      type(bounds_type), intent(in) :: bounds
      !
      ! !LOCAL VARIABLES:
      integer :: p,l                        ! indices
      !-----------------------------------------------------------------------
  
    end subroutine InitCold

   !-----------------------------------------------------------------------

   subroutine CN_soil_veg_exchange (bounds, cnveg_nitrogenstate_inst, leaf_prof_patch, froot_prof_patch, croot_prof_patch, soilbiogeochem_nitrogenflux_inst)

   ! !DESCRIPTION:
   
   ! Step 1:
   ! Calculate Nstress | from CTSM?
   ! Plant tissues have a fixed C:N ratio, but vary by PFT
   ! If N is limited = biomass is limited to avaliable N and C is left in pool


   ! !USES:
      use clm_varcon, only: pct_to_frac

   
   !
   ! !ARGUMENTS:
   type(bounds_type)                       , intent(in)    :: bounds 
   type(cnveg_nitrogenstate_type)          , intent(in)    :: cnveg_nitrogenstate_inst
   type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
   
   real(r8)                             , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: croot_prof_patch(bounds%begp:,1:) 
   !
   ! !LOCAL VARIABLES
   real(r8)                             , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: croot_prof_patch(bounds%begp:,1:) 

   real(r8), parameter :: N_stress_max = 2.0_r8       ! Maximum N demand of plant, based on current N amount in plant []
   real(r8)    :: N_stress                            ! N demand of plant, based on current N amount in plant []
   real(r8)    :: C_transfer                          ! Carbon allocated to N acquisition, higher N_stress leads to higher C_transfer
   real(r8)    :: N_to_plant_scav
   real(r8)    :: N_to_plant_mine
   real(r8)    :: N_to_plant_fix
   real(r8)    :: t_soi_degC
   real(r8)    :: no3_passiv     ! Passive NO3 (nitrate) uptake
   real(r8)    :: nh4_passiv     ! Passive NH4 (ammonium) uptake
   

   associate(                                                   &
   sulman_cn_m          => params_inst%sulman_cn_m            , &   !Soil microbial C:N ratio
   sulman_v_nh4         => params_inst%sulman_v_nh4           , &   !Maximum NH4+ immobilization rate [year-1]
   sulman_v_no3         => params_inst%sulman_v_no3           , &   !Maximum NO3- immobilization rate [year-1]
   sulman_vmax_denit    => params_inst%sulman_vmax_denit      , &   !Maximum denitrification decomposition rate at reference temperature [year-1]
   sulman_fden          => params_inst%sulman_fden            , &   !Maximum denitrification decomposition rate at reference temperature [unitless]
   sulman_kdenit        => params_inst%sulman_kdenit          , &   !Half-saturation constant for nitrate concentration in denitrification [kg NO3-N kg NO3-N demand-1 year-1]
   sulman_root_no3      => params_inst%sulman_root_no3        , &   !Maximum root active nitrate uptake rate [kg N m-3 year-1]
   sulman_root_nh4      => params_inst%sulman_root_nh4        , &   !Maximum root active ammonium uptake rate [kg N m-3 year-1]
   sulman_km_no3        => params_inst%sulman_km_no3          , &   !Half-saturation nitrate concentration for root active uptake [kg N m-3]
   sulman_km_nh4        => params_inst%sulman_km_nh4          , &   !Half-saturation nitrate concentration for root active uptake [kg N m-3]
   sulman_r_rhiz        => params_inst%sulman_r_rhiz          , &   !Radius of the rhizosphere [m]
   sulman_v_scav        => params_inst%sulman_v_scav          , &   !Maximum N uptake rate by scavenging mycorrhizae [kg N m-3 year-1]
   sulman_k_scav_Ninorg => params_inst%sulman_k_scav_Ninorg   , &   !Half-saturation inorganic N concentration for mycorrhizal uptake [kg N m-3]
   sulman_k_scav        => params_inst%sulman_k_scav          , &   !Half-saturation mycorrhizal biomass concentration for scavenging [kg C m-3]
   sulman_km_mine       => params_inst%sulman_km_mine         , &   !Half-saturation mycorrhizal biomass concentration for mining [kg C m-3]
   sulman_cue_mine      => params_inst%sulman_cue_mine        , &   !Carbon use efficiency of mycorrhizal mining [fraction]
   sulman_nue_mine      => params_inst%sulman_nue_mine        , &   !Nitrogen use efficiency of mycorrhizal mining [fraction]
   sulman_vmax_ref_mine => params_inst%sulman_vmax_ref_mine   , &   !Maximum decomposition rate at reference temperature for mycorrhizal mining [year-1]
   sulman_rfix          => params_inst%sulman_rfix            , &   !N fixation rate per unit symbiotic biomass [kg N kg biomass C-1 year-1]
   sulman_kgrowth       => params_inst%sulman_kgrowth         , &   !Half-saturation of intermediate C pool for symbiotic growth [kg C m -2]
   sulman_rgrowth       => params_inst%sulman_rgrowth         , &   !Maximum symbiont growth rate [kg C m-2 year-1]
   sulman_tau_sym       => params_inst%sulman_tau_sym         , &   !Fraction of symbiotic biomass turnover not used for maintenance respiration [fraction]
   sulman_growth_scav   => params_inst%sulman_growth_scav     , &   !N scavenger growth efficiency [unitless]
   sulman_growth_mine   => params_inst%sulman_growth_mine     , &   !N miner growth efficiency [unitless]
   sulman_growth_fix    => params_inst%sulman_growth_fix      , &   !N fixer growth efficiency [unitless]
   sulman_tau_scav      => params_inst%sulman_tau_scav        , &   !N scavenger turnover time [year-1]
   sulman_tau_mine      => params_inst%sulman_tau_mine        , &   !N miner turnover time [year-1]
   sulman_tau_fix       => params_inst%sulman_tau_fix         , &   !N fixer turnover time [year-1]
   sulman_cn_scav       => params_inst%sulman_cn_scav         , &   !N scavenger C:N [unitless]
   sulman_cn_mine       => params_inst%sulman_cn_mine         , &   !N miner C:N [unitless]
   sulman_cn_fix        => params_inst%sulman_cn_fix          , &   !N fixer C:N [unitless]
   sulman_tau_int       => params_inst%sulman_tau_int         , &   !Turnover time of intermediate C pool [year-1]
   sulman_rup_veg       => params_inst%sulman_rup_veg         , &   !Vegetation N uptake rate from intermediate N pool [year-1]
   sulman_fnalloc       => params_inst%sulman_fnalloc         , & !Fraction of NPP allocated to N uptake per unit N stress [fraction] 

   leafn                => cnveg_nitrogenstate_inst%leafn_patch                  , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
   leafn_storage        => cnveg_nitrogenstate_inst%leafn_storage_patch          , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
   frootn               => cnveg_nitrogenstate_inst%frootn_patch                 , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N                               
   frootn_storage       => cnveg_nitrogenstate_inst%frootn_storage_patch         , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
   livecrootn           => cnveg_nitrogenstate_inst%livecrootn_patch             , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
   livecrootn_storage   => cnveg_nitrogenstate_inst%livecrootn_storage_patch     , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
   h2osoi_liq           => waterstatebulk_inst%h2osoi_liq_col                    , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)  
   t_soisno             => temperature_inst%t_soisno_col                         , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       
   availc               => cnveg_carbonflux_inst%availc_patch                    & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
   )
   
   sulman_vmax_denit_fast = params_inst%sulman_vmax_denit(1)
   sulman_vmax_denit_slow = params_inst%sulman_vmax_denit(2)
   sulman_vmax_denit_necr = params_inst%sulman_vmax_denit(3)

   sulman_cue_mine_fast = params_inst%sulman_cue_mine(1)
   sulman_cue_mine_slow = params_inst%sulman_cue_mine(2)
   sulman_cue_mine_necr = params_inst%sulman_cue_mine(3)

   sulman_nue_mine_fast = params_inst%sulman_nue_mine(1)
   sulman_nue_mine_slow = params_inst%sulman_nue_mine(2)
   sulman_nue_mine_necr = params_inst%sulman_nue_mine(3)

   sulman_vmax_ref_mine_fast = params_inst%sulman_vmax_ref_mine(1)
   sulman_vmax_ref_mine_slow = params_inst%sulman_vmax_ref_mine(2)
   sulman_vmax_ref_mine_necr = params_inst%sulman_vmax_ref_mine(3)
   
   !-----------------------------------------------------------------------

   ! Calculation of N_stress: N stress indicates the nitrogen content in a plant and is used to calculate how much carbon 
   ! the plant allocates belowground = the higher N stress, the higher is the carbon allocation belowground
   N_stress = (N_stress_max*(leafn + frootn + livecrootn) - (leafn_storage + frootn_storage + livecrootn_storage)) / &
               (leafn + frootn + livecrootn)

   ! Calculate the amount of C transferent to N aquisation
   ! Dynamic allocation of a fraction of NPP to root exudation
   C_transfer = max(availc(p),0) * sulman_fnalloc * N_stress

   
   ! Return of investment function calculates plant nitrogen aquisition per unit of symbiont biomass 
   ! multiplied by the efficency of symbiont biomass production and turnover time

   ! For symbionts miner & scavenger
   
   ! Calculate plant N uptake per symbiont pathway (Nj)
   ! Nj = Nint,pool (from previous timestep) * sulman_rup_veg  

   N_to_plant_scav = scav_N_reservoir * sulman_rup_veg * dt
   N_to_plant_mine = mine_N_reservoir * sulman_rup_veg * dt 
   N_to_plant_fix = fix_N_reservoir  * sulman_rup_veg * dt

   !-----------------------------------------------------------------------

   ! Plant uptake strategies:

   ! Passive nitrogen uptake from the rhizosphere by roots with soil water uptake:
   
   ! smin_no3_vr is per column and has to be per patch
   ! if water in layer:

   t_soi_degC = t_soisno_col(c,j) - tfrz     ! Soil temperature in degrees Celcius

   if (t_soi_degC < 0.01) then
   no3_passiv = waterflux_type%qflx_tran_veg_patch * (smin_no3_vr / waterstatebulk_inst%h2osoi_liq_col)
   nh4_passiv = waterflux_type%qflx_tran_veg_patch * (smin_nh4_vr / waterstatebulk_inst%h2osoi_liq_col)
   end if



   ! Root uptake of nitrogen
   call root_N_uptake()

   ! Scavenging (AM-style)
   call myc_scavenger_N_uptake()

   ! Mycorrhizal N mining (ECM-style)
   call myc_miner_N_uptake()

   ! Symbiotic N2 Fixation
   N_fixation = fix_C_biomass*sulman_rfix*dt ! sulman_rfix different values

   !-----------------------------------------------------------------------
   
  
   ! Mycorrhizal scavengers
   scav_growth = sulman_rgrowth * scav_C_reservoir / (scav_C_reservoir + sulman_kgrowth) * sulman_growth_scav * dt
   maint_resp = min(scav_C_biomass/sulman_tau_scav * (1.0 - sulman_tau_sym) * dt, scav_growth)
   ! Nitrogen limitation
   if (scav_growth - maint_resp > sulman_cn_scav * scav_N_reservoir * 0.9) then
      ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
      scav_growth = sulman_cn_scav * scav_N_reservoir * 0.9 + maint_resp
   end if

   total_myc_resp = total_myc_resp + scav_growth / sulman_growth_scav * (1.0 - sulman_growth_scav)

   ! Intermediate pool scavenger
   scav_N_reservoir = scav_N_reservoir + scav_N_biomass * (1 - sulman_tau_sym / sulman_tau_scav *dt)
   scav_C_biomass = scav_C_biomass + scav_growth - scav_C_biomass / sulman_tau_scav * dt 
   ! scav_N_biomass = scav_N_biomass + scav_growth - maint_resp) / sulman_cn_scav - scav_N_biomass / sulman_tau_scav * sulman_tau_sym * dt
   scav_N_biomass = scav_C_biomass / sulman_cn_scav
   scav_C_reservoir = scav_C_reservoir - scav_growth / sulman_growth_scav
   ! scav_N_reservoir = scav_N_reservoir - (scav_growth - maint_resp) / sulman_cn_scav
   scav_N_reservoir = scav_N_reservoir - scav_N_biomass


   ! Mycorrhizal miners
   mine_growth = sulman_rgrowth * mine_C_reservoir / (mine_C_reservoir + sulman_kgrowth) * sulman_growth_mine * dt
   maint_resp = min(mine_C_biomass/sulman_tau_mine * (1.0 - sulman_tau_sym) * dt, mine_growth)
   ! Nitrogen limitation
   if (mine_growth - maint_resp > sulman_cn_mine * mine_N_reservoir * 0.9) then
      ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
      mine_growth = sulman_cn_mine * mine_N_reservoir * 0.9 + maint_resp
   end if

   total_myc_resp = total_myc_resp + mine_growth / sulman_growth_mine * (1.0 - sulman_growth_mine)

   ! Intermediate pool miner
   mine_N_reservoir = mine_N_reservoir + mine_N_biomass * (1 - sulman_tau_sym / sulman_tau_mine *dt)
   mine_C_biomass = mine_C_biomass + mine_growth - mine_C_biomass / sulman_tau_mine * dt 
   ! mine_N_biomass = mine_N_biomass + mine_growth - maint_resp) / sulman_cn_mine - mine_N_biomass / sulman_tau_mine * sulman_tau_sym * dt
   mine_N_biomass = mine_C_biomass / sulman_cn_mine
   mine_C_reservoir = mine_C_reservoir - mine_growth / sulman_growth_mine
   ! mine_N_reservoir = mine_N_reservoir - (mine_growth - maint_resp) / sulman_cn_mine
   mine_N_reservoir = mine_N_reservoir - mine_N_biomass


   ! Nitrogen Fixation
   fix_growth = sulman_rgrowth * fix_C_reservoir / (fix_C_reservoir + sulman_kgrowth) * sulman_growth_fix * dt
   maint_resp = min(fix_C_biomass / sulman_tau_fix * (1.0 - sulman_tau_sym) * dt, fix_growth)
   ! if (fix_growth > sulman_cn_fix * fix_N_reservoir * 0.9) then
      ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
      ! fix_growth = sulman_cn_fix * fix_N_reservoir * 0.9
   ! end if

   total_myc_resp = total_myc_resp + fix_growth / sulman_growth_fix * (1.0 - sulman_growth_fix)

   ! NOT SURE WHAT EXACTLY N_FIXATION IS
   N_fixation = N_fixation - fix_N_biomass * (1 - sulman_tau_sym / sulman_tau_fix * dt)
   fix_C_biomass = fix_C_biomass + fix_growth - fix_C_biomass / sulman_tau_fix * dt
   ! fix_N_biomass = fix_N_biomass + (fix_growth - maint_resp) / sulman_cn_fix - fix_N_biomass / sulman_tau_fix * sulman_tau_sym * dt
   fix_N_biomass = fix_N_biomass / sulman_cn_fix
   fix_C_reservoir = fix_C_reservoir - (fix_growth) / sulman_growth_fix
   ! fix_N_reservoir = fix_N_reservoir - fix_growth / sulman_cn_mine ! WHY C:N MINER?
   ! N fixers just make all the N they need for their biomass
   ! N_fixation = N_fixation + (fix_growth - maint_resp) / sulman_cn_mine ! WHY C:N MINER?
   N_fixation = N_fixation + fix_N_biomass


   !----------------------------------------------------------------

   ! Calculate N released to plants per uptake pathway - Return Of Investment (line 553 in vegn_dynamics)
   ! multiplied with secspday * days_per_year to get values per second (parameter is in per year)

   days_per_year = get_average_days_per_year()     ! to get average number of days per year 

   ! Scavengers
   if () ! IF SCAV PATHWAY IS ACTIVE DO...
   scav_N_to_plant = scav_N_reservoir * sulman_rup_veg * secspday * days_per_year
   if (scav_C_biomass < 0.0_r8) then 
      scav_roi = (max(0.0_r8, scav_N_to_plant)dt) / (scav_C_biomass / sulman_growth_scav / sulman_tau_scav)
   else 
      ! scav_efficiency is calculated in myc_scavenger_N_uptake under myc_efficiency
      scav_roi = scav_efficiency / (dt * sulman_growth_scav * sulman_tau_scav) 
   end if 
   else
      scav_N_to_plant = 0.0_r8 ; scav_roi = 0.0_r8
   end if 

   scav_N_reservoir = scav_N_reservoir - scav_N_to_plant


   ! Miners
   if () ! IF MINE PATHWAY IS ACTIVE DO...
   mine_N_to_plant = mine_N_reservoir * sulman_rup_veg * secspday * days_per_year
   if (mine_C_biomass < 0.0_r8) then 
      mine_roi = (max(0.0_r8, mine_N_to_plant)dt) / (mine_C_biomass / sulman_growth_mine / sulman_tau_mine)
   else 
      ! mine is calculated in one of the mining routines under myc_efficiency
      mine_roi = mine_efficiency / (dt * sulman_growth_mine * sulman_tau_mine) 
   end if 
   else
      mine_N_to_plant = 0.0_r8 ; mine_roi = 0.0_r8
   end if 

   mine_N_reservoir = mine_N_reservoir - mine_N_to_plant


   !NOT SUPER SURE WITH THIS ONE
   ! Root uptake of Nitrogen
   if (C_allocation_to_N_acq > 0.0_r8) then                                ! also called: C_alloc_to_N_acq, calculatet in vegn_carbon_int_lm3
      root_roi = max(0.001, root_N_uptake / dt) / C_allocation_to_N_acq    ! C_allocation_to_N_acq & root_N_uptake ?
   else
      root_roi = (scav_roi + mine_roi) * 0.25
   end if 


   ! Nitrogen Fixers
   if () ! IF FIXER PATHWAY IS ACTIVE DO...
   fix_N_to_plant = fix_N_reservoir * sulman_rup_veg * secspday * days_per_year
   if (fix_C_biomass < 0.0_r8) then 
      fix_roi = (fix_N_to_plant / dt) / (fix_C_biomass / sulman_growth_fix / sulman_tau_fix)
   else 
      ! mine is calculated in one of the mining routines under myc_efficiency
      fix_roi =  sulman_rfix * sulman_growth_fix * sulman_tau_fix
   end if 
   else
      fix_N_to_plant = 0.0_r8 ; fix_roi = 0.0_r8
   end if 

   fix_N_reservoir = fix_N_reservoir - fix_N_to_plant


   end associate

   end subroutine CN_soil_veg_exchange

   !----------------------------------------

   ! PLANT UPTAKE STRATEGIES

    !----------------------------------------

   function Vmax_myc(soil_T)
      real, intent(in)  :: soil_T               ! Soil temperature in Kelvin
      real, parameter   :: Tref=293.15
      real(r8)          :: alpha
      real(r8)          :: Vmax_myc
      real, parameter   :: Ea=(37e3)              ! Activation energy (kJ/mol)
      real, public, parameter :: R_gas = 8.314472 ! universal gas constant, J K-1 mol-1

      alpha = params_inst%sulman_vmax_ref_mine / exp(-Ea /(R_gas*Tref))
      Vmax_myc = alpha * exp(-Ea / R_gas * soil_T)
   end function Vmax_myc

   !
   function theta 
    

   function resp_myc(Cavail, myc_biomass, soil_T, wliq, wair)
      ! DO STH WITH CAVAIL

      real, intent(in) :: myc_biomass                ! mycorrhizal biomass [kgC/m2]
      real, intent(in) :: soil_T                     ! Soil temperature in Kelvin
      real, intent(in) :: wliq                       ! water liquid                        
      real, intent(in) :: wair                       ! air in soil

      real, parameter  :: enzyme_frac=1.0            ! Relative amount of enzymes produced by microbes
      real(r8)         :: enzymes
      real(r8)         :: Cavail

      enzymes = myc_biomass * enzyme_frac

      ! Pre check: if ther eis no carbon / no water content in soil / no enzymes = no mycorrhizal repiration
      if (sum(Cavail).eq.0.0_r8 .or. wliq .eq. 0.0_r8 .or. enzymes .eq. 0.0_r8) then
         resp_myc = 0.0_r8
         return
      endif

      ! If there is carbon avaliable, calculate mycorrhizal repiration
      if (Cavail < 0.0_r8) then 
         resp_myc = Vmax_myc(soil_T) * Cavail * enzymes / (sum(Cavail) * params_inst%sulman_km_mine + enzymes) * &
         theta_func(wliq, wair)
      else 
         resp_myc = 0.0_r8
      end if 

   end function resp_myc

   !----------------------------------------
  
   ! PLANT UPTAKE STRATEGIES
  
   ! Nitrogen uptake from the rhizosphere by roots (active transport across root-soil interface)
   ! Mineral nitrogen is taken up from the rhizosphere only
   subroutine active_root_N_uptake(bounds, froot_carbon, soilbiogeochem_nitrogenstate_inst, soilstate_inst, &
      no3_uptake, nh4_uptake) ! not done

      use clm_varcon        , only : rpi
      use decompMod         , only : bounds_type

      type(bounds_type)      , intent(in)    :: bounds
      type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
      type(soilstate_type)            , intent(in)    :: soilstate_inst

    ! !LOCAL VARIABLES:
      integer :: begp, endp
      integer :: begc, endc
   
      begp = bounds%begp; endp= bounds%endp
      begc = bounds%begc; endc= bounds%endc


      real,intent(in)    :: froot_carbon( bounds%begp: )    ! fine root carbon (gC/m2) [pft]   
      real,intent(inout) :: no3_uptake
      real,intent(inout) :: nh4_uptake

      real(r8) :: root_cross_sec_area
      real(r8) :: root_length_density
      real(r8) :: root_biomass_density
      real(r8) :: rhizosphere_frac
      real(r8) :: no3_uptake
      real(r8) :: nh4_uptake
      real(r8), parameter :: c_to_b = 2.0_r8           !(g biomass /g C)

      associate(                                                                           &
         sminn_vr               => soilbiogeochem_nitrogenstate_inst%sminn_vr_col        , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N                
         smin_nh4_vr            => soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col     , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NH4              
         smin_no3_vr            => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col     , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NO3              
         ivt                    => patch%itype                                           , & ! Input:  [integer  (:)   ]  patch vegetation type
         rootfr                 => soilstate_inst%rootfr_patch                           , & ! Input:   [real(r8) (:,:)]
         dz                     => col%dz                                                , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         root_radius  => pftcon%root_radius                  , & ! Input: 0.29e-03_r8 !(m) 
         root_density => pftcon%root_density                  & ! Input: 0.31e06_r8 !(g biomass / m3 root) 
         )
      
     

         no3_uptake = 0.0_r8
         nh4_uptake = 0.0_r8
         total_active_root_N_uptake = 0.0_r8

      do fp = 1,num_soilp        ! PFT Starts
         p = filter_soilp(fp)
         c = patch%column(p)
         do j = 1, nlevdecomp

           ! calculate conversion from conductivity to conductance
           root_biomass_density = c_to_b * froot_carbon(p) * rootfr(p,j) / dz(c,j)
           ! ensure minimum root biomass (using 1gC/m2)
           root_biomass_density = max(c_to_b*1._r8,root_biomass_density)
           ! Root length density: m root per m3 soil 
           root_cross_sec_area = rpi*params_inst%sulman_r_rhiz(ivt(p))**2
           root_length_density = root_biomass_density / (root_density(ivt(p)) * root_cross_sec_area)
     
           rhizosphere_frac = min(rpi*((sulman_r_rhiz+r_r)**2-r_r**2)*root_length_density,1.0_r8) 
           
           no3_uptake = rhizosphere_frac * params_inst%sulman_root_no3 * (smin_no3_vr(c,j) / (smin_no3_vr(c,j) + params_inst%sulman_km_no3))
           nh4_uptake = rhizosphere_frac * params_inst%sulman_root_nh4 * (smin_nh4_vr(c,j) / (smin_nh4_vr(c,j) + params_inst%sulman_km_nh4))
  
           ! add some if statements that if plant wants to take up more no3 /nh4 than in soil sth happens
  
           total_active_root_N_uptake(p) = total_active_root_N_uptake(p) + no3_uptake + nh4_uptake
         end do
      end do

      end associate
     
   end subroutine active_root_N_uptake


   
   subroutine myc_scavenger_N_uptake(X, X, X, myc_efficiency, dt, X)

      real,intent(in)     :: myc_biomass             ! (kgC/m2)
      real,intent(in)     :: layer_thickness         ! (m)
      real,intent(inout)  :: no3_uptake              ! (kgN/m2/year)
      real,intent(inout)  :: nh4_uptake              ! (kgN/m2/year)
      real, intent(out)  :: myc_efficiency           ! units: kgN/kg myc biomass C. Should give N uptake efficiency even when myc biomass is zero

   
      associate(  
      crootfr           => soilstate_inst%crootfr_patch                          , & ! Input:   [real(r8) (:,:)] fraction of roots for carbon in each soil layer  (nlevgrnd)
      frootc            => cnveg_carbonstate_inst%frootc_patch                   , & ! Input:   [real(r8) (:)]  (gC/m2) fine root C
      smin_nh4_vr       => soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col     , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NH4              
      smin_no3_vr       => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col     , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NO3              
      sulman_k_scav_Ninorg => params_inst%sulman_k_scav_Ninorg                   , & ! Half-saturation inorganic N concentration for mycorrhizal uptake [kg N m-3]
      sulman_k_scav        => params_inst%sulman_k_scav                            & ! Half-saturation mycorrhizal biomass concentration for scavenging [kg C m-3]
      )

      ! Calculationg a root profile
      ! https://escomp.github.io/ctsm-docs/versions/master/html/tech_note/Plant_Hydraulics/CLM50_Tech_Note_Plant_Hydraulics.html?highlight=root
      root_dens_sum = sum(crootfr(p,:) * frootc(p,:))
      do j = 1:nlev
         root_dens(p,j) = (crootfr(p,j) * frootc(p,j)) / root_dens_sum
      end do 

      ! Calculating mycorrhizal biomass per soil layer
      myc%C_pool(p)
      myc_c_layer = myc%C_pool(p) * root_dens(p,j)

      ! Check if there is  mycorrhizal biomass in soil layer
      do j = 1:nlev
         if (myc_c_layer > 0) then 

            ! If there is mycorrhizal biomass in soil layer, calculate N uptake
            do j = 1, nlevdecomp
               ! sulman_v_scav is different in Sulman code
               no3_uptake = sulman_v_scav * (smin_no3_vr / dz) / ((smin_no3_vr / dz) + sulman_k_scav_Ninorg) * &
                           (scav_C_biomass / dz) / ((scav_C_biomass / dz) + sulman_k_scav)

               nh4_uptake = sulman_v_scav * (smin_nh4_vr / dz) / ((smin_nh4_vr / dz) + sulman_k_scav_Ninorg) * &
                           (scav_C_biomass / dz) / ((scav_C_biomass / dz) + sulman_k_scav)

               !no3_uptake = min() 
               !nh4_uptake = min()
               total_n_uptake = total_n_uptake + (no3_uptake + nh4_uptake) * dt
            end do
         end if
      end do 

      end associate
   end subroutine myc_scavenger_N_uptake

   !----------------------------------------

   ! Miner
   subroutine myc_miner_N_uptake (temperature_inst)

      ! USES
      use TemperatureType                 , only : temperature_type
      use clm_varcon                      , only : tfrz
      use clm_varcon                      , only : denh2o, denice
       

      ! ARGUMENTS
      type(temperature_type)              , intent(in)    :: temperature_inst


      real(r8)         :: t_soi_degC              ! Soil temperature            (degrees Celsius)
      real, intent(in) :: myc_biomass             ! (kgC/m2)
      
      real, intent(out) :: myc_resp
      real, intent(out) :: N_uptake

      real(r8), intent(in)    :: h2osoi_liq  ! liquid water content kg/m2
      real(r8), intent(in)    :: watsat      ! porosity m3/m3
      real(r8), intent(in)    :: h2osoi_ice  ! ice content kg/m2 
      real(r8), intent(in)    :: dz          ! soil layer thickness 
      
      real(r8):: wliq ! water liquid
      real(r8):: wice ! water ice
      real(r8):: wair ! air in soil

     

      associate(
         t_soisno       => temperature_inst%t_soisno_col        , & ! Input:  [real(r8) (:,:)]soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         watsat         => soilstate_inst%watsat_col            , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity)  
         h2osoi_liq     => waterstatebulk_inst%h2osoi_liq_col   , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)             
         h2osoi_ice     => waterstatebulk_inst%h2osoi_ice_col   , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2) (new)
         t_soisno       => temperature_inst%t_soisno_col          & ! Input:   [real(r8) (:,:)] soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         
         )

      ! Calculate soil temperature for each soil layer
      do j = 1, nlevdecomp
         t_soi_degC         = t_soisno(c,j)  -   tfrz  ! tfrz = 273.15
      end do 
      

      ! Calculating water, ice and air content in soil (equivalent to air_filled porosity, theta, theta sat in Sulman)
      wliq = h2osoi_liq / dz * denh2o
      wice = h2osoi_ice / dz * denice
      wliq = min(1.0_r8, wliq/watsat)            ! fraction of liquid water-filled pore space (0.0 - 1.0)
      wice = min(1.0_r8, wice/watsat)            ! fraction of frozen water-filled pore space (0.0 - 1.0)
      wair = max(0.0_r8, 1.0_r8 - wliq- wice)    ! fraction of air-filled pore space (0.0 - 1.0)
 
      ! Calculationg a root profile
      ! https://escomp.github.io/ctsm-docs/versions/master/html/tech_note/Plant_Hydraulics/CLM50_Tech_Note_Plant_Hydraulics.html?highlight=root
      root_dens_sum = sum(crootfr(p,:) * frootc(p,:))
      do j = 1:nlev
         root_dens(p,j) = (crootfr(p,j) * frootc(p,j)) / root_dens_sum
      end do 

      ! Calculating mycorrhizal biomass per soil layer
      myc%C_pool(p)
      myc_c_layer = myc%C_pool(p) * root_dens(p,j)


      do j = 1, nlevdecomp
         call miner_decomposition(XX, XX, t_soisno(c,j), wliq(c,j), wair(c,j), XX, XX)
      end do

     ! Organic N Mining by Mycorrhizal Fungi
      end associate

   end subroutine myc_miner_N_uptake


   subroutine miner_decomposition()

      real, intent(in):: wliq ! water liquid
      real, intent(in):: wice ! water ice
      real, intent(in):: wair ! air in soil
      real, intent(in):: t_soi_degC ! soil temperature
      real, intent(in):: myc_biomass

      real, intent(out):: N_uptake
      real, intent(out):: myc_resp
      

      myc_resp    = 0.0_r8
      N_uptake    = 0.0_r8

      


      ! Equation 35
      real(r8) :: N_mine ! Flux from SOM pools to intermediate miner pool

      ! Vmax Michaelis Menten Kinetics
      ! A function? that gives Vmax based on a reference temperature

      ! Get the moisture function von MIMICS
      moist_mod = r_moist(h2osoi_liq(c,j),watsat(c,j), h2osoi_ice(c,j), col%dz(c,j)) 


      N_mine = params_inst%sulman_vmax_ref_mine * moist_mod



      myc_resp = myc_resp + sum() * (1 - params_inst%sulman_cue_mine)
      N_uptake = N_uptake + sum() * params_inst%sulman_nue_mine


     

      end associate
   end subroutine miner_decomposition
  
        
   !----------------------------------------

   ! GROWTH AND TURNOVER OF SYMBIOTIC BIOMASS
   ! Intermediate pools:
   ! line 553 - 604 in veg_dynamics_mod

   subroutine intermediate_pools(fixer_inst, mine_inst, scav_inst )

   real(r8) :: scav_growth             ! Symbiotic biomass growth rate for scavangers 
   real(r8) :: mine_growth             ! Symbiotic biomass growth rate for miners     
   real(r8) :: fix_growth              ! Symbiotic biomass growth rate for fixers     
   real(r8) :: maint_resp              ! Respiration, calculated for scav, mine, fix and updated to total_myc_resp
   real(r8) :: total_myc_resp          ! Total repiration of scav, mine and fix
   
   real(r8) :: scav_C_biomass          ! Scavenger C biomass
   real(r8) :: scav_N_biomass          ! Scavenger N biomass
   real(r8) :: mine_C_biomass          ! Miner C biomass
   real(r8) :: mine_N_biomass          ! Miner N biomass
   real(r8) :: fix_C_biomass           ! Fixer C biomass
   real(r8) :: fix_N_biomass           ! Fixer N biomass

   real(r8) :: scav_C_reservoir        ! Carbon reservoir from previous timestep for scavangers
   real(r8) :: scav_N_reservoir        ! Nitrogen reservoir from previous timestep for scavangers
   real(r8) :: mine_C_reservoir        ! Carbon reservoir from previous timestep for miners
   real(r8) :: mine_N_reservoir        ! Nitrogen reservoir from previous timestep for miners
   real(r8) :: fix_C_reservoir         ! Carbon reservoir from previous timestep for fixers   
   real(r8) :: fix_N_reservoir         ! Nitrogen reservoir from previous timestep for fixers
   real(r8) :: N_fixation              !


   ! Carbon intermediate pools
   !fixer_inst%C_inter_biomass
   !mine_inst%C_inter_biomass
   !scav_inst%C_inter_biomass

   ! Nitrogen intermediate pools
   !fixer_inst%N_inter_biomass
   !miner_inst%N_inter_biomass
   !scav_inst%N_inter_biomass

   type(symbiont_type)   , intent(inout) :: fixer_inst, mine_inst, scav_inst

   !----------------------------------------------------------------

   
   ! WARNING FLAGS

   end subroutine intermediate_pools
           



end module CNSoilVegMIMICSplus
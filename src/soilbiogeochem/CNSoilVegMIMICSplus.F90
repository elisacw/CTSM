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
  use SoilStateType                      , only : soilstate_type

  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: init_mimicsplus_veg_myc         !
  private :: CNallocation                    !
  private :: calc_roi                        !
 
  type, public :: symbiont_type

  real(r8), pointer, private :: C_biomass             (:) ! Carbon biomass
  real(r8), pointer, private :: N_biomass             (:) ! Nitrogen biomass
  real(r8), pointer, private :: C_inter_biomass       (:) ! Carbon intermediate pool biomass
  real(r8), pointer, private :: N_inter_biomass       (:) ! Nitrogen intermediate pool biomass
  
   contains

     ! Public procedures
     procedure, public  :: Init
     procedure, public  :: Restart

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
   real(r8)    :: t_soi_degC
   real(r8)    :: no3_passiv     ! Passive NO3 (nitrate) uptake
   real(r8)    :: nh4_passiv     ! Passive NH4 (ammonium) uptake
   
      
   
!ECW get parameters from the mimics module
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
   sulman_k_scav_Ninorg => params_inst%sulman_k_scav_Ninorg  , &   !Half-saturation inorganic N concentration for mycorrhizal uptake [kg N m-3]
   sulman_k_scav        => params_inst%sulman_k_scav          , &   !Half-saturation mycorrhizal biomass concentration for scavenging [kg C m-3]
   sulman_km_mine       => params_inst%sulman_km_mine         , &   !Half-saturation mycorrhizal biomass concentration for scavenging [kg C m-3]
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
   t_soisno             => temperature_inst%t_soisno_col               , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       
            
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
   C_transfer = max(NPP,0) * sulman_fnalloc * N_stress

   
   ! Return of investment function calculates plant nitrogen aquisition per unit of symbiont biomass 
   ! multiplied by the efficency of symbiont biomass production and turnover time

   ! For symbionts miner & scavenger
   call calc_roi(plantN_uptake, sulman_tau_scav, sulman_growth_scav, scav_C_biomass, roi, f_alloc)
   call calc_roi(plantN_uptake, sulman_tau_mine, sulman_growth_mine, mine_C_biomass, roi, f_alloc)

   !-----------------------------------------------------------------------

   ! Plant uptake strategies:

   ! Passive nitrogen uptake from the rhizosphere by roots with soil water uptake:
   
   ! smin_no3_vr is per column and has to be per patch
   ! if water in layer:
   
   t_soi_degC = t_soisno_col(c,j) - tfrz

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


   end associate
   end subroutine CN_soil_veg_exchange


   ! Calculate fraction of rhizosphere for each soil layer
   ! r_rhiz, 3686 in soil.f90

   !---------------------------------

   subroutine calc_roi(plantN_uptake, tau_symb, growth_eff, biomass_symb, roi, f_alloc)
      !
      ! DESCRIPTION:
      ! Calculates ROI (Return Of Investment) of nitrogen acquisition strategy 
      !
      ! !USES:
      use clm_varcon       , only : secspday, secsphr, tfrz, spval
            
      !
      ! !ARGUMENTS:
      real(r8), intent(in)     :: plantN_uptake     ! Plant N aquisition
      real(r8), intent(in)     :: tau_symb          ! Turnover time of symbionet (sulman_tau_sym, sulman_tau_fix)
      real(r8), intent(in)     :: growth_eff        ! Growth efficiency of symbionts
      real(r8), intent(in)     :: biomass_symb      ! Symbiont biomass
      real(r8), intent(inout)  :: roi               ! RoI nitrogen per carbon invested from vegetation to mycorrhiza [gN/gC]
      real(r8), intent(inout)  :: f_alloc           ! Fraction allocated to uptake strategy

      ! This might need to be limited if there is little N, but not sure
      roi = (plantN_uptake * tau_symb * growth_eff) / biomass_symb

   end subroutine calc_roi

   
   !----------------------------------------

   ! Include Nitrogen limitation from FUN

   !----------------------------------------

   ! PLANT UPTAKE STRATEGIES
  
  
   ! Nitrogen uptake from the rhizosphere by roots (active transport across root-soil interface)
   ! Mineral nitrogen is taken up from the rhizosphere only
   subroutine active_root_N_uptake(bounds, froot_carbon, soilbiogeochem_nitrogenstate_inst, soilstate_inst, &
      no3_uptake, nh4_uptake) ! not done

      use clm_varcon        , only : rpi
      use decompMod         , only : bounds_type, 

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
     
   end subroutine active_root_N_uptake


   
   subroutine myc_scavenger_N_uptake ()

      ! Check if there is mycorrhizal bimass in the soil layers - with a loop
      ! If there is no mycorrhizal biomass, still do calculation so efficiency can be estimated
      ! This allows plants to grow some biomass from zero
      ! Check if there is any living biomass in the soil - if not uptake is 0 

      ! You need to know how much scavenger biomass is in the soil layers
      ! So you need to calculate a soil profile

      ! Once you have that you can call the things below (need to turn into subroution) and calculate the no3 and nh4 uptake

      ! Then you somehow scale the N uptake? 

      ! And then I am not sure what happens

   end subroutine myc_scavenger_N_uptake


   ! Inorganic N Scaveging by Mycorrhizal Fungi
   subroutine scavenger_mineral_N_uptake_rate (myc_biomass,layer_thickness,no3_uptake,nh4_uptake)
   
   real,intent(in)     :: myc_biomass             ! (kgC/m2)
   real,intent(in)     :: layer_thickness         ! (m)
   real,intent(inout)  :: no3_uptake              ! (kgN/m2/year)
   real,intent(inout)  :: nh4_uptake              ! (kgN/m2/year)
   
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
   
   end subroutine scavenger_mineral_N_uptake_rate 



   subroutine myc_miner_N_uptake ()

     ! Organic N Mining by Mycorrhizal Fungi

   end subroutine myc_miner_N_uptake
  
        
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
   fixer_inst%C_inter_biomass
   mine_inst%C_inter_biomass
   scav_inst%C_inter_biomass

   ! Nitrogen intermediate pools
   fixer_inst%N_inter_biomass
   miner_inst%N_inter_biomass
   scav_inst%N_inter_biomass

   type(symbiont_type)   , intent(inout) :: fixer_inst, mine_inst, scav_inst

   !----------------------------------------------------------------

   ! Mycorrhizal scavengers
   scav_growth = sulman_rgrowth * scav_C_reservoir / (scav_C_reservoir + sulman_kgrowth) * sulman_growth_scav * dt
   maint_resp = min(scav_C_biomass/sulman_tau_scav * (1.0 - sulman_tau_sym) * dt, scav_growth)
   ! Nitrogen limitation
   if (scav_growth - maint_resp > sulman_cn_scav * scav_N_reservoir * 0.9) then
      ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
      scav_growth = sulman_cn_scav * scav_N_reservoir * 0.9 + maint_resp
   end if

   total_myc_resp = total_myc_resp + scav_growth / sulman_growth_scav * (1.0 - sulman_growth_scav)

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

   ! WARNING FLAGS

   end subroutine intermediate_pools
           



end module CNSoilVegMIMICSplus
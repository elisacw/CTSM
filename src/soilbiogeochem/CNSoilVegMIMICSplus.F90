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
     procedure, private :: InitCold

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
    call this%InitCold     (bounds)

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

   !------------------------------------------------------------------------
   subroutine Clean(this)
      !
      ! !ARGUMENTS:
      class(symbiont_type) :: this
      !
      ! !LOCAL VARIABLES:
      !------------------------------------------------------------------------
  
      call params_inst%cleanParams()
      ! deallocate(this%ac_phs_patch      )
   
   end subroutine Clean

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

   !-----------------------------------------------------------------------
   subroutine InitCold(this, bounds)
     !
     ! !ARGUMENTS:
     class(symbiont_type) :: this
     type(bounds_type), intent(in) :: bounds
     !
     ! !LOCAL VARIABLES:
     integer :: p,l                        ! indices
     !-----------------------------------------------------------------------
  
     do p = bounds%begp,bounds%endp
        l = patch%landunit(p)
  
        this%alphapsnsun_patch(p) = spval
        this%alphapsnsha_patch(p) = spval
  
        if (lun%ifspecial(l)) then
           this%psnsun_patch(p) = 0._r8
           this%psnsha_patch(p) = 0._r8
           if ( use_c13 ) then
              this%c13_psnsun_patch(p) = 0._r8
              this%c13_psnsha_patch(p) = 0._r8
           endif
           if ( use_c14 ) then
              this%c14_psnsun_patch(p) = 0._r8
              this%c14_psnsha_patch(p) = 0._r8
           endif
        end if
     end do

   end subroutine InitCold

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
   !-----------------------------------------------------------------------

   subroutine init_mimicsplus_veg_myc(bounds)

   ! !DESCRIPTION:

   ! !USES:
    use clm_varcon, only: pct_to_frac
   !
   ! !ARGUMENTS:
   type(bounds_type)               , intent(in)    :: bounds 
   !
   
   ! !LOCAL VARIABLES

   

   !ECW get parameters from the mimics module
      associate(      
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
      sulman_k_scav_Ninorg  => params_inst%sulman_k_scav_Ninorg    , &   !Half-saturation inorganic N concentration for mycorrhizal uptake [kg N m-3]
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
      sulman_rup_veg       => params_inst%sulman_rup_veg           &   !Vegetation N uptake rate from intermediate N pool [year-1]

      )
      end associate


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


   !-------------------  list of pools and their attributes  ------------

   ! Miner 
   ! Intermediate miner 

   ! Scavanger 
   ! Intermediate scavanger 

   ! Root uptake 

   ! Fixer 
   ! Intermediate fixer 

   !----------------  list of transitions and their time-independent coefficients  ---------------!

   

   end subroutine init_mimicsplus_veg_myc

   ! Step 1:
   ! Calculate Nstress | from CTSM?
   ! Plant tissues have a fixed C:N ratio, but vary by PFT
   ! If N is limited = biomass is limited to avaliable N and C is left in pool

   subroutine CN_soil_veg_exchange (cnveg_nitrogenstate_inst, leaf_prof_patch, froot_prof_patch, croot_prof_patch)
   
   ! !DESCRIPTION:
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   type(cnveg_nitrogenstate_type)       , intent(in)    :: cnveg_nitrogenstate_inst
   real(r8) :: N_stress                               ! N demand of plant, based on current N amount in plant []
   real(r8), parameter :: N_stress_max = 2.0_r8       ! Maximum N demand of plant, based on current N amount in plant []

   real(r8)                             , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: croot_prof_patch(bounds%begp:,1:) 
   
   associate(  
   leafn                               => cnveg_nitrogenstate_inst%leafn_patch                  , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
   leafn_storage                       => cnveg_nitrogenstate_inst%leafn_storage_patch          , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
   frootn                              => cnveg_nitrogenstate_inst%frootn_patch                 , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N                               
   frootn_storage                      => cnveg_nitrogenstate_inst%frootn_storage_patch         , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
   livecrootn                          => cnveg_nitrogenstate_inst%livecrootn_patch             , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
   livecrootn_storage                  => cnveg_nitrogenstate_inst%livecrootn_storage_patch     , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
   mimicsplus_fnalloc                  => params_inst%mimicsplus_fnalloc                          & !Fraction of NPP allocated to N uptake per unit N stress [fraction]
      
   )

   !N_stress = (2(Nleaf + Nroot)-Nstorage) / (Nleaf + Nroot)
   N_stress = (N_stress_max*(leafn + frootn + livecrootn) - (leafn_storage + frootn_storage + livecrootn_storage)) / &
               (leafn + frootn + livecrootn)

   
   ! Calculate the amount of C transferent to Naquisation
   ! Ctransfer = max(NPP,0)fNallocNstress
   real(r8) :: C_transfer                                   ! Carbon allocated to N acquisition, higher N_stress leads to higher C_transfer

   ! Dynamic allocation of a fraction of NPP to root exudation
   C_transfer = max(NPP,0) * mimicsplus_fnalloc * N_stress
      

   end associate
   end subroutine CN_soil_veg_exchange

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
      real(r8), intent(in)     :: tau_symb          ! Turnover time of symbionet (mimicsplus_tau_sym, mimicsplus_tau_fix)
      real(r8), intent(in)     :: growth_eff        ! Growth efficiency of symbionts
      real(r8), intent(in)     :: biomass_symb      ! Symbiont biomass
      real(r8), intent(inout)  :: roi               ! RoI nitrogen per carbon invested from vegetation to mycorrhiza [gN/gC]
      real(r8), intent(inout)  :: f_alloc           ! Fraction allocated to uptake strategy

      ! This might need to be limited if there is little N, but not sure
      roi = (plantN_uptake * tau_symb * growth_eff) / biomass_symb

   end subroutine calc_roi

   ! For symbionts mimer & scavenger
   call calc_roi(plantN_uptake, mimicsplus_tau_sym, mimicsplus_mge_ecm, biomass_symb, roi, f_alloc)   ! might need renaming
   call calc_roi(plantN_uptake, mimicsplus_tau_sym, mimicsplus_mge_am, biomass_symb, roi, f_alloc)

   ! For fixers
   call calc_roi(plantN_uptake, mimicsplus_tau_fix, mimicsplus_growth_fix, biomass_symb, roi, f_alloc)
   
   !----------------------------------------

   ! Include Nitrogen limitation from FUN

   !----------------------------------------

   ! PLANT UPTAKE STRATEGIES

   ! make functions for Nj & Bj to calculate ROI

   ! Direct N Uptake by Roots
   real(r8), pointer :: Npassive_patch                            (:)     ! N acquired by passive uptake      (gN/m2/s)
                                                                          ! in Sulman this is in (kgN/m2/s)

   do j = 1, nlevdecomp
      Npassive_patch = root water uptake flux * 
   end do 
         

   ! Inorganic N Scaveging by Mycorrhizal Fungi

   ! Organic N Mining by Mycorrhizal Fungi
        
   ! Symbiotic N2 Fixation

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
           


           




    

   ! Carbon intermediate pools
   fixer_inst%C_inter_biomass
   mine_inst%C_inter_biomass
   scav_inst%C_inter_biomass

   ! Nitrogen intermediate pools
   fixer_inst%N_inter_biomass
   miner_inst%N_inter_biomass
   scav_inst%N_inter_biomass

   end associate

end module CNSoilVegMIMICSplus
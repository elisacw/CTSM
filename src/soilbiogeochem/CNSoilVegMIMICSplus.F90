module CNSoilVegMIMICSplus

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module connects the Soil decomposition module MIMICS+ (Aas et al. 2023) with the vegetation through the 
  ! symbiosis between mycorrhizal fungi and plants.
  ! Coupling follows Sulman et al. (2019)
  
  ! !USES:
  use shr_kind_mod                        , only : r8 => shr_kind_r8
  use shr_infnan_mod                      , only : isnan => shr_infnan_isnan
  use clm_time_manager                    , only : get_step_size_real
  use clm_varpar                          , only : nlevdecomp
  use spmdMod                             , only : masterproc
  use abortutils                          , only : endrun
  use clm_varctl                          , only : iulog
  use clm_varpar                          , only : i_litr_min, i_litr_max, i_cwd
  use clm_varpar                          , only : i_met_lit, i_str_lit, i_phys_som, i_chem_som, i_avl_som
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
  use SoilBiogeochemDecompCascadeMIMICSMod, only : params_inst
  use WaterStateType                      , only : waterstate_type
  use SoilStateType                       , only : soilstate_type
  use WaterStateBulkType                  , only : waterstatebulk_type
  use WaterFluxBulkType                   , only : waterfluxbulk_type
  use TemperatureType                     , only : temperature_type
  use pftconMod                           , only : pftcon, noveg
  use abortutils                          , only : endrun
  use shr_log_mod                         , only : errMsg => shr_log_errMsg
  use ColumnType                          , only : col
  use GridcellType                        , only : grc
 
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  ! !PRIVATE MEMBER FUNCTIONS:
  public  :: CN_soil_veg_exchange         ! Calculates plant nitrogen stress, carbon transfer, growth & turnover of symbionts, passive root uptake and calls uptake routions
  private :: active_root_N_uptake         ! Nitrogen uptake strategy for active root uptake
  private :: myc_scavenger_N_uptake       ! Nitrogen uptake strategy for scavengers
  private :: myc_miner_N_uptake           ! Nitrogen uptake strategy for miners
  private :: potential_mined_n            ! Helper function for miner_nuptake and leftover_n_mining
  private :: miner_nuptake                ! Function to calculate mining from soil organic matter pools
  private :: leftover_n_mining            ! Function to calculate leftover N after applying NUE (this N is send to SOM pool)
  public  :: roi_symbionts                ! Calculates Return Of Investement, N to plant, C allocation to symbionts

  real(r8) :: dt                          ! decomp timestep [s]

  ! !FUNCTIONS:
  private :: resp_myc                     ! Respiration of mycorrhiza
  private :: Vmax_myc                     ! Michaelis Menten Kinetics for mycorrhiza [s-1]

  ! ! PUBLIC DATA
  integer, public :: i_fixer = 1 
  integer, public :: i_scav  = 2 
  integer, public :: i_miner = 3
  integer, public :: n_symb  = 3

  integer, public :: som_c   = 1 
  integer, public :: som_a   = 2 
  integer, public :: som_p   = 3
  integer, public :: n_som   = 3
  
  type, public :: symbiont_type

  real(r8), pointer           :: C_biomass             (:,:) ! [patch,n_symb] Carbon biomass                      [gC/m2]
  real(r8), pointer           :: N_biomass             (:,:) ! [patch,n_symb] Nitrogen biomass                    [gN/m2]
  real(r8), pointer           :: C_reservoir           (:,:) ! [patch,n_symb] Carbon intermediate pool biomass    [gC/m2]
  real(r8), pointer           :: N_reservoir           (:,:) ! [patch,n_symb] Nitrogen intermediate pool biomass  [gN/m2]

  real(r8), pointer           :: N_stress                (:) ! N demand of plant, based on current N amount in plant  [-]
  real(r8), pointer           :: C_allocation_to_N_acq   (:) ! [patch] C allocated to N acquisition                     [gC/m2/s]
  
  real(r8), pointer           :: symb_eff              (:,:) ! [patch,n_symb] Symbiont efficiency when biomass 0  [gC/gN]
  real(r8), pointer           :: symb_growth           (:,:) ! [patch,n_symb] Symbiotic biomass growth rate       [gC/m2/s]
  real(r8), pointer           :: N_symb_up             (:,:) ! Symbiont nitrogen uptake                    [gN/m3/s], fixers: [gN/m2/s]
  real(r8), pointer           :: N_to_plant            (:,:) ! Nitrogen send to plant                      [gN/m2/s]
  real(r8), pointer           :: C_alloc               (:,:) ! Carbon allocation from plant to symbiont    [gC/m2/s]
  real(r8), pointer           :: C_mortality           (:,:) ! [col,nlevdecomp]Turnover of symbionts per layer and column [gC/m3/s]
  real(r8), pointer           :: N_mortality           (:,:) ! Turnover of symbionts per layer and column                 [gN/m3/s]

  real(r8), pointer           :: total_symbiont_turnover_C  (:,:) ! Turnover of symbionts (necromass and maintanance respiration) [gC/m2/s]
  real(r8), pointer           :: total_symbiont_turnover_N  (:,:) ! Turnover of symbionts (necromass and maintanance respiration) [gN/m2/s]

  real(r8), pointer           :: no3_passiv_up           (:,:) ! [patch,nlevdecomp] Passive root NO3 (nitrate) uptake  [gN/m3/s]
  real(r8), pointer           :: nh4_passiv_up           (:,:) ! [patch,nlevdecomp] Passive root NH4 (ammonium) uptake [gN/m3/s]
  real(r8), pointer           :: no3_active_up           (:,:) ! [patch,nlevdecomp] Active root NO3 (nitrate) uptake   [gN/m3/s]
  real(r8), pointer           :: nh4_active_up           (:,:) ! [patch,nlevdecomp] Active root NH4 (ammonium) uptake  [gN/m3/s]
  real(r8), pointer           :: no3_scav_up             (:,:) ! [patch,nlevdecomp] Scavenger NO3 (nitrate) uptake     [gN/m3/s]
  real(r8), pointer           :: nh4_scav_up             (:,:) ! [patch,nlevdecomp] Scavenger NH4 (ammonium) uptake    [gN/m3/s]

  real(r8), pointer           :: N_mine_somc2soma_col  (:,:) ! Leftover part of co-mineralized N, not taken up by miners  [gN/m3/s]
  real(r8), pointer           :: N_mine_somp2soma_col  (:,:) ! Leftover part of co-mineralized N, not taken up by miners  [gN/m3/s]
  real(r8), pointer           :: somc_nuptake_col      (:,:) ! Nitrogen uptake from SOMc pool per column  [gN/m3/s]
  real(r8), pointer           :: somp_nuptake_col      (:,:) ! Nitrogen uptake from SOMp pool per column  [gN/m3/s]
  real(r8), pointer           :: somc_cuptake_col      (:,:) ! Carbon uptake from SOMc pool per column  [gC/m3/s]
  real(r8), pointer           :: somp_cuptake_col      (:,:) ! Carbon uptake from SOMp pool per column  [gC/m3/s]
  real(r8), pointer           :: root_exudate_C_col    (:,:) ! Leftover C from allocation to symbionts    [gC/m3/s]
  logical, pointer            :: is_active             (:,:) ! If the symbiont type is active             [-]
  character(len=10), pointer  :: symb_name             (:)   ! Symbiont name                              [-]
  character(len=10), pointer  :: symb_hist_name        (:)   ! Symbiont name on history tapes             [-]
   
contains

   ! Public procedures
   procedure, public  :: Init
   procedure, public  :: Restart
   procedure, public  :: InitCold
   procedure, public  :: Summary
   procedure, public  :: InitAllocate
   procedure, public  :: InitHistory

  end type symbiont_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(symbiont_type)          :: this
    type(bounds_type), intent(in) :: bounds

   ! !LOCAL variables
    integer :: iveg
    integer :: p,i                        ! indices

    logical :: local_active(bounds%begp:bounds%endp, 1:n_symb)   ! Indicates active pathway, based on symbiont type of PFT
   
    local_active(:,:) = .false. 

    allocate(this%symb_name(n_symb))
    allocate(this%symb_hist_name(n_symb))

    this%symb_name(i_fixer)      = 'fixer'
    this%symb_name(i_scav)       = 'scavanger'
    this%symb_name(i_miner)      = 'miner'
    this%symb_hist_name(i_fixer) = 'FIX'
    this%symb_hist_name(i_scav)  = 'SCAV'
    this%symb_hist_name(i_miner) = 'MINE'
    
    do p = bounds%begp,bounds%endp
      iveg = patch%itype(p)
      if (iveg /= noveg) then 
        if (pftcon%FUN_fracfixers(iveg) > 0.0_r8) then  ! PFTs have a "fixer fraction" that determines the fraction of C that can be used for fixation.     
         local_active(p,i_fixer) = .true.               ! FUN_fracfixers is a placeholder
        endif
        if (pftcon%myc_symbiont(iveg) == 1.0_r8) then
         local_active(p,i_miner) = .true.
         local_active(p,i_scav) = .false.
        elseif (pftcon%myc_symbiont(iveg) == 0.0_r8) then
         local_active(p,i_miner) = .false.
         local_active(p,i_scav) = .true.
        else
         local_active(p,i_miner) = .true.
         local_active(p,i_scav) = .true.
        endif
      endif
   enddo 
   allocate(this%is_active(bounds%begp:bounds%endp,1:n_symb)) ; this%is_active(bounds%begp:bounds%endp,1:n_symb)=.false.
  
   this%is_active(bounds%begp:bounds%endp,1:n_symb)=local_active(bounds%begp:bounds%endp,1:n_symb)
    call this%InitAllocate (bounds)
    call this%InitHistory (bounds)
    call this%InitCold (bounds)

  end subroutine Init

  !------------------------------------------------------------------------

  subroutine InitAllocate(this, bounds)
   ! ! USES:
   !
   ! !ARGUMENTS:
   class(symbiont_type)          :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:
   integer :: begp, endp, begc, endc
  
   begp = bounds%begp; endp= bounds%endp
   begc = bounds%begc; endc= bounds%endc

    allocate(this%C_biomass(begp:endp,1:n_symb)) ; this%C_biomass(begp:endp,1:n_symb)     = 0.0_r8
    allocate(this%N_biomass(begp:endp,1:n_symb)) ; this%N_biomass(begp:endp,1:n_symb)     = 0.0_r8
    allocate(this%C_reservoir(begp:endp,1:n_symb)) ; this%C_reservoir(begp:endp,1:n_symb) = 0.0_r8
    allocate(this%N_reservoir(begp:endp,1:n_symb)) ; this%N_reservoir(begp:endp,1:n_symb) = 0.0_r8

    allocate(this%N_stress(begp:endp))             ; this%N_stress(begp:endp)             = 0.0_r8
    allocate(this%C_allocation_to_N_acq(begp:endp)) ; this%C_allocation_to_N_acq(begp:endp) = 0.0_r8

    allocate(this%symb_eff(begp:endp,1:n_symb)) ; this%symb_eff(begp:endp,1:n_symb)       = 0.0_r8
    allocate(this%symb_growth(begp:endp,1:n_symb)) ; this%symb_growth(begp:endp,1:n_symb) = 0.0_r8
    allocate(this%N_symb_up(begp:endp,1:n_symb)) ; this%N_symb_up(begp:endp,1:n_symb)     = 0.0_r8
    allocate(this%N_to_plant(begp:endp,1:n_symb)) ; this%N_to_plant(begp:endp,1:n_symb)   = 0.0_r8
    allocate(this%C_alloc(begp:endp,1:n_symb))   ; this%C_alloc(begp:endp,1:n_symb)       = 0.0_r8

    allocate(this%total_symbiont_turnover_C(begp:endp,1:n_symb))   ; this%total_symbiont_turnover_C(begp:endp,1:n_symb)       = 0.0_r8
    allocate(this%total_symbiont_turnover_N(begp:endp,1:n_symb))   ; this%total_symbiont_turnover_N(begp:endp,1:n_symb)       = 0.0_r8

    allocate(this%root_exudate_C_col(begc:endc,1:nlevdecomp)) ; this%root_exudate_C_col(begc:endc,1:nlevdecomp) = 0.0_r8

    allocate(this%C_mortality(begc:endc,1:nlevdecomp)) ; this%C_mortality(begc:endc,1:nlevdecomp) = 0.0_r8
    allocate(this%N_mortality(begc:endc,1:nlevdecomp)) ; this%N_mortality(begc:endc,1:nlevdecomp) = 0.0_r8

    allocate(this%no3_passiv_up(begp:endp,1:nlevdecomp))    ; this%no3_passiv_up(begp:endp,1:nlevdecomp)  = 0.0_r8
    allocate(this%nh4_passiv_up(begp:endp,1:nlevdecomp))    ; this%nh4_passiv_up(begp:endp,1:nlevdecomp)     = 0.0_r8
    allocate(this%no3_active_up(begp:endp,1:nlevdecomp))    ; this%no3_active_up(begp:endp,1:nlevdecomp)     = 0.0_r8
    allocate(this%nh4_active_up(begp:endp,1:nlevdecomp))    ; this%nh4_active_up(begp:endp,1:nlevdecomp)     = 0.0_r8
    allocate(this%no3_scav_up(begp:endp,1:nlevdecomp))      ; this%no3_scav_up(begp:endp,1:nlevdecomp)       = 0.0_r8
    allocate(this%nh4_scav_up(begp:endp,1:nlevdecomp))      ; this%nh4_scav_up(begp:endp,1:nlevdecomp)       = 0.0_r8

    allocate(this%N_mine_somc2soma_col(begc:endc,1:nlevdecomp)) ; this%N_mine_somc2soma_col(begc:endc,1:nlevdecomp) = 0.0_r8
    allocate(this%N_mine_somp2soma_col(begc:endc,1:nlevdecomp)) ; this%N_mine_somp2soma_col(begc:endc,1:nlevdecomp) = 0.0_r8
    allocate(this%somc_nuptake_col(begc:endc,1:nlevdecomp)) ; this%somc_nuptake_col(begc:endc,1:nlevdecomp) = 0.0_r8
    allocate(this%somp_nuptake_col(begc:endc,1:nlevdecomp)) ; this%somp_nuptake_col(begc:endc,1:nlevdecomp) = 0.0_r8
    allocate(this%somc_cuptake_col(begc:endc,1:nlevdecomp)) ; this%somc_cuptake_col(begc:endc,1:nlevdecomp) = 0.0_r8
    allocate(this%somp_cuptake_col(begc:endc,1:nlevdecomp)) ; this%somp_cuptake_col(begc:endc,1:nlevdecomp) = 0.0_r8
   
   end subroutine InitAllocate

   !-----------------------------------------------------------------------

   subroutine InitHistory(this, bounds)
      !
      ! !USES:
      use histFileMod                , only: hist_addfld1d, hist_addfld2d
      use clm_varcon                 , only: spval
      !
      ! !ARGUMENTS:
      class(symbiont_type)          :: this
      type(bounds_type), intent(in) :: bounds
      !
      ! !LOCAL VARIABLES:
      integer :: begp, endp, begc, endc
      integer :: i, j
      real(r8), pointer :: data1dptr (:)
      real(r8), pointer :: data2dptr(:,:)
     
      begp = bounds%begp; endp= bounds%endp
      begc = bounds%begc; endc= bounds%endc

      do i = 1,n_symb 
         this%C_biomass(begp:endp,i) = spval
         data1dptr => this%C_biomass(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_C', units='gC/m2', &
         avgflag='A', long_name=('Carbon pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='active')
   
         this%C_reservoir(begp:endp,i) = spval
         data1dptr => this%C_reservoir(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_C_INTER', units='gC/m2', &
         avgflag='A', long_name=('Carbon intermediate pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%N_biomass(begp:endp,i) = spval
         data1dptr => this%N_biomass(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_N', units='gN/m2', &
         avgflag='A', long_name=('Nitrogen pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='active')
   
         this%N_reservoir(begp:endp,i) = spval
         data1dptr => this%N_reservoir(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_N_INTER', units='gN/m2', &
         avgflag='A', long_name=('Nitrogen intermediate pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%symb_growth(begp:endp,i) = spval
         data1dptr => this%symb_growth(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_C_GROWTH', units='gC/m2/s', &
         avgflag='A', long_name=('Symbiont growth of '//this%symb_name(i)), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%N_symb_up(begp:endp,i) = spval
         data1dptr => this%N_symb_up(:,i)
         call hist_addfld1d (fname=trim('N_SYM_'//this%symb_hist_name(i))//'_UPTAKE', units='gN/m3/s', &
         avgflag='A', long_name=('Nitrogen uptake of symbiont '//this%symb_name(i)), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%N_to_plant(begp:endp,i) = spval
         data1dptr => this%N_to_plant(:,i)
         call hist_addfld1d (fname=trim('N_FROM_'//this%symb_hist_name(i))//'_TO_PLANT', units='gN/m2/s', &
         avgflag='A', long_name=('Nitrogen from symbiont send to plant '//this%symb_name(i)), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%C_alloc(begp:endp,i) = spval
         data1dptr => this%C_alloc(:,i)
         call hist_addfld1d (fname=trim('C_ALLOC_TO_'//this%symb_hist_name(i))//'_FROM_PLANT', units='gC/m2/s', &
         avgflag='A', long_name=('Plant carbon allocation to symbiont '//this%symb_name(i)), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%total_symbiont_turnover_C(begp:endp,i) = spval
         data1dptr => this%total_symbiont_turnover_C(:,i)
         call hist_addfld1d (fname=trim('C_'//this%symb_hist_name(i))//'_TURNOVER', units='gC/m2/s', &
         avgflag='A', long_name=('Turnover (necromass & maintainance respiration) of symbionts '//this%symb_name(i)), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%total_symbiont_turnover_N(begp:endp,i) = spval
         data1dptr => this%total_symbiont_turnover_N(:,i)
         call hist_addfld1d (fname=trim('N_'//this%symb_hist_name(i))//'_TURNOVER', units='gN/m2/s', &
         avgflag='A', long_name=('Turnover (necromass & maintainance respiration) of symbionts '//this%symb_name(i)), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

      end do

         this%C_mortality(begc:endc,1:nlevdecomp) = spval
         data2dptr => this%C_mortality
         call hist_addfld2d (fname='C_MORTALITY_SYMB_MIMICSPLUS', units='gC/m^3/s', type2d='levsoi', &
         avgflag='A', long_name='Symbiotic C turnover per soil layer and column', &
         ptr_col= this%C_mortality, set_spec=spval, default='inactive')
   
         this%N_mortality(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MORTALITY_SYMB_MIMICSPLUS', units='gN/m^3/s', type2d='levsoi', &
         avgflag='A', long_name='Symbiotic N turnover per soil layer and column', &
         ptr_col=this%N_mortality, set_spec=spval, default='inactive')

         this%no3_passiv_up(begp:endp,1:nlevdecomp) = spval
         call hist_addfld2d (fname='NO3_PASSIVE_UP', units='gN/m^3/s', type2d='levsoi', &
         avgflag='A', long_name='Passive root NO3 (nitrate) uptake', &
         ptr_col=this%no3_passiv_up, set_spec=spval, default='inactive')

         this%nh4_passiv_up(begp:endp,1:nlevdecomp) = spval
         call hist_addfld2d (fname='NH4_PASSIVE_UP', units='gN/m^3/s', type2d='levsoi', &
         avgflag='A', long_name='Passive root NH4 (ammonium) uptake', &
         ptr_col=this%nh4_passiv_up, set_spec=spval, default='inactive')
       
         this%no3_active_up(begp:endp,1:nlevdecomp) = spval
         call hist_addfld2d (fname='NO3_ACTIVE_UP', units='gN/m^3/s', type2d='levsoi', &
         avgflag='A', long_name='Active root NO3 (nitrate) uptake', &
         ptr_col=this%no3_active_up, set_spec=spval, default='inactive')
         
         this%nh4_active_up(begp:endp,1:nlevdecomp) = spval
         call hist_addfld2d (fname='NH4_ACTIVE_UP', units='gN/m^3/s', type2d='levsoi', &
         avgflag='A', long_name='Active root NH4 (ammonium) uptake', &
         ptr_col=this%nh4_active_up, set_spec=spval, default='inactive')
       
         this%no3_scav_up  (begp:endp,1:nlevdecomp) = spval
         call hist_addfld2d (fname='NO3_SCAV_UP', units='gN/m^3/s', type2d='levsoi', &
         avgflag='A', long_name='Scavenger NO3 (nitrate) uptake', &
         ptr_col=this%no3_scav_up, set_spec=spval, default='inactive')
       
         this%nh4_scav_up  (begp:endp,1:nlevdecomp) = spval
         call hist_addfld2d (fname='NH4_SCAV_UP', units='gN/m^3/s', type2d='levsoi', &
         avgflag='A', long_name='Scavenger NH4 (ammonium) uptake', &
         ptr_col=this%nh4_scav_up, set_spec=spval, default='inactive')
       
         this%N_mine_somc2soma_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MINE_SOMC_TO_SOMA_MIMICSPLUS', units='gN/m^3/s', type2d='levsoi',&
         avgflag='A', long_name='Leftover N from SOMc mineralization (not taken up by miners)', &
         ptr_col=this%N_mine_somc2soma_col, set_spec=spval, default='inactive')
   
         this%N_mine_somp2soma_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MINE_SOMP_TO_SOMA_MIMICSPLUS', units='gN/m3/s', type2d='levsoi', &
         avgflag='A', long_name='Leftover N from SOMp mineralization (not taken up by miners)', &
         ptr_col=this%N_mine_somp2soma_col, set_spec=spval, default='inactive')
   
         this%somc_nuptake_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MINE_UPTAKE_SOMC_MIMICSPLUS', units='gN/m3/s', type2d='levsoi', &
         avgflag='A', long_name='N uptake from SOMc pool via mining', &
         ptr_col=this%somc_nuptake_col, set_spec=spval, default='inactive')
   
         this%somp_nuptake_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MINE_UPTAKE_SOMP_MIMICSPLUS', units='gN/m3/s', type2d='levsoi',&
         avgflag='A', long_name='N uptake from SOMp pool via mining', &
         ptr_col=this%somp_nuptake_col, set_spec=spval, default='inactive')
   
         this%somc_cuptake_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='C_MINE_UPTAKE_SOMC_MIMICSPLUS', units='gC/m3/s',type2d='levsoi', &
         avgflag='A', long_name='co-decomposed C from SOMc pool during mining', &
         ptr_col=this%somc_cuptake_col, set_spec=spval, default='inactive')
   
         this%somp_cuptake_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='C_MINE_UPTAKE_SOMP_MIMICSPLUS', units='gC/m3/s', type2d='levsoi', &
         avgflag='A', long_name='co-decomposed C from SOMp pool during mining', &
         ptr_col=this%somp_cuptake_col, set_spec=spval, default='inactive')
   
         this%root_exudate_C_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='ROOT_EXUDATE_C_MIMICSPLUS', units='gC/m3/s', type2d='levsoi', &
         avgflag='A', long_name='Leftover root exudate C from symbiont allocation', &
         ptr_col=this%root_exudate_C_col, set_spec=spval, default='inactive')

         this%N_stress(begp:endp) = spval
         call hist_addfld1d (fname='N_STRESS', units='-', &
         avgflag='A', long_name='N stress of plant', &
         ptr_patch=this%N_stress)

         this%C_allocation_to_N_acq(begp:endp) = spval
         call hist_addfld1d (fname='C_ALLOC_TO_N_ACQ', units='-', &
         avgflag='A', long_name='C allocated from plant to recive N', &
         ptr_patch=this%C_allocation_to_N_acq)


   end subroutine InitHistory

   !-----------------------------------------------------------------------

   subroutine InitCold(this, bounds)
      ! !USES:
      !
      ! !ARGUMENTS:
      class(symbiont_type)          :: this
      type(bounds_type), intent(in) :: bounds
      !
      ! !LOCAL VARIABLES:
      character(len=128)            :: varname   ! temporary
      integer :: p,i,j                           ! indices

      ! Initial carbon stock in symbiont pools for coldstart (from parameter file ) are separated into:
      ! 90 % = initial C, 10% = initial intermediate C
      do i = 1,n_symb
         do p = bounds%begp,bounds%endp
            this%symb_eff(p,i) = 0.0001_r8 ! start with an arbitrary number
            if (i == i_fixer) then
               this%symb_eff(p,i) = params_inst%sulman_rfix
            end if
            if (i == i_scav) then
               this%symb_eff(p,i) = 1.0_r8
            endif
            if (i == i_miner) then
               this%symb_eff(p,i) = 1.0_r8
            endif
           
            if (this%is_active(p,i)) then
               this%C_biomass(p,i) = params_inst%sulman_initial_C_stocks(i)
               this%C_reservoir(p,i) = params_inst%sulman_initial_C_stocks(i)
               if (params_inst%sulman_cn_symbionts(i) == 0.0_r8) then
                  this%N_biomass(p,i) = 0.0_r8
                  this%N_reservoir(p,i) = 0.0_r8
               else
                  this%N_biomass(p,i) = this%C_biomass(p,i) / params_inst%sulman_cn_symbionts(i)
                  this%N_reservoir(p,i) = this%C_reservoir(p,i) / params_inst%sulman_cn_symbionts(i)
               endif
            else
               this%C_biomass(p,i)     = 0.0_r8
               this%C_reservoir(p,i)   = 0.0_r8
               this%N_biomass(p,i)     = 0.0_r8
               this%N_reservoir(p,i)   = 0.0_r8
            endif
         enddo
      enddo

      this%C_mortality(bounds%begc:bounds%endc,1:nlevdecomp)            = 0.0_r8
      this%N_mortality(bounds%begc:bounds%endc,1:nlevdecomp)            = 0.0_r8
      this%N_mine_somc2soma_col(bounds%begc:bounds%endc,1:nlevdecomp)   = 0.0_r8
      this%N_mine_somp2soma_col(bounds%begc:bounds%endc,1:nlevdecomp)   = 0.0_r8
      this%somc_nuptake_col(bounds%begc:bounds%endc,1:nlevdecomp)       = 0.0_r8
      this%somp_nuptake_col(bounds%begc:bounds%endc,1:nlevdecomp)       = 0.0_r8
      this%somc_cuptake_col(bounds%begc:bounds%endc,1:nlevdecomp)       = 0.0_r8
      this%somp_cuptake_col(bounds%begc:bounds%endc,1:nlevdecomp)       = 0.0_r8
      this%root_exudate_C_col(bounds%begc:bounds%endc,1:nlevdecomp)     = 0.0_r8

    end subroutine InitCold

   !------------------------------------------------------------------------

   subroutine Restart(this, bounds, ncid, flag)

     !
     ! !USES:
     use restUtilMod
     use ncdio_pio                     , only: file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
     !
     ! !ARGUMENTS:
     class(symbiont_type)              :: this
     type(bounds_type), intent(in)     :: bounds
     type(file_desc_t), intent(inout)  :: ncid         ! netcdf id
     character(len=*) , intent(in)     :: flag         ! 'read' or 'write'
     !
     ! !LOCAL VARIABLES:
     integer                           :: i            ! indices
     logical                           :: readvar      ! determine if variable is on initial file
     character(len=128)                :: varname      ! temporary
     real(r8), pointer                 :: data1dptr (:)
     real(r8), pointer                 :: data2dptr (:,:)
     
     !-----------------------------------------------------------------------
   do i = 1,n_symb
      data1dptr => this%C_biomass(:,i)
      call restartvar(ncid=ncid, flag=flag, varname=trim('c_pool_'//this%symb_name(i)), xtype=ncd_double,  &
      dim1name='pft', long_name=('Carbon pool of '//this%symb_name(i)//' symbionts'), units='g/m2', &
      interpinic_flag='interp', readvar=readvar, data=data1dptr)
      if (flag == 'read' .and. (.not. readvar)) then
         call endrun(msg = "ERROR: Symbiont "//this%symb_name(i)// " is not on the restart file."// & 
             errMsg(sourcefile, __LINE__))
      endif
      data1dptr => this%C_reservoir(:,i)
      call restartvar(ncid=ncid, flag=flag, varname=trim('c_int_pool_'//this%symb_name(i)), xtype=ncd_double,  &
      dim1name='pft', long_name=('Carbon intermediate pool of '//this%symb_name(i)//' symbionts'), units='g/m2', &
      interpinic_flag='interp', readvar=readvar, data=data1dptr)
      if (flag == 'read' .and. (.not. readvar)) then
         call endrun(msg = "ERROR: Symbiont "//this%symb_name(i)// " is not on the restart file."// & 
             errMsg(sourcefile, __LINE__))
      endif
      data1dptr => this%N_biomass(:,i)
      call restartvar(ncid=ncid, flag=flag, varname=trim('n_pool_'//this%symb_name(i)), xtype=ncd_double,  &
      dim1name='pft', long_name=('Nitrogen pool of '//this%symb_name(i)//' symbionts'), units='g/m2', &
      interpinic_flag='interp', readvar=readvar, data=data1dptr)
      if (flag == 'read' .and. (.not. readvar)) then
         call endrun(msg = "ERROR: Symbiont "//this%symb_name(i)// " is not on the restart file."// & 
             errMsg(sourcefile, __LINE__))
      endif
      data1dptr => this%N_reservoir(:,i)
      call restartvar(ncid=ncid, flag=flag, varname=trim('n_int_pool_'//this%symb_name(i)), xtype=ncd_double,  &
      dim1name='pft', long_name=('Nitrogen intermediate pool of '//this%symb_name(i)//' symbionts'), units='g/m2', &
      interpinic_flag='interp', readvar=readvar, data=data1dptr)
      if (flag == 'read' .and. (.not. readvar)) then
         call endrun(msg = "ERROR: Symbiont "//this%symb_name(i)// " is not on the restart file."// & 
             errMsg(sourcefile, __LINE__))
      endif

   end do 

   end subroutine Restart

   
   subroutine Summary(this, bounds, filter_soilp, num_soilp, filter_bgc_soilc, num_bgc_soilc, comp, total)
      ! Summing up all symbiont and intermediate pools for CN Balance check

      ! !USES:
      use subgridAveMod                      , only: p2c
      !
      ! !ARGUMENTS:
      class(symbiont_type)             :: this
      type(bounds_type), intent(in)    :: bounds
      integer                                , intent(in)    :: filter_soilp(:)     ! filter for soil patches
      integer                                , intent(in)    :: num_soilp           ! number of soil patches in filter
      integer                                , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
      integer                                , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
      character(len=*), intent(in)     :: comp     ! what to sum
      real(r8) ,intent(out)            :: total(bounds%begc:bounds%endc)    ! column level total output
      ! LOCAL VARIABLES
      integer :: p, fp, c, fc, j, k, l, s, i, g
      real(r8) :: total_patch(bounds%begp:bounds%endp) 
      real(r8) :: total_col(bounds%begc:bounds%endc)

      total_patch(bounds%begp:bounds%endp)   = 0.0_r8
      total_col(bounds%begc:bounds%endc)     = 0.0_r8


      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)

         select case( trim(comp) )
         case( 'c12' ) 
            total_patch(p) = this%C_biomass(p,i_scav) + this%C_biomass(p,i_miner) + this%C_biomass(p,i_fixer) &
                         + this%C_reservoir(p,i_scav) + this%C_reservoir(p,i_miner) + this%C_reservoir(p,i_fixer)          
         case( 'n' ) 
            total_patch(p) = this%N_biomass(p,i_scav) + this%N_biomass(p,i_miner) + this%N_biomass(p,i_fixer) &
                         + this%N_reservoir(p,i_scav) + this%N_reservoir(p,i_miner) + this%N_reservoir(p,i_fixer)

         case default
            call endrun('Bad compound name ='//comp )
         end select
       end do

      
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
       total_patch(bounds%begp:bounds%endp), &
       total_col(bounds%begc:bounds%endc))
       
      total(bounds%begc:bounds%endc)  = total_col(bounds%begc:bounds%endc) 
     
   end subroutine Summary
   !-----------------------------------------------------------------------

   subroutine CN_soil_veg_exchange (filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds, symbiont_inst, &
      cnveg_nitrogenstate_inst, waterstatebulk_inst, temperature_inst, cnveg_carbonflux_inst, &
      soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst, cnveg_state_inst, &
      waterfluxbulk_inst, soilstate_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst, cnveg_nitrogenflux_inst, soilbiogeochem_state_inst)

   ! !DESCRIPTION:
   !
   ! Calculates plant nitrogen stress, carbon transfer, growth & turnover of symbionts, passive root uptake and calls uptake routions
   !
   ! !USES:
   use clm_varcon                         , only: pct_to_frac
   use WaterFluxType                      , only: waterflux_type
   use SoilStateType                      , only: soilstate_type
   use TemperatureType                    , only: temperature_type
   use CNVegCarbonStateType               , only: cnveg_carbonstate_type
   use SoilBiogeochemCarbonStateType      , only: soilbiogeochem_carbonstate_type
   use SoilBiogeochemNitrogenStateType    , only: soilbiogeochem_nitrogenstate_type   
   use subgridAveMod                      , only: p2c
   use CNVegnitrogenfluxType              , only : cnveg_nitrogenflux_type
   use clm_varcon                         , only : smallValue
   use pftconMod                          , only : noveg
   !
   ! !ARGUMENTS:
   type(symbiont_type)                    , intent(inout) :: symbiont_inst
   type(bounds_type)                      , intent(in)    :: bounds 
   type(cnveg_nitrogenstate_type)         , intent(in)    :: cnveg_nitrogenstate_inst
   type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
   type(waterstatebulk_type)              , intent(in)    :: waterstatebulk_inst
   type(temperature_type)                 , intent(in)    :: temperature_inst
   type(cnveg_carbonflux_type)            , intent(in)    :: cnveg_carbonflux_inst
   type(soilbiogeochem_nitrogenstate_type), intent(in)    :: soilbiogeochem_nitrogenstate_inst
   type(waterfluxbulk_type)               , intent(in)    :: waterfluxbulk_inst
   type(soilstate_type)                   , intent(in)    :: soilstate_inst
   type(cnveg_carbonstate_type)           , intent(in)    :: cnveg_carbonstate_inst
   type(soilbiogeochem_carbonstate_type)  , intent(in)    :: soilbiogeochem_carbonstate_inst
   type(cnveg_nitrogenflux_type)          , intent(inout) :: cnveg_nitrogenflux_inst
   type(cnveg_state_type)                 , intent(in)    :: cnveg_state_inst
   type(soilbiogeochem_state_type)        , intent(in)    :: soilbiogeochem_state_inst

   integer                                , intent(in)    :: num_soilp           ! number of soil patches in filter
   integer                                , intent(in)    :: filter_soilp(:)     ! filter for soil patches
   integer                                , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
   integer                                , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
   
   ! !LOCAL VARIABLES
   integer :: p, fp, c, fc, j, k, l, s, i, g
   integer :: begp, endp, begc, endc
   
   real(r8), parameter :: N_stress_max = 2.0_r8             ! Maximum N demand of plant, based on current N amount in plant []
   real(r8), parameter :: N_stress_min = 0.05_r8            ! Miminmum value of N_stress
   real(r8), parameter :: sulman_fnalloc = 0.05       ! Fraction of NPP allocated to N uptake per unit N stress [fraction] 

   real(r8) :: root_dens_sum                                            ! Fine root C per layer             [gC/m2]
   real(r8) :: t_soi_degC                                               ! Soil temperature                  [degrees Celcius]
   real(r8) :: availc_alloc(bounds%begp:bounds%endp)                    ! The avaible C pool for allocation [gC/m2/s]
   real(r8) :: plantCN(bounds%begp:bounds%endp)                         ! plant C:N
   real(r8) :: npp_growth_potential(bounds%begp:bounds%endp) 
   real(r8) :: root_dens_frac(bounds%begp:bounds%endp,1:nlevdecomp)     ! Fraction of root density          [-]

   ! VARIABLES FOR LOCAL BALANCE CHECK
   real(r8) :: errbalc, errbaln                                         ! balance error
   real(r8) :: N_to_plant(bounds%begp:bounds%endp,1:n_symb)             ! Nitrogen allocation to plant                [gN/m2/s]
   real(r8) :: symb_growth_gross(bounds%begp:bounds%endp,1:n_symb)      ! Gross symbiotic growth (without  CUE)       [gC/m2/s]

   real(r8) :: old_C_biomass(bounds%begp:bounds%endp,1:n_symb)
   real(r8) :: old_N_biomass(bounds%begp:bounds%endp,1:n_symb)
   real(r8) :: old_C_reservoir(bounds%begp:bounds%endp,1:n_symb)
   real(r8) :: old_N_reservoir(bounds%begp:bounds%endp,1:n_symb)

   ! Nitrogen uptake variables for pathways into intermediated pools
   real(r8) :: smin_no3_avail(bounds%begp:bounds%endp, 1:nlevdecomp)     ! no3 available for uptake per soil layer     [gN/m3/s]
   real(r8) :: smin_nh4_avail(bounds%begp:bounds%endp, 1:nlevdecomp)     ! nh4 available for uptake per soil layer     [gN/m3/s]
   real(r8) :: smin_no3_avail_col(bounds%begc:bounds%endc, 1:nlevdecomp) ! col no3 available for uptake per soil layer [gN/m3/s]
   real(r8) :: smin_nh4_avail_col(bounds%begc:bounds%endc, 1:nlevdecomp) ! col nh4 available for uptake per soil layer [gN/m3/s]
      
   real(r8) :: sum_no3_up(bounds%begp:bounds%endp,1:nlevdecomp)      ! Scavenger NO3 (nitrate) uptake    [gN/m3/s]
   real(r8) :: sum_nh4_up(bounds%begp:bounds%endp,1:nlevdecomp)      ! Scavenger NH4 (ammonium) uptake    [gN/m3/s]
   
   real(r8) :: N_fixation(bounds%begp:bounds%endp)                   ! Nitrogen uptake from fixation      [gN/m2/s]
  
   real(r8) :: somc_nuptake(bounds%begp:bounds%endp, 1:nlevdecomp)   ! Nitrogen uptake from SOMc pool by miners   [gN/m3/s]
   real(r8) :: somp_nuptake(bounds%begp:bounds%endp, 1:nlevdecomp)   ! Nitrogen uptake from SOMp pool by miners   [gN/m3/s]
   real(r8) :: somc_cuptake(bounds%begp:bounds%endp, 1:nlevdecomp)   ! Carbon uptake from SOMc pool by miners     [gC/m3/s]
   real(r8) :: somp_cuptake(bounds%begp:bounds%endp, 1:nlevdecomp)   ! Carbon uptake from SOMp pool by miners     [gC/m3/s]
     
   real(r8) :: maint_resp_N(bounds%begp:bounds%endp, 1:n_symb)       ! N remaining in N pool from C maintainance respiration [gN/m2/s]
   real(r8) :: maint_resp                                            ! Maintainace respiration, C from symbiont pool used to sustain existing biomass [gC/m2/s]
   real(r8) :: growth_resp(bounds%begp:bounds%endp,1:n_symb)         ! Growth respiration [gC/m2/s]

   real(r8) :: total_symbiont_turnover_C(bounds%begp:bounds%endp, 1:n_symb)       ! Part of symbiont turnover going into SOM [gC/m2/s]
   real(r8) :: total_symbiont_turnover_N(bounds%begp:bounds%endp, 1:n_symb)       ! Part of symbiont turnover going into SOM [gN/m2/s]
   real(r8) :: symbiont_turnover_C_to_som(bounds%begp:bounds%endp, 1:n_symb)      ! Part of symbiont turnover going into SOM [gC/m2/s]
   real(r8) :: symbiont_turnover_N_to_som(bounds%begp:bounds%endp, 1:n_symb)      ! Part of symbiont turnover going into SOM [gN/m2/s]
   real(r8) :: symb_turnover_layer_C(bounds%begp:bounds%endp, 1:nlevdecomp)       ! Symbiotic turnover per soil layer        [gC/m3/s]
   real(r8) :: symb_turnover_layer_N(bounds%begp:bounds%endp, 1:nlevdecomp)       ! Symbiotic turnover per soil layer        [gN/m3/s]
   real(r8) :: N_mine_somc2soma(bounds%begp:bounds%endp, 1:nlevdecomp)            ! Leftover co-mineralized N, not taken up by miners  [gN/m3/s]
   real(r8) :: N_mine_somp2soma(bounds%begp:bounds%endp, 1:nlevdecomp)            ! Leftover co-mineralized N, not taken up by miners  [gN/m3/s]

   real(r8) :: root_exudate_C(bounds%begp:bounds%endp)                            ! Leftover C from allocation to symbionts    [gC/m2/s]
   real(r8) :: root_exudate_C_layer(bounds%begp:bounds%endp,1:nlevdecomp)         ! Leftover C from allocation to symbionts    [gC/m3/s]
   real(r8) :: root_N_active_uptake(bounds%begp:bounds%endp)                      ! active root N uptake [gN/m2/s]
   real(r8) :: root_N_to_plant(bounds%begp:bounds%endp)                           ! toatl (active + passive) root N uptake [gN/m2/s]
   real(r8) :: scale_N_to_plant(bounds%begp:bounds%endp)       ! Scale factor to scale N uptake to plant if it is bigger that the uptake capazitiy of plant

   real(r8) :: norm_froot_prof(bounds%begp:bounds%endp,1:nlevdecomp)    ! normalized fine root profile

   ! Nstress
   real(r8) :: N_demand(bounds%begp:bounds%endp)
   real(r8) :: total_N(bounds%begp:bounds%endp)
   real(r8) :: potential_stored_N(bounds%begp:bounds%endp)
   
   begp = bounds%begp; endp= bounds%endp
   dt   = get_step_size_real()
   
   associate(                                                   &
   sulman_root_no3      => params_inst%sulman_root_no3        , &   ! Maximum root active nitrate uptake rate                        [gN/m3/s]
   sulman_root_nh4      => params_inst%sulman_root_nh4        , &   ! Maximum root active ammonium uptake rate                       [gN/m3/s]
   sulman_km_no3        => params_inst%sulman_km_no3          , &   ! Half-saturation nitrate concentration for root active uptake     [gN/m3]
   sulman_km_nh4        => params_inst%sulman_km_nh4          , &   ! Half-saturation nitrate concentration for root active uptake     [gN/m3]
   sulman_r_rhiz        => params_inst%sulman_r_rhiz          , &   ! Radius of the rhizosphere                                            [m]
   sulman_v_scav        => params_inst%sulman_v_scav          , &   ! Maximum N uptake rate by scavenging mycorrhizae                [gN/m3/s]
   sulman_k_scav_Ninorg => params_inst%sulman_k_scav_Ninorg   , &   ! Half-saturation inorganic N concentration for mycorrhizal uptake [gN/m3]
   sulman_k_scav        => params_inst%sulman_k_scav          , &   ! Half-saturation mycorrhizal biomass concentration for scavenging [gC/m3]
   sulman_km_mine       => params_inst%sulman_km_mine         , &   ! Half-saturation mycorrhizal biomass concentration for mining     [gC miners/gC substrate]
   sulman_nue_mine      => params_inst%sulman_nue_mine        , &   ! Nitrogen use efficiency of mycorrhizal mining                        [-]
   sulman_vmax_ref_mine => params_inst%sulman_vmax_ref_mine   , &   ! Maximum decomposition rate at reference temp for mycorrhizal mining  [s]
   sulman_rfix          => params_inst%sulman_rfix            , &   ! N fixation rate per unit symbiotic biomass                     [gN/gC/s]
   sulman_kgrowth       => params_inst%sulman_kgrowth         , &   ! Half-saturation of intermediate C pool for symbiotic growth      [gC/m2]
   sulman_max_symb_growth => params_inst%sulman_max_symb_growth, &  ! Maximum symbiont growth rate                                   [gC/m2/s]
   symbiont_necromass     => params_inst%symbiont_necromass  , &    ! Fraction of symbiotic biomass turnover into SOM as necromass [-]
   symbiont_mr            => params_inst%symbiont_mr         , &    ! Fraction of symbiotic biomass turnover used for maintenance respiration [-]
  
   symbiont_CUE         => params_inst%symbiont_CUE           , &   ! Symbiont growth efficiency / CUE [-]
   symbiont_tau         => params_inst%symbiont_tau           , &   ! Turnover of symbionts         [s-1]
   sulman_cn_symbionts  => params_inst%sulman_cn_symbionts    , &   ! N symbiont C:N               [-]
   sulman_tau_int       => params_inst%sulman_tau_int         , &   ! Turnover time of intermediate C pool                    [s-1]
   sulman_rup_veg       => params_inst%sulman_rup_veg         , &   ! Vegetation N uptake rate from intermediate N pool       [s-1]
      
   symbiont_tau_som        => params_inst%symbiont_tau_som    , &   ! Fraction symbiont necromass into SOM pools              [-]

   froot_prof           => soilbiogeochem_state_inst%froot_prof_patch        , & ! Input: (:,:)   (1/m) profile of fine roots                                 
   crootfr              => soilstate_inst%crootfr_patch                      , & ! Input: (:,:)     (-) patch fraction of roots for carbon in each soil layer (nlevgrnd)
   frootc               => cnveg_carbonstate_inst%frootc_patch               , & ! Input:   (:) (gC/m2) fine root C
   leafn                => cnveg_nitrogenstate_inst%leafn_patch              , & ! Input:   (:) (gN/m2) leaf N                                    
   leafn_storage        => cnveg_nitrogenstate_inst%leafn_storage_patch      , & ! Input:   (:) (gN/m2) leaf N storage                            
   frootn               => cnveg_nitrogenstate_inst%frootn_patch             , & ! Input:   (:) (gN/m2) fine root N                               
   frootn_storage       => cnveg_nitrogenstate_inst%frootn_storage_patch     , & ! Input:   (:) (gN/m2) fine root N storage                       
   livecrootn           => cnveg_nitrogenstate_inst%livecrootn_patch         , & ! Input:   (:) (gN/m2) live coarse root N                        
   livecrootn_storage   => cnveg_nitrogenstate_inst%livecrootn_storage_patch , & ! Input:   (:) (gN/m2) live coarse root N storage                
   livestemn            => cnveg_nitrogenstate_inst%livestemn_patch          , & ! Input:   [real(r8)  (:)] (gN/m2) live stem N
   livestemn_storage    => cnveg_nitrogenstate_inst%livestemn_storage_patch  , & ! Input:  [real(r8) (:)     ]  (gN/m2) live stem N storage
   h2osoi_liq           => waterstatebulk_inst%h2osoi_liq_col                , & ! Output:(:,:) (kg/m2) liquid water (new)  
   t_soisno             => temperature_inst%t_soisno_col                     , & ! Input: (:,:) (Kelvin) soil temperature     
   availc               => cnveg_carbonflux_inst%availc_patch                , & ! Output:  (:) (gC/m2/s) C flux available for allocation 
   npp_growth           => cnveg_carbonflux_inst%npp_growth_patch            , & ! Output:  (:) (gC/m2/s) Total C u for growth in FUN / MIMICSplus
   c_allometry          => cnveg_state_inst%c_allometry_patch                , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)
   n_allometry          => cnveg_state_inst%n_allometry_patch                , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)
  
   sminn_to_symbiont_vr     => cnveg_nitrogenflux_inst%sminn_to_symbiont_mimicsplus_vr_patch    , & ! Output: (:,:) (gN/m2/s) Total layer soil N uptake of MIMICSplus 
   smin_no3_to_symbiont_vr  => cnveg_nitrogenflux_inst%smin_no3_to_symbiont_mimicsplus_vr_patch , & ! Output: (:,:) (gN/m2/s) Total layer soil NO3 uptake of MIMICSplus 
   smin_nh4_to_symbiont_vr  => cnveg_nitrogenflux_inst%smin_nh4_to_symbiont_mimicsplus_vr_patch , & ! Output: (:,:) (gN/m2/s) Total layer soil NH4 uptake of MIMICSplus
   
   smin_no3_to_plant_vr => soilbiogeochem_nitrogenflux_inst%smin_no3_to_plant_vr_col , & ! Input:  (:,:) (gN/m3/s) col vertically-resolved plant uptake of soil NO3 
   smin_nh4_to_plant_vr => soilbiogeochem_nitrogenflux_inst%smin_nh4_to_plant_vr_col , & ! Input:  (:,:) (gN/m3/s) col vertically-resolved plant uptake of soil NH4 
   
   decomp_cpools_vr     => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col   , &  ! Input: (:,:,:) (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) C pools
   decomp_npools_vr     => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col , &  ! Input: (:,:,:) (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   totsymbc             => soilbiogeochem_carbonstate_inst%totsymbc_col         , &  ! Total C of symbiont pools [gC/m2]
   totsymbn             => soilbiogeochem_nitrogenstate_inst%totsymbn_col       , &  ! Total N of symbiont pools [gC/m2]
 
   ! Symbiont variables  
   is_active            => symbiont_inst%is_active      , &     ! If symbiont uptake pathway is active for patch   [-]
   myc_symbiont         => pftcon%myc_symbiont                , &     ! The fraction of ECM-associated PFT               [-]
   symb_eff             => symbiont_inst%symb_eff       , &     ! Symbiont efficiency in nitrogen uptake       [gC/gN]
   symb_name            => symbiont_inst%symb_name      , &     ! 
  
   C_reservoir          => symbiont_inst%C_reservoir    , &     ! Carbon reservoir in intermediate pools       [gC/m2]
   N_reservoir          => symbiont_inst%N_reservoir    , &     ! Nitrogen reservoir in intermediate pools     [gN/m2]
   C_biomass            => symbiont_inst%C_biomass      , &     ! Carbon biomass of symbiont                   [gC/m2]
   N_biomass            => symbiont_inst%N_biomass      , &     ! Nitrogen biomass of symbiont                 [gN/m2]

   C_alloc              => symbiont_inst%C_alloc        , &     ! Carbon allocation to symbionts based on ROI [gC/m2/s]
   
   total_symbiont_turnover_C              => symbiont_inst%total_symbiont_turnover_C        , &  
   total_symbiont_turnover_N              => symbiont_inst%total_symbiont_turnover_N        , &  

   symb_growth          => symbiont_inst%symb_growth    , &                ! Symbiotic biomass growth rate                [gC/m2]
   symbiont_gr_patch    => cnveg_carbonflux_inst%symbiont_gr_patch     , & ! Total C loss of symbionts that is repired (growth)    [gC/m2/s] 
   symbiont_maint_patch => cnveg_carbonflux_inst%symbiont_maint_patch  , & ! Total C loss of symbionts repired (maintainence)      [gC/m2/s] 
   
   C_mortality          => symbiont_inst%C_mortality    , &     ! Symbiotic turnover per soil layer and column [gC/m3/s]
   N_mortality          => symbiont_inst%N_mortality    , &     ! Symbiotic turnover per soil layer and column [gN/m3/s]
   
   no3_passiv_up  => symbiont_inst%no3_passiv_up  , &   ! Passive root NO3 (nitrate) uptake     [gN/m3/s]
   nh4_passiv_up  => symbiont_inst%nh4_passiv_up  , &   ! Passive root NH4 (ammonium) uptake    [gN/m3/s]
   no3_active_up  => symbiont_inst%no3_active_up  , &   ! Active root NO3 (nitrate) uptake      [gN/m3/s]
   nh4_active_up  => symbiont_inst%nh4_active_up  , &   ! Active root NH4 (ammonium) uptake     [gN/m3/s]
   no3_scav_up    => symbiont_inst%no3_scav_up    , &   ! Scavenger NO3 (nitrate) uptake        [gN/m3/s]
   nh4_scav_up    => symbiont_inst%nh4_scav_up    , &   ! Scavenger NH4 (ammonium) uptake       [gN/m3/s]

   somc_nuptake_col     => symbiont_inst%somc_nuptake_col , &   ! Nitrogen uptake from SOMc via mining         [gN/m3/s]
   somp_nuptake_col     => symbiont_inst%somp_nuptake_col , &   ! Nitrogen uptake from SOMp via mining         [gN/m3/s]
   n_mine_somc2soma_col    => symbiont_inst%n_mine_somc2soma_col, &   ! Leftover part of co-mineralized N, not taken up by miners  [gN/m3/s]
   n_mine_somp2soma_col    => symbiont_inst%n_mine_somp2soma_col, &   ! Leftover part of co-mineralized N, not taken up by miners  [gN/m3/s]

   N_stress    => symbiont_inst%N_stress, &  
   C_allocation_to_N_acq    => symbiont_inst%C_allocation_to_N_acq, &  

   n_to_plant_mimicsplus => cnveg_nitrogenflux_inst%n_to_plant_mimicsplus_patch, & ! Output:[real(r8) (:)]  nitrogen sent to plant from symbionts (gN/m2/s)
   N_fixation            => cnveg_nitrogenflux_inst%Nfix_patch                   , & ! Output:  [real(r8) (:) ]  Symbiotic BNF (gN/m2/s)
   nfix_to_sminn_mimicsplus  => soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_mimicsplus_col   , & ! Output:  [real(r8) (:)] symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
  
   somc_cuptake_col     => symbiont_inst%somc_cuptake_col   , &   ! Nitrogen uptake from SOMc via mining       [gC/m3/s]
   somp_cuptake_col     => symbiont_inst%somp_cuptake_col   , &   ! Nitrogen uptake from SOMp via mining       [gC/m3/s]
   root_exudate_C_col   => symbiont_inst%root_exudate_C_col , &   ! Leftover C from allocation to symbionts    [gC/m3/s]
   N_symb_up            => symbiont_inst%N_symb_up          , &   ! Symbiont nitrogen uptake                   [gN/m3/s]
   N_to_plant           => symbiont_inst%N_to_plant           &
    )

   !-----------------------------------------------------------------------
   
   ! Calculationg a root profile
   ! https://escomp.github.io/ctsm-docs/versions/master/html/tech_note/Plant_Hydraulics/CLM50_Tech_Note_Plant_Hydraulics.html?highlight=root
   availc_alloc(bounds%begp:bounds%endp)                          = 0.0_r8
   C_allocation_to_N_acq(bounds%begp:bounds%endp)                 = 0.0_r8
   root_dens_frac(bounds%begp:bounds%endp, 1:nlevdecomp)          = 0.0_r8
   root_exudate_C(bounds%begp:bounds%endp)                        = 0.0_r8
   root_exudate_C_layer(bounds%begp:bounds%endp,1:nlevdecomp)     = 0.0_r8
   smin_nh4_avail(bounds%begp:bounds%endp, 1:nlevdecomp)          = 0.0_r8
   smin_no3_avail(bounds%begp:bounds%endp, 1:nlevdecomp)          = 0.0_r8
   symb_turnover_layer_C(bounds%begp:bounds%endp, 1:nlevdecomp)   = 0.0_r8 
   symb_turnover_layer_N(bounds%begp:bounds%endp, 1:nlevdecomp)   = 0.0_r8 
   N_mine_somc2soma(bounds%begp:bounds%endp, 1:nlevdecomp)        = 0.0_r8
   N_mine_somp2soma(bounds%begp:bounds%endp, 1:nlevdecomp)        = 0.0_r8
   symbiont_gr_patch(bounds%begp:bounds%endp)                     = 0.0_r8
   symbiont_maint_patch(bounds%begp:bounds%endp)                  = 0.0_r8

  
   do i = 1,n_symb
      old_C_biomass(bounds%begp:bounds%endp,i)     = C_biomass(bounds%begp:bounds%endp,i)
      old_C_reservoir(bounds%begp:bounds%endp,i)   = C_reservoir(bounds%begp:bounds%endp,i)
      old_N_biomass(bounds%begp:bounds%endp,i)     = N_biomass(bounds%begp:bounds%endp,i)
      old_N_reservoir(bounds%begp:bounds%endp,i)   = N_reservoir(bounds%begp:bounds%endp,i)
   end do

   do fp = 1,num_soilp        
      p = filter_soilp(fp)
      c = patch%column(p)
      g = patch%gridcell(p)
      
      norm_froot_prof(p,1:nlevdecomp) = 0.0_r8

        do j = 1,nlevdecomp
           norm_froot_prof(p,j) = froot_prof(p,j) * col%dz(c,j)
        end do
      
        if (sum(norm_froot_prof(p,1:nlevdecomp)) > 0._r8) then
           do j = 1,nlevdecomp
              norm_froot_prof(p,j) = norm_froot_prof(p,j) / sum(norm_froot_prof(p,1:nlevdecomp))
           end do
        else
           norm_froot_prof(p,1:nlevdecomp) = 0._r8
        endif

    ! Inorganic soil nitrogen flux avaliable for plant uptake
      do j = 1, nlevdecomp
         smin_nh4_avail(p,j) = smin_nh4_to_plant_vr(c,j)  
         smin_no3_avail(p,j) = smin_no3_to_plant_vr(c,j) 
      end do
   end do

   !-----------------------------------------------------------------------

   maint_resp                                                     = 0.0_r8
   symbiont_turnover_C_to_som(bounds%begp:bounds%endp, 1:n_symb)  = 0.0_r8
   symbiont_turnover_N_to_som(bounds%begp:bounds%endp, 1:n_symb)  = 0.0_r8
   npp_growth(bounds%begp:bounds%endp)                            = 0.0_r8
   plantCN(bounds%begp:bounds%endp)                               = 0.0_r8
   N_demand(bounds%begp:bounds%endp)                              = 0.0_r8
   total_N(bounds%begp:bounds%endp)                               = 0.0_r8
   potential_stored_N(bounds%begp:bounds%endp)                    = 0.0_r8
   npp_growth_potential(bounds%begp:bounds%endp)                  = 0.0_r8
  
   somc_nuptake(bounds%begp:bounds%endp, 1:nlevdecomp)      = 0.0_r8
   somp_nuptake(bounds%begp:bounds%endp, 1:nlevdecomp)      = 0.0_r8
   somc_cuptake(bounds%begp:bounds%endp, 1:nlevdecomp)      = 0.0_r8
   somp_cuptake(bounds%begp:bounds%endp, 1:nlevdecomp)      = 0.0_r8
   N_mine_somc2soma(bounds%begp:bounds%endp, 1:nlevdecomp)  = 0.0_r8
   N_mine_somp2soma(bounds%begp:bounds%endp, 1:nlevdecomp)  = 0.0_r8
   no3_scav_up(begp:endp,1:nlevdecomp)                      = 0.0_r8
   nh4_scav_up(begp:endp,1:nlevdecomp)                      = 0.0_r8
   no3_active_up(begp:endp,1:nlevdecomp)                    = 0.0_r8
   nh4_active_up(begp:endp,1:nlevdecomp)                    = 0.0_r8
   sum_no3_up(begp:endp,1:nlevdecomp)                       = 0.0_r8
   sum_nh4_up(begp:endp,1:nlevdecomp)                       = 0.0_r8
   N_fixation(begp:endp)                                    = 0.0_r8
   root_N_active_uptake(begp:endp)                          = 0.0_r8
   root_N_to_plant(begp:endp)                               = 0.0_r8
   n_to_plant_mimicsplus(bounds%begp:bounds%endp)           = 0.0_r8
   

   ! Scavenging (AM-style)
   call myc_scavenger_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                                 bounds, symbiont_inst, norm_froot_prof(begp:endp,1:nlevdecomp), &
                                 smin_no3_avail(begp:endp,1:nlevdecomp), smin_nh4_avail(begp:endp,1:nlevdecomp), &
                                 no3_scav_up(begp:endp,1:nlevdecomp), nh4_scav_up(begp:endp,1:nlevdecomp)) 

   ! Mycorrhizal N mining (ECM-style)
   call myc_miner_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                             bounds, symbiont_inst, temperature_inst, soilstate_inst, waterstatebulk_inst, &
                             soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst, &
                             norm_froot_prof(begp:endp,1:nlevdecomp), &
                             somc_nuptake(begp:endp,1:nlevdecomp), somp_nuptake(begp:endp,1:nlevdecomp), &
                             somc_cuptake(begp:endp,1:nlevdecomp), somp_cuptake(begp:endp,1:nlevdecomp), &
                             N_mine_somc2soma(begp:endp,1:nlevdecomp), N_mine_somp2soma(begp:endp,1:nlevdecomp))
                           
   ! Active root uptake 
   call active_root_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds, &
                                soilstate_inst, norm_froot_prof(begp:endp,1:nlevdecomp), frootc(begp:endp), &
                                smin_no3_avail(begp:endp,1:nlevdecomp), smin_nh4_avail(begp:endp,1:nlevdecomp), &
                                no3_active_up(begp:endp,1:nlevdecomp), nh4_active_up(begp:endp,1:nlevdecomp))

   ! Passive root uptake

   nh4_passiv_up(bounds%begp:bounds%endp, 1:nlevdecomp) = 0.0_r8
   no3_passiv_up(bounds%begp:bounds%endp, 1:nlevdecomp) = 0.0_r8

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
      
      ! Passive nitrogen uptake from the rhizosphere by roots with soil water uptake:
      ! smin_no3_vr is per column and has to be per patch
      ! if water in layer:
      do j = 1,nlevdecomp
         t_soi_degC = t_soisno(c,j) - tfrz     ! Soil temperature in degrees Celcius
         if (t_soi_degC > 0.01_r8 .and. h2osoi_liq(c,j) > 0.01_r8) then
            no3_passiv_up(p,j) = (waterfluxbulk_inst%qflx_tran_veg_patch(p) * dt) * (smin_no3_avail(p,j) / h2osoi_liq(c,j))
            nh4_passiv_up(p,j) = (waterfluxbulk_inst%qflx_tran_veg_patch(p) * dt) * (smin_nh4_avail(p,j) / h2osoi_liq(c,j))
         else
            nh4_passiv_up(p,j) = 0.0_r8
            no3_passiv_up(p,j) = 0.0_r8
         end if
      enddo 
   enddo

   !----------------------------------------------------------------------
 
   ! Fluxes: Inorganic nitrogen pool -> Symbionts 
   ! Nitrogen uptake by symbionts is limited to not deplete inorganic N pool
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
      
      root_N_active_uptake(p)     = 0.0_r8
      root_N_to_plant(p)          = 0.0_r8
      do j = 1, nlevdecomp
         sum_no3_up(p,j) = (no3_passiv_up(p,j) + no3_active_up(p,j) + no3_scav_up(p,j))
         sum_nh4_up(p,j) = (nh4_passiv_up(p,j) + nh4_active_up(p,j) + nh4_scav_up(p,j))
         
        if (smin_no3_avail(p,j) <= 0.0_r8 .and. sum_no3_up(p,j) > 0.0_r8) then
            write(iulog,*) 'Warning: NO3 uptake attempted from layer with zero availability.'
            no3_passiv_up(p,j) = 0.0_r8
            no3_active_up(p,j) = 0.0_r8
            no3_scav_up(p,j)   = 0.0_r8
        endif

        if (smin_nh4_avail(p,j) <= 0.0_r8 .and. sum_nh4_up(p,j) > 0.0_r8) then
            write(iulog,*) 'Warning: NH4 uptake attempted from layer with zero availability.'
            nh4_passiv_up(p,j) = 0.0_r8
            nh4_active_up(p,j) = 0.0_r8
            nh4_scav_up(p,j)   = 0.0_r8
        endif

        
         ! If nitrogen uptake exceeds avaliable nitrogen, scale each uptake pathway down
         ! Without multipling by 0.9, I scale to N uptake down, but still allow to take up all avaliable N from soil (maybe not so good)
         ! Therefore I multiply with 0.9 to leave 10% in soil
         if ( (sum_no3_up(p,j) > smin_no3_avail(p,j)) .and. &
              (smin_no3_avail(p,j) > 0.0_r8) ) then
           write(iulog,*)'NO3 uptake by passive / active / scavenger pathway exceeds soil N uptake and was scaled down, leaving 10% N in soil'
           no3_passiv_up(p,j)  = no3_passiv_up(p,j)  * ((smin_no3_avail(p,j) / sum_no3_up(p,j)) * 0.9_r8)
           no3_active_up(p,j)  = no3_active_up(p,j)  * ((smin_no3_avail(p,j) / sum_no3_up(p,j)) * 0.9_r8)
           no3_scav_up(p,j)    = no3_scav_up(p,j)    * ((smin_no3_avail(p,j) / sum_no3_up(p,j)) * 0.9_r8)
         endif

         if ( (sum_nh4_up(p,j) > smin_nh4_avail(p,j)) .and. &
              (smin_nh4_avail(p,j) > 0.0_r8) ) then
            write(iulog,*)'NH4 uptake by passive / active / scavenger pathway exceeds soil N uptake and was scaled down, leaving 10% N in soil'
            nh4_passiv_up(p,j)  = nh4_passiv_up(p,j)  * ((smin_nh4_avail(p,j) / sum_nh4_up(p,j)) * 0.9_r8)
            nh4_active_up(p,j)  = nh4_active_up(p,j)  * ((smin_nh4_avail(p,j) / sum_nh4_up(p,j)) * 0.9_r8)
            nh4_scav_up(p,j)    = nh4_scav_up(p,j)    * ((smin_nh4_avail(p,j) / sum_nh4_up(p,j)) * 0.9_r8)
         endif

         root_N_active_uptake(p) = root_N_active_uptake(p) + (no3_active_up(p,j) + nh4_active_up(p,j)) * col%dz(c,j)

         root_N_to_plant(p) = root_N_to_plant(p) + (no3_active_up(p,j) + nh4_active_up(p,j) +  no3_passiv_up(p,j) + nh4_passiv_up(p,j)) * col%dz(c,j)

         ! Total NO3 and NH4 soil uptake that needs to be substracted from inorganic N soil pool
         smin_no3_to_symbiont_vr(p,j) = (no3_passiv_up(p,j) + no3_active_up(p,j) + no3_scav_up(p,j))
         smin_nh4_to_symbiont_vr(p,j) = (nh4_passiv_up(p,j) + nh4_active_up(p,j) + nh4_scav_up(p,j))

         sminn_to_symbiont_vr(p,j) = smin_nh4_to_symbiont_vr(p,j) + smin_no3_to_symbiont_vr(p,j)

      end do
   end do

  
   !--------------------------
   ! UPDATEING RESERVOIRS 
 
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
      ! Amount of nitrogen fixed by fixer biomass
      N_fixation(p) = C_biomass(p,i_fixer) * sulman_rfix ! units: gN/m2/s

      do j = 1,nlevdecomp
         N_reservoir(p,i_miner) =  N_reservoir(p,i_miner) + ((somc_nuptake(p,j) * dt + somp_nuptake(p,j) * dt) * col%dz(c,j))
         C_reservoir(p,i_miner) =  C_reservoir(p,i_miner) + ((somc_cuptake(p,j) * dt + somp_cuptake(p,j) * dt) * col%dz(c,j))
         N_reservoir(p,i_scav) =  N_reservoir(p,i_scav) + ((no3_scav_up(p,j) * dt + nh4_scav_up(p,j) * dt) * col%dz(c,j))
         ! T dependence of N fixation from Houlton et al. (2008) Nature paper (normalized to peak at 1.0)
         ! a=-3.62, b=0.27, c=25.15, T effect = exp(-0.5*b*c+b*Ts*(1-0.5*Ts/c))
         ! Could be used as if statement: if(N_fix_Tdep_Houlton) 
         N_fixation(p) = N_fixation(p) * exp(-0.5*0.27*25.15 + 0.27*(t_soisno(c,j)-273.15)*(1.0-0.5*(t_soisno(c,j)-273.15)/25.15))
      enddo
      N_reservoir(p,i_fixer) = N_reservoir(p,i_fixer) + (N_fixation(p) * dt)
   enddo

   !------------------------


   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)

       ! GROWTH AND TURNOVER
      
       ! gross symbiotic growth (without applying CUE) [gC/m2/s]
       symb_growth_gross(p,i_scav) = ((sulman_max_symb_growth) * C_reservoir(p,i_scav) / (C_reservoir(p,i_scav) + sulman_kgrowth))

       ! Growth respiration: C loss during growth from C reservoir to C biomass pool due to CUE [gC/m2/s]
       growth_resp(p,i_scav) = symb_growth_gross(p,i_scav) * (1.0 - symbiont_CUE(i_scav))

       ! Net symbiotic growth [gC/m2/s]
       symb_growth(p,i_scav) =  symb_growth_gross(p,i_scav) *  symbiont_CUE(i_scav)

       ! Maintainance respiration [gC/m2]
       maint_resp = min(C_biomass(p,i_scav) * (symbiont_tau(i_scav) * dt) * symbiont_mr, symb_growth(p,i_scav) * dt)
       !Nitrogen limitation
       if ((symb_growth(p,i_scav) * dt) - maint_resp > sulman_cn_symbionts(i_scav) * N_reservoir(p,i_scav) * 0.9_r8)  then
          ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
        symb_growth(p,i_scav) = (sulman_cn_symbionts(i_scav) * N_reservoir(p,i_scav) * 0.9_r8 + maint_resp) / dt
        ! Growth respiration updated in case of N limitation [gC/m2/s]
        growth_resp(p,i_scav) = symb_growth(p,i_scav) * (1.0 - symbiont_CUE(i_scav)) / symbiont_CUE(i_scav)
       end if

       ! Fraction of N from maintainace respiration stays in N reservoir while C is respiered
       ! N_reservoir(p,i_scav) = N_reservoir(p,i_scav) + ((N_biomass(p,i_scav) * symbiont_tau(i_scav) * dt) * symbiont_mr))
       
       ! Total symbiont Turnover (including necromass and maintanance respiration) [gC/m2/s]
       total_symbiont_turnover_C(p,i_scav)  = C_biomass(p,i_scav)  * symbiont_tau(i_scav) 
       total_symbiont_turnover_N(p,i_scav)  = N_biomass(p,i_scav)  * symbiont_tau(i_scav) 

       ! Fraction of N from maintainace respiration stays in N reservoir while C is respiered
       maint_resp_N(p,i_scav) = total_symbiont_turnover_N(p,i_scav) * symbiont_mr
      
       ! C biomass plus growth flux from reservoir minus the turnover (including maint. respiration) [gC/m2]
       C_biomass(p,i_scav) = (C_biomass(p,i_scav) + symb_growth(p,i_scav) * dt) - total_symbiont_turnover_C(p,i_scav) * dt

       ! N_reservoir(p,i_scav) = N_reservoir(p,i_scav) + N_biomass(p,i_scav) !MAT

       ! N_biomass [gN/m2]
       N_biomass(p,i_scav) = (C_biomass(p,i_scav) / sulman_cn_symbionts(i_scav))
        
       ! C resevoir minus growth and growth respiration [gC/m2]
       C_reservoir(p,i_scav) = C_reservoir(p,i_scav) - ((symb_growth(p,i_scav) + growth_resp(p,i_scav)) * dt)
       !if (C_reservoir(p,i_scav) < 0._r8) then
       !   write(iulog,*) "WARNING: C_reservoir scav went negative"
       !   C_reservoir(p,i_scav) = 0._r8
       !end if

       ! N reservoir doesn't have growth respiration so only minus growth
       N_reservoir(p,i_scav) = N_reservoir(p,i_scav) - ((symb_growth(p,i_scav) * dt) / sulman_cn_symbionts(i_scav)) +  (maint_resp_N(p,i_scav) * dt)  
      ! if (N_reservoir(p,i_scav) < 0._r8) then
      !    write(iulog,*) "WARNING: N_reservoir scav went negative"
      !    N_reservoir(p,i_scav) = 0._r8
      ! end if

       ! Mycorrhizal miners
       symb_growth_gross(p,i_miner) = (sulman_max_symb_growth * C_reservoir(p,i_miner) / (C_reservoir(p,i_miner) + sulman_kgrowth))
       growth_resp(p,i_miner) = symb_growth_gross(p,i_miner) * (1.0 - symbiont_CUE(i_miner))
       symb_growth(p,i_miner) = symb_growth_gross(p,i_miner) * symbiont_CUE(i_miner)
       
       maint_resp = min(C_biomass(p,i_miner) * (symbiont_tau(i_miner) * dt) * symbiont_mr, symb_growth(p,i_miner) * dt)

       if ((symb_growth(p,i_miner) * dt) - maint_resp > sulman_cn_symbionts(i_miner) * N_reservoir(p,i_miner) * 0.9_r8) then
           symb_growth(p,i_miner) = (sulman_cn_symbionts(i_miner) * N_reservoir(p,i_miner) * 0.9_r8 + maint_resp) / dt
           growth_resp(p,i_miner) = symb_growth(p,i_miner) * (1.0 - symbiont_CUE(i_miner)) / symbiont_CUE(i_miner)
       end if

       total_symbiont_turnover_C(p,i_miner) = C_biomass(p,i_miner) * symbiont_tau(i_miner)
       total_symbiont_turnover_N(p,i_miner) = N_biomass(p,i_miner) * symbiont_tau(i_miner)

       !N_reservoir(p,i_miner) = N_reservoir(p,i_miner) + ((total_symbiont_turnover_N(p,i_miner) * dt) * symbiont_mr)
       maint_resp_N(p,i_miner) = total_symbiont_turnover_N(p,i_miner) * symbiont_mr
      
       C_biomass(p,i_miner) = (C_biomass(p,i_miner) + symb_growth(p,i_miner) * dt) - total_symbiont_turnover_C(p,i_miner) * dt
       
       N_biomass(p,i_miner) = C_biomass(p,i_miner) / sulman_cn_symbionts(i_miner)
       C_reservoir(p,i_miner) = C_reservoir(p,i_miner) - ((symb_growth(p,i_miner) + growth_resp(p,i_miner)) * dt)
       N_reservoir(p,i_miner) = N_reservoir(p,i_miner) - (symb_growth(p,i_miner) * dt / sulman_cn_symbionts(i_miner)) + (maint_resp_N(p,i_miner) * dt)
       

       ! Nitrogen Fixation
       symb_growth_gross(p,i_fixer) = (sulman_max_symb_growth * C_reservoir(p,i_fixer) / (C_reservoir(p,i_fixer) + sulman_kgrowth))
       growth_resp(p,i_fixer) = symb_growth_gross(p,i_fixer) * (1.0 - symbiont_CUE(i_fixer))
       symb_growth(p,i_fixer) = symb_growth_gross(p,i_fixer) * symbiont_CUE(i_fixer) !C30
   
       ! Fixation has to be done at the biomas update, since it is reduced by the growth

       total_symbiont_turnover_C(p,i_fixer) = C_biomass(p,i_fixer) * symbiont_tau(i_fixer)
       total_symbiont_turnover_N(p,i_fixer) = N_biomass(p,i_fixer) * symbiont_tau(i_fixer)

       maint_resp_N(p,i_fixer) = total_symbiont_turnover_N(p,i_fixer) * symbiont_mr
      
       C_biomass(p,i_fixer) = (C_biomass(p,i_fixer) + symb_growth(p,i_fixer) * dt) - total_symbiont_turnover_C(p,i_fixer) * dt

       N_biomass(p,i_fixer) = C_biomass(p,i_fixer) / sulman_cn_symbionts(i_fixer)  !
       C_reservoir(p,i_fixer) = C_reservoir(p,i_fixer) - (symb_growth_gross(p,i_fixer)) * dt


       
       ! N reservoir grows by the amount of N that was fixed and by the N that is left after maintainance respiration 
       ! symb_growth is not substracted here, because N fixers just make all the N they need for their biomass out of the air
       N_reservoir(p,i_fixer) = N_reservoir(p,i_fixer) + (maint_resp_N(p,i_fixer) * dt)

       !---Terje
       !N_reservoir(p,i_fixer) = N_reservoir(p,i_fixer) - symb_growth(p,i_fixer) / sulman_cn_symbionts(i_fixer)
       !---
       
       ! N fixers just make all the N they need for their biomass
       N_fixation(p) = N_fixation(p) + (symb_growth(p,i_fixer) / sulman_cn_symbionts(i_fixer))

       do i = 1,n_symb
          if (is_active(p,i)) then
             if (C_reservoir(p,i) .eq. 0._r8) then
                write(iulog,*) 'WARNING: C_res_'//trim(symb_name(i))//' is 0'
             else if (C_reservoir(p,i) < 0._r8) then 
                call endrun(msg = 'ERROR: C_res_'//trim(symb_name(i))//' is negative' // &
                    errMsg(sourcefile, __LINE__))
             endif
       
             if (C_biomass(p,i) .eq. 0._r8) then
                write(iulog,*) 'WARNING: C_bio_'//trim(symb_name(i))//' is 0'
             else if (C_biomass(p,i) < 0._r8) then 
                call endrun(msg = 'ERROR: C_bio_'//trim(symb_name(i))//' is negative' // &
                    errMsg(sourcefile, __LINE__))

             if (N_reservoir(p,i) .eq. 0._r8) then
                write(iulog,*) 'WARNING: N_res_'//trim(symb_name(i))//' is 0'
             else if (C_reservoir(p,i) < 0._r8) then 
                call endrun(msg = 'ERROR: N_res_'//trim(symb_name(i))//' is negative' // &
                    errMsg(sourcefile, __LINE__))
             endif
       
             if (C_biomass(p,i) .eq. 0._r8) then
                write(iulog,*) 'WARNING: N_bio_'//trim(symb_name(i))//' is 0'
             else if (C_biomass(p,i) < 0._r8) then 
                call endrun(msg = 'ERROR: N_bio_'//trim(symb_name(i))//' is negative' // &
                    errMsg(sourcefile, __LINE__))
             endif
             endif
          endif
       enddo

      
       !----------------------------------------------------------------------------------------------------------------------------

       ! Respiration during symbiont growth 
       symbiont_gr_patch(p) = growth_resp(p,i_scav) + growth_resp(p,i_miner) + growth_resp(p,i_fixer)
    
       ! Symbiotic necromass
       symbiont_turnover_C_to_som(p,i_miner) = total_symbiont_turnover_C(p,i_miner) * params_inst%symbiont_necromass
       symbiont_turnover_C_to_som(p,i_scav)  = total_symbiont_turnover_C(p,i_scav)  * params_inst%symbiont_necromass 
       symbiont_turnover_C_to_som(p,i_fixer) = total_symbiont_turnover_C(p,i_fixer) * params_inst%symbiont_necromass
 
       symbiont_turnover_N_to_som(p,i_miner) = total_symbiont_turnover_N(p,i_miner) * params_inst%symbiont_necromass
       symbiont_turnover_N_to_som(p,i_scav)  = total_symbiont_turnover_N(p,i_scav)  * params_inst%symbiont_necromass
       symbiont_turnover_N_to_som(p,i_fixer) = total_symbiont_turnover_N(p,i_fixer) * params_inst%symbiont_necromass
 
       ! Maintainance respiration as fraction of turnover [gC/m2/s]             
       symbiont_maint_patch(p) = params_inst%symbiont_mr * &
                                 (total_symbiont_turnover_C(p,i_miner) + total_symbiont_turnover_C(p,i_scav) + total_symbiont_turnover_C(p,i_fixer))
      ! Scavengers
   if (is_active(p,i_scav)) then 
      N_to_plant(p,i_scav) = N_reservoir(p,i_scav) * params_inst%sulman_rup_veg
   else
      N_to_plant(p,i_scav) = 0.0_r8
   endif
   if (is_active(p,i_miner)) then 
      N_to_plant(p,i_miner) = N_reservoir(p,i_miner) * params_inst%sulman_rup_veg
   else
      N_to_plant(p,i_miner) = 0.0_r8
   endif
      if (is_active(p,i_fixer)) then 
      N_to_plant(p,i_fixer) = N_reservoir(p,i_fixer) * params_inst%sulman_rup_veg
   else
      N_to_plant(p,i_fixer) = 0.0_r8
   endif

   ! Total nitrogen uptake by plant from intermediate symbiont pools gN/m2/s
   n_to_plant_mimicsplus(p) = N_to_plant(p,i_scav) + N_to_plant(p,i_miner) + N_to_plant(p,i_fixer) + root_N_to_plant(p)
   

   ! Calculating Plant-Microbe C-N exchange
   plantCN(p) = max(1._r8, c_allometry(p) / max(n_allometry(p), (1.e-12_r8)))
      
   npp_growth_potential(p) = n_to_plant_mimicsplus(p) * plantCN(p)
      
   ! Limit growth by available C, taking the smaller value
   npp_growth(p) = min(availc(p), npp_growth_potential(p))
      
   C_allocation_to_N_acq(p) = availc(p) - npp_growth(p)
   ! C_allocation_to_N_acq(p) = min(availc(p) *0.1, availc(p) - npp_growth(p))
   ! npp_growth(p) = availc - C_allocation_to_N_acq(p) 

   
   
   ! Limit C allocation to be only 10 % of avaliable C, the rest goes to growth
   C_allocation_to_N_acq(p) =  C_allocation_to_N_acq(p) * 0.1_r8

   npp_growth(p) = npp_growth(p) + C_allocation_to_N_acq(p) * 0.9_r8

   
   
   end do
      !----------------------------------------------------------------------------------------------------------------------------

    call roi_symbionts(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                        bounds, symbiont_inst, root_N_active_uptake, root_N_to_plant,  &
                        root_exudate_C(bounds%begp:bounds%endp))

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
   
      ! RESERVOIR UPDATES
      
      ! Updating C reservoirs with allocated plant C (calculated with ROI)
       C_reservoir(p,i_scav) = C_reservoir(p,i_scav) + (C_alloc(p,i_scav) * dt)
       C_reservoir(p,i_miner) = C_reservoir(p,i_miner) + (C_alloc(p,i_miner) * dt)
       C_reservoir(p,i_fixer) = C_reservoir(p,i_fixer) + (C_alloc(p,i_fixer) * dt)

      ! Updating N reservoirs
       N_reservoir(p,i_scav) = N_reservoir(p,i_scav) - (N_to_plant(p,i_scav) * dt)
       N_reservoir(p,i_miner) = N_reservoir(p,i_miner) - (N_to_plant(p,i_miner) * dt)
       N_reservoir(p,i_fixer) = N_reservoir(p,i_fixer) - (N_to_plant(p,i_fixer) * dt)
        
     
      ! THIS IS NOT OKE
      ! if (N_to_plant_mimicsplus(p) + root_N_to_plant(p) > C_allocation_to_N_acq(p) / plantCN(p)) then
      !    scale_N_to_plant(p) = (C_allocation_to_N_acq(p) / plantCN(p) - root_N_to_plant(p)) /  (N_to_plant_mimicsplus(p))
      ! else
      !   scale_N_to_plant(p) = 1.0_r8
      ! endif

      ! Nitrogen uptake to plant
      !N_to_plant(p,i_scav) = N_to_plant(p,i_scav) * scale_N_to_plant(p)
      !N_to_plant(p,i_miner) = N_to_plant(p,i_miner) * scale_N_to_plant(p)
      !N_to_plant(p,i_fixer)  = N_to_plant(p,i_fixer) * scale_N_to_plant(p)
        
      ! Scale uptake and return leftovers to reservoirs
      !N_reservoir(p,i_scav)  = N_reservoir(p,i_scav)  + (N_to_plant(p,i_scav) * dt) * (1.0_r8 - scale_N_to_plant(p))
      !N_reservoir(p,i_miner) = N_reservoir(p,i_miner) + (N_to_plant(p,i_miner) * dt) * (1.0_r8 - scale_N_to_plant(p))
      !N_reservoir(p,i_fixer) = N_reservoir(p,i_fixer) + (N_to_plant(p,i_fixer) * dt)  * (1.0_r8 - scale_N_to_plant(p))
       
      ! Since this variable is what plant actually gets, add root nitrogen here
      !N_to_plant_mimicsplus(p) = N_to_plant_mimicsplus(p) * scale_N_to_plant(p) + root_N_to_plant(p)

      !----------------------------------------------------------------------------------------------------------------------------
       N_symb_up(p,i_scav)  = 0._r8
       N_symb_up(p,i_miner) = 0._r8
       N_symb_up(p,i_fixer) = 0._r8

       do j = 1, nlevdecomp
            N_symb_up(p,i_scav) = N_symb_up(p,i_scav) + (no3_scav_up(p,j) + nh4_scav_up(p,j)) * col%dz(c,j)
            N_symb_up(p,i_miner) = N_symb_up(p,i_miner) + (somc_nuptake(p,j) + somp_nuptake(p,j)) * col%dz(c,j)
       end do

       N_symb_up(p,i_fixer) = N_fixation(p)
   
          
   ! LOCAL CARBON BALANCE CHECK

      do i = 1,n_symb
          errbalc =  C_reservoir(p,i) + C_biomass(p,i)                                         &   ! after updates
                  - (old_C_reservoir(p,i) + old_C_biomass(p,i))                                &   ! before updates
                  - ((C_alloc(p,i) * dt)                                                       &   ! incoming C
                  - growth_resp(p,i) * dt                                                      &   ! growth resp
                  - (total_symbiont_turnover_C(p,i) * dt * params_inst%symbiont_mr) &   ! maint. resp
                  - (total_symbiont_turnover_C(p,i) * dt * params_inst%symbiont_necromass))            ! to soil
          errbaln =  N_reservoir(p,i) + N_biomass(p,i)                                         &   ! after updates
                   - (old_N_reservoir(p,i)  + old_N_biomass(p,i))                              &   ! before updates

                   - ((N_symb_up(p,i) * dt)                                                    &   ! input
                   - (N_to_plant(p,i) * dt)                                                    &   ! outgoing N
                   - (total_symbiont_turnover_N(p,i) * dt * params_inst%symbiont_necromass))           ! to soil
          if (abs(errbalc) > 1.0e-10_r8 .or. abs(errbaln) > 1.0e-10_r8) then
          !if (abs(errbalc) > 1.0e-10_r8 .or. abs(errbaln) > 0.0_r8) then
            ! write(iulog,*), '------------------------------------------------'
            ! write(iulog,*), 'ECW: NOT BALANCED', trim(symb_name(i))
            ! write(iulog,*), 'ECW: C/N BALANCE ERRORS: ', errbalc, errbaln
            ! write(iulog,*), 'ECW:C Reservoir '         , C_reservoir(p,i)   , old_C_reservoir(p,i)
            ! write(iulog,*), 'ECW:C Biomass '           , C_biomass(p,i)     , old_C_biomass(p,i)
            ! write(iulog,*), 'ECW:C alloc'              , C_alloc(p,i)     * dt
            ! write(iulog,*), 'ECW:C Growth'             , symb_growth(p,i) * dt
            ! write(iulog,*), 'ECW:C Growth resp'        , growth_resp(p,i)
            ! write(iulog,*), 'ECW:Maint resp'           , total_symbiont_turnover_C(p,i) * params_inst%symbiont_mr
            ! write(iulog,*), 'ECW:C Turnover'           , total_symbiont_turnover_C(p,i) * params_inst%symbiont_necromass
            ! write(iulog,*), '------------------------------------------------'
            ! write(iulog,*), 'ECW:N Reservoir '                  , N_reservoir(p,i)   , old_N_reservoir(p,i)
            ! write(iulog,*), 'ECW:N Biomass '                    , N_biomass(p,i)     , old_N_biomass(p,i)
            ! write(iulog,*), 'ECW:N from symbiont uptake '       , N_symb_up(p,i)      * dt
            ! write(iulog,*), 'ECW:N send to plant '              , N_to_plant(p,i)     * dt
            ! write(iulog,*), 'ECW:N Maint resp. to reservoir '   , maint_resp_N(p,i)   * dt
            ! write(iulog,*), 'ECW:N Turnover '                   , total_symbiont_turnover_N(p,i) * params_inst%symbiont_necromass * dt
          
             if (i == i_miner .and. C_biomass(p,i) > 0.0_r8) then 
             do j = 1,nlevdecomp
               if (norm_froot_prof(p,j) > 0.0_r8) then
                !output N in somp loop and fluxes
                 !write(iulog,*) 'Layer: ', j
                 !write(iulog,*), 'norm_froot_prof'     , norm_froot_prof(p,j)
                 !write(iulog,*), 'somc_nuptake'        , somc_nuptake(p,j)
                 !write(iulog,*), 'N_mine_somc2soma'    , N_mine_somc2soma(p,j)
                 !write(iulog,*), 'N SOMc'              , decomp_npools_vr(c,j,i_chem_som)
                 !write(iulog,*), ' t_soisno'           , t_soisno(c,j)
               end if
             end do 
             end if
          end if
       end do
   end do



    

   ! Update the symb_turnover_C and _N to make them per layer with root_dens_frac
   do j = 1, nlevdecomp
      do p = bounds%begp,bounds%endp
         c = patch%column(p)
         symb_turnover_layer_C(p,j) = (symbiont_turnover_C_to_som(p,i_miner) + symbiont_turnover_C_to_som(p,i_scav)  & 
                                       + symbiont_turnover_C_to_som(p,i_fixer)) / col%dz(c,j) * norm_froot_prof (p,j)
         symb_turnover_layer_N(p,j) = (symbiont_turnover_N_to_som(p,i_miner) + symbiont_turnover_N_to_som(p,i_scav) &
                                       + symbiont_turnover_N_to_som(p,i_fixer)) / col%dz(c,j) * norm_froot_prof (p,j)
         root_exudate_C_layer(p,j) =  (root_exudate_C(p) / col%dz(c,j)) * norm_froot_prof(p,j)
      end do
   
      ! All variables going into soil have to be per column and per layer
      ! Respiration can stay per patch
      ! subroutine p2c_2d_filter performs patch to column averaging for multi level patch arrays

      ! symb_turnover_layer_C is the patch array, C_mortality is the column array 
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      symb_turnover_layer_C(bounds%begp:bounds%endp,j), &
      C_mortality(bounds%begc:bounds%endc,j))
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      symb_turnover_layer_N(bounds%begp:bounds%endp,j), &
      N_mortality(bounds%begc:bounds%endc,j))

      ! Leftover part of co-mineralized N, into SOMa
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      N_mine_somc2soma(bounds%begp:bounds%endp,j), &
      N_mine_somc2soma_col(bounds%begc:bounds%endc,j))
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      N_mine_somp2soma(bounds%begp:bounds%endp,j), &
      N_mine_somp2soma_col(bounds%begc:bounds%endc,j))

      ! Transforming nuptake from patch,layer to column,layer
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      somc_nuptake(bounds%begp:bounds%endp,j), &
      somc_nuptake_col(bounds%begc:bounds%endc,j))
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      somp_nuptake(bounds%begp:bounds%endp,j), &
      somp_nuptake_col(bounds%begc:bounds%endc,j))

      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      somc_cuptake(bounds%begp:bounds%endp,j), &
      somc_cuptake_col(bounds%begc:bounds%endc,j))
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      somp_cuptake(bounds%begp:bounds%endp,j), &
      somp_cuptake_col(bounds%begc:bounds%endc,j))
      
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
      root_exudate_C_layer(bounds%begp:bounds%endp,j), &
      root_exudate_C_col(bounds%begc:bounds%endc,j))  

   end do

    ! This needs to be here, bc N_fixation gets a different value later
    call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
    N_fixation(bounds%begp:bounds%endp), &
    nfix_to_sminn_mimicsplus(bounds%begc:bounds%endc))
   
   end associate
   end subroutine CN_soil_veg_exchange

   !----------------------------------------

   ! FUNCTIONS
   ! Calculates the maximum enzymatic activity rate (Vmax) of mycorrhizal organisms as a function of soil temperature
   ! Determines how temperature regulates enzyme activity, which is a key control on how fast mycorrhizae can process organic matter
   function Vmax_myc(soil_T)
      real(r8), intent(in)   :: soil_T                 ! Soil temperature in Kelvin
      real(r8), parameter    :: Tref=293.15            ! Reference Temperature in Kelvin
      !real(r8), parameter    :: Ea=37000_r8           ! Activation energy (kJ/mol) Sulman et al. (2019)
      real(r8), parameter    :: Ea=54000_r8            ! Activation energy (J/mol) ELIN
      real(r8), parameter    :: R_gas = 8.314472       ! Universal gas constant, J/mol*K
      real(r8)               :: alpha                  ! Scaling factor that normalizes the exponential temperature response to match a specified reference value. [s]
      real(r8)               :: Vmax_myc               ! [s-1]

      ! exp(-Ea / (R_gas * Tref)) is the Arrhenius term evaluated at the reference temperature
      alpha = params_inst%sulman_vmax_ref_mine / exp(-Ea /(R_gas*Tref))
      Vmax_myc = alpha * exp(-Ea / (R_gas * soil_T))
   end function Vmax_myc

   
   function resp_myc(soil_carbon, myc_biomass_layer, soil_T, wliq, wair)
      ! This is the rate of C removed from soil pool as respiration(not actual respiration, rename)
      ! Respiration is driven by mycorrhizae, and depends on how much mycorrhizal biomass there is and how many enzymes they produce 
      ! also limited by environmental conditions (soil moisture & temperature)

      ! USES
      use decompMod         , only : bounds_type
      !
      real(r8), intent(in) :: soil_carbon                      ! Soil carbon stocks, vertically resolved            [gC/m3]
      
      real(r8), intent(in) :: myc_biomass_layer                ! Mycorrhyzal biomass per soillayer                  [gC/m3]
      real(r8), intent(in) :: soil_T                           ! Soil temperature                                       [K]
      real(r8), intent(in) :: wliq                             ! Fraction of liquid water-filled pore space (0.0 - 1.0) [-]  
      real(r8), intent(in) :: wair                             ! Fraction of air-filled pore space (0.0 - 1.0)          [-]
      
      real(r8), parameter  :: enzyme_frac=0.1_r8                  ! Relative amount of enzymes produced by microbes     [-]
      real(r8), parameter  :: substrate_diffusion_exp = 3.0_r8    ! Exponent for theta dependence at low theta. See Davison et al DAMM model paper
      real(r8), parameter  :: gas_diffusion_exp = 2.5_r8          ! Exponent for gas diffusion power law dependence on theta See Meslin et al 2010, SSAJ
      real(r8), parameter  :: min_anaerobic_resp_factor = 0.0_r8  ! Minimum for high soil moisture Resp limitation CHECK [-]
      real(r8), parameter  :: min_dry_resp_factor       = 0.0_r8  ! Minimum for low soil moisture Resp limitation CHECK  [-]

      real(r8) :: enzymes                                      ! Enzymes released by mycorrhiza based on biomass     [gC/m3]

      ! LOCAL VARIABLES:
      real(r8) :: theta_resp_max                               ! normalization factors for soil moisture aerobic respiration depencence
      real(r8) :: aerobic_max                                  ! normalization factors for soil moisture aerobic respiration depencence
      real(r8) :: theta_func                                   !
      real(r8) :: resp_myc                                     ! decomposed carbon during mining process [gC/m3/s]
     
      ! initialize normalization factor for aerobic respiration soil moisture function
      ! From solving theta dependence for maximum:
      theta_resp_max = substrate_diffusion_exp/(gas_diffusion_exp*(1.0_r8 +substrate_diffusion_exp/gas_diffusion_exp))

      aerobic_max=theta_resp_max**substrate_diffusion_exp*(1.0_r8 - theta_resp_max)**gas_diffusion_exp !ECW

      ! Functional dependence on soil moisture, normalized so max is 1
      theta_func=(wliq**substrate_diffusion_exp)*(wair**gas_diffusion_exp)/aerobic_max
      
      ! On the wet side of the function, make sure it does not go below min_anaerobic_resp_factor
      if(wliq>theta_resp_max .and. theta_func<min_anaerobic_resp_factor) theta_func=min_anaerobic_resp_factor
      ! On the dry side of the function, make sure it does not go below min_dry_resp_factor
      if(wliq<theta_resp_max .and. theta_func<min_dry_resp_factor) theta_func=min_dry_resp_factor

      enzymes = myc_biomass_layer * enzyme_frac      
     
      ! If there is carbon avaliable, calculate mycorrhizal repiration     
      if (soil_carbon > 0.0_r8 .and. wliq > 0.0_r8) then 
         resp_myc = Vmax_myc(soil_T) * soil_carbon * enzymes / (soil_carbon * params_inst%sulman_km_mine + enzymes) * theta_func
         resp_myc = Vmax_myc(soil_T) * enzymes * (soil_carbon / (soil_carbon + params_inst%sulman_km_mine)) * theta_func
      else 
         resp_myc = 0.0_r8
      end if 

   end function resp_myc


   !----------------------------------------
  
   ! PLANT UPTAKE STRATEGIES

   ! Nitrogen uptake from the rhizosphere by roots (active transport across root-soil interface)
   ! Mineral nitrogen is taken up from the rhizosphere only

   subroutine active_root_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds, &
                                    soilstate_inst, norm_froot_prof, froot_carbon, no3_soil, nh4_soil, no3_uptake, nh4_uptake)
      ! ! USES:
      use clm_varcon        , only: rpi
      use decompMod         , only: bounds_type
      use SoilStateType     , only: soilstate_type

      ! ! LOCAL VARIABLES:
      integer :: p, fp, c, fc, j, k, l, s  ! indices'
      integer :: begp, endp, begc, endc

      ! ! ARGUMENTS
      integer, intent(in)    :: num_soilp           ! number of soil patches in filter
      integer, intent(in)    :: filter_soilp(:)     ! filter for soil patches
      integer, intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
      integer, intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
   
      type(bounds_type)      , intent(in)    :: bounds
      type(soilstate_type)   , intent(in)    :: soilstate_inst

      real(r8), intent(in)   :: froot_carbon(bounds%begp:bounds%endp)               ! fine root carbon           [gC/m2]
      real(r8), intent(in)   :: norm_froot_prof(bounds%begp:bounds%endp,1:nlevdecomp)! Fraction of root density   [-] 
      real(r8), intent(in)   :: no3_soil(bounds%begp:bounds%endp,1:nlevdecomp)      ! Avaliable soil mineral NO3 [gN/m3/s]
      real(r8), intent(in)   :: nh4_soil(bounds%begp:bounds%endp,1:nlevdecomp)      ! Avaliable soil mineral NH4 [gN/m3/s]
      real(r8), intent(inout):: no3_uptake(bounds%begp:bounds%endp,1:nlevdecomp)    ! NO3 uptake from soil       [gN/m3/s]
      real(r8), intent(inout):: nh4_uptake(bounds%begp:bounds%endp,1:nlevdecomp)    ! NH4 uptake from soil       [gN/m3/s]

      real(r8) :: root_biomass_density                                              ! Root biomass density       [g/m3]
      real(r8) :: root_cross_sec_area                                               ! Root cross sectional area  [m2]
      real(r8) :: root_length_density                                               ! Root length density        [m/m3]
      real(r8) :: rhizosphere_frac                                                  ! Fraction of rihzosphere    [-] 
                                                                                    ! sulman_r_rhiz              [m] 
      real(r8), parameter :: root_radius = 0.29e-03_r8                              ! Root radius                [m]
      real(r8), parameter :: c_to_b = 2.0_r8                                        !                            [g biomass /g C]

      associate(                                                         &
         sulman_r_rhiz          => params_inst%sulman_r_rhiz           , & ! Radius of the rhizosphere           [m]
         ivt                    => patch%itype                         , & ! Input: (:) patch vegetation type    [-]
         rootfr                 => soilstate_inst%rootfr_patch         , & ! Input: (:,:)                        [-]
         root_radius            => pftcon%root_radius                  , & ! Input: 0.00029                      [m] 
         root_density           => pftcon%root_density                   & ! Input: 0.31e06_r8  [g biomass/m3 root]
         )
      
         no3_uptake(bounds%begp:bounds%endp,1:nlevdecomp) = 0.0_r8
         nh4_uptake(bounds%begp:bounds%endp,1:nlevdecomp) = 0.0_r8
         rhizosphere_frac = 0.0_r8

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)
         do j = 1, nlevdecomp
            rhizosphere_frac = 0.0_r8
         if  (norm_froot_prof(p,j) > 0.0_r8) then

           ! Calculate Nitrogen concentration in soil layers
           ! smin_nh4_vr_col and smin_no3_vr_col should be already per layer and tell how much N is there           

           ! Root calculations (Sulman calculated root surface, nor sure how different that is)
           ! calculate conversion from conductivity to conductance
           root_biomass_density = c_to_b * froot_carbon(p) * rootfr(p,j) / col%dz(c,j)
           ! ensure minimum root biomass (using 1gC/m2)
           root_biomass_density = max(c_to_b*1._r8,root_biomass_density)
           ! Root length density: m root per m3 soil 
           root_cross_sec_area = rpi*root_radius(ivt(p))**2
           if (root_density(ivt(p)) > 0.0_r8 .and. root_cross_sec_area > 0.0_r8) then
            root_length_density = root_biomass_density / (root_density(ivt(p)) * root_cross_sec_area)
           else
              root_length_density = 0.0_r8
           endif
            rhizosphere_frac = min(rpi*((params_inst%sulman_r_rhiz+root_radius(ivt(p)))**2-root_radius(ivt(p))**2)*root_length_density,1.0_r8)
           
           
            ! Calculate Nitrogen uptake by roots
            if (no3_soil (p,j) > 0.0_r8) then
             no3_uptake(p,j) = rhizosphere_frac * ((params_inst%sulman_root_no3 * dt) * (no3_soil(p,j) / (no3_soil(p,j) + params_inst%sulman_km_no3)))
            else
             no3_uptake(p,j) = 0.0_r8
            end if 

            if (nh4_soil (p,j) > 0.0_r8) then
             nh4_uptake(p,j) = rhizosphere_frac * ((params_inst%sulman_root_nh4 * dt) * (nh4_soil(p,j) / (nh4_soil(p,j) + params_inst%sulman_km_nh4)))
            else
             nh4_uptake(p,j) = 0.0_r8
            end if 
         else 
            no3_uptake(p,j) = 0.0_r8
            nh4_uptake(p,j) = 0.0_r8
         end if 
         end do
      end do
      end associate

   end subroutine active_root_N_uptake

   !----------------------------------------

   subroutine myc_scavenger_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                                       bounds, symbiont_inst,  & 
                                       norm_froot_prof, no3_soil, nh4_soil, no3_uptake, nh4_uptake)
      !
      ! !USES:
      !
      ! ! ARGUMENTS
      type(bounds_type)           , intent(in)     :: bounds
      type(symbiont_type)         , intent(inout)  :: symbiont_inst 

      integer                     , intent(in)     :: num_soilp           ! number of soil patches in filter
      integer                     , intent(in)     :: filter_soilp(:)     ! filter for soil patches
      integer                     , intent(in)     :: num_bgc_soilc       ! number of soil columns in filter
      integer                     , intent(in)     :: filter_bgc_soilc(:) ! filter for soil columns
      !
      ! ! LOCAL VARIABLES:
      integer :: p, fp, c, fc, j, k, l, s  ! indices
      integer :: begp, endp, begc, endc

      real(r8),intent(in)     :: norm_froot_prof(bounds%begp:bounds%endp,1:nlevdecomp)    ! Fraction of root desity     [-]
      real(r8),intent(inout)  :: no3_soil(bounds%begp:bounds%endp,1:nlevdecomp)           ! avaliable soil mineral NO3  [gN/m3/s]
      real(r8),intent(inout)  :: nh4_soil(bounds%begp:bounds%endp,1:nlevdecomp)           ! avaliable soil mineral NH4  [gN/m3/s]
      real(r8),intent(inout)  :: no3_uptake(bounds%begp:bounds%endp,1:nlevdecomp)         ! NO3 uptake from soil        [gN/m3/s]
      real(r8),intent(inout)  :: nh4_uptake(bounds%begp:bounds%endp,1:nlevdecomp)         ! NH4 uptake from soil        [gN/m3/s]

      real(r8) :: myc_biomass_layer(bounds%begp:bounds%endp, 1:nlevdecomp)                ! Mycorrhyzal biomass per soil layer [gC/m3]
      real(r8) :: total_scav_nuptake(bounds%begp:bounds%endp)                             ! Total n uptake by scavengers [gN/m3/s]

      associate(                                                                  &
      sulman_k_scav_Ninorg => params_inst%sulman_k_scav_Ninorg                  , & ! Half-saturation inorganic N concentration for mycorrhizal uptake [gN/m3]
      sulman_k_scav        => params_inst%sulman_k_scav                         , & ! Half-saturation mycorrhizal biomass concentration for scavenging [gC/m3]
      sulman_v_scav        => params_inst%sulman_v_scav                         , & ! Maximum N uptake rate by scavenging mycorrhizae                [gN/m3/s]
      C_biomass            => symbiont_inst%C_biomass                           , & ! Carbon biomass of symbiont                                       [gC/m2]
      N_biomass            => symbiont_inst%N_biomass                           , & ! Nitrogen biomass of symbiont                                     [gN/m2]
      symb_efficiency      => symbiont_inst%symb_eff                              & ! Symbiont efficiency in nitrogen uptake                           [gN/gC]
      )

      myc_biomass_layer(bounds%begp:bounds%endp, 1:nlevdecomp) = 0.0_r8

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = patch%column(p)
               
            ! Check if there is mycorrhizal biomass per patch and per layer
            ! Calculating mycorrhizal biomass per soil layer
            ! Used for inorganic N (scavangers and root uptake)
            do j = 1, nlevdecomp
               myc_biomass_layer(p,j) = C_biomass(p,i_scav) * norm_froot_prof(p,j)  / col%dz(c,j)
            end do

            ! Check if there is  mycorrhizal biomass in soil layer
            do j = 1, nlevdecomp
               no3_uptake(p,j) = 0.0_r8
               nh4_uptake(p,j) = 0.0_r8
               if (myc_biomass_layer(p,j) > 0) then 
            
                  ! If there is mycorrhizal biomass in the soil layer, calculate N uptake
                  no3_uptake(p,j) = sulman_v_scav * ((no3_soil(p,j) * dt) / ((no3_soil(p,j) * dt) + sulman_k_scav_Ninorg)) * &
                      myc_biomass_layer(p,j) / ((myc_biomass_layer(p,j) + sulman_k_scav))

                  nh4_uptake(p,j) = sulman_v_scav * (nh4_soil(p,j) * dt) / ((nh4_soil(p,j) * dt) + sulman_k_scav_Ninorg) * &
                      myc_biomass_layer(p,j) / (myc_biomass_layer(p,j) + sulman_k_scav)
                  
                  total_scav_nuptake(p) = (no3_uptake(p,j) + nh4_uptake(p,j)) * col%dz(c,j)

               else 
                  no3_uptake(p,j) = 0.0_r8
                  nh4_uptake(p,j) = 0.0_r8
                  total_scav_nuptake(p) = 0.0_r8
               end if
            end do

            ! Compute efficiency, ensuring no division by zero
           if (C_biomass(p,i_scav) > 0.0_r8) then
              symb_efficiency(p,i_scav) = (total_scav_nuptake(p) * dt) / C_biomass(p,i_scav)
           else 
              symb_efficiency(p,i_scav) = 0.0_r8
           endif
         enddo
      end associate

   end subroutine myc_scavenger_N_uptake

   !----------------------------------------

   subroutine myc_miner_N_uptake (filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                                 bounds, symbiont_inst, temperature_inst, soilstate_inst, waterstatebulk_inst, &
                                 soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst, &
                                 norm_froot_prof, somc_nuptake,somp_nuptake, somc_cuptake, somp_cuptake, N_mine_somc2soma, N_mine_somp2soma)

      ! USES
      use TemperatureType                       , only: temperature_type
      use clm_varcon                            , only: tfrz
      use clm_varcon                            , only: denh2o, denice
      use SoilStateType                         , only: soilstate_type
      use WaterStateBulkType                    , only: waterstatebulk_type
      use SoilBiogeochemDecompCascadeConType    , only: decomp_method, mimics_decomp, mimicsplus_decomp
      use SoilBiogeochemDecompCascadeMIMICSMod  , only: decomp_rates_mimics
      use SoilBiogeochemCarbonStateType         , only: soilbiogeochem_carbonstate_type 
      use SoilBiogeochemNitrogenStateType       , only: soilbiogeochem_nitrogenstate_type 

      ! ARGUMENTS
      type(bounds_type)                         , intent(in)    :: bounds
      type(symbiont_type)                       , intent(inout) :: symbiont_inst 
      type(temperature_type)                    , intent(in)    :: temperature_inst
      type(soilstate_type)                      , intent(in)    :: soilstate_inst
      type(waterstatebulk_type)                 , intent(in)    :: waterstatebulk_inst
      type(soilbiogeochem_carbonstate_type)     , intent(in)    :: soilbiogeochem_carbonstate_inst
      type(soilbiogeochem_nitrogenstate_type)   , intent(in)    :: soilbiogeochem_nitrogenstate_inst

      integer                                   , intent(in)     :: num_soilp           ! number of soil patches in filter
      integer                                   , intent(in)     :: filter_soilp(:)     ! filter for soil patches
      integer                                   , intent(in)     :: num_bgc_soilc       ! number of soil columns in filter
      integer                                   , intent(in)     :: filter_bgc_soilc(:) ! filter for soil columns

      ! LOCAL VARIABLES:
      integer :: p, fp, c, fc, j, k, l, s ! indices
      integer :: begp, endp, begc, endc

      real(r8), intent(in)    :: norm_froot_prof(bounds%begp:bounds%endp, 1:nlevdecomp) ! Fraction of root desity [-]
      real(r8), intent(inout) :: somc_nuptake(bounds%begp:bounds%endp,1:nlevdecomp)     ! Nitrogen removed from SOMc pool into int miner pool [gN/m3/s]
      real(r8), intent(inout) :: somp_nuptake(bounds%begp:bounds%endp,1:nlevdecomp)     ! Nitrogen removed from SOMp pool into int miner pool [gN/m3/s]
      real(r8), intent(inout) :: somc_cuptake(bounds%begp:bounds%endp,1:nlevdecomp)     ! Carbon removed from SOMc pool into int miner pool   [gC/m3/s]
      real(r8), intent(inout) :: somp_cuptake(bounds%begp:bounds%endp,1:nlevdecomp)     ! Carbon removed from SOMp pool into int miner pool   [gC/m3/s]
      real(r8), intent(inout) :: N_mine_somc2soma(bounds%begp:bounds%endp,1:nlevdecomp) ! Leftover N due to NUE after mining from SOMc to SOMa [gN/m3/s]
      real(r8), intent(inout) :: N_mine_somp2soma(bounds%begp:bounds%endp,1:nlevdecomp) ! Leftover N due to NUE after mining from SOMp to SOMa [gN/m3/s]

      real(r8) :: wliq                                                     ! Fraction of liquid water-filled pore space (0.0 - 1.0)
      real(r8) :: wice                                                     ! Fraction of frozen water-filled pore space (0.0 - 1.0)
      real(r8) :: wair                                                     ! Fraction of air-filled pore space (0.0 - 1.0)
      real(r8) :: myc_biomass_layer(bounds%begp:bounds%endp, 1:nlevdecomp) ! Mycorrhyzal biomass in soil                   [gC/m3]
      real(r8) :: total_org_nuptake(bounds%begp:bounds%endp)               ! Total N uptake from all soil layers           [gN/m3/s]
      
      associate(                                                     &
         t_soisno          => temperature_inst%t_soisno_col        , &     ! Input:  [real(r8) (:,:)] soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         watsat            => soilstate_inst%watsat_col            , &     ! Input:  [real(r8) (:,:)]  volumetric soil water at saturation (porosity)  
         h2osoi_liq        => waterstatebulk_inst%h2osoi_liq_col   , &     ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)             
         h2osoi_ice        => waterstatebulk_inst%h2osoi_ice_col   , &     ! Output: [real(r8) (:,:) ] ice lens (kg/m2) (new)
         decomp_cpools_vr  => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col   , &  ! Input: [real(r8) (:,:,:) ] (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) C pools
         decomp_npools_vr  => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col , &  ! Input: [real(r8) (:,:,:) ] (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
         C_reservoir       => symbiont_inst%C_reservoir            , &     ! Carbon reservoir in intermediate pools   [gC/m2]
         N_reservoir       => symbiont_inst%N_reservoir            , &     ! Nitrogen reservoir in intermediate pools [gN/m2]
         C_biomass         => symbiont_inst%C_biomass              , &     ! Carbon biomass of symbiont               [gC/m2]
         N_biomass         => symbiont_inst%N_biomass              , &     ! Nitrogen biomass of symbiont             [gN/m2]
         symb_efficiency   => symbiont_inst%symb_eff                 &     ! Symbiont efficiency N uptake per unit of mycorrhizal C biomass  [gN/gC]
         )
    
      begp = bounds%begp; endp= bounds%endp

      myc_biomass_layer(bounds%begp:bounds%endp, 1:nlevdecomp) = 0.0_r8

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = patch%column(p)
            
            ! Calculating water, ice and air content in soil (equivalent to air_filled porosity, theta, theta sat in Sulman)
      
            total_org_nuptake(p) = 0.0_r8
            do j = 1, nlevdecomp
               ! Miners per layer
               myc_biomass_layer(p,j) = C_biomass(p,i_miner) * norm_froot_prof(p,j) / col%dz(c,j)

               ! this is necessary for miners when they somehow all die.
               !if (C_biomass(p,i_miner) <= 0.0_r8) then ! ==
               !   myc_biomass_layer(p,j) = 0.01_r8 * norm_froot_prof(p,j) / col%dz(c,j)
               !endif 

               if (myc_biomass_layer(p,j) > 0.0_r8) then               

                  wliq = h2osoi_liq(c,j) / col%dz(c,j) * denh2o
                  wice = h2osoi_ice(c,j) / col%dz(c,j) * denice
                  wliq = min(1.0_r8, wliq/watsat(c,j))            ! fraction of liquid water-filled pore space (0.0 - 1.0)
                  wice = min(1.0_r8, wice/watsat(c,j))            ! fraction of frozen water-filled pore space (0.0 - 1.0)
                  wair = max(0.0_r8, 1.0_r8 - wliq- wice)         ! fraction of air-filled pore space (0.0 - 1.0)
            
                  ! N uptake by mineres from SOMc & SOMp
                  somc_nuptake(p,j) = miner_nuptake(decomp_cpools_vr(c,j,i_chem_som), decomp_npools_vr(c,j,i_chem_som), &
                                       myc_biomass_layer(p,j), t_soisno(c,j), wliq, wair)
                  
                  somp_nuptake(p,j) = miner_nuptake(decomp_cpools_vr(c,j,i_phys_som), decomp_npools_vr(c,j,i_phys_som), &
                                        myc_biomass_layer(p,j), t_soisno(c,j), wliq, wair)
                 
                  total_org_nuptake(p) = total_org_nuptake(p) + (somc_nuptake(p,j) + somp_nuptake(p,j)) * col%dz(c,j)
                  
                  ! Co-mineralized carbon during mining, is send to SOMa pool (could also be respiered)
                  somc_cuptake(p,j) = resp_myc(decomp_cpools_vr(c,j,i_chem_som), myc_biomass_layer(p,j), t_soisno(c,j), wliq, wair)
                  
                  somp_cuptake(p,j) = resp_myc(decomp_cpools_vr(c,j,i_phys_som), myc_biomass_layer(p,j), t_soisno(c,j), wliq, wair)
                  
                  
                  ! Leftover N in soil after applying  NUE added to SOMa pool
                  N_mine_somc2soma(p,j) = leftover_n_mining(decomp_cpools_vr(c,j,i_chem_som), decomp_npools_vr(c,j,i_chem_som), &   
                                          myc_biomass_layer(p,j), t_soisno(c,j), wliq, wair)
   
                  N_mine_somp2soma(p,j) = leftover_n_mining(decomp_cpools_vr(c,j,i_phys_som), decomp_npools_vr(c,j,i_phys_som), &
                                           myc_biomass_layer(p,j), t_soisno(c,j), wliq, wair)
                                          
               else 
                  somc_nuptake(p,j)    = 0.0_r8
                  somp_nuptake(p,j)    = 0.0_r8
                  total_org_nuptake(p) = 0.0_r8
                  somc_cuptake(p,j)    = 0.0_r8
                  somp_cuptake(p,j)    = 0.0_r8
               end if
            end do
      
            if (C_biomass(p,i_miner) > 0._r8) then
               symb_efficiency(p,i_miner) = (total_org_nuptake(p) * dt) / C_biomass(p,i_miner)
            else
               symb_efficiency(p,i_miner) = 0._r8
            end if
            
         end do 
      end associate

   end subroutine myc_miner_N_uptake

   !------------------------------------------------------------------------------------------------
   ! MINER FUNCTIONS

   function potential_mined_n(soil_carbon, soil_nitrogen, myc_biomass_layer, soil_T, soil_water, soil_air)
     ! DESCRIPTION
     ! This helper returns the co-mineralized N before NUE is applied

     ! ARGUMENTS
      real(r8), intent(in) :: soil_carbon        ! Soil carbon stocks, vertically resolved   [gC/m3]
      real(r8), intent(in) :: soil_nitrogen      ! Soil nitrogen stocks, vertically resolved [gN/m3]
      real(r8), intent(in) :: myc_biomass_layer  ! Mycorrhyzal biomass                       [gC/m3]
      real(r8), intent(in) :: soil_T             ! Soil temperature                              [K]
      real(r8), intent(in) :: soil_water         ! Fraction of liquid water-filled pore space    [-]
      real(r8), intent(in) :: soil_air           ! Fraction of air-filled pore space             [-]
      
      ! LOCAL VARIABLES:
      real(r8) :: potential_mined_n              ! potential mined N [gN/m3/s]
      real(r8) :: potential_tempResp             ! decomposed carbon during mining process [gC/m3/s]
    
      potential_mined_n = 0.0_r8

      ! Call resp_myc(...) to get the soil carbon respiration by mycorrhiza
      potential_tempResp=resp_myc(soil_carbon, myc_biomass_layer, soil_T, soil_water, soil_air)

      ! Don't exceed avaliable C
      if(dt*potential_tempResp > soil_carbon) then
         potential_tempResp = soil_carbon / dt
      end if

      if(soil_carbon > 0) then
         ! Use C:N ratio to estimate how much nitrogen is co-mineralized.
         potential_mined_n = potential_tempResp * (soil_nitrogen / soil_carbon)
      else 
         potential_mined_n=0.0_r8
      end if

   end function potential_mined_n
   
  
   function miner_nuptake(soil_carbon, soil_nitrogen, myc_biomass_layer, soil_T, soil_water, soil_air)
   ! This calculates the amount of nitrogen taken up from SOM pools by mycorrhizal mining.
   ! Converts carbon processing (decomposition) into nitrogen uptake, which is the key biogeochemical role of mining.
     
      ! ARGUMENTS
      real(r8), intent(in) :: soil_carbon        ! Soil carbon stocks, vertically resolved   [gC/m3]
      real(r8), intent(in) :: soil_nitrogen      ! Soil nitrogen stocks, vertically resolved [gN/m3]
      real(r8), intent(in) :: myc_biomass_layer  ! Mycorrhyzal biomass                       [gC/m3]
      real(r8), intent(in) :: soil_T             ! Soil temperature                              [K]
      real(r8), intent(in) :: soil_water         ! Fraction of liquid water-filled pore space    [-]
      real(r8), intent(in) :: soil_air           ! Fraction of air-filled pore space             [-]
      real(r8) :: miner_nuptake                  ! Mycorrhizal nitrogen uptake from soil organic matter pools via mining process [gN/m3/s]
      ! ! LOCAL VARIABLES:
      real(r8) :: pot_tempN_decomposed           ! 

      miner_nuptake    = 0.0_r8

      !Call potential_mined_n to get the potential mined N before NUE
      pot_tempN_decomposed = potential_mined_n(soil_carbon, soil_nitrogen, myc_biomass_layer, soil_T, soil_water, soil_air)

     ! Apply nitrogen use efficiency (sulman_nue_mine) to scale the actual uptake.
      miner_nuptake  = pot_tempN_decomposed*params_inst%sulman_nue_mine

   end function miner_nuptake


   function leftover_n_mining(soil_carbon, soil_nitrogen, myc_biomass_layer, soil_T, soil_water, soil_air)

       ! ARGUMENTS
      real(r8), intent(in) :: soil_carbon        ! Soil carbon stocks, vertically resolved   [gC/m3]
      real(r8), intent(in) :: soil_nitrogen      ! Soil nitrogen stocks, vertically resolved [gN/m3]
      real(r8), intent(in) :: myc_biomass_layer  ! Mycorrhyzal biomass                       [gC/m3]
      real(r8), intent(in) :: soil_T             ! Soil temperature                              [K]
      real(r8), intent(in) :: soil_water         ! Fraction of liquid water-filled pore space    [-]
      real(r8), intent(in) :: soil_air           ! Fraction of air-filled pore space             [-]
      real(r8) :: leftover_n_mining              ! Mycorrhizal nitrogen uptake from soil organic matter pools via mining process [gN/m3/s]

      ! ! LOCAL VARIABLES:
      real(r8) :: n_som_to_miner                  ! Mycorrhizal nitrogen uptake from soil organic matter pools via mining process [gN/m3/s]
      real(r8) :: pot_tempN_decomposed            ! 
   
      leftover_n_mining = 0.0_r8

       !Call potential_mined_n to get the potential mined N before NUE
      pot_tempN_decomposed = potential_mined_n(soil_carbon, soil_nitrogen, myc_biomass_layer, soil_T, soil_water, soil_air)

     ! Apply nitrogen use efficiency (sulman_nue_mine) to scale the actual uptake.
      n_som_to_miner  = pot_tempN_decomposed*params_inst%sulman_nue_mine

      ! Remaining N that is not taken up
      leftover_n_mining = pot_tempN_decomposed - n_som_to_miner

   end function leftover_n_mining
   
  !------------------------------------------------------------------------------------------------
   
  subroutine roi_symbionts(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc,  &
                              bounds, symbiont_inst, root_N_active_uptake, root_N_to_plant,  &
                              root_exudate_C)
  
   ! ARGUMENTS
   integer                        , intent(in)     :: filter_soilp(:)     ! filter for soil patches
   integer                        , intent(in)     :: filter_bgc_soilc(:) ! filter for soil columns   
   integer                        , intent(in)     :: num_soilp           ! number of soil patches in filter
   integer                        , intent(in)     :: num_bgc_soilc       ! number of soil columns in filter
   type(bounds_type)              , intent(in)     :: bounds              
   type(symbiont_type)            , intent(inout)  :: symbiont_inst       
  
   real(r8), intent(inout) :: root_exudate_C(bounds%begp:bounds%endp)        ! Leftover C from allocation to symbionts     [gC/m2/s]
   real(r8), intent(in)    :: root_N_active_uptake(bounds%begp:bounds%endp)  ! active root N uptake                        [gN/m2/s]
   real(r8), intent(in)    :: root_N_to_plant(bounds%begp:bounds%endp)       ! total (active + passive) root N uptake      [gN/m2/s]
   !
   ! ! LOCAL VARIABLES:
   integer :: p, fp, c, fc, j, k, l, s  ! indices
   integer :: begp, endp, begc, endc
   
   real(r8) :: local_active            ! Indicates which pathway is active, based on symbiont type of PFT
   
   real(r8) :: roi(bounds%begp:bounds%endp)                      ! Return of investment for all pathways [gN/gC]
   real(r8) :: mine_roi(bounds%begp:bounds%endp)                 ! Nitrogen return of carbon investment  [gN/gC]
   real(r8) :: scav_roi(bounds%begp:bounds%endp)                 ! Nitrogen return of carbon investment  [gN/gC]
   real(r8) :: fix_roi(bounds%begp:bounds%endp)                  ! Nitrogen return of carbon investment  [gN/gC]
   real(r8) :: root_roi(bounds%begp:bounds%endp)                 ! 
   real(r8) :: scav_roi_frac                                     ! Scavenger fraction of ROI [-]
   real(r8) :: mine_roi_frac                                     ! Miner fraction of ROI [-]
   real(r8) :: fix_roi_frac                                      ! Fixer fraction of ROI [-]
   real(r8) :: root_roi_frac                                     ! 
   real(r8) :: fix_alloc_accum(bounds%begp:bounds%endp)          ! Accumulated carbon allocation from plant to fixer pool     [gC/m2]
   real(r8) :: mine_alloc_accum(bounds%begp:bounds%endp)         ! Accumulated carbon allocation from plant to scavenger pool [gC/m2]
   real(r8) :: scav_alloc_accum(bounds%begp:bounds%endp)         ! Accumulated carbon allocation from plant to miner pool     [gC/m2]
  
   associate(                                                                    &
   is_active              => symbiont_inst%is_active                           , & ! Input: [logical (:,:)] if symbiont uptake pathway is active for patch
   N_reservoir            => symbiont_inst%N_reservoir                         , & ! Nitrogen reservoir in intermediate pools[gN/m2/s]
   C_biomass              => symbiont_inst%C_biomass                           , & ! Carbon biomass of symbiont [gC/m2]
   symb_eff               => symbiont_inst%symb_eff                            , & ! Symbiont efficiency N uptake per unit of mycorrhizal C biomass  [gN/gC]
   C_alloc                => symbiont_inst%C_alloc                             , & ! Carbon allocation to symbionts based on ROI [gC/m2/s]
   N_to_plant             => symbiont_inst%N_to_plant                          , &
   C_allocation_to_N_acq  => symbiont_inst%C_allocation_to_N_acq                 &
   
   )

   !--------------------------------------------------------------------------------------------------------------------------------

   ! RETURN OF INVESTMENT
  
  do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)

   ! Scavengers
   if (is_active(p,i_scav)) then 
      !N_to_plant(p,i_scav) = N_reservoir(p,i_scav) * params_inst%sulman_rup_veg
      if (C_biomass(p,i_scav) > 0.0_r8) then  ! or (C_biomass(p,i_scav) < 0.0_r8)
          scav_roi(p) = (max(0.0_r8, N_to_plant(p,i_scav)))  / (C_biomass(p,i_scav) * params_inst%symbiont_CUE(i_scav) * (params_inst%symbiont_tau(i_scav)))
      else 
      ! scav_efficiency is calculated in myc_scavenger_N_uptake under myc_efficiency
      scav_roi(p) = symb_eff(p,i_scav) / (params_inst%symbiont_CUE(i_scav) * (params_inst%symbiont_tau(i_scav) * dt))
      end if 
   else
      scav_roi = 0.0_r8
   end if 


   ! Miners

   ! N_to_plant needs to be in a new routione / in CNveg routine
   ! ROI stays in roi routien

   if (is_active(p,i_miner)) then
      !N_to_plant(p,i_miner) = N_reservoir(p,i_miner) * params_inst%sulman_rup_veg 
      if (C_biomass(p,i_miner) > 0.0_r8) then 
         mine_roi(p) = ((max(0.0_r8, N_to_plant(p,i_miner))) / (C_biomass(p,i_miner))) * params_inst%symbiont_CUE(i_miner) / (params_inst%symbiont_tau(i_miner))
      else 
         ! mine is calculated in one of the mining routines under myc_efficiency
         mine_roi(p) = symb_eff(p,i_miner) / (params_inst%symbiont_CUE(i_miner) * (params_inst%symbiont_tau(i_miner) * dt))
      end if 
   else
       mine_roi(p) = 0.0_r8
   end if 

   if (C_allocation_to_N_acq(p) > 0.0_r8) then
      root_roi(p) = max(0.0000001,(root_N_active_uptake(p))/C_allocation_to_N_acq(p))
   else
      root_roi(p) = (mine_roi(p) + scav_roi(p))*0.25_r8 ! make as above root roi should be akways super small 0.0000000000001
   endif


   ! Nitrogen Fixers
   if (is_active(p,i_fixer)) then
      !N_to_plant(p,i_fixer) = N_reservoir(p,i_fixer) * params_inst%sulman_rup_veg
      if (C_biomass(p,i_fixer) > 0.0_r8) then 
         fix_roi(p) = ((N_to_plant(p,i_fixer)) / (C_biomass(p,i_fixer))) * params_inst%symbiont_CUE(i_fixer) / (params_inst%symbiont_tau(i_fixer))
      else 
         fix_roi(p) =  params_inst%sulman_rfix / params_inst%symbiont_CUE(i_fixer) * params_inst%symbiont_tau(i_fixer)
      end if 
   else
     fix_roi = 0.0_r8
   end if 

   
   ! Calculate relative fractions
   if (scav_roi(p) + mine_roi(p) + fix_roi(p) + root_roi(p) > 0.0_r8) then
      roi(p) = scav_roi(p) + mine_roi(p) + fix_roi(p) + root_roi(p)

      scav_roi_frac = scav_roi(p) / roi(p)
      mine_roi_frac = mine_roi(p) / roi(p)
      fix_roi_frac = fix_roi(p) / roi(p)
      root_roi_frac = root_roi(p) / roi(p)
   else !should add up to 1
      scav_roi_frac = 0.4_r8 * 0.7_r8
      mine_roi_frac = 0.3_r8 * 0.7_r8
      fix_roi_frac = 0.3_r8 * 0.7_r8
      root_roi_frac = 0.3_r8
   endif 

   ! Calculate carbon allocation to pathways (without smoothing filters)
    C_alloc(p,i_fixer) = C_allocation_to_N_acq(p) * fix_roi_frac
    C_alloc(p,i_miner) = C_allocation_to_N_acq(p) * mine_roi_frac
    C_alloc(p,i_scav)  = C_allocation_to_N_acq(p) * scav_roi_frac

    ! Carbon that wasn't spend on scav, miner or fixer (including root)
    root_exudate_C(p) = C_allocation_to_N_acq(p) - C_alloc(p,i_scav) - C_alloc(p,i_miner) - C_alloc(p,i_fixer)

   end do
   end associate
   end subroutine roi_symbionts
end module CNSoilVegMIMICSplus
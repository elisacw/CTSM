module CNSoilVegMIMICSplus

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module connects the Soil decomposition module MIMICS+ (Aas et al. 2023) with the vegetation through the 
  ! symbiosis between mycorrhizal fungi and plants.
  ! Coupling follows Sulman et al. (2019)
  
  ! !USES:
  use shr_kind_mod                        , only : r8 => shr_kind_r8
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
  use SoilBiogeochemDecompCascadeMIMICSMod, only : params_inst
  use WaterStateType                      , only : waterstate_type
  use SoilStateType                       , only : soilstate_type
  use WaterStateBulkType                  , only : waterstatebulk_type
  use TemperatureType                     , only : temperature_type
  use pftconMod                           , only : pftcon, noveg
  use abortutils                          , only : endrun
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:

  ! !PRIVATE MEMBER FUNCTIONS:
  public  :: CN_soil_veg_exchange         !
  private :: active_root_N_uptake         !
  private :: myc_scavenger_N_uptake       !
  private :: myc_miner_N_uptake           !
  private :: miner_decomposition          !

  real(r8) :: dt              ! decomp timestep (seconds)

  ! !FUNCTIONS:
  private :: resp_myc                     ! Respiration of mycorrhiza
  private :: Vmax_myc                     ! Michaelis Menten Kinetics for mycorrhiza


  ! ! PUBLIC DATA
  integer, public :: i_fixer = 1 
  integer, public :: i_scav  = 2
  integer, public :: i_miner = 3
  integer, public :: n_symb  = 3

  type, public :: symbiont_type

  real(r8), pointer           :: C_biomass             (:,:) ! Carbon biomass
  real(r8), pointer           :: N_biomass             (:,:) ! Nitrogen biomass
  real(r8), pointer           :: C_reservoir           (:,:) ! Carbon intermediate pool biomass
  real(r8), pointer           :: N_reservoir           (:,:) ! Nitrogen intermediate pool biomass
  logical, pointer            :: is_active             (:,:) ! if the symbiont type is active
  character(len=10), pointer  :: symb_name             (:)   ! symbiont name
  character(len=10), pointer  :: symb_hist_name        (:)   ! symbiont name on history tapes
   
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

    class(symbiont_type)          :: this
    type(bounds_type), intent(in) :: bounds

   ! !LOCAL variables
    integer :: iveg
    integer :: p,i                        ! indices

    logical(r8) :: local_active(bounds%begp:bounds%endp, 1:n_symb)
   
    local_active(:,:) = .false. 
   
    
    
    this%symb_name(i_fixer) = 'fixer'
    this%symb_name(i_scav)  = 'scavanger'
    this%symb_name(i_miner) = 'miner'
    this%symb_hist_name(i_fixer) = 'FIX'
    this%symb_hist_name(i_scav)  = 'SCAV'
    this%symb_hist_name(i_miner) = 'MINE'
    
    do p = bounds%begp,bounds%endp
      iveg = patch%itype(p)
      if (iveg /= noveg) then 
      ! check for fixers 
      ! use fun fraction as place holder
        if (pftcon%FUN_fracfixers(iveg) > 0.0_r8) then
          local_active(p,i_fixer) = .true.
        endif
        if (pftcon%perecm(iveg) == 1.0_r8) then
         local_active(p,i_miner) = .true.
         local_active(p,i_scav) = .false.
        elseif (pftcon%perecm(iveg) == 0.0_r8) then
         local_active(p,i_miner) = .false.
         local_active(p,i_scav) = .true.
        else
         local_active(p,i_miner) = .true.
         local_active(p,i_scav) = .true.
        endif
      endif
   enddo 
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
   integer :: begp, endp
  
   begp = bounds%begp; endp= bounds%endp

    allocate(this%C_biomass(begp:endp,1:n_symb)) ; this%C_biomass(begp:endp,1:n_symb) = 0.0_r8
    allocate(this%N_biomass(begp:endp,1:n_symb)) ; this%N_biomass(begp:endp,1:n_symb) = 0.0_r8
    allocate(this%C_reservoir(begp:endp,1:n_symb)) ; this%C_reservoir(begp:endp,1:n_symb) = 0.0_r8
    allocate(this%N_reservoir(begp:endp,1:n_symb)) ; this%N_reservoir(begp:endp,1:n_symb) = 0.0_r8

   end subroutine InitAllocate

   !-----------------------------------------------------------------------

   subroutine InitHistory(this, bounds)
      !
      ! !USES:
      use histFileMod      , only: hist_addfld1d, hist_addfld2d
      use clm_varcon       , only : spval
      !
      ! !ARGUMENTS:
      class(symbiont_type)          :: this
      type(bounds_type), intent(in) :: bounds
      !
      ! !LOCAL VARIABLES:
      integer :: begp, endp
      integer :: i
      real(r8), pointer :: data1dptr (:)
     
      begp = bounds%begp; endp= bounds%endp

      do i = 1,n_symb 
         this%C_biomass(begp:endp,i) = spval
         data1dptr => this%C_biomass(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_C', units='g/m2', &
         avgflag='A', long_name=('Carbon pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='active')
   
         this%C_reservoir(begp:endp,i) = spval
         data1dptr => this%C_reservoir(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_C_inter', units='g/m2', &
         avgflag='A', long_name=('Carbon intermediate pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%N_biomass(begp:endp,i) = spval
         data1dptr => this%N_biomass(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_N', units='g/m2', &
         avgflag='A', long_name=('Nitrogen pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='active')
   
         this%N_reservoir(begp:endp,i) = spval
         data1dptr => this%N_reservoir(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_N_inter', units='g/m2', &
         avgflag='A', long_name=('Nitrogen intermediate pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')
      end do
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
      integer :: p,i                        ! indices
      integer :: begp, endp

      ! Initial carbon stock in symbiont pools for coldstart (from parameter file ) are separated into:
      ! 90 % = initial C, 10% = initial intermediate C
      do i = 1,n_symb
         do p = bounds%begp,bounds%endp
            if (this%is_active(p,i)) then
               this%C_biomass(p,i) = params_inst%sulman_initial_C_stocks(i) * 0.9_r8
               this%C_reservoir(p,i) = params_inst%sulman_initial_C_stocks(i) * 0.1_r8
               if (params_inst%sulman_cn_symbionts(i) == 0.0) then
                  this%N_biomass(p,i) = 0.0_r8
                  this%N_reservoir(p,i) = 0.0_r8
               else
                  this%N_biomass(p,i) = this%C_biomass(p,i) / params_inst%sulman_cn_symbionts(i)
                  this%N_reservoir(p,i) = this%C_reservoir(p,i) / params_inst%sulman_cn_symbionts(i)
               endif
            else
               this%C_biomass(p,i) = 0.0_r8
               this%C_reservoir(p,i) = 0.0_r8
               this%N_biomass(p,i) = 0.0_r8
               this%N_reservoir(p,i) = 0.0_r8
            endif
         enddo
      enddo
    end subroutine InitCold

   !------------------------------------------------------------------------

   subroutine Restart(this, bounds, ncid, flag)
     !
     ! !USES:
      use restUtilMod
      use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
     
     
          !
     ! !ARGUMENTS:
     class(symbiont_type)             :: this
     type(bounds_type), intent(in)    :: bounds
     type(file_desc_t), intent(inout) :: ncid   ! netcdf id
     character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
     !
     ! !LOCAL VARIABLES:
     integer                  :: i            ! indices
     logical                  :: readvar      ! determine if variable is on initial file
     character(len=128)       :: varname   ! temporary
     real(r8), pointer :: data1dptr (:)
     
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
      enddo
   end subroutine Restart

   !-----------------------------------------------------------------------

   subroutine CN_soil_veg_exchange (filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds, symbiont_inst, &
      cnveg_nitrogenstate_inst, waterstatebulk_inst, temperature_inst, cnveg_carbonflux_inst, &
      soilbiogeochem_nitrogenstate_inst, leaf_prof_patch, froot_prof_patch, croot_prof_patch, soilbiogeochem_nitrogenflux_inst, &
      waterflux_inst, soilstate_inst, cnveg_carbonstate_inst)

   ! !DESCRIPTION:

   ! !USES:
      use clm_varcon          , only: pct_to_frac
      use clm_time_manager    , only: get_step_size_real
      use WaterFluxType       , only: waterflux_type
      use SoilStateType       , only: soilstate_type
      use TemperatureType     , only: temperature_type
      use CNVegCarbonStateType, only: cnveg_carbonstate_type
      
   !
   ! !ARGUMENTS:
   type(symbiont_type)                     , intent(inout) :: symbiont_inst
   type(bounds_type)                       , intent(in)    :: bounds 
   type(cnveg_nitrogenstate_type)          , intent(in)    :: cnveg_nitrogenstate_inst
   type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
   type(waterstatebulk_type)               , intent(in)    :: waterstatebulk_inst
   type(temperature_type)                  , intent(in)    :: temperature_inst
   type(cnveg_carbonflux_type)             , intent(in)    :: cnveg_carbonflux_inst
   type(soilbiogeochem_nitrogenstate_type) , intent(in)    :: soilbiogeochem_nitrogenstate_inst
   type(waterflux_type)                    , intent(in)    :: waterflux_inst
   type(soilstate_type)                    , intent(in)    :: soilstate_inst
   type(cnveg_carbonstate_type)            , intent(in)    :: cnveg_carbonstate_inst

   
   integer                                 , intent(in)    :: num_soilp           ! number of soil patches in filter
   integer                                 , intent(in)    :: filter_soilp(:)     ! filter for soil patches
   integer                                 , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
   integer                                 , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
   
   ! !LOCAL VARIABLES
   integer :: p, fp, c, fc, j, k, l, s  ! indices
   integer :: begp, endp, begc,endc
   
   real(r8)                             , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: croot_prof_patch(bounds%begp:,1:) 
   
   real(r8)            :: t_soi_degC
   real(r8), parameter :: N_stress_max = 2.0_r8       ! Maximum N demand of plant, based on current N amount in plant []
   real(r8)    :: N_stress(bounds%begp:bounds%endp)   ! N demand of plant, based on current N amount in plant []
   real(r8)    :: C_transfer(bounds%begp:bounds%endp) ! Carbon allocated to N acquisition, higher N_stress leads to higher C_transfer
   real(r8)    :: root_dens_frac(bounds%begp:bounds%endp,1:nlevdecomp)
   real(r8)    :: root_dens_sum
   
   real(r8)    :: scav_N_reservoir
   real(r8)    :: mine_N_reservoir
   real(r8)    :: fix_N_reservoir

   ! Uptake variables for pathways
   real(r8)    :: no3_passiv_up(bounds%begp:bounds%endp,1:nlevdecomp)     ! Passive root NO3 (nitrate) uptake
   real(r8)    :: nh4_passiv_up(bounds%begp:bounds%endp,1:nlevdecomp)     ! Passive root NH4 (ammonium) uptake
   real(r8)    :: no3_active_up(bounds%begp:bounds%endp,1:nlevdecomp)     ! Active root NO3 (nitrate) uptake
   real(r8)    :: nh4_active_up(bounds%begp:bounds%endp,1:nlevdecomp)     ! Active root NH4 (ammonium) uptake
   real(r8)    :: no3_scav_up(bounds%begp:bounds%endp,1:nlevdecomp)       ! Scavenger NO3 (nitrate) uptake
   real(r8)    :: nh4_scav_up(bounds%begp:bounds%endp,1:nlevdecomp)       ! Scavenger NH4 (ammonium) uptake
   real(r8)    :: n_somc_miner_up(bounds%begp:bounds%endp,1:nlevdecomp)   ! Miner nitrogen uptake from SOMc pool
   real(r8)    :: n_somp_miner_up(bounds%begp:bounds%endp,1:nlevdecomp)   ! Minier nitrogen uptake from SOMp pool

   ! N_to_plant(bounds%begp:bounds%endp,1:n_symb)
   real(r8)    :: N_to_plant_scav(bounds%begp:bounds%endp)  ! N_to_plant_scav(p) = N_to_plant(p,i_scav)
   real(r8)    :: N_to_plant_mine(bounds%begp:bounds%endp)
   real(r8)    :: N_to_plant_fix(bounds%begp:bounds%endp)

   real(r8) :: scav_C_biomass(bounds%begp:bounds%endp)          ! Scavenger C biomass
   real(r8) :: scav_N_biomass(bounds%begp:bounds%endp)          ! Scavenger N biomass
   real(r8) :: mine_C_biomass(bounds%begp:bounds%endp)          ! Miner C biomass
   real(r8) :: mine_N_biomass(bounds%begp:bounds%endp)          ! Miner N biomass
   real(r8) :: fix_C_biomass(bounds%begp:bounds%endp)           ! Fixer C biomass
   real(r8) :: fix_N_biomass(bounds%begp:bounds%endp)           ! Fixer N biomass

   dt           = get_step_size_real()

   
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
   sulman_fnalloc       => params_inst%sulman_fnalloc         , &   !Fraction of NPP allocated to N uptake per unit N stress [fraction] 
   sulman_initial_C_stocks => params_inst%sulman_initial_C_stocks , &   ! Initial carbon stocks in fixer, miner, scavenger pools (only for coldstart)
   sulman_cn_symbionts     => params_inst%sulman_cn_symbionts     , &   ! C:N ratios of fixer, miner, scavengers as array

   crootfr              => soilstate_inst%crootfr_patch                      , & ! Input:   [real(r8) (:,:)] fraction of roots for carbon in each soil layer  (nlevgrnd)
   frootc               => cnveg_carbonstate_inst%frootc_patch               , & ! Input:   [real(r8) (:)]  (gC/m2) fine root C
   leafc                => cnveg_carbonstate_inst%leafc_patch                , & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
   leafn                => cnveg_nitrogenstate_inst%leafn_patch              , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
   leafn_storage        => cnveg_nitrogenstate_inst%leafn_storage_patch      , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
   frootn               => cnveg_nitrogenstate_inst%frootn_patch             , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N                               
   frootn_storage       => cnveg_nitrogenstate_inst%frootn_storage_patch     , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
   livecrootn           => cnveg_nitrogenstate_inst%livecrootn_patch         , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
   livecrootn_storage   => cnveg_nitrogenstate_inst%livecrootn_storage_patch , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
   h2osoi_liq           => waterstatebulk_inst%h2osoi_liq_col                , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)  
   t_soisno             => temperature_inst%t_soisno_col                     , & ! Input:  [real(r8) (:,:)   ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)       
   availc               => cnveg_carbonflux_inst%availc_patch                , & ! Output: [real(r8) (:)   ]  C flux available for allocation (gC/m2/s)
   
   sminn_vr               => soilbiogeochem_nitrogenstate_inst%sminn_vr_col        , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N                
   smin_nh4_vr            => soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col     , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NH4              
   smin_no3_vr            => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col     , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral NO3              
   smin_no3_to_plant_vr   => soilbiogeochem_nitrogenflux_inst%smin_no3_to_plant_vr_col     , & ! Input:  [real(r8) (:,:) ]
   smin_nh4_to_plant_vr   => soilbiogeochem_nitrogenflux_inst%smin_nh4_to_plant_vr_col     , & ! Input:  [real(r8) (:,:) ]
   
   ! Symbiont variables  
   is_active              => symbiont_inst%is_active                                       , & ! Input: [logical (:,:)] if symbiont uptake pathway is active for patch
   perecm                 => pftcon%perecm                                                 , & ! Input:   The fraction of ECM-associated PFT 
   C_reservoir            => symbiont_inst%C_reservoir                                     , &
   C_biomass              => symbiont_inst%C_biomass                                       , & 
   N_biomass              => symbiont_inst%N_biomass                                         &
   )
   
   !-----------------------------------------------------------------------
   
   begp = bounds%begp
   endp = bounds%endp
   begc = bounds%begc
   endc = bounds%endc
   
   !scav_c_biomass(:) = 0.0_r8
   ! do p.... 
   !scav_c_biomass(p) = symbiont_inst%C_biomass(p,i_scav)
   !....
   !enddo

   scav_C_biomass(begp:endp)  = 0.0_r8
   mine_C_biomass(begp:endp)  = 0.0_r8
   fix_C_biomass(begp:endp)   = 0.0_r8
   scav_N_biomass(begp:endp)  = 0.0_r8
   mine_N_biomass(begp:endp)  = 0.0_r8
   fix_N_biomass(begp:endp)   = 0.0_r8

   do fp = 1,num_soilp        
      p = filter_soilp(fp)
      c = patch%column(p)

      scav_C_biomass(p)    = symbiont_inst%C_biomass(p,i_scav)
      mine_C_biomass(p)    = symbiont_inst%C_biomass(p,i_miner)
      fix_C_biomass(p)     = symbiont_inst%C_biomass(p,i_fixer)
      scav_N_biomass(p)    = symbiont_inst%N_biomass(p,i_scav)
      mine_N_biomass(p)    = symbiont_inst%N_biomass(p,i_miner)
      fix_N_biomass(p)     = symbiont_inst%N_biomass(p,i_fixer)
   enddo 
   
   !-----------------------------------------------------------------------
   
   ! Calculationg a root profile
   ! https://escomp.github.io/ctsm-docs/versions/master/html/tech_note/Plant_Hydraulics/CLM50_Tech_Note_Plant_Hydraulics.html?highlight=root
   root_dens_frac(begp:endp,1:nlevdecomp) = 0.0_r8
   do fp = 1,num_soilp        
      p = filter_soilp(fp)
      c = patch%column(p)

      root_dens_sum = 0.0_r8

      do j = 1, nlevdecomp
         root_dens_sum = root_dens_sum + crootfr(p,j) * frootc(p,j)
      end do
      do j = 1, nlevdecomp
         root_dens_frac(p,j) = (crootfr(p,j) * frootc(p,j)) / root_dens_sum
      end do
   end do

   !-----------------------------------------------------------------------

   ! Calculation of N_stress: N stress indicates the nitrogen content in a plant and is used to calculate how much carbon 
   ! the plant allocates belowground = the higher N stress, the higher is the carbon allocation belowground 
   do fp = 1,num_soilp        
      p = filter_soilp(fp)
      c = patch%column(p)

      N_stress(p) = (N_stress_max*(leafn + frootn + livecrootn) - (leafn_storage + frootn_storage + livecrootn_storage)) / &
                    (leafn + frootn + livecrootn)
      
      ! Calculate the amount of C transferent to N aquisation
      ! Dynamic allocation of a fraction of NPP to root exudation
      C_transfer(p) = max(availc(p),0) * sulman_fnalloc * N_stress

      ! Plant N uptake per symbiont pathway
      N_to_plant_scav = scav_N_reservoir * sulman_rup_veg * dt
      N_to_plant_mine = mine_N_reservoir * sulman_rup_veg * dt 
      N_to_plant_fix = fix_N_reservoir  * sulman_rup_veg * dt
   end do

      ! Passive nitrogen uptake from the rhizosphere by roots with soil water uptake:
      
      ! smin_no3_vr is per column and has to be per patch
      ! if water in layer:
   
      do j = 1,nlevdecomp
         t_soi_degC = t_soisno(c,j) - tfrz     ! Soil temperature in degrees Celcius
         if (t_soi_degC < 0.01 .and. h2osoi_liq(c,j) > 0.0_r8) then
            no3_passiv_up(p,j)= waterflux_inst%qflx_tran_veg_patch(p) * (smin_no3_vr(c,j) / h2osoi_liq(c,j))
            nh4_passiv_up(p,j) = waterflux_inst%qflx_tran_veg_patch(p) * (smin_nh4_vr(c,j) / h2osoi_liqc,j))
         end if
      enddo 
   
   
   
   end associate
   end subroutine CN_soil_veg_exchange

   !----------------------------------------

   ! FUNCTIONS

   function Vmax_myc(soil_T)
      real, intent(in)   :: soil_T                 ! Soil temperature in Kelvin
      real, parameter    :: Tref=293.15            ! Reference Temperature in Kelvin
      real, dimension(:) :: alpha
      real, dimension(:) :: Vmax_myc
      real, parameter    :: Ea=(37e3)              ! Activation energy (kJ/mol)
      real, parameter    :: R_gas = 8.314472       ! Universal gas constant, J K-1 mol-1

      alpha = params_inst%sulman_vmax_ref_mine / exp(-Ea /(R_gas*Tref))
      Vmax_myc = alpha * exp(-Ea / R_gas * soil_T)
   end function Vmax_myc

   
   function resp_myc(Cavail, myc_biomass, soil_T, wliq, wair)
      ! DO STH WITH CAVAIL
      real, intent(in) :: myc_biomass                ! mycorrhizal biomass [kgC/m2]
      real, intent(in) :: soil_T                     ! Soil temperature in Kelvin
      real, intent(in) :: wliq                       ! water liquid                        
      real, intent(in) :: wair                       ! air in soil

      real, parameter  :: enzyme_frac=1.0            ! Relative amount of enzymes produced by microbes
      real, parameter  :: substrate_diffusion_exp = 3.0    ! Exponent for theta dependence at low theta. See Davison et al DAMM model paper
      real, parameter  :: gas_diffusion_exp = 2.5          ! Exponent for gas diffusion power law dependence on theta See Meslin et al 2010, SSAJ
      real, parameter  :: min_anaerobic_resp_factor = 0.0  ! Minimum for high soil moisture Resp limitation
      real, parameter  :: min_dry_resp_factor       = 0.0  ! Minimum for low soil moisture Resp limitation

      real(r8)         :: resp_myc
      real(r8)         :: enzymes
      real(r8)         :: Cavail
      real(r8)         :: theta_func

      ! LOCAL VARIABLES:
      real(r8)         :: theta_resp_max
      real(r8)         :: aerobic_max
      
      ! From solving theta dependence for maximum:
      theta_resp_max = substrate_diffusion_exp/(gas_diffusion_exp*(1.0+substrate_diffusion_exp/gas_diffusion_exp))
      aerobic_max=theta_resp_max**substrate_diffusion_exp*(1.0-theta_resp_max)**gas_diffusion_exp

      ! Functional dependence on soil moisture, normalized so max is 1
      theta_func=(wliq**substrate_diffusion_exp)*(wair**gas_diffusion_exp)/aerobic_max
      ! On the wet side of the function, make sure it does not go below min_anaerobic_resp_factor
      if(wliq>theta_resp_max .and. theta_func<min_anaerobic_resp_factor) theta_func=min_anaerobic_resp_factor
      ! On the dry side of the function, make sure it does not go below min_dry_resp_factor
      if(wliq<theta_resp_max .and. theta_func<min_dry_resp_factor) theta_func=min_dry_resp_factor

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
   subroutine active_root_N_uptake( filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds, &
                                    soilstate_inst, froot_carbon, no3_soil, nh4_soil, no3_uptake, nh4_uptake)
      ! ! USES:
      use clm_varcon        , only : rpi
      use decompMod         , only : bounds_type
      use SoilStateType     , only: soilstate_type

      ! ! ARGUMENTS
      integer, intent(in)    :: num_soilp           ! number of soil patches in filter
      integer, intent(in)    :: filter_soilp(:)     ! filter for soil patches
      integer, intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
      integer, intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
   
      type(bounds_type)      , intent(in)    :: bounds
      type(soilstate_type)   , intent(in)    :: soilstate_inst

      real,intent(in)    :: froot_carbon(bounds%begp:bounds%endp) !,1:nlevdecomp)    ! fine root carbon (gC/m2)   [pft]   
      real,intent(in)    :: no3_soil(bounds%begc:bounds%endc,1:nlevdecomp)                                             ! avaliable soil mineral NO3 [gN/m3]
      real,intent(in)    :: nh4_soil(bounds%begc:bounds%endc,1:nlevdecomp)                                             ! avaliable soil mineral NH4 [gN/m3]
      real,intent(inout) :: no3_uptake(bounds%begp:bounds%endp,1:nlevdecomp)
      real,intent(inout) :: nh4_uptake(bounds%begp:bounds%endp,1:nlevdecomp)

      real(r8) :: root_biomass_density    ! root biomass density [g/m3]
      real(r8) :: root_cross_sec_area     ! root cross sectional area [m2]
      real(r8) :: root_length_density     ! root length density [m/m3]
      real(r8) :: rhizosphere_frac        ! fraction of rihzosphere [-] 
                                          ! sulman_r_rhiz [m] 

      real(r8), parameter :: root_radius = 0.29e-03_r8      !(m)
      real(r8), parameter :: c_to_b = 2.0_r8                ![g biomass /g C]

      associate(                                                         &
         ivt                    => patch%itype                         , & ! Input: [integer (:)] patch vegetation type
         rootfr                 => soilstate_inst%rootfr_patch         , & ! Input: [real(r8) (:,:)]
         dz                     => col%dz                              , & ! Input: [real(r8) (:,:)] layer thickness (m)
         root_radius            => pftcon%root_radius                  , & ! Input: 0.29e-03_r8 (m) 
         root_density           => pftcon%root_density                   & ! Input: 0.31e06_r8 (g biomass / m3 root) 
         )
      
         no3_uptake = 0.0_r8
         nh4_uptake = 0.0_r8
         total_active_root_N_uptake = 0.0_r8

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = patch%column(p)
         do j = 1, nlevdecomp

           ! Calculate Nitrogen concentration in soil layers
           ! smin_nh4_vr_col and smin_no3_vr_col should be already per layer and tell how much N is there           

           ! Root calculations (Sulman calculated root surface, nor sure how different that is)
           ! calculate conversion from conductivity to conductance
           root_biomass_density = c_to_b * froot_carbon(p) * rootfr(p,j) / dz(c,j)
           ! ensure minimum root biomass (using 1gC/m2)
           root_biomass_density = max(c_to_b*1._r8,root_biomass_density)
           ! Root length density: m root per m3 soil 
           root_cross_sec_area = rpi*params_inst%sulman_r_rhiz(ivt(p))**2
           root_length_density = root_biomass_density / (root_density(ivt(p)) * root_cross_sec_area)
     
           rhizosphere_frac = min(rpi*((sulman_r_rhiz+root_radius)**2-root_radius**2)*root_length_density,1.0_r8)  ! Kinda sus
           
           ! Calculate Nitrocen uptake by roots
           no3_uptake(p,j) = rhizosphere_frac * params_inst%sulman_root_no3 * (no3_soil(c,j) / (no3_soil(c,j) + params_inst%sulman_km_no3))
           nh4_uptake(p,j) = rhizosphere_frac * params_inst%sulman_root_nh4 * (nh4_soil(c,j) / (nh4_soil(c,j) + params_inst%sulman_km_nh4))
  
           ! add some if statements that if plant wants to take up more no3 /nh4 than in soil sth happens
  
           !total_active_root_N_uptake(p) = total_active_root_N_uptake(p) + no3_uptake + nh4_uptake
         end do
      end do
      end associate
   end subroutine active_root_N_uptake

   !----------------------------------------
  
   subroutine myc_scavenger_N_uptake(myc_biomass,root_dens_frac, no3_soil, nh4_soil, dz, myc_efficiency, no3_uptake, nh4_uptake)

      real,intent(in)     :: myc_biomass            ! (kgC/m2)
      real,intent(in)     :: dz                     ! (m)
      real,intent(in)     :: root_dens_frac
      real,intent(in)     :: no3_soil               ! avaliable soil mineral NO3 [gN/m3]
      real,intent(in)     :: nh4_soil               ! avaliable soil mineral NH4 [gN/m3]
      
      real,intent(inout)  :: no3_uptake             ! (kgN/m2/year)
      real,intent(inout)  :: nh4_uptake             ! (kgN/m2/year)
      real, intent(out)   :: myc_efficiency         ! [kgN/kgC per myc biomass] Nitrogen uptake efficiency of symbiont type (for ROI calculations)

      associate(                                                                  &
      sulman_k_scav_Ninorg => params_inst%sulman_k_scav_Ninorg                  , & ! Half-saturation inorganic N concentration for mycorrhizal uptake [kg N m-3]
      sulman_k_scav        => params_inst%sulman_k_scav                         , & ! Half-saturation mycorrhizal biomass concentration for scavenging [kg C m-3]
      sulman_v_scav        => params_inst%sulman_v_scav                         , & ! Maximum N uptake rate by scavenging mycorrhizae [kg N m-3 year-1]
      crootfr              => soilstate_inst%crootfr_patch                      , & ! Input:   [real(r8) (:,:)] fraction of roots for carbon in each soil layer  (nlevgrnd)
      frootc               => cnveg_carbonstate_inst%frootc_patch               , & ! Input:   [real(r8) (:)]  (gC/m2) fine root C
      leafc                => cnveg_carbonstate_inst%leafc_patch                  & ! Input:  [real(r8) (:)     ]  (gC/m2) leaf C                                    
      )

      ! Check if there is mycorrhizal biomass per patch and per layer
      ! Calculating mycorrhizal biomass per soil layer
      !myc%C_pool(p)
      !myc_c_layer = myc%C_pool(p) * root_dens(p,j)
      do j = 1, nlevdecomp
         myc_biomass_layer(p, j) = myc_biomass(p) * root_dens(p, j)
      end do

      ! In Sulman code calculations are still done if there is no biomass in the soil. They multiplu by 0.0001
      ! and creat biomass out of nowher (so that mycorrhiza can grow)
      ! THIS DOES NOT WORK IN CTSM
      ! It probably wont pass balance checks (Rosie). Instead:
      ! Make sure that mycorrhizal pool never goes to 0

      ! Check if there is  mycorrhizal biomass in soil layer
      do j = 1, nlevdecomp
         if (myc_biomass_layer >= 0) then 

            ! If there is mycorrhizal biomass in soil layer, calculate N uptake
            do j = 1, nlevdecomp
               ! sulman_v_scav is different in Sulman code
               no3_uptake = sulman_v_scav * no3_soil(c,j) / (no3_soil(c,j) + sulman_k_scav_Ninorg) * &
                           scav_C_biomass(c,j) / (scav_C_biomass + sulman_k_scav)

               nh4_uptake = sulman_v_scav * nh4_soil(c,j) / (nh4_soil(c,j) + sulman_k_scav_Ninorg) * &
                           scav_C_biomass(c,j) / (scav_C_biomass(c,j) + sulman_k_scav)

               ! NO3 and NH4 uptake depends on how much N is avaliable in soil
               no3_uptake = min(no3_uptake(c,j), no3_soil(c,j)) 
               nh4_uptake = min(nh4_uptake(c,j), nh4_soil(c,j))
               
               myc_efficiency = (no3_uptake(c,j) + nh4_uptake(c,j)) / scav_C_biomass(c,j)

               ! Total uptake is scaled 
            end do
         end if 

         if ( (leafc + frootc) == 0) then ! leaves, fine roots, and sapwood biomass find a variable for living biomass in soil
            no3_uptake     = 0.0_r8
            nh4_uptake     = 0.0_r8
            myc_efficiency = 0.0_r8
         end if
      end do
      end associate

   end subroutine myc_scavenger_N_uptake

   !----------------------------------------

   ! Miner
   subroutine myc_miner_N_uptake (temperature_inst, root_dens_frac, soil_carbon, myc_biomass, dz, total_N_uptake, total_myc_resp, dt)

      ! USES
      use TemperatureType                 , only : temperature_type
      use clm_varcon                      , only : tfrz
      use clm_varcon                      , only : denh2o, denice
       

      ! ARGUMENTS
      type(temperature_type)              , intent(in)    :: temperature_inst

      real(r8)                :: t_soi_degC              ! Soil temperature            (degrees Celsius)
      real,intent(in)         :: root_dens_frac
      real, intent(in)        :: myc_biomass             ! (kgC/m2)
      real, intent(in)        :: soil_carbon             ! 
      real(r8), intent(in)    :: dz                      ! soil layer thickness 
      real(r8), intent(in)    :: dt                      ! decomposition time step
      
      real(r8), intent(inout) :: total_N_uptake
      real(r8), intent(inout) :: total_myc_resp

      ! LOCAL VARIABLES:
      real(r8) :: myc_biomass_layer ! Mycorrhizal biomass per layer
      real(r8) :: wliq ! water liquid
      real(r8) :: wice ! water ice
      real(r8) :: wair ! air in soil

      associate(                                                  &
         t_soisno       => temperature_inst%t_soisno_col        , & ! Input:  [real(r8) (:,:)]soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         watsat         => soilstate_inst%watsat_col            , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity)  
         h2osoi_liq     => waterstatebulk_inst%h2osoi_liq_col   , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)             
         h2osoi_ice     => waterstatebulk_inst%h2osoi_ice_col   , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2) (new)
         t_soisno       => temperature_inst%t_soisno_col          & ! Input:   [real(r8) (:,:)] soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         )

         ! Check if there is mycorrhizal biomass per patch and per layer   
         !myc%C_pool(p)
         !myc_c_layer = myc%C_pool(p) * root_dens(p,j)
         do j = 1, nlevdecomp
            myc_biomass_layer(p, j) = myc_biomass(p) * root_dens(p, j)
         end do

         ! Check if there is  mycorrhizal biomass in soil layer

         ! In Sulman code calculations are still done if there is no biomass in the soil. They multiplu by 0.0001
         ! and creat biomass out of nowher (so that mycorrhiza can grow)
         ! THIS DOES NOT WORK IN CTSM
         ! It probably wont pass balance checks (Rosie). Instead:
         ! Make sure that mycorrhizal pool never goes to 0

         do j = 1,nlevdecomp
            if (myc_biomass_layer >= 0) then 
               mine_biomass = myc_biomass_layer    ! DOES THIS HAVE STH TO DO WITH THE mine_C_biomass calculations above?
            else 
            N_uptake       = 0.0_r8
            mine_resp      = 0.0_r8
            myc_efficiency = 0.0_r8
            end if 
         end do
      
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
   
         ! Make sure fluxes are zero before the loop
         total_N_uptake       = 0.0_r8
         total_myc_resp       = 0.0_r8
   
         do j = 1, nlevdecomp
            ! SOMc to miner
            call miner_decomposition(decomp_cpools_vr(c,j,i_chem_som), &
            (decomp_npools_vr(c,j,i_chem_som) + decomp_npools_vr(c,j,i_phys)), &
            myc_biomass, t_soisno(c,j), wliq(c,j), wair(c,j), N_uptake, mine_resp)

            ! SOMp to miner
            call miner_decomposition(decomp_cpools_vr(c,j,i_phys), &
            (decomp_npools_vr(c,j,i_chem_som) + decomp_npools_vr(c,j,i_phys)), &
            myc_biomass, t_soisno(c,j), wliq(c,j), wair(c,j), N_uptake, mine_resp)
   
            total_myc_resp = total_myc_resp + myc_resp ! NOT SURE
         end do
   
         myc_efficiency = N_uptake(c,j) / mine_C_biomass(c,j)
   
   
        ! Organic N Mining by Mycorrhizal Fungi
      end associate

   end subroutine myc_miner_N_uptake

   !----------------------------------------
  
   subroutine miner_decomposition(soil_carbon, soil_nitrogen, myc_biomass, soil_T, soil_water, soil_air, N_uptake, myc_resp)

      real, intent(in) :: soil_carbon
      real, intent(in) :: soil_nitrogen
      real, intent(in) :: myc_biomass
      real, intent(in) :: soil_T
      real, intent(in) :: soil_water ! water liquid
      real, intent(in) :: soil_air   ! air in soil
      
      real, intent(inout) :: N_uptake  ! Nitrogen taken up by plant
      real, intent(inout) :: myc_resp  ! Mycorrhizal repiration (CO2) during symbiosis

      ! LOCAL VARIABLES:
      real(r8) :: wliq ! water liquid
      real(r8) :: wice ! water ice
      real(r8) :: wair ! air in soil

      associate(                                                  &
         t_soisno       => temperature_inst%t_soisno_col        , & ! Input:  [real(r8) (:,:)]soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         watsat         => soilstate_inst%watsat_col            , & ! Input:  [real(r8) (:,:)  ]  volumetric soil water at saturation (porosity)  
         h2osoi_liq     => waterstatebulk_inst%h2osoi_liq_col   , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2) (new)             
         h2osoi_ice     => waterstatebulk_inst%h2osoi_ice_col   , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2) (new)
         t_soisno       => temperature_inst%t_soisno_col          & ! Input:   [real(r8) (:,:)] soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         )

      myc_resp    = 0.0_r8
      N_uptake    = 0.0_r8

      ! Calculating water, ice and air content in soil (equivalent to air_filled porosity, theta, theta sat in Sulman)
      wliq = h2osoi_liq / dz * denh2o
      wice = h2osoi_ice / dz * denice
      wliq = min(1.0_r8, wliq/watsat)            ! fraction of liquid water-filled pore space (0.0 - 1.0)
      wice = min(1.0_r8, wice/watsat)            ! fraction of frozen water-filled pore space (0.0 - 1.0)
      wair = max(0.0_r8, 1.0_r8 - wliq- wice)    ! fraction of air-filled pore space (0.0 - 1.0)

      potential_tempResp=Resp_myc(soil_carbon,myc_biomass,T,wliq,wair)

      ! Don't exceed avaliable C
      if(dt*potential_tempResp > soil_carbon) then
         potential_tempResp = soil_carbon / dt
      end if

      if(soil_carbon > 0) then
         pot_tempN_decomposed = potential_tempResp * soil_carbon / soil_nitrogen
      else 
         pot_tempN_decomposed=0.0
      end if

      N_uptake = N_uptake + sum(pot_tempN_decomposed*params_inst%sulman_nue_mine)*dt
      myc_resp = CO2prod + sum(potential_tempResp*(1 - params_inst%sulman_cue_mine))*dt

    end associate
   end subroutine miner_decomposition
   
   !----------------------------------------

   ! GROWTH AND TURNOVER OF SYMBIOTIC BIOMASS
   ! Intermediate pools:
   ! line 553 - 604 in veg_dynamics_mod

   subroutine update_symbionts(this)
   
   class(symbiont_type)                                    :: this

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

   real(r8) :: scav_C_reservoir        ! Carbon reservoir for scavangers
   real(r8) :: scav_N_reservoir        ! Nitrogen reservoir for scavangers
   real(r8) :: mine_C_reservoir        ! Carbon reservoir for miners
   real(r8) :: mine_N_reservoir        ! Nitrogen reservoir for miners
   real(r8) :: fix_C_reservoir         ! Carbon reservoir for fixers   
   real(r8) :: fix_N_reservoir         ! Nitrogen reservoir for fixers
   real(r8) :: N_fixation              !

   real(r8) :: reservoir_C_leakage

   real(r8) :: symbiont_turnover_C
   real(r8) :: symbiont_turnover_N
   real(r8) :: fix_C_alloc
   real(r8) :: scav_C_alloc
   real(r8) :: mine_C_alloc
   real(r8) :: fix_N_alloc
   real(r8) :: scav_N_alloc
   real(r8) :: mine_N_alloc
   real(r8) :: root_exudate_N_frac
   real(r8) :: total_plant_N_uptake
   real(r8) :: plant_N_storage         ! sum of plant N store from (leafn_storage + frootn + frootn_storage)

   associate(                                                                  &
   leafn_storage        => cnveg_nitrogenstate_inst%leafn_storage_patch      , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
   frootn               => cnveg_nitrogenstate_inst%frootn_patch             , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N                               
   frootn_storage       => cnveg_nitrogenstate_inst%frootn_storage_patch       & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
   
   )

   ! Root uptake of nitrogen
   ! patch loop whithin
   ! loops are inside routines
   do p = bounds%begp,bounds%endp

      ! Nitrogen limitation:
         ! Nitrogen is limited in individual uptake pathways, but also needs to be limited as a total
         ! Start with scavengers taking up N from the inorganic N pool
         ! Then active root uptake takes up what is left in inorganic N pool
         ! Then passive root uptake takes up what is left in inorganic N pool
   
      call active_root_N_uptake(froot_carbon(begp:endp), smin_no3_to_plant_vr(begc:endc,1:nlevdecomp), &
                                 smin_nh4_to_plant_vr(begc:endc,1:nlevdecomp), &
                                 no3_active_up(begc:endc,1:nlevdecomp), nh4_active_up(begc:endc,1:nlevdecomp))
   
      ! Scavenging (AM-style)
      ! C_reservoir(p,i_scav)
      call myc_scavenger_N_uptake(C_biomass(p,i_scav), smin_no3_to_plant_vr(begc:endc,1:nlevdecomp), smin_nh4_to_plant_vr(begc:endc,1:nlevdecomp), &
                                    myc_efficiency, no3_scav_up(begc:endc,1:nlevdecomp), nh4_scav_up(begc:endc,1:nlevdecomp))
      
      ! Mycorrhizal N mining (ECM-style)
      call myc_miner_N_uptake(temperature_inst, decomp_cpools_vr(c,j,i_chem_som), mine_C_biomass, &
                              dz, n_somc_miner_up(begc:endc,1:nlevdecomp), total_myc_resp, dt)
      call myc_miner_N_uptake(temperature_inst, decomp_cpools_vr(c,j,i_phys_som), mine_C_biomass, &
                              dz, n_somp_miner_up(begc:endc,1:nlevdecomp), total_myc_resp, dt)
      
      ! Symbiotic N2 Fixation
      N_fixation = fix_C_biomass*sulman_rfix*dt ! sulman_rfix different values
   
      enddo
    
      !----------------------------------------------------------------
   

   ! Add mycorrhizal and fixer turnover (necromass) in SOM pools:
      symbiont_turnover_C = (mine_C_biomass / params_inst%sulman_tau_mine + scav_C_biomass / params_inst%sulman_tau_scav + &
                              fix_C_biomass / params_inst%sulman_tau_fix) * dt * params_inst%sulman_tau_sym
      symbiont_turnover_N = (mine_N_biomass / params_inst%sulman_tau_mine + scav_N_biomass / params_inst%sulman_tau_scav + &
                              fix_N_biomass / params_inst%sulman_tau_fix) * dt * params_inst%sulman_tau_sym
   
      myc_resp = myc_resp + (1.0 - params_inst%sulman_tau_sym) * symbiont_turnover_C / sulman_tau_sym
   

   ! Update biomass of symbionts
   ! GROWTH AND TURNOVER
 
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
    ! Here I could make myc stop from dying if the are very little myc they don't die!!!!!!!!!!!!!!!!!!!
    ! To adress problem from line 597
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
 
    ! ADD SOME ERRORS IN CASE BIOMASS IS ZERO

 

    ! Excess C leaks out of reservoir into root exudates at a time scale of one day
    reservoir_C_leakage = 0.0_r8
    reservoir_C_leakage = reservoir_C_leakage + (scav_C_reservoir + mine_C_reservoir + fix_C_reservoir)*dt*365
    scav_C_reservoir = scav_C_reservoir - scav_C_reservoir*dt*365
    mine_C_reservoir = mine_C_reservoir - mine_C_reservoir*dt*365
    nfix_C_reservoir = nfix_C_reservoir - nfix_C_reservoir*dt*365


  
    ! RETURN OF INVESTMENT
   ! Calculate N released to plants per uptake pathway - Return Of Investment (line 553 in vegn_dynamics)
   ! multiplied with secspday * days_per_year to get values per second (parameter is in per year)

   days_per_year = get_average_days_per_year()     ! to get average number of days per year 

   ! Scavengers
   if (is_active(p,i_scav)) then 
      scav_N_to_plant = scav_N_reservoir * sulman_rup_veg * secspday * days_per_year
      if (scav_C_biomass < 0.0_r8) then  ! or (C_biomass(p,i_scav) < 0.0_r8)
         scav_roi = (max(0.0_r8, scav_N_to_plant) * dt) / (scav_C_biomass / sulman_growth_scav / sulman_tau_scav)
      else 
      ! scav_efficiency is calculated in myc_scavenger_N_uptake under myc_efficiency
      scav_roi = scav_efficiency / (dt * sulman_growth_scav * sulman_tau_scav) 
      end if 
   else
      scav_N_to_plant = 0.0_r8 ; scav_roi = 0.0_r8
   end if 

   scav_N_reservoir = scav_N_reservoir - scav_N_to_plant

   ! Miners
   if (is_active(p,i_miner)) then
      mine_N_to_plant = mine_N_reservoir * sulman_rup_veg * secspday * days_per_year
      if (mine_C_biomass < 0.0_r8) then 
         mine_roi = (max(0.0_r8, mine_N_to_plant) * dt) / (mine_C_biomass / sulman_growth_mine / sulman_tau_mine)
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
   if (local_active(p,i_fixer)) then
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


    ! Adding Smoothing filters to avoid abrupt N changes


   ! Calculate relative fractions
   if (scav_roi + mine_roi + fix_roi + root_roi > 0) then
      roi = scav_roi + mine_roi + fix_roi + root_roi

      scav_roi_frac = scav_roi / roi
      mine_roi_frac = mine_roi / roi
      fix_roi_frac = fix_roi / roi
      root_roi_frac = root_roi / roi
   else
      scav_roi_frac = 0.4*0.7
      mine_roi_frac = 0.3*0.7
      fix_roi_frac = 0.3*0.7
      root_roi_frac = 0.3
   endif 

   ! Calculate carbon allocation to pathways (without smoothing filters)
    fix_C_alloc = C_allocation_to_N_acq * fix_roi_frac * dt
    mine_C_alloc = C_allocation_to_N_acq * mine_roi_frac * dt
    scav_C_alloc = C_allocation_to_N_acq * scav_roi_frac * dt

   ! We are accumulating potential (not limited/smoothed) allocation so max allocation can increase over time
    fix_alloc_accum = Nfix_alloc_accum + C_allocation_to_N_acq * fix_roi_frac * dt
    mine_alloc_accum = mine_alloc_accum + C_allocation_to_N_acq * mine_roi_frac * dt
    scav_alloc_accum = scav_alloc_accum + C_allocation_to_N_acq * scav_roi_frac * dt

    ! real    :: root_exudate_N_frac = 0.0 ! N fraction of root exudates. See e.g. Drake et al 2013
    
    fix_N_alloc = 0.0
    mine_N_alloc = mine_C_alloc*sp%root_exudate_N_frac
    scav_N_alloc = scav_C_alloc*sp%root_exudate_N_frac

    plant_N_storage = leafn_storage + frootn_storage + livecrootn_storage

    ! Make sure N allocation from plant to symbionts does not completely deplete stored N
    if (plant_N_storage <= 0.0) then 
      fix_N_alloc=0.0
      mine_N_alloc=0.0
      scav_N_alloc=0.0
    elseif (fix_N_alloc + mine_N_alloc + scav_N_alloc > 0.0 .AND. &
      fix_N_alloc + mine_N_alloc + scav_N_alloc > plant_N_storage * 0.1) then

      lim_factor = plant_N_storage * 0.1 / (fix_N_alloc + mine_N_alloc + scav_N_alloc)

      fix_N_alloc=nfix_N_alloc*lim_factor
      mine_N_alloc=mine_N_alloc*lim_factor
      scav_N_alloc=scav_N_alloc*lim_factor
    endif

    scav_C_reservoir = scav_C_reservoir + scav_C_alloc
    scav_N_reservoir = scav_N_reservoir + scav_N_alloc
    mine_C_reservoir = mine_C_reservoir + mine_C_alloc
    mine_N_reservoir = mine_N_reservoir + mine_N_alloc
    fix_C_reservoir =  fix_C_reservoir + fix_C_alloc
    fix_N_reservoir =  fix_N_reservoir + fix_N_alloc

    total_plant_N_uptake = scav_N_to_plant + mine_N_to_plant + fix_N_to_plant + root_N_uptake
    plant_N_storage = plant_N_storage + total_plant_N_uptake
    stored_N_loss = max(0.0,min(plant_N_storage), C_allocation_to_N_acq * dt *root_exudate_N_frac)
    plant_N_storage = plant_N_storage - stored_N_loss
    root_exudate_N = stored_N_loss - scav_N_alloc - fix_N_alloc - mine_N_alloc

    root_exudate_C = C_allocation_to_N_acq * dt - scav_C_alloc - mine_C_alloc -  fix_C_alloc + reservoir_C_leakage

   end associate

   end subroutine update_symbionts

end module CNSoilVegMIMICSplus
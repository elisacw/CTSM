module CNSoilVegMIMICSplus

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module connects the Soil decomposition module MIMICS+ (Aas et al. 2023) with the vegetation through the 
  ! symbiosis between mycorrhizal fungi and plants.
  ! Coupling follows Sulman et al. (2019)
  
  ! !USES:
  use shr_kind_mod                        , only : r8 => shr_kind_r8
  use shr_infnan_mod                      , only :  isnan => shr_infnan_isnan
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
  public  :: update_symbionts             !

  real(r8) :: dt                          ! decomp timestep [s]

  ! !FUNCTIONS:
  private :: resp_myc                     ! Respiration of mycorrhiza
  private :: Vmax_myc                     ! Michaelis Menten Kinetics for mycorrhiza

  ! ! PUBLIC DATA
  integer, public :: i_fixer = 1 
  integer, public :: i_scav  = 2 
  integer, public :: i_miner = 3
  integer, public :: n_symb  = 3
  
  type, public :: symbiont_type

  real(r8), pointer           :: C_biomass             (:,:) ! [patch,n_symb] Carbon biomass    [gC/m2]
  real(r8), pointer           :: N_biomass             (:,:) ! [patch,n_symb] Nitrogen biomass  [gN/m2]
  real(r8), pointer           :: C_reservoir           (:,:) ! [patch,n_symb] Carbon intermediate pool biomass    [gC/m2]
  real(r8), pointer           :: N_reservoir           (:,:) ! [patch,n_symb] Nitrogen intermediate pool biomass  [gN/m2]
  real(r8), pointer           :: symb_eff              (:,:) ! [patch,n_symb] Symbiont efficiency when its biomass is 0
  real(r8), pointer           :: symb_growth           (:,:) ! [patch,n_symb] Symbiotic biomass growth rate [gC/m2]
  real(r8), pointer           :: C_mortality           (:,:) ! [col,nlevdecomp]Turnover of symbionts per layer and column [gC/m3/s]
  real(r8), pointer           :: N_mortality           (:,:) ! Turnover of symbionts per layer and column                 [gN/m3/s]
  real(r8), pointer           :: N_mine_somc2soma_col  (:,:) ! Leftover part of co-mineralized N, not taken up by miners  [gN/m3/s]
  real(r8), pointer           :: N_mine_somp2soma_col  (:,:) ! Leftover part of co-mineralized N, not taken up by miners  [gN/m3/s]
  real(r8), pointer           :: somc_nuptake_col      (:,:) ! Nitrogen uptake from SOMc pool per column  [gN/m3/s]
  real(r8), pointer           :: somp_nuptake_col      (:,:) ! Nitrogen uptake from SOMp pool per column  [gN/m3/s]
  real(r8), pointer           :: somc_cuptake_col      (:,:) ! Nitrogen uptake from SOMc pool per column  [gC/m3/s]
  real(r8), pointer           :: somp_cuptake_col      (:,:) ! Nitrogen uptake from SOMp pool per column  [gC/m3/s]
  real(r8), pointer           :: root_exudate_C_col    (:,:) ! Leftover C from allocation to symbionts    [gC/m2/s]
  logical, pointer            :: is_active             (:,:) ! If the symbiont type is active
  character(len=10), pointer  :: symb_name             (:)   ! Symbiont name
  character(len=10), pointer  :: symb_hist_name        (:)   ! Symbiont name on history tapes
   
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

    logical(r8) :: local_active(bounds%begp:bounds%endp, 1:n_symb)   ! Indicates active pathway, based on symbiont type of PFT
   
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

    allocate(this%symb_eff(begp:endp,1:n_symb)) ; this%symb_eff(begp:endp,1:n_symb)       = 0.0_r8
    allocate(this%symb_growth(begp:endp,1:n_symb)) ; this%symb_growth(begp:endp,1:n_symb) = 0.0_r8
    allocate(this%root_exudate_C_col(begc:endc,1:nlevdecomp)) ; this%root_exudate_C_col(begc:endc,1:nlevdecomp) = 0.0_r8

    allocate(this%C_mortality(begc:endc,1:nlevdecomp)) ; this%C_mortality(begc:endc,1:nlevdecomp) = 0.0_r8
    allocate(this%N_mortality(begc:endc,1:nlevdecomp)) ; this%N_mortality(begc:endc,1:nlevdecomp) = 0.0_r8

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
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_C_inter', units='gC/m2', &
         avgflag='A', long_name=('Carbon intermediate pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%N_biomass(begp:endp,i) = spval
         data1dptr => this%N_biomass(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_N', units='gN/m2', &
         avgflag='A', long_name=('Nitrogen pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='active')
   
         this%N_reservoir(begp:endp,i) = spval
         data1dptr => this%N_reservoir(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_N_inter', units='gN/m2', &
         avgflag='A', long_name=('Nitrogen intermediate pool of '//this%symb_name(i)//' symbionts'), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')

         this%symb_growth(begp:endp,i) = spval
         data1dptr => this%symb_growth(:,i)
         call hist_addfld1d (fname=trim('SYM_'//this%symb_hist_name(i))//'_C_GROWTH', units='gC/m2', &
         avgflag='A', long_name=('Symbiont growth of '//this%symb_name(i)), &
         ptr_patch=data1dptr, set_spec=spval, default='inactive')
         
      end do

         ! CHANGE NAME & CHECK UNIT
         this%C_mortality(begc:endc,1:nlevdecomp) = spval
         data2dptr => this%C_mortality
         call hist_addfld2d (fname='C_MORTALITY_SYMB_MIMICSPLUS', units='gC/m^2', type2d='levsoi', &
         avgflag='A', long_name='Symbiotic C turnover per soil layer and column', &
         ptr_col= this%C_mortality, set_spec=spval, default='inactive')
   
         this%N_mortality(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MORTALITY_SYMB_MIMICSPLUS', units='gN/m^2', type2d='levsoi', &
         avgflag='A', long_name='Symbiotic N turnover per soil layer and column', &
         ptr_col=this%N_mortality, set_spec=spval, default='inactive')
   
         this%N_mine_somc2soma_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MINE_SOMC_TO_SOMA_MIMICSPLUS', units='gN/m^2', type2d='levsoi',&
         avgflag='A', long_name='Leftover N from SOMc mineralization (not taken up by miners)', &
         ptr_col=this%N_mine_somc2soma_col, set_spec=spval, default='inactive')
   
         this%N_mine_somp2soma_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MINE_SOMP_TO_SOMA_MIMICSPLUS', units='gN/m^2', type2d='levsoi', &
         avgflag='A', long_name='Leftover N from SOMp mineralization (not taken up by miners)', &
         ptr_col=this%N_mine_somp2soma_col, set_spec=spval, default='inactive')
   
         this%somc_nuptake_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MINE_UPTAKE_SOMC_MIMICSPLUS', units='gN/m^2', type2d='levsoi', &
         avgflag='A', long_name='N uptake from SOMc pool via mining', &
         ptr_col=this%somc_nuptake_col, set_spec=spval, default='inactive')
   
         this%somp_nuptake_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='N_MINE_UPTAKE_SOMP_MIMICSPLUS', units='gN/m^2', type2d='levsoi',&
         avgflag='A', long_name='N uptake from SOMp pool via mining', &
         ptr_col=this%somp_nuptake_col, set_spec=spval, default='inactive')
   
         this%somc_cuptake_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='C_MINE_UPTAKE_SOMC_MIMICSPLUS', units='gC/m^2',type2d='levsoi', &
         avgflag='A', long_name='co-decomposed C from SOMc pool during mining', &
         ptr_col=this%somc_cuptake_col, set_spec=spval, default='inactive')
   
         this%somp_cuptake_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='C_MINE_UPTAKE_SOMP_MIMICSPLUS', units='gC/m^2', type2d='levsoi', &
         avgflag='A', long_name='co-decomposed C from SOMp pool during mining', &
         ptr_col=this%somp_cuptake_col, set_spec=spval, default='inactive')
   
         this%root_exudate_C_col(begc:endc,1:nlevdecomp) = spval
         call hist_addfld2d (fname='ROOT_EXUDATE_C_MIMICSPLUS', units='gC/m^2/s', type2d='levsoi', &
         avgflag='A', long_name='Leftover root exudate C from symbiont allocation', &
         ptr_col=this%root_exudate_C_col, set_spec=spval, default='inactive')

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
               this%C_biomass(p,i) = params_inst%sulman_initial_C_stocks(i) * 0.9_r8
               this%C_reservoir(p,i) = params_inst%sulman_initial_C_stocks(i) * 0.1_r8
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

      !ECW DOES THIS HAVE TO BE = spval
      this%C_mortality(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8
      this%N_mortality(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8
      this%N_mine_somc2soma_col(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8
      this%N_mine_somp2soma_col(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8
      this%somc_nuptake_col(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8
      this%somp_nuptake_col(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8
      this%somc_cuptake_col(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8
      this%somp_cuptake_col(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8
      this%root_exudate_C_col(bounds%begc:bounds%endc,1:nlevdecomp) = 0.0_r8

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

      data1dptr => this%symb_eff(:,i)
      call restartvar(ncid=ncid, flag=flag, varname=trim('_efficiency'//this%symb_name(i)), xtype=ncd_double,  &
      dim1name='pft', long_name=('Symbiont'//this%symb_name(i)//' efficiency'), units='g/m2', &
      interpinic_flag='interp', readvar=readvar, data=data1dptr)
      if (flag == 'read' .and. (.not. readvar)) then
         call endrun(msg = "ERROR: Efficiency "//this%symb_name(i)// " is not on the restart file."// & 
             errMsg(sourcefile, __LINE__))
      endif
   enddo
      
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

       write(iulog,*) 'Biomass: ', this%C_biomass(p,i_scav)
       write(iulog,*) 'total patch: ',  total_patch(p)


      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
       total_patch(bounds%begp:bounds%endp), &
       total_col(bounds%begc:bounds%endc))
       
      total(bounds%begc:bounds%endc)  = total_col(bounds%begc:bounds%endc) 
     
   end subroutine Summary
   !-----------------------------------------------------------------------

   subroutine CN_soil_veg_exchange (filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds, symbiont_inst, &
      cnveg_nitrogenstate_inst, waterstatebulk_inst, temperature_inst, cnveg_carbonflux_inst, &
      soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst, cnveg_state_inst, &
      waterfluxbulk_inst, soilstate_inst, cnveg_carbonstate_inst, soilbiogeochem_carbonstate_inst, cnveg_nitrogenflux_inst)

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

   integer                                , intent(in)    :: num_soilp           ! number of soil patches in filter
   integer                                , intent(in)    :: filter_soilp(:)     ! filter for soil patches
   integer                                , intent(in)    :: num_bgc_soilc       ! number of soil columns in filter
   integer                                , intent(in)    :: filter_bgc_soilc(:) ! filter for soil columns
   
   ! !LOCAL VARIABLES
   integer :: p, fp, c, fc, j, k, l, s, i, g
   integer :: begp, endp, begc, endc
   
   real(r8) :: root_dens_sum                                            ! Fine root C per layer             [gC/m2]
   real(r8) :: t_soi_degC                                               ! Soil temperature                  [degrees Celcius]
   real(r8) :: availc_alloc(bounds%begp:bounds%endp)                    ! The avaible C pool for allocation [gC/m2/s]
   real(r8) :: C_allocation_to_N_acq(bounds%begp:bounds%endp)           ! C allocated to N acquisition      [gC/m2/s]
   real(r8) :: root_dens_frac(bounds%begp:bounds%endp,1:nlevdecomp)     ! Fraction of root density          [-]
   
   ! Nitrogen uptake variables for pathways into intermediated pools
   real(r8) :: smin_no3_avail(bounds%begp:bounds%endp, 1:nlevdecomp)    ! no3 available for uptake per soil layer      [gN/m2]
   real(r8) :: smin_nh4_avail(bounds%begp:bounds%endp, 1:nlevdecomp)    ! nh4 available for uptake per soil layer      [gN/m2]
   real(r8) :: smin_no3_avail_col(bounds%begc:bounds%endc, 1:nlevdecomp) ! col no3 available for uptake per soil layer [gN/m2]
   real(r8) :: smin_nh4_avail_col(bounds%begc:bounds%endc, 1:nlevdecomp) ! col nh4 available for uptake per soil layer [gN/m2]
      
   real(r8) :: no3_passiv_up(bounds%begp:bounds%endp,1:nlevdecomp)   ! Passive root NO3 (nitrate) uptake  [gN/m2/s]
   real(r8) :: nh4_passiv_up(bounds%begp:bounds%endp,1:nlevdecomp)   ! Passive root NH4 (ammonium) uptake [gN/m2/s]
   real(r8) :: no3_active_up(bounds%begp:bounds%endp,1:nlevdecomp)   ! Active root NO3 (nitrate) uptake   [gN/m2/s]
   real(r8) :: nh4_active_up(bounds%begp:bounds%endp,1:nlevdecomp)   ! Active root NH4 (ammonium) uptake  [gN/m2/s]
   real(r8) :: no3_scav_up(bounds%begp:bounds%endp,1:nlevdecomp)     ! Scavenger NO3 (nitrate) uptake     [gN/m2/s]
   real(r8) :: nh4_scav_up(bounds%begp:bounds%endp,1:nlevdecomp)     ! Scavenger NH4 (ammonium) uptake    [gN/m2/s]
   real(r8) :: N_fixation(bounds%begp:bounds%endp)                   ! Nitrogen uptake from fixation      [gN/m2/s]

   real(r8) :: somc_nuptake(bounds%begp:bounds%endp, 1:nlevdecomp)   ! Nitrogen uptake from SOMc pool by miners   [gN/m3/s]
   real(r8) :: somp_nuptake(bounds%begp:bounds%endp, 1:nlevdecomp)   ! Nitrogen uptake from SOMp pool by miners   [gN/m3/s]
   real(r8) :: somc_cuptake(bounds%begp:bounds%endp, 1:nlevdecomp)   ! Carbon uptake from SOMc pool by miners     [gC/m3/s]
   real(r8) :: somp_cuptake(bounds%begp:bounds%endp, 1:nlevdecomp)   ! Carbon uptake from SOMp pool by miners     [gC/m3/s]
     
   real(r8) :: maint_resp                              ! Maintainace respiration, C from symbiont pool used to sustain existing biomass [gC/m2/s]
   real(r8) :: symb_CO2_prod(bounds%begp:bounds%endp)  ! Hetrotrophic respiration, C loss during growth of symbionts            [gC/m2/s]
   
   real(r8) :: myc_biomass_layer(bounds%begp:bounds%endp, 1:nlevdecomp)           ! Mycorrhyzal biomass per soil layer       [gC/m2]

   real(r8) :: total_symbiont_turnover_C(bounds%begp:bounds%endp, 1:n_symb)       ! Part of symbiont turnover going into SOM [gC/m2/s]
   real(r8) :: symbiont_turnover_C(bounds%begp:bounds%endp, 1:n_symb)             ! Part of symbiont turnover going into SOM [gC/m3/s]
   real(r8) :: symbiont_turnover_N(bounds%begp:bounds%endp, 1:n_symb)             ! Part of symbiont turnover going into SOM [gN/m3/s]
   real(r8) :: symb_turnover_layer_C(bounds%begp:bounds%endp, 1:nlevdecomp)       ! Symbiotic turnover per soil layer        [gC/m3/s]
   real(r8) :: symb_turnover_layer_N(bounds%begp:bounds%endp, 1:nlevdecomp)       ! Symbiotic turnover per soil layer        [gN/m3/s]
   real(r8) :: N_mine_somc2soma(bounds%begp:bounds%endp, 1:nlevdecomp)            ! Leftover co-mineralized N, not taken up by miners  [gN/m3/s]
   real(r8) :: N_mine_somp2soma(bounds%begp:bounds%endp, 1:nlevdecomp)            ! Leftover co-mineralized N, not taken up by miners  [gN/m3/s]

   real(r8) :: root_exudate_C(bounds%begp:bounds%endp)                            ! Leftover C from allocation to symbionts    [gC/m2/s]
   real(r8) :: root_exudate_C_layer(bounds%begp:bounds%endp,1:nlevdecomp)         ! Leftover C from allocation to symbionts    [gC/m2/s]
   real(r8) :: N_biomass_old(bounds%begp:bounds%endp)                             ! temporary variable, to calculate delta N bimass for fixers [gN/m2]
   real(r8) :: root_N_uptake(bounds%begp:bounds%endp)                             ! active root N uptake [gN/m2/s]
   real(r8) :: root_N_to_plant(bounds%begp:bounds%endp)                           ! toatl (active + passive) root N uptake [gN/m2/s]
   

   ! VARIABLES FOR LOCAL BALANCE CHECK
   real(r8) :: old_C_biomass_scav(bounds%begp:bounds%endp)
   real(r8) :: old_C_biomass_mine(bounds%begp:bounds%endp)
   real(r8) :: old_C_biomass_fix(bounds%begp:bounds%endp)
   real(r8) :: old_C_reservoir_scav(bounds%begp:bounds%endp)
   real(r8) :: old_C_reservoir_mine(bounds%begp:bounds%endp)
   real(r8) :: old_C_reservoir_fix(bounds%begp:bounds%endp)


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
   sulman_rfix          => params_inst%sulman_rfix            , &   ! N fixation rate per unit symbiotic biomass             [gN per gC per s]
   sulman_kgrowth       => params_inst%sulman_kgrowth         , &   ! Half-saturation of intermediate C pool for symbiotic growth      [gC/m2]
   sulman_max_symb_growth => params_inst%sulman_max_symb_growth, &  ! Maximum symbiont growth rate                                   [gC/m2/s]
   sulman_tau_sym       => params_inst%sulman_tau_sym         , &   ! Fraction of symbiotic biomass turnover not used for maint resp       [-]
   sulman_growth_scav   => params_inst%sulman_growth_scav     , &   ! N scavenger growth efficiency [-]
   sulman_growth_mine   => params_inst%sulman_growth_mine     , &   ! N miner growth efficiency     [-]
   sulman_growth_fix    => params_inst%sulman_growth_fix      , &   ! N fixer growth efficiency     [-]
   sulman_tau_scav      => params_inst%sulman_tau_scav        , &   ! N scavenger turnover time     [s]
   sulman_tau_mine      => params_inst%sulman_tau_mine        , &   ! N miner turnover time         [s]
   sulman_tau_fix       => params_inst%sulman_tau_fix         , &   ! N fixer turnover time         [s]
   sulman_cn_scav       => params_inst%sulman_cn_scav         , &   ! N scavenger C:N               [-]
   sulman_cn_mine       => params_inst%sulman_cn_mine         , &   ! N miner C:N                   [-]
   sulman_cn_fix        => params_inst%sulman_cn_fix          , &   ! N fixer C:N                   [-]
   sulman_tau_int       => params_inst%sulman_tau_int         , &   ! Turnover time of intermediate C pool                    [s]
   sulman_rup_veg       => params_inst%sulman_rup_veg         , &   ! Vegetation N uptake rate from intermediate N pool       [s]
   sulman_fnalloc       => params_inst%sulman_fnalloc         , &   ! Fraction of NPP allocated to N uptake per unit N stress [-] 
   symb_tau_somc        => params_inst%symb_tau_somc          , &   ! Fractiom pf symbiont turnover into SOMc                 [-]
   symb_tau_soma        => params_inst%symb_tau_soma          , &   ! Fractiom pf symbiont turnover into SOMa                 [-]
   symb_tau_somp        => params_inst%symb_tau_somp          , &   ! Fractiom pf symbiont turnover into SOMp                 [-]
      
   sulman_initial_C_stocks => params_inst%sulman_initial_C_stocks, &! Initial carbon stocks of symbionts (coldstart)      [gC/m2]
   sulman_cn_symbionts     => params_inst%sulman_cn_symbionts , &   ! C:N ratios of fixer, miner, scavengers as array         [-]
   symb_tau_som            => params_inst%symb_tau_som        , &   ! Fraction symbiont necromass into SOM pools              [-]

   crootfr              => soilstate_inst%crootfr_patch                      , & ! Input: (:,:)     (-) patch fraction of roots for carbon in each soil layer (nlevgrnd)
   frootc               => cnveg_carbonstate_inst%frootc_patch               , & ! Input:   (:) (gC/m2) fine root C
   leafn                => cnveg_nitrogenstate_inst%leafn_patch              , & ! Input:   (:) (gN/m2) leaf N                                    
   leafn_storage        => cnveg_nitrogenstate_inst%leafn_storage_patch      , & ! Input:   (:) (gN/m2) leaf N storage                            
   frootn               => cnveg_nitrogenstate_inst%frootn_patch             , & ! Input:   (:) (gN/m2) fine root N                               
   frootn_storage       => cnveg_nitrogenstate_inst%frootn_storage_patch     , & ! Input:   (:) (gN/m2) fine root N storage                       
   livecrootn           => cnveg_nitrogenstate_inst%livecrootn_patch         , & ! Input:   (:) (gN/m2) live coarse root N                        
   livecrootn_storage   => cnveg_nitrogenstate_inst%livecrootn_storage_patch , & ! Input:   (:) (gN/m2) live coarse root N storage                
   h2osoi_liq           => waterstatebulk_inst%h2osoi_liq_col                , & ! Output:(:,:) (kg/m2) liquid water (new)  
   t_soisno             => temperature_inst%t_soisno_col                     , & ! Input: (:,:) (Kelvin) soil temperature     
   plantCN              => cnveg_state_inst%plantCN_patch                    , & ! Input:   (:)           Plant C:N used by FUN
   availc               => cnveg_carbonflux_inst%availc_patch                , & ! Output:  (:) (gC/m2/s) C flux available for allocation 
   npp_growth           => cnveg_carbonflux_inst%npp_growth_patch            , & ! Output:  (:) (gC/m2/s) Total C u for growth in FUN / MIMICSplus
  
   smin_nh4_vr          => soilbiogeochem_nitrogenstate_inst%smin_nh4_vr_col         , & ! Input:  (:,:) (gN/m3) col soil mineral NH4              
   smin_no3_vr          => soilbiogeochem_nitrogenstate_inst%smin_no3_vr_col         , & ! Input:  (:,:) (gN/m3) col soil mineral NO3              
   smin_no3_to_plant_vr => soilbiogeochem_nitrogenflux_inst%smin_no3_to_plant_vr_col , & ! Input:  (:,:) (gN/m3/s) col vertically-resolved plant uptake of soil NO3 
   smin_nh4_to_plant_vr => soilbiogeochem_nitrogenflux_inst%smin_nh4_to_plant_vr_col , & ! Input:  (:,:) (gN/m3/s) col vertically-resolved plant uptake of soil NH4 
   
   total_inorgN_uptake_vr  => cnveg_nitrogenflux_inst%sminn_to_symbiont_mimicsplus_vr_patch     , & ! Output: (:,:) (gN/m2/s) Total layer soil N uptake of MIMICSplus 
   total_inorg_no3_uptake  => cnveg_nitrogenflux_inst%sminn_to_symbiont_mimicsplus_no3_vr_patch , & ! Output: (:,:) (gN/m2/s) Total layer soil NO3 uptake of MIMICSplus 
   total_inorg_nh4_uptake  => cnveg_nitrogenflux_inst%sminn_to_symbiont_mimicsplus_nh4_vr_patch , & ! Output: (:,:) (gN/m2/s) Total layer soil NH4 uptake of MIMICSplus

   decomp_cpools_vr     => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col   , &  ! Input: (:,:,:) (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) C pools
   decomp_npools_vr     => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col , &  ! Input: (:,:,:) (gN/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   totsymbc             =>   soilbiogeochem_carbonstate_inst%totsymbc_col         , &  ! Total C of symbiont pools [gC/m2]
   totsymbn             =>   soilbiogeochem_nitrogenstate_inst%totsymbn_col       , &  ! Total N of symbiont pools [gC/m2]
 
   ! Symbiont variables  
   is_active            => symbiont_inst%is_active      , &     ! If symbiont uptake pathway is active for patch   [-]
   perecm               => pftcon%perecm                , &     ! The fraction of ECM-associated PFT               [-]
   symb_eff             => symbiont_inst%symb_eff       , &     ! Symbiont efficiency in nitrogen uptake           [-]
   C_reservoir          => symbiont_inst%C_reservoir    , &     ! Carbon reservoir in intermediate pools       [gC/m2]
   N_reservoir          => symbiont_inst%N_reservoir    , &     ! Nitrogen reservoir in intermediate pools     [gN/m2]
   C_biomass            => symbiont_inst%C_biomass      , &     ! Carbon biomass of symbiont                   [gC/m2]
   N_biomass            => symbiont_inst%N_biomass      , &     ! Nitrogen biomass of symbiont                 [gN/m2]
   symb_growth          => symbiont_inst%symb_growth    , &     ! Symbiotic biomass growth rate                [gC/m2]
   C_mortality          => symbiont_inst%C_mortality    , &     ! Symbiotic turnover per soil layer and column [gC/m3/s]
   N_mortality          => symbiont_inst%N_mortality    , &     ! Symbiotic turnover per soil layer and column [gN/m3/s]
   somc_nuptake_col     => symbiont_inst%somc_nuptake_col , &   ! Nitrogen uptake from SOMc via mining         [gN/m3/s]
   somp_nuptake_col     => symbiont_inst%somp_nuptake_col , &   ! Nitrogen uptake from SOMp via mining         [gN/m3/s]
   N_mine_somc2soma_col    => symbiont_inst%N_mine_somc2soma_col, &   ! Leftover part of co-mineralized N, not taken up by miners  [gN/m3/s]
   N_mine_somp2soma_col    => symbiont_inst%N_mine_somp2soma_col, &   ! Leftover part of co-mineralized N, not taken up by miners  [gN/m3/s]
   symbiont_gr_patch    => cnveg_carbonflux_inst%symbiont_gr_patch     , & ! Total C loss of symbionts that is repired (growth)    [gC/m2/s] 
   symbiont_maint_patch => cnveg_carbonflux_inst%symbiont_maint_patch  , & ! Total C loss of symbionts repired (maintainence)      [gC/m2/s] 

   plant_ndemand             => cnveg_nitrogenflux_inst%plant_ndemand_patch            , & ! Output:[real(r8) (:)]  N flux required to support initial GPP (gN/m2/s)
   sminn_to_plant_mimicsplus => cnveg_nitrogenflux_inst%sminn_to_plant_mimicsplus_patch, & ! Output:[real(r8) (:)]  nitrogen sent to plant from symbionts (gN/m2/s)
   
   nfix_to_sminn          => soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col   , & ! Output:  [real(r8) (:)] symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
   N_fixation             => cnveg_nitrogenflux_inst%Nfix_patch                   , & ! Output:  [real(r8) (:) ]  Symbiotic BNF (gN/m2/s)
       
   somc_cuptake_col     => symbiont_inst%somc_cuptake_col , &   ! Nitrogen uptake from SOMc via mining       [gC/m3/s]
   somp_cuptake_col     => symbiont_inst%somp_cuptake_col , &   ! Nitrogen uptake from SOMp via mining       [gC/m3/s]
   root_exudate_C_col   => symbiont_inst%root_exudate_C_col &   ! Leftover C from allocation to symbionts    [gC/m2/s]
   )

   !-----------------------------------------------------------------------
   
   ! Calculationg a root profile
   ! https://escomp.github.io/ctsm-docs/versions/master/html/tech_note/Plant_Hydraulics/CLM50_Tech_Note_Plant_Hydraulics.html?highlight=root
   availc_alloc(bounds%begp:bounds%endp)                          = 0.0_r8
   C_allocation_to_N_acq(bounds%begp:bounds%endp)                 = 0.0_r8
   myc_biomass_layer(bounds%begp:bounds%endp, 1:nlevdecomp)       = 0.0_r8
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
   
   old_C_biomass_scav(bounds%begp:bounds%endp)                    = C_biomass(p,i_scav)
   old_C_biomass_mine(bounds%begp:bounds%endp)                    = C_biomass(p,i_miner)
   old_C_biomass_fix(bounds%begp:bounds%endp)                     = C_biomass(p,i_fixer)
   old_C_reservoir_scav(bounds%begp:bounds%endp)                  = C_reservoir(p,i_scav)
   old_C_reservoir_mine(bounds%begp:bounds%endp)                  = C_reservoir(p,i_miner)
   old_C_reservoir_fix(bounds%begp:bounds%endp)                   = C_reservoir(p,i_fixer)
   
   do fp = 1,num_soilp        
      p = filter_soilp(fp)
      c = patch%column(p)
      g = patch%gridcell(p)
      root_dens_sum = 0.0_r8

      do j = 1, nlevdecomp
         root_dens_sum = root_dens_sum + crootfr(p,j) * frootc(p)
         if (isnan(decomp_cpools_vr(c,j,i_avl_som))) then
            write(iulog,*)'ECW, AVL IS NAN ', grc%latdeg(g),grc%londeg(g)
         endif
      end do

      do j = 1, nlevdecomp
         if (root_dens_sum > 0.0_r8) then
            root_dens_frac(p,j) = (crootfr(p,j) * frootc(p)) / root_dens_sum
         else
            root_dens_frac(p,j) = 0.0_r8
         end if         
      
         smin_nh4_avail(p,j) = smin_nh4_to_plant_vr(c,j)  
         smin_no3_avail(p,j) = smin_no3_to_plant_vr(c,j) 

      end do
   end do

   !-----------------------------------------------------------------------

   C_allocation_to_N_acq(bounds%begp:bounds%endp)                 = 0.0_r8
   maint_resp                                                     = 0.0_r8
   symb_CO2_prod(bounds%begp:bounds%endp)                         = 0.0_r8
   total_symbiont_turnover_C(bounds%begp:bounds%endp, 1:n_symb)   = 0.0_r8
   symbiont_turnover_C(bounds%begp:bounds%endp, 1:n_symb)         = 0.0_r8
   symbiont_turnover_N(bounds%begp:bounds%endp, 1:n_symb)         = 0.0_r8

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)

      ! Carbon provided, to be used for either growth or Nitrogen uptake
      availc_alloc(p)           = availc(p)        *  0.5_r8

      C_allocation_to_N_acq(p)  = availc_alloc(p)
      npp_growth(p) = availc(p)        *  0.5_r8

   end do

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
   N_fixation(begp:endp)                                    = 0.0_r8
   root_N_uptake(begp:endp)                                 = 0.0_r8
   root_N_to_plant(begp:endp)                               = 0.0_r8
   
   ! Mycorrhizal N mining (ECM-style)
   call myc_miner_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                             bounds, symbiont_inst, temperature_inst, soilstate_inst, waterstatebulk_inst, &
                             soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst, &
                             root_dens_frac(begp:endp,1:nlevdecomp), &
                             somc_nuptake(begp:endp,1:nlevdecomp), somp_nuptake(begp:endp,1:nlevdecomp), &
                             somc_cuptake(begp:endp,1:nlevdecomp), somp_cuptake(begp:endp,1:nlevdecomp), &
                             N_mine_somc2soma(begp:endp,1:nlevdecomp), N_mine_somp2soma(begp:endp,1:nlevdecomp))

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
      do j = 1,nlevdecomp
         N_reservoir(p,i_miner) =  N_reservoir(p,i_miner) + somc_nuptake(p,j) + somp_nuptake(p,j)
      enddo
   enddo
                           
   ! Scavenging (AM-style)
   call myc_scavenger_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                                 bounds, symbiont_inst, root_dens_frac(begp:endp,1:nlevdecomp), &
                                 smin_no3_avail(begp:endp,1:nlevdecomp), smin_nh4_avail(begp:endp,1:nlevdecomp), &
                                 no3_scav_up(begp:endp,1:nlevdecomp), nh4_scav_up(begp:endp,1:nlevdecomp)) 
   
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
      do j = 1,nlevdecomp
         N_reservoir(p,i_scav) =  N_reservoir(p,i_scav) + no3_scav_up(p,j) + nh4_scav_up(p,j)
      enddo
   enddo

 
   ! Active root uptake 
   call active_root_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds, &
                                soilstate_inst, root_dens_frac(begp:endp,1:nlevdecomp), frootc(begp:endp), &
                                smin_no3_avail(begp:endp,1:nlevdecomp), smin_nh4_avail(begp:endp,1:nlevdecomp), &
                                no3_active_up(begp:endp,1:nlevdecomp), nh4_active_up(begp:endp,1:nlevdecomp))

  
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)

      nh4_passiv_up(bounds%begp:bounds%endp, 1:nlevdecomp) = 0.0_r8
      no3_passiv_up(bounds%begp:bounds%endp, 1:nlevdecomp) = 0.0_r8
      ! Symbiotic N2 Fixation
      
      ! Passive nitrogen uptake from the rhizosphere by roots with soil water uptake:
      
      ! smin_no3_vr is per column and has to be per patch
      ! if water in layer:
      do j = 1,nlevdecomp
         t_soi_degC = t_soisno(c,j) - tfrz     ! Soil temperature in degrees Celcius
         if (t_soi_degC > 0.01_r8 .and. h2osoi_liq(c,j) > 0.01_r8) then
            no3_passiv_up(p,j) = waterfluxbulk_inst%qflx_tran_veg_patch(p) * (smin_no3_avail(p,j) / h2osoi_liq(c,j)) !per patch?
            nh4_passiv_up(p,j) = waterfluxbulk_inst%qflx_tran_veg_patch(p) * (smin_nh4_avail(p,j) / h2osoi_liq(c,j))
         else
            nh4_passiv_up(p,j) = 0.0_r8
            no3_passiv_up(p,j) = 0.0_r8
         end if
         ! NO3 and NH4 uptake depends on how much N is available in soil
         nh4_passiv_up(p,j) = min(nh4_passiv_up(p,j), smin_nh4_avail(p,j)) !ECW not sure, check
         no3_passiv_up(p,j) = min(no3_passiv_up(p,j), smin_no3_avail(p,j))
      enddo 
   enddo

   !----------------------------------------------------------------------
 
   ! Fluxes: Inorganic nitrogen pool -> Symbionts 
   ! Nitrogen uptake by symbionts is limited to not deplete inorganic N pool
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
      
      !This might be unnecessary, it's the C that could be respired while mining, but I think it should go to SOMa pool
      !miner_n_patch(p) = 0.0_r8
      !do j = 1, nlevdecomp
      !   miner_n_patch(p) = miner_n_patch(p) + somc_cuptake(p,j) + somp_cuptake(p,j)
      !end do
      root_N_uptake(p) = 0.0_r8
      root_N_to_plant(p) = 0.0_r8
      do j = 1, nlevdecomp
         total_inorg_no3_uptake(p,j) = (no3_passiv_up(p,j) + no3_active_up(p,j) + no3_scav_up(p,j)) * dt
         total_inorg_nh4_uptake(p,j) = (nh4_passiv_up(p,j) + nh4_active_up(p,j) + nh4_scav_up(p,j)) * dt
         if (total_inorg_nh4_uptake(p,j) > 0.0_r8) then
            write(iulog,*)' '
         endif
         if(smin_nh4_avail(p,j) > 0.0_r8) then
            write(iulog,*)' '
         endif
         if (total_inorg_no3_uptake(p,j) > 0.0_r8) then
            write(iulog,*)' '
         endif
         if(smin_no3_avail(p,j) > 0.0_r8) then
            write(iulog,*)' '
         endif

         if (smin_no3_avail(p,j) <= 0.0_r8 .and. total_inorg_no3_uptake(p,j) > 0.0_r8) then
            write(iulog,*) 'Warning: NO3 uptake attempted from layer with zero availability.'
            no3_passiv_up(p,j) = 0.0_r8
            no3_active_up(p,j) = 0.0_r8
            no3_scav_up(p,j)   = 0.0_r8
        endif

        if (smin_nh4_avail(p,j) <= 0.0_r8 .and. total_inorg_nh4_uptake(p,j) > 0.0_r8) then
            write(iulog,*) 'Warning: NH4 uptake attempted from layer with zero availability.'
            nh4_passiv_up(p,j) = 0.0_r8
            nh4_active_up(p,j) = 0.0_r8
            nh4_scav_up(p,j)   = 0.0_r8
        endif

        
         ! If nitrogen uptake exceeds avaliable nitrogen, scale each uptake pathway down
         ! Without multipling by 0.9, I scale to N uptake down, but still allow to take up all avaliable N from soil (maybe not so good)
         ! Therefore I multiply with 0.9 to leave 10% in soil
         if ( (total_inorg_no3_uptake(p,j) > smin_no3_avail(p,j)) .and. &
              (smin_no3_avail(p,j) > 0.0_r8) ) then
           write(iulog,*)'NO3 uptake by passive / active / scavenger pathway exceeds soil N uptake and was scaled down, leaving 10% N in soil'
           no3_passiv_up(p,j)  = no3_passiv_up(p,j)  * (smin_no3_avail(p,j) / total_inorg_no3_uptake(p,j)) * 0.9_r8
           no3_active_up(p,j)  = no3_active_up(p,j)  * (smin_no3_avail(p,j) / total_inorg_no3_uptake(p,j)) * 0.9_r8
           no3_scav_up(p,j)    = no3_scav_up(p,j)    * (smin_no3_avail(p,j) / total_inorg_no3_uptake(p,j)) * 0.9_r8
         endif

         if ( (total_inorg_nh4_uptake(p,j) > smin_nh4_avail(p,j)) .and. &
              (smin_nh4_avail(p,j) > 0.0_r8) ) then
            write(iulog,*)'NH4 uptake by passive / active / scavenger pathway exceeds soil N uptake and was scaled down, leaving 10% N in soil'
            nh4_passiv_up(p,j)  = nh4_passiv_up(p,j)  * (smin_nh4_avail(p,j) / total_inorg_nh4_uptake(p,j)) * 0.9_r8
            nh4_active_up(p,j)  = nh4_active_up(p,j)  * (smin_nh4_avail(p,j) / total_inorg_nh4_uptake(p,j)) * 0.9_r8
            nh4_scav_up(p,j)    = nh4_scav_up(p,j)    * (smin_nh4_avail(p,j) / total_inorg_nh4_uptake(p,j)) * 0.9_r8
         endif

         root_N_uptake(p) = root_N_uptake(p) + no3_active_up(p,j) + nh4_active_up(p,j)

         root_N_to_plant(p) = root_N_to_plant(p) + no3_active_up(p,j) + nh4_active_up(p,j) +  no3_passiv_up(p,j) + nh4_passiv_up(p,j)

         total_inorg_no3_uptake(p,j) = (no3_passiv_up(p,j) + no3_active_up(p,j) + no3_scav_up(p,j)) * dt !ECW is this correct updated?
         total_inorg_nh4_uptake(p,j) = (nh4_passiv_up(p,j) + nh4_active_up(p,j) + nh4_scav_up(p,j)) * dt

         total_inorgN_uptake_vr(p,j) = total_inorg_nh4_uptake(p,j) + total_inorg_no3_uptake(p,j)

      end do
   end do

   ! smin_avail is the inorganic N avaliable for plant uptake and should equal the sum of all symbiont uptakes
   do j = 1, nlevdecomp
      do p = bounds%begp,bounds%endp
       smin_no3_avail(p,j) = no3_passiv_up(p,j) + no3_active_up(p,j) + no3_scav_up(p,j)
       smin_nh4_avail(p,j) = nh4_passiv_up(p,j) + nh4_active_up(p,j) + nh4_scav_up(p,j) 
      end do
      ! make patch to column
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, smin_no3_avail(bounds%begp:bounds%endp,j), smin_no3_avail_col(bounds%begc:bounds%endc,j))
      call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, smin_nh4_avail(bounds%begp:bounds%endp,j), smin_nh4_avail_col(bounds%begc:bounds%endc,j))
   end do 

 
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)

      ! GROWTH AND TURNOVER
      
       ! Mycorrhizal scavengers
       symb_growth(p,i_scav) = sulman_max_symb_growth * C_reservoir(p,i_scav) / (C_reservoir(p,i_scav) + sulman_kgrowth) * sulman_growth_scav * dt
       maint_resp = min(C_biomass(p,i_scav)*sulman_tau_scav * (1.0_r8 - sulman_tau_sym) * dt, symb_growth(p,i_scav))
       ! Nitrogen limitation
       if (symb_growth(p,i_scav)  > sulman_cn_scav * N_reservoir(p,i_scav) * 0.9_r8 + maint_resp)  then
          ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
          symb_growth(p,i_scav) = sulman_cn_scav * N_reservoir(p,i_scav) * 0.9_r8 + maint_resp
       end if
    
       ! C loss during growth from C reservoir to C biomass pool, due to CUE efficiency
       symb_CO2_prod(p) = symb_CO2_prod(p) + symb_growth(p,i_scav) / sulman_growth_scav * (1.0 - sulman_growth_scav)
    
       ! Fraction of N from maintainace respiration stays in N reservoir while C is respiered
       N_reservoir(p,i_scav) = N_reservoir(p,i_scav) + N_biomass(p,i_scav) * ((1 - sulman_tau_sym) * sulman_tau_scav *dt)
       
       ! Total symbiont Turnover (including necromass and maintanance respiration)
       total_symbiont_turnover_C(p,i_scav)  = C_biomass(p,i_scav)  * params_inst%sulman_tau_scav
      
       ! C biomass plus growth flux from reservoir minus the turnover (including maint. respiration)
       C_biomass(p,i_scav) = (C_biomass(p,i_scav) + symb_growth(p,i_scav)) - (C_biomass(p,i_scav) * (sulman_tau_scav * dt))
       ! N_biomass(p,i_scav) = N_biomass(p,i_scav) + (symb_growth(p,i_scav) - maint_resp) / sulman_cn_scav - N_biomass(p,i_scav) / sulman_tau_scav * sulman_tau_sym * dt
       N_biomass(p,i_scav) = C_biomass(p,i_scav) / sulman_cn_scav
       C_reservoir(p,i_scav) = C_reservoir(p,i_scav) - symb_growth(p,i_scav) / sulman_growth_scav
       ! N poool doesn't have growth respiration so I don't need to account for it
       N_reservoir(p,i_scav) = N_reservoir(p,i_scav) - (symb_growth(p,i_scav)) / sulman_cn_scav
      
      
       ! Mycorrhizal miners
       symb_growth(p,i_miner) = sulman_max_symb_growth * C_reservoir(p,i_miner) / (C_reservoir(p,i_miner) + sulman_kgrowth) * sulman_growth_mine * dt
       maint_resp = min(C_biomass(p,i_miner) * (sulman_tau_mine * (1.0_r8 - sulman_tau_sym)) * dt, symb_growth(p,i_miner))
       ! Nitrogen limitation
       if (symb_growth(p,i_miner) > sulman_cn_mine * N_reservoir(p,i_miner) * 0.9_r8 + maint_resp) then
          symb_growth(p,i_miner) = sulman_cn_mine * N_reservoir(p,i_miner) * 0.9_r8 + maint_resp
       end if
    
       symb_CO2_prod(p) = symb_CO2_prod(p) + symb_growth(p,i_miner) / sulman_growth_mine * (1.0 - sulman_growth_mine)
    
       N_reservoir(p,i_miner) = N_reservoir(p,i_miner) + N_biomass(p,i_miner) * ((1.0_r8 - sulman_tau_sym) * sulman_tau_mine * dt)
      
       total_symbiont_turnover_C(p,i_miner) = C_biomass(p,i_miner) * params_inst%sulman_tau_mine

       C_biomass(p,i_miner) = (C_biomass(p,i_miner) + symb_growth(p,i_miner)) - (C_biomass(p,i_miner) * (sulman_tau_mine * dt))
       ! N_biomass(p,i_miner) = N_biomass(p,i_miner) + (symb_growth(p,i_miner) - maint_resp) / sulman_cn_mine - N_biomass(p,i_miner) / sulman_tau_mine * sulman_tau_sym * dt
       N_biomass(p,i_miner) = C_biomass(p,i_miner) / sulman_cn_mine
       C_reservoir(p,i_miner) = (C_reservoir(p,i_miner) - symb_growth(p,i_miner) / sulman_growth_mine)
       N_reservoir(p,i_miner) = N_reservoir(p,i_miner) - (symb_growth(p,i_miner)) / sulman_cn_mine
      

       ! Nitrogen Fixation
       symb_growth(p,i_fixer) = sulman_max_symb_growth * C_reservoir(p,i_fixer) / (C_reservoir(p,i_fixer) + sulman_kgrowth) * sulman_growth_fix * dt
       maint_resp = min(C_biomass(p,i_fixer) * sulman_tau_fix * (1.0_r8 - sulman_tau_sym) * dt, symb_growth(p,i_fixer))
       ! if (symb_growth(p,i_fixer) > sulman_cn_fix * N_reservoir(p,i_fixer) * 0.9) then
          ! Not enough nitrogen to support growth. Limit to available N, and leave a little bit left over for plant
          ! symb_growth(p,i_fixer) = sulman_cn_fix * N_reservoir(p,i_fixer) * 0.9
       ! end if
    
       symb_CO2_prod(p) = symb_CO2_prod(p) + symb_growth(p,i_fixer) / sulman_growth_fix * (1.0 - sulman_growth_fix)
       
       ! Fixation has to be done at the biomas update, since it is reduced by the growth
       ! will use later N_fixation=N_fixation*exp(-0.5*0.27*25.15 + 0.27*(soilT-273.15)*(1.0-0.5*(soilT-273.15)/25.15))

       ! Amount of nitrogen fixed by fixer biomass
       N_fixation(p) = C_biomass(p,i_fixer) * sulman_rfix * dt

       ! N reservoir grows by the amount of N that was fixed and by the N that is left after maintainance respiration
       N_reservoir(p,i_fixer) = N_reservoir(p,i_fixer) + N_fixation(p) + N_biomass(p,i_fixer) * ((1 - sulman_tau_sym) * sulman_tau_fix *dt)

       ! Substract the maintainence respiration and turnover from N fixation
       !N_fixation(p) = N_fixation(p) - N_biomass(p,i_fixer) * sulman_tau_fix * dt * sulman_tau_sym ! (1 - sulman_tau_sym)


       ! MOVE TOTAL STUFF HERE
       total_symbiont_turnover_C(p,i_fixer) = C_biomass(p,i_fixer) * params_inst%sulman_tau_fix

       C_biomass(p,i_fixer) = (C_biomass(p,i_fixer) + symb_growth(p,i_fixer)) - (C_biomass(p,i_fixer) * (sulman_tau_fix * dt))
       ! C_biomass(p,i_fixer) = max(C_biomass(p,i_fixer), 0.0)
      
       ! N Biomass from last timestep
       N_biomass_old(p) = N_biomass(p,i_fixer)
       !N_biomass(p,i_fixer) = N_biomass(p,i_fixer) + (symb_growth(p,i_fixer) - maint_resp) / sulman_cn_fix - N_biomass(p,i_fixer) / sulman_tau_fix * sulman_tau_sym * dt

       N_biomass(p,i_fixer) = C_biomass(p,i_fixer) / sulman_cn_fix
       C_reservoir(p,i_fixer) = C_reservoir(p,i_fixer) - (symb_growth(p,i_fixer)) / sulman_growth_fix
       
       ! N fixers just make all the N they need for their biomass
       ! Don't understant, but this should be the variable going into nfix? 
       N_fixation(p) = N_fixation(p) + N_biomass(p,i_fixer) - N_biomass_old(p) 


      !----------------------------------------------------------------------------------------------------------------------------

      ! Respiration during symbiont growth 
      symbiont_gr_patch(p) = symb_CO2_prod(p) / dt 

      ! Symbiont Turnover
      ! gC/m2/s
     
      

      symbiont_turnover_C(p,i_miner) = total_symbiont_turnover_C(p,i_miner) * params_inst%sulman_tau_sym
      symbiont_turnover_C(p,i_scav)  = total_symbiont_turnover_C(p,i_scav) * params_inst%sulman_tau_sym 
      symbiont_turnover_C(p,i_fixer) = total_symbiont_turnover_C(p,i_fixer) * params_inst%sulman_tau_sym

      !gN/m3/s
      symbiont_turnover_N(p,i_miner) = (N_biomass(p,i_miner) * params_inst%sulman_tau_mine) * params_inst%sulman_tau_sym
      symbiont_turnover_N(p,i_scav)  = (N_biomass(p,i_scav)  * params_inst%sulman_tau_scav) * params_inst%sulman_tau_sym
      symbiont_turnover_N(p,i_fixer) = (N_biomass(p,i_fixer) * params_inst%sulman_tau_fix)  * params_inst%sulman_tau_sym

      ! Maintainance respiration as fraction of turnover 
      !gC/m2/s               
      symbiont_maint_patch(p) = (1.0 - params_inst%sulman_tau_sym) * &
                       (total_symbiont_turnover_C(p,i_miner) + total_symbiont_turnover_C(p,i_scav) + total_symbiont_turnover_C(p,i_fixer))
   
   end do
      !----------------------------------------------------------------------------------------------------------------------------

   call update_symbionts(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                        bounds, symbiont_inst, root_N_uptake, root_N_to_plant, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, cnveg_carbonflux_inst, &
                        cnveg_state_inst, C_allocation_to_N_acq(bounds%begp:bounds%endp), root_exudate_C(bounds%begp:bounds%endp))

   ! Add nitrogen that was taken through the roots straight to the plant

      ! Maybe I need to make my own variables instead of using these:  
   call p2c(bounds, num_bgc_soilc, filter_bgc_soilc, &
           N_fixation(bounds%begp:bounds%endp), &
           nfix_to_sminn(bounds%begc:bounds%endc))
     
     
                        !Check if reservoirs are not zero
     ! if (C_reservoir(p,i_scav) <= 0._r8 .or. N_reservoir(p,i_scav) <= 0._r8) then 
      !   write(iulog,*), 'C_reservoir_scav =', C_reservoir(p,i_scav), 'N_reservoir_scav =', N_reservoir(p,i_scav)
      !   write(iulog,*), 'no3_scav_up=', no3_scav_up(begp:endp,1:nlevdecomp), 'nh4_scav_up=', nh4_scav_up(begp:endp,1:nlevdecomp)
   
      !call endrun(msg = "ERROR: Scavenger symbiont C or N reservoirs are zero or negative." // &
      !  errMsg(sourcefile, __LINE__))
      !end if

      !if (C_reservoir(p,i_miner) <= 0._r8 .or. N_reservoir(p,i_miner) <= 0._r8) then 
      !  write(iulog,*), 'C_reservoir_miner =', C_reservoir(p,i_miner)
      !  write(iulog,*), 'N_reservoir_miner =', N_reservoir(p,i_miner) 
      !  write(iulog,*), ' somp_nuptake =', somp_nuptake(begp:endp,1:nlevdecomp)
      !  write(iulog,*), ' somc_cuptake =', somc_cuptake(begp:endp,1:nlevdecomp)
      !  write(iulog,*), ' somp_cuptake =', somp_cuptake(begp:endp,1:nlevdecomp)
      !  write(iulog,*), ' N_mine_somc2soma =', N_mine_somc2soma(begp:endp,1:nlevdecomp)
      !  write(iulog,*), ' N_mine_somp2soma =', N_mine_somp2soma(begp:endp,1:nlevdecomp)
        
      !call endrun(msg = "ERROR: Miner symbiont C or N reservoirs are zero or negative." // &
      !   errMsg(sourcefile, __LINE__))
      !end if
     
     !if (C_reservoir(p,i_fixer) <= 0._r8 .or. N_reservoir(p,i_fixer) <= 0._r8) then 
     !   write(iulog,*), 'C_reservoir_fixer =', C_reservoir(p,i_fixer), 'N_reservoir_fixer =', N_reservoir(p,i_fixer)
     ! call endrun(msg = "ERROR: Fixer symbiont C or N reservoirs are zero or negative." // &
     !    errMsg(sourcefile, __LINE__))
     ! end if
      
   !npp_growth(bounds%begp:bounds%endp) = sminn_to_plant_mimicsplus(p) * plantCN + root_exudate_C(p)

   ! Update the symb_turnover_C and _N to make them per layer with root_dens_frac
   do j = 1, nlevdecomp
      do p = bounds%begp,bounds%endp
         symb_turnover_layer_C(p,j) = (symbiont_turnover_C(p,i_miner) + symbiont_turnover_C(p,i_scav)  & 
                                       + symbiont_turnover_C(p,i_fixer)) * root_dens_frac (p,j) / col%dz(c,j)
         symb_turnover_layer_N(p,j) = (symbiont_turnover_N(p,i_miner) + symbiont_turnover_N(p,i_scav) &
                                       + symbiont_turnover_N(p,i_fixer)) * root_dens_frac (p,j) / col%dz(c,j)
         root_exudate_C_layer(p,j) =  (root_exudate_C(p) / col%dz(c,j))  * root_dens_frac(p,j)
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
  
  !CHANGER p to c
!   call symbiont_inst%Summary(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds,'c12',totsymbc(bounds%begc:bounds%endc))
!   call symbiont_inst%Summary(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds,'n',totsymbn(bounds%begc:bounds%endc))


   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
      do j = 1, nlevdecomp
         if (decomp_cpools_vr(c,j,i_phys_som) - somp_cuptake_col(c,j) * dt < 0.0_r8 ) then
         write(iulog,*), ' somp diag ', col%z(c,j), decomp_cpools_vr(c,j,i_phys_som),somp_cuptake_col(c,j)
         endif
         if (decomp_cpools_vr(c,j,i_chem_som) - somc_cuptake_col(c,j) * dt < 0.0_r8 ) then
         write(iulog,*), ' somc diag ', col%z(c,j), decomp_cpools_vr(c,j,i_chem_som),somc_cuptake_col(c,j)
         end if
         if(abs(somc_cuptake(p,j))< 0.0_r8 .or. &
         abs(somp_cuptake(p,j))< 0.0_r8 ) then
         write(iulog,*),'ERROR : som fluxes'
         write(iulog,*),'somc,somp',somc_cuptake(p,j),somp_cuptake(p,j)
         call endrun(msg=errMsg(sourcefile, __LINE__))
         endif

      enddo
   enddo

   end associate
   end subroutine CN_soil_veg_exchange

   !----------------------------------------

   ! FUNCTIONS
   ! Calculates the maximum enzymatic activity rate (Vmax) of mycorrhizal organisms as a function of soil temperature
   ! Determines how temperature regulates enzyme activity, which is a key control on how fast mycorrhizae can process organic matter
   function Vmax_myc(soil_T)
      real(r8), intent(in)   :: soil_T                 ! Soil temperature in Kelvin
      real(r8), parameter    :: Tref=293.15            ! Reference Temperature in Kelvin
      !real(r8), parameter    :: Ea=37000_r8            ! Activation energy (kJ/mol) Sulman et al. (2019)
      real(r8), parameter    :: Ea=54000_r8            ! Activation energy (J/mol) ELIN
      real(r8), parameter    :: R_gas = 8.314472       ! Universal gas constant, J/mol*K
      real(r8)               :: alpha                  ! Scaling factor that normalizes the exponential temperature response to match a specified reference value. [s]
      real(r8)               :: Vmax_myc               ! [s]

      ! exp(-Ea / (R_gas * Tref)) is the Arrhenius term evaluated at the reference temperature
      alpha = params_inst%sulman_vmax_ref_mine / exp(-Ea /(R_gas*Tref))
      Vmax_myc = alpha * exp(-Ea / (R_gas * soil_T))
   end function Vmax_myc

   
   function resp_myc(soil_carbon, myc_biomass, soil_T, wliq, wair)
      ! This is the rate of C removed from soil pool as respiration(not actual respiration, rename)
      ! Respiration is driven by mycorrhizae, and depends on how much mycorrhizal biomass there is and how many enzymes they produce 
      ! also limited by environmental conditions (soil moisture & temperature)

      ! USES
      use decompMod         , only : bounds_type
      !
      real(r8), intent(in) :: soil_carbon                      ! Soil carbon stocks, vertically resolved            [gC/m3]
      !ECW MYC_BIOMASS SHOULD BE M3
      real(r8), intent(in) :: myc_biomass                      ! Mycorrhyzal biomass per soillayer                  [gC/m3]
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


      enzymes = myc_biomass * enzyme_frac      
     

      ! If there is carbon avaliable, calculate mycorrhizal repiration

      ! ECW change when enzymes go into SOM pools
      ! enzymes can be calculated in resp_myc
      ! This needs to be connected to soil carbon in CNCStateUpdate
      ! Doubel check this in Sulman code

      !CHECK UNITS
      !CHECK UNITS C UPTAKE
      if (soil_carbon > 0.0_r8 .and. wliq > 0.0_r8) then 
         resp_myc = Vmax_myc(soil_T) * soil_carbon * enzymes / (soil_carbon * params_inst%sulman_km_mine + enzymes) * theta_func
      else 
         resp_myc = 0.0_r8
      end if 

   end function resp_myc

   !----------------------------------------
  
   ! PLANT UPTAKE STRATEGIES

   ! Nitrogen uptake from the rhizosphere by roots (active transport across root-soil interface)
   ! Mineral nitrogen is taken up from the rhizosphere only

   subroutine active_root_N_uptake(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, bounds, &
                                    soilstate_inst, root_dens_frac, froot_carbon, no3_soil, nh4_soil, no3_uptake, nh4_uptake)
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
      real(r8), intent(in)   :: root_dens_frac(bounds%begp:bounds%endp,1:nlevdecomp)! Fraction of root density   [-] 
      real(r8), intent(in)   :: no3_soil(bounds%begp:bounds%endp,1:nlevdecomp)      ! Avaliable soil mineral NO3 [gN/m2]
      real(r8), intent(in)   :: nh4_soil(bounds%begp:bounds%endp,1:nlevdecomp)      ! Avaliable soil mineral NH4 [gN/m2]
      real(r8), intent(inout):: no3_uptake(bounds%begp:bounds%endp,1:nlevdecomp)    ! NO3 uptake from soil       [gN/m2/s]
      real(r8), intent(inout):: nh4_uptake(bounds%begp:bounds%endp,1:nlevdecomp)    ! NH4 uptake from soil       [gN/m2/s]

      real(r8) :: root_biomass_density                                              ! Root biomass density       [g/m3]
      real(r8) :: root_cross_sec_area                                               ! Root cross sectional area  [m2]
      real(r8) :: root_length_density                                               ! Root length density        [m/m3]
      real(r8) :: rhizosphere_frac                         ! Fraction of rihzosphere    [-] 
                                                                                    ! sulman_r_rhiz              [m] 
      real(r8), parameter :: root_radius = 0.29e-03_r8                              ! Root radius                [m]
      real(r8), parameter :: c_to_b = 2.0_r8                                        !                            [g biomass /g C]

      associate(                                                         &
         sulman_r_rhiz          => params_inst%sulman_r_rhiz           , & ! Radius of the rhizosphere          [m]
         ivt                    => patch%itype                         , & ! Input: (:) patch vegetation type   [-]
         rootfr                 => soilstate_inst%rootfr_patch         , & ! Input: (:,:)                       [-]
         root_radius            => pftcon%root_radius                  , & ! Input: 0.00029                     [m] 
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
         if  (root_dens_frac(p,j) > 0.0_r8) then

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
             no3_uptake(p,j) = rhizosphere_frac * (params_inst%sulman_root_no3 * col%dz(c,j)) * (no3_soil(p,j) / (no3_soil(p,j) + params_inst%sulman_km_no3 * col%dz(c,j)))
            else
             no3_uptake(p,j) = 0.0_r8
            end if 

            if (nh4_soil (p,j) > 0.0_r8) then
             nh4_uptake(p,j) = rhizosphere_frac * (params_inst%sulman_root_nh4 * col%dz(c,j)) * (nh4_soil(p,j) / (nh4_soil(p,j) + params_inst%sulman_km_nh4 * col%dz(c,j)))
            else
             nh4_uptake(p,j) = 0.0_r8
            end if 

            ! NO3 and NH4 uptake depends on how much N is available in soil
            no3_uptake(p,j) = min(no3_uptake(p,j), no3_soil(p,j)) !ECW not sure, check
            nh4_uptake(p,j) = min(nh4_uptake(p,j), nh4_soil(p,j))
  
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
                                       root_dens_frac, no3_soil, nh4_soil, no3_uptake, nh4_uptake)
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

      real(r8),intent(in)     :: root_dens_frac(bounds%begp:bounds%endp,1:nlevdecomp)     ! Fraction of root desity        [-]
      real(r8),intent(inout)  :: no3_soil(bounds%begp:bounds%endp,1:nlevdecomp)           ! avaliable soil mineral NO3 [gN/m2]
      real(r8),intent(inout)  :: nh4_soil(bounds%begp:bounds%endp,1:nlevdecomp)           ! avaliable soil mineral NH4 [gN/m2]
      real(r8),intent(inout)  :: no3_uptake(bounds%begp:bounds%endp,1:nlevdecomp)         ! NO3 uptake from soil     [gN/m2/s]
      real(r8),intent(inout)  :: nh4_uptake(bounds%begp:bounds%endp,1:nlevdecomp)         ! NH4 uptake from soil     [gN/m2/s]

      real(r8) :: myc_biomass_layer(bounds%begp:bounds%endp, 1:nlevdecomp)                ! Mycorrhyzal biomass per soil layer [gC/m2]

      associate(                                                                  &
      sulman_k_scav_Ninorg => params_inst%sulman_k_scav_Ninorg                  , & ! Half-saturation inorganic N concentration for mycorrhizal uptake [gN/m3]
      sulman_k_scav        => params_inst%sulman_k_scav                         , & ! Half-saturation mycorrhizal biomass concentration for scavenging [gC/m3]
      sulman_v_scav        => params_inst%sulman_v_scav                         , & ! Maximum N uptake rate by scavenging mycorrhizae                [gN/m3/s]
      C_biomass            => symbiont_inst%C_biomass                           , & ! Carbon biomass of symbiont                                       [gC/m2]
      N_biomass            => symbiont_inst%N_biomass                           , & ! Nitrogen biomass of symbiont                                     [gN/m2]
      symb_efficiency      => symbiont_inst%symb_eff                              & ! Symbiont efficiency in nitrogen uptake                         [gN/gC/s]
      )

      myc_biomass_layer(bounds%begp:bounds%endp, 1:nlevdecomp) = 0.0_r8

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = patch%column(p)
               
            ! Check if there is mycorrhizal biomass per patch and per layer
            ! Calculating mycorrhizal biomass per soil layer
            ! Used for inorganic N (scavangers and root uptake)
            do j = 1, nlevdecomp
               myc_biomass_layer(p,j) = C_biomass(p,i_scav) * root_dens_frac(p,j)
            end do

            ! It probably wont pass balance checks (Rosie). Instead:
            ! Make sure that mycorrhizal pool never goes to 0
   
            ! Check if there is  mycorrhizal biomass in soil layer
            do j = 1, nlevdecomp
               no3_uptake(p,j) = 0.0_r8
               nh4_uptake(p,j) = 0.0_r8
               if (myc_biomass_layer(p,j) >= 0) then 
            
                  ! If there is mycorrhizal biomass in the soil layer, calculate N uptake
                  no3_uptake(p,j) = (sulman_v_scav * col%dz(c,j)) * no3_soil(p,j) / (no3_soil(p,j) + sulman_k_scav_Ninorg * col%dz(c,j)) * &
                      myc_biomass_layer(p,j) / (myc_biomass_layer(p,j) + (sulman_k_scav * col%dz(c,j)))

                  nh4_uptake(p,j) = (sulman_v_scav * col%dz(c,j)) * nh4_soil(p,j) / (nh4_soil(p,j) + sulman_k_scav_Ninorg * col%dz(c,j)) * &
                      myc_biomass_layer(p,j) / (myc_biomass_layer(p,j) + (sulman_k_scav * col%dz(c,j)))

                  
                  ! NO3 and NH4 uptake depends on how much N is available in soil
                  !no3_uptake(p,j) = min(no3_uptake(p,j), no3_soil(p,j))
                  !nh4_uptake(p,j) = min(nh4_uptake(p,j), nh4_soil(p,j))
            
                  ! Compute efficiency, ensuring no division by zero
                  if (myc_biomass_layer(p,j) > 0.0_r8) then
                    ! symb_efficiency(p,i_scav) = symb_efficiency(p,i_scav) + (no3_uptake(p,j) + nh4_uptake(p,j)) / C_biomass(p,i_scav) !ECW 
                     symb_efficiency(p,i_scav) = (no3_uptake(p,j) + nh4_uptake(p,j)) / C_biomass(p,i_scav)
                  endif
                  no3_soil(p,j) = max(0.0_r8,no3_soil(p,j) - no3_uptake(p,j)) !ECW
                  nh4_soil(p,j) = max(0.0_r8,nh4_soil(p,j) - nh4_uptake(p,j))
               end if
            end do
         enddo
      end associate

   end subroutine myc_scavenger_N_uptake

   !----------------------------------------

   subroutine myc_miner_N_uptake (filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                                 bounds, symbiont_inst, temperature_inst, soilstate_inst, waterstatebulk_inst, &
                                 soilbiogeochem_carbonstate_inst, soilbiogeochem_nitrogenstate_inst, &
                                 root_dens_frac, somc_nuptake,somp_nuptake, somc_cuptake, somp_cuptake, N_mine_somc2soma, N_mine_somp2soma)

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

      real(r8), intent(in)    :: root_dens_frac(bounds%begp:bounds%endp, 1:nlevdecomp)  ! Fraction of root desity [-]
      real(r8), intent(inout) :: somc_nuptake(bounds%begp:bounds%endp,1:nlevdecomp)     ! Nitrogen removed from SOMc pool into int miner pool [gN/m3/s]
      real(r8), intent(inout) :: somp_nuptake(bounds%begp:bounds%endp,1:nlevdecomp)     ! Nitrogen removed from SOMp pool into int miner pool [gN/m3/s]
      real(r8), intent(inout) :: somc_cuptake(bounds%begp:bounds%endp,1:nlevdecomp)     ! Carbon removed from SOMc pool into int miner pool [gC/m3/s]
      real(r8), intent(inout) :: somp_cuptake(bounds%begp:bounds%endp,1:nlevdecomp)     ! Carbon removed from SOMp pool into int miner pool [gC/m3/s]
      real(r8), intent(inout) :: N_mine_somc2soma(bounds%begp:bounds%endp,1:nlevdecomp) ! Leftover N due to NUE after mining from SOMc to SOMa [gN/m3/s]
      real(r8), intent(inout) :: N_mine_somp2soma(bounds%begp:bounds%endp,1:nlevdecomp) ! Leftover N due to NUE after mining from SOMp to SOMa [gN/m3/s]

      real(r8) :: wliq                                                     ! Fraction of liquid water-filled pore space (0.0 - 1.0)
      real(r8) :: wice                                                     ! Fraction of frozen water-filled pore space (0.0 - 1.0)
      real(r8) :: wair                                                     ! Fraction of air-filled pore space (0.0 - 1.0)
      real(r8) :: myc_biomass(bounds%begp:bounds%endp, 1:nlevdecomp)       ! Mycorrhyzal biomass in soil  [gC/m3]
      real(r8) :: total_org_nuptake                                        ! Total N uptake from all soil layers []
      
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
         symb_efficiency   => symbiont_inst%symb_eff                 &     ! Symbiont efficiency N uptake per unit of mycorrhizal C biomass  [gN/gC/s]
         )
    
      begp = bounds%begp; endp= bounds%endp

      myc_biomass(bounds%begp:bounds%endp, 1:nlevdecomp) = 0.0_r8

         do fp = 1,num_soilp
            p = filter_soilp(fp)
            c = patch%column(p)
            
            ! Calculating water, ice and air content in soil (equivalent to air_filled porosity, theta, theta sat in Sulman)
      
            ! Make sure fluxes are zero before the loop
            total_org_nuptake = 0.0_r8
            do j = 1, nlevdecomp
               ! Miners per layer
               myc_biomass(p,j) = C_biomass(p,i_miner) * root_dens_frac(p,j) / col%dz(c,j)


               ! this is necessary for miners when they somehow all die.
               if (C_biomass(p,i_miner) == 0.0_r8) then 
                  myc_biomass(p,j) = 0.01_r8 * root_dens_frac(p,j) / col%dz(c,j)
               endif 
               wliq = h2osoi_liq(c,j) / col%dz(c,j) * denh2o
               wice = h2osoi_ice(c,j) / col%dz(c,j) * denice
               wliq = min(1.0_r8, wliq/watsat(c,j))            ! fraction of liquid water-filled pore space (0.0 - 1.0)
               wice = min(1.0_r8, wice/watsat(c,j))            ! fraction of frozen water-filled pore space (0.0 - 1.0)
               wair = max(0.0_r8, 1.0_r8 - wliq- wice)         ! fraction of air-filled pore space (0.0 - 1.0)
         
         
               ! SOMc to miner
               somc_nuptake(p,j) = miner_nuptake(decomp_cpools_vr(c,j,i_chem_som), decomp_npools_vr(c,j,i_chem_som), &            ! N uptake from miners
                                    myc_biomass(p,j), t_soisno(c,j), wliq, wair)
              
               somc_cuptake(p,j) = resp_myc(decomp_cpools_vr(c,j,i_chem_som), myc_biomass(p,j), t_soisno(c,j), wliq, wair)  ! Respired soil C
               
               ! SOMp to miner
               somp_nuptake(p,j) = miner_nuptake(decomp_cpools_vr(c,j,i_phys_som), decomp_npools_vr(c,j,i_phys_som), &            ! N uptake from miners
                                     myc_biomass(p,j), t_soisno(c,j), wliq, wair)
               
               somp_cuptake(p,j) = resp_myc(decomp_cpools_vr(c,j,i_phys_som), myc_biomass(p,j), t_soisno(c,j), wliq, wair)  ! Respired soil C
               
               total_org_nuptake = total_org_nuptake + somc_nuptake(p,j) + somp_nuptake(p,j)
               
               
               ! Leftover N in soil after applying  NUE
               N_mine_somc2soma(p,j) = leftover_n_mining(decomp_cpools_vr(c,j,i_chem_som), decomp_npools_vr(c,j,i_chem_som), &   ! N left in soil during mining needs to be added to Npool
                                       myc_biomass(p,j), t_soisno(c,j), wliq, wair)

               N_mine_somp2soma(p,j) = leftover_n_mining(decomp_cpools_vr(c,j,i_phys_som), decomp_npools_vr(c,j,i_phys_som), &    ! N left in soil during mining
                                        myc_biomass(p,j), t_soisno(c,j), wliq, wair)
   

            end do
      
            if (C_biomass(p, i_miner) > 0._r8) then
               symb_efficiency(p, i_miner) = total_org_nuptake / C_biomass(p, i_miner)
            else
               symb_efficiency(p, i_miner) = 0._r8
            end if
            
         end do 
      end associate

   end subroutine myc_miner_N_uptake

   !------------------------------------------------------------------------------------------------
   ! MINER FUNCTIONS

   function potential_mined_n(soil_carbon, soil_nitrogen, myc_biomass, soil_T, soil_water, soil_air)
     ! DESCRIPTION
     ! This helper returns the co-mineralized N before NUE is applied

     ! ARGUMENTS
      real(r8), intent(in) :: soil_carbon        ! Soil carbon stocks, vertically resolved   [gC/m3]
      real(r8), intent(in) :: soil_nitrogen      ! Soil nitrogen stocks, vertically resolved [gN/m3]
      real(r8), intent(in) :: myc_biomass        ! Mycorrhyzal biomass                       [gC/m2]
      real(r8), intent(in) :: soil_T             ! Soil temperature                              [K]
      real(r8), intent(in) :: soil_water         ! Fraction of liquid water-filled pore space    [-]
      real(r8), intent(in) :: soil_air           ! Fraction of air-filled pore space             [-]
      
      ! LOCAL VARIABLES:
      real(r8) :: potential_mined_n              ! potential mined N [gN/m3/s]
      real(r8) :: potential_tempResp             ! decomposed carbon during mining process [gC/m3/s]
    
      potential_mined_n = 0.0_r8

      ! Call resp_myc(...) to get the soil carbon respiration by mycorrhiza
      potential_tempResp=resp_myc(soil_carbon, myc_biomass, soil_T, soil_water, soil_air)

      ! Don't exceed avaliable C
      if(dt*potential_tempResp > soil_carbon) then
         potential_tempResp = soil_carbon / dt
      end if

      if(soil_carbon > 0) then
         ! Use C:N ratio to estimate how much nitrogen is co-mineralized.
         potential_mined_n = potential_tempResp * (soil_nitrogen / soil_carbon)
      else 
         potential_mined_n=0.0
      end if
   
   end function potential_mined_n
   
  
   function miner_nuptake(soil_carbon, soil_nitrogen, myc_biomass, soil_T, soil_water, soil_air)
   ! This calculates the amount of nitrogen taken up from SOM pools by mycorrhizal mining.
   ! Converts carbon processing (decomposition) into nitrogen uptake, which is the key biogeochemical role of mining.
     
      ! ARGUMENTS
      real(r8), intent(in) :: soil_carbon        ! Soil carbon stocks, vertically resolved   [gC/m3]
      real(r8), intent(in) :: soil_nitrogen      ! Soil nitrogen stocks, vertically resolved [gN/m3]
      real(r8), intent(in) :: myc_biomass        ! Mycorrhyzal biomass                       [gC/m2]
      real(r8), intent(in) :: soil_T             ! Soil temperature                              [K]
      real(r8), intent(in) :: soil_water         ! Fraction of liquid water-filled pore space    [-]
      real(r8), intent(in) :: soil_air           ! Fraction of air-filled pore space             [-]
      real(r8) :: miner_nuptake                  ! Mycorrhizal nitrogen uptake from soil organic matter pools via mining process [gN/m3/s]
      ! ! LOCAL VARIABLES:
      real(r8) :: pot_tempN_decomposed           ! 

      miner_nuptake    = 0.0_r8

      !Call potential_mined_n to get the potential mined N before NUE
      pot_tempN_decomposed = potential_mined_n(soil_carbon, soil_nitrogen, myc_biomass, soil_T, soil_water, soil_air)

     ! Apply nitrogen use efficiency (sulman_nue_mine) to scale the actual uptake.
      miner_nuptake  = pot_tempN_decomposed*params_inst%sulman_nue_mine

   end function miner_nuptake


   function leftover_n_mining(soil_carbon, soil_nitrogen, myc_biomass, soil_T, soil_water, soil_air)

       ! ARGUMENTS
      real(r8), intent(in) :: soil_carbon        ! Soil carbon stocks, vertically resolved   [gC/m3]
      real(r8), intent(in) :: soil_nitrogen      ! Soil nitrogen stocks, vertically resolved [gN/m3]
      real(r8), intent(in) :: myc_biomass        ! Mycorrhyzal biomass                       [gC/m2]
      real(r8), intent(in) :: soil_T             ! Soil temperature                              [K]
      real(r8), intent(in) :: soil_water         ! Fraction of liquid water-filled pore space    [-]
      real(r8), intent(in) :: soil_air           ! Fraction of air-filled pore space             [-]
      real(r8) :: leftover_n_mining              ! Mycorrhizal nitrogen uptake from soil organic matter pools via mining process [gN/m3/s]

      ! ! LOCAL VARIABLES:
      real(r8) :: miner_nuptake                  ! Mycorrhizal nitrogen uptake from soil organic matter pools via mining process [gN/m3/s]
      real(r8) :: pot_tempN_decomposed           ! 
   
      leftover_n_mining = 0.0_r8

       !Call potential_mined_n to get the potential mined N before NUE
      pot_tempN_decomposed = potential_mined_n(soil_carbon, soil_nitrogen, myc_biomass, soil_T, soil_water, soil_air)

     ! Apply nitrogen use efficiency (sulman_nue_mine) to scale the actual uptake.
      miner_nuptake  = pot_tempN_decomposed*params_inst%sulman_nue_mine

      ! Remaining N that is not taken up
      leftover_n_mining = pot_tempN_decomposed - miner_nuptake

   end function leftover_n_mining
   
  !------------------------------------------------------------------------------------------------
   
  subroutine update_symbionts(filter_soilp, filter_bgc_soilc, num_soilp, num_bgc_soilc, &
                              bounds, symbiont_inst, root_N_uptake, root_N_to_plant, cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, &
                              cnveg_carbonflux_inst, cnveg_state_inst, C_allocation_to_N_acq, root_exudate_C)

   ! ! USES:
   use CNVegNitrogenStateType          , only: cnveg_nitrogenstate_type
   use CNVegCarbonFluxType             , only: cnveg_carbonflux_type
   use CNVegStateType                  , only : cnveg_state_type

   ! ARGUMENTS
   integer                        , intent(in)     :: filter_soilp(:)     ! filter for soil patches
   integer                        , intent(in)     :: filter_bgc_soilc(:) ! filter for soil columns   
   integer                        , intent(in)     :: num_soilp           ! number of soil patches in filter
   integer                        , intent(in)     :: num_bgc_soilc       ! number of soil columns in filter
   type(bounds_type)              , intent(in)     :: bounds              
   type(symbiont_type)            , intent(inout)  :: symbiont_inst       
   type(cnveg_nitrogenstate_type) , intent(in)     :: cnveg_nitrogenstate_inst
   type(cnveg_nitrogenflux_type)  , intent(in)     :: cnveg_nitrogenflux_inst
   type(cnveg_carbonflux_type)    , intent(in)     :: cnveg_carbonflux_inst
   type(cnveg_state_type)         , intent(in)     :: cnveg_state_inst

   integer :: begp, endp, begc, endc
  
   real(r8), intent(in)    :: C_allocation_to_N_acq(bounds%begp:bounds%endp) ! Carbon allocated to nitrogen acquisition [gC/m2/s)]
   real(r8), intent(in)    :: root_N_uptake(bounds%begp:bounds%endp)         !
   real(r8), intent(in)    :: root_N_to_plant(bounds%begp:bounds%endp)         !
   real(r8), intent(inout) :: root_exudate_C(bounds%begp:bounds%endp)        ! Leftover C from allocation to symbionts    [gC/m2/s]

   
   !
   ! ! LOCAL VARIABLES:
   integer :: p, fp, c, fc, j, k, l, s  ! indices
   real(r8):: days_per_year
   real(r8) :: total_symbiont_turnover_C(bounds%begp:bounds%endp, 1:n_symb)       ! Total symbiont turnover                  [gC/m2]

   ! Nitrogen uptake by plant from intermediate pools 
   ! Rewrite as N_to_plant(bounds%begp:bounds%endp,1:n_symb) N_to_plant_scav(p) = N_to_plant(p,i_scav)
   real(r8) :: N_to_plant_scav(bounds%begp:bounds%endp)        ! Plant nitrogen uptake from scavengers [gN/m2/s]
   real(r8) :: N_to_plant_mine(bounds%begp:bounds%endp)        ! Plant nitrogen uptake from miners     [gN/m2/s]
   real(r8) :: N_to_plant_fix(bounds%begp:bounds%endp)         ! Plant nitrogen uptake from fixers     [gN/m2/s]
   real(r8) :: scale_N_to_plant(bounds%begp:bounds%endp)       ! Scale factor to scale N uptake to plant if it is bigger that the uptake capazitiy of plant

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
   real(r8) :: fix_C_alloc(bounds%begp:bounds%endp)              ! Carbon allocation from plant to fixer pool     [gC/m2]
   real(r8) :: scav_C_alloc(bounds%begp:bounds%endp)             ! Carbon allocation from plant to scavenger pool [gC/m2]
   real(r8) :: mine_C_alloc(bounds%begp:bounds%endp)             ! Carbon allocation from plant to miner pool     [gC/m2]
   real(r8) :: fix_alloc_accum(bounds%begp:bounds%endp)          ! Accumulated carbon allocation from plant to fixer pool     [gC/m2]
   real(r8) :: mine_alloc_accum(bounds%begp:bounds%endp)         ! Accumulated carbon allocation from plant to scavenger pool [gC/m2]
   real(r8) :: scav_alloc_accum(bounds%begp:bounds%endp)         ! Accumulated carbon allocation from plant to miner pool     [gC/m2]
  
   real(r8),parameter :: root_exudate_frac = 0.05                ! Fraction of NPP that goes into root exudates, same as sulman_fn_alloc
  
   associate(                                                                    &
   leafn_storage          => cnveg_nitrogenstate_inst%leafn_storage_patch      , & ! Input: [real(r8) (:)] (gN/m2) leaf N storage                            
   frootn                 => cnveg_nitrogenstate_inst%frootn_patch             , & ! Input: [real(r8) (:)] (gN/m2) fine root N                               
   frootn_storage         => cnveg_nitrogenstate_inst%frootn_storage_patch     , & ! Input: [real(r8) (:)] (gN/m2) fine root N storage                       
   livecrootn_storage     => cnveg_nitrogenstate_inst%livecrootn_storage_patch , & ! Input:   (:) (gN/m2) live coarse root N storage                
   plantCN                => cnveg_state_inst%plantCN_patch                    , & ! Input:  [real(r8)  (:)]  Plant C:N used by FUN
   availc                 => cnveg_carbonflux_inst%availc_patch                , & ! Output:  (:) (gC/m2/s) C flux available for allocation 
   is_active              => symbiont_inst%is_active                           , & ! Input: [logical (:,:)] if symbiont uptake pathway is active for patch
   perecm                 => pftcon%perecm                                     , & ! Input: The fraction of ECM-associated PFT 
   C_reservoir            => symbiont_inst%C_reservoir                         , & ! Carbon reservoir in intermediate pools[gC/m2/s]
   N_reservoir            => symbiont_inst%N_reservoir                         , & ! Nitrogen reservoir in intermediate pools[gN/m2/s]
   C_biomass              => symbiont_inst%C_biomass                           , & ! Carbon biomass of symbiont [gC/m2]
   N_biomass              => symbiont_inst%N_biomass                           , & ! Nitrogen biomass of symbiont [gN/m2]
   symb_eff               => symbiont_inst%symb_eff                            , & ! Symbiont efficiency N uptake per unit of mycorrhizal C biomass  [gN/gC/s]
   C_mortality            => symbiont_inst%C_mortality                         , & ! Symbiotic turnover per soil layer and column [gC/m3/s]
   N_mortality            => symbiont_inst%N_mortality                         , & ! Symbiotic turnover per soil layer and column [gN/m3/s]
   sminn_to_plant_mimicsplus => cnveg_nitrogenflux_inst%sminn_to_plant_mimicsplus_patch & ! Output:[real(r8) (:) ]  nitrogen sent to plant (gN/m2/s)
   
   )

   !--------------------------------------------------------------------------------------------------------------------------------

   ! RETURN OF INVESTMENT

   ! Calculate N released to plants per uptake pathway - Return Of Investment (line 553 in vegn_dynamics)
   ! multiplied with secspday * days_per_year to get values per second (parameter is in per year)

   days_per_year = get_average_days_per_year()     ! to get average number of days per year 

  
  do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)

   ! Scavengers
   if (is_active(p,i_scav)) then 
      N_to_plant_scav(p) = N_reservoir(p,i_scav) * params_inst%sulman_rup_veg
      if (C_biomass(p,i_scav) > 0.0_r8) then  ! or (C_biomass(p,i_scav) < 0.0_r8)
        ! scav_roi(p) = ((max(0.0_r8, N_to_plant_scav(p)) * dt) / (C_biomass(p,i_scav))) * params_inst%sulman_growth_scav / (params_inst%sulman_tau_scav * dt)
          scav_roi(p) = (max(0.0_r8, N_to_plant_scav(p)))  / (C_biomass(p,i_scav) * params_inst%sulman_growth_scav * (params_inst%sulman_tau_scav))
      else 
      ! scav_efficiency is calculated in myc_scavenger_N_uptake under myc_efficiency
      scav_roi(p) = symb_eff(p,i_scav) / (params_inst%sulman_growth_scav * (params_inst%sulman_tau_scav))
      end if 
   else
      N_to_plant_scav(p) = 0.0_r8 ; scav_roi = 0.0_r8
   end if 

   ! Miners
   if (is_active(p,i_miner)) then
      N_to_plant_mine(p) = N_reservoir(p,i_miner) * params_inst%sulman_rup_veg !* secspday * days_per_year
      if (C_biomass(p,i_miner) > 0.0_r8) then 
         mine_roi(p) = ((max(0.0_r8, N_to_plant_mine(p))) / (C_biomass(p,i_miner))) * params_inst%sulman_growth_mine / (params_inst%sulman_tau_mine)
      else 
         ! mine is calculated in one of the mining routines under myc_efficiency
         mine_roi(p) = symb_eff(p,i_miner) / (params_inst%sulman_growth_mine * params_inst%sulman_tau_mine)
      end if 
   else
      N_to_plant_mine(p) = 0.0_r8 ; mine_roi(p) = 0.0_r8
   end if 

   if (C_allocation_to_N_acq(p) > 0.0_r8) then
      root_roi(p) = max(0.001,(root_N_uptake(p)/dt)/C_allocation_to_N_acq(p))
   else
      root_roi(p) = (mine_roi(p) + scav_roi(p))*0.25_r8
   endif


   ! Nitrogen Fixers
   if (is_active(p,i_fixer)) then
      N_to_plant_fix(p) = N_reservoir(p,i_fixer) * params_inst%sulman_rup_veg !* secspday * days_per_year
      if (C_biomass(p,i_fixer) > 0.0_r8) then 
         fix_roi(p) = ((N_to_plant_fix(p)) / (C_biomass(p,i_fixer))) * params_inst%sulman_growth_fix / (params_inst%sulman_tau_fix)
      else 
         ! 
         fix_roi(p) =  params_inst%sulman_rfix / params_inst%sulman_growth_fix * params_inst%sulman_tau_fix
      end if 
   else
      N_to_plant_fix(p) = 0.0_r8 ; fix_roi = 0.0_r8
   end if 

   N_reservoir(p,i_scav) = N_reservoir(p,i_scav) - (N_to_plant_scav(p) * dt)
   N_reservoir(p,i_miner) = N_reservoir(p,i_miner) - (N_to_plant_mine(p) * dt)
   N_reservoir(p,i_fixer) = N_reservoir(p,i_fixer) - (N_to_plant_fix(p) * dt)
   
   !ECW THIS might cause problems with plant growth

    ! Total nitrogen uptake by plant from intermediate symbiont pools
   sminn_to_plant_mimicsplus(p) = N_to_plant_scav(p) + N_to_plant_mine(p) + N_to_plant_fix(p) 

   !  Scale N uptake to plant if it is bigger that the uptake capazitiy of plant
   ! make 15 a parameter, comes from leafcn_max = leafcn(ivt(p)) + 15.0_r8
   if (sminn_to_plant_mimicsplus(p) + root_N_to_plant(p) > C_allocation_to_N_acq(p) / plantCN(p)) then
      scale_N_to_plant(p) = (C_allocation_to_N_acq(p) / plantCN(p) - root_N_to_plant(p)) /  (sminn_to_plant_mimicsplus(p))
   else
      scale_N_to_plant(p) = 1.0_r8
   endif
      
   
   ! Scale uptake and return leftovers to reservoirs
      N_reservoir(p,i_scav)  = N_reservoir(p,i_scav)  + (N_to_plant_scav(p) * dt) * (1.0_r8 - scale_N_to_plant(p))
      N_reservoir(p,i_miner) = N_reservoir(p,i_miner) + (N_to_plant_mine(p) * dt) * (1.0_r8 - scale_N_to_plant(p))
      N_reservoir(p,i_fixer) = N_reservoir(p,i_fixer) + (N_to_plant_fix(p) * dt)  * (1.0_r8 - scale_N_to_plant(p))
      ! Nitrogen uptake to plant
      N_to_plant_scav(p) = N_to_plant_scav(p) * scale_N_to_plant(p)
      N_to_plant_mine(p) = N_to_plant_mine(p) * scale_N_to_plant(p)
      N_to_plant_fix(p)  = N_to_plant_fix(p)  * scale_N_to_plant(p)
       
      ! Since this variable is what plant actually gets, add root nitrogen here
       sminn_to_plant_mimicsplus(p) = sminn_to_plant_mimicsplus(p) * scale_N_to_plant(p) + root_N_to_plant(p)

      
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
    fix_C_alloc(p) = C_allocation_to_N_acq(p) * fix_roi_frac
    mine_C_alloc(p) = C_allocation_to_N_acq(p) * mine_roi_frac
    scav_C_alloc(p) = C_allocation_to_N_acq(p) * scav_roi_frac

   !Fraction of the plant to miner flux is send to SOMa, representing enzyme flux
   !MAYBE ADD HERE

   ! We are accumulating potential (not limited/smoothed) allocation so max allocation can increase over time
    !WHY are we doing this?
   ! fix_alloc_accum(p) = fix_alloc_accum(p) + C_allocation_to_N_acq(p) * fix_roi_frac
    !mine_alloc_accum(p) = mine_alloc_accum(p) + C_allocation_to_N_acq(p) * mine_roi_frac
    !scav_alloc_accum(p) = scav_alloc_accum(p) + C_allocation_to_N_acq(p) * scav_roi_frac

    C_reservoir(p,i_scav) = C_reservoir(p,i_scav) + (scav_C_alloc(p) * dt)
    C_reservoir(p,i_miner) = C_reservoir(p,i_miner) + (mine_C_alloc(p) * dt) 
    C_reservoir(p,i_fixer) = C_reservoir(p,i_fixer) + (fix_C_alloc(p) * dt)

      
    ! Carbon that wasn't spend on scav, miner or fixer (including root)
    root_exudate_C(p) = C_allocation_to_N_acq(p) - scav_C_alloc(p) - mine_C_alloc(p) -  fix_C_alloc(p)


   end do
   end associate

   end subroutine update_symbionts

end module CNSoilVegMIMICSplus
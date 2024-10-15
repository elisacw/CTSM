module CNFUNMIMICSplusMod

#include "shr_assert.h"

  !--------------------------------------------------------------------
  ! ! DESCRIPTION
  ! ! Updated version of theThe FUN model developed by Fisher et al. 2010
  ! ! for coupling with MIMICSplus (Aas et al., 2023) decomposition model
  !
  ! !USES: 
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use shr_log_mod                     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod                  , only : nan => shr_infnan_nan, assignment(=)
  use clm_varctl                      , only : iulog
  use decompMod                       , only : subgrid_level_gridcell, subgrid_level_column, subgrid_level_patch
  use clm_varpar                      , only : i_phys_som, i_chem_som, i_avl_som, i_ecm_myc, i_am_myc
  use PatchType                       , only : patch
  use ColumnType                      , only : col
  use pftconMod                       , only : pftcon, npcropmin
  use decompMod                       , only : bounds_type
  use clm_varctl                      , only : use_nitrif_denitrif,use_flexiblecn
  use CNSharedParamsMod               , only : use_matrixcn
  use abortutils                      , only : endrun
  use CNVegstateType                  , only : cnveg_state_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegnitrogenstateType          , only : cnveg_nitrogenstate_type
  use CNVegnitrogenfluxType           , only : cnveg_nitrogenflux_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemNitrogenStateType , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemStateType         , only : soilbiogeochem_state_type
  use SoilBiogeochemDecompCascadeConType, only : mimicsplus_decomp, decomp_method
  use WaterStateBulkType              , only : waterstatebulk_type
  use WaterFluxBulkType               , only : waterfluxbulk_type
  use TemperatureType                 , only : temperature_type
  use SoilStateType                   , only : soilstate_type
  use CanopyStateType                 , only : canopystate_type
  use perf_mod                        , only : t_startf, t_stopf
  use clm_varpar                      , only : nlevdecomp
  use clm_varcon                      , only : secspday,  spval
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! TYPES

  ! !PUBLIC MEMBER FUNCTIONS:
  public:: readParams            ! Read in parameters needed for FUNMimicplus
  public:: CNFUNMIMICSplusInit   ! FUNMIMICSplus calculation initialization
  public:: CNFUNMIMICSplus       ! FUNMIMICSplus calculation itself
  public:: updateCNFUNMIMICSplus 
  private :: fun_cost_nonmyc
  private :: fun_cost_fix
  private :: fun_retranslocation

  ! DATA STRUCTURES
  integer,  parameter :: ipecm    = 1
  integer,  parameter :: ipam     = 2
  integer,  parameter :: ipnmno3    = 3             ! Process number for active nonmyc uptake of NO3.
  integer,  parameter :: ipnmnh4    = 4
  integer,  parameter :: ipfix      = 5            ! Process number for fixing.                          orginal version in FUN: [icost...]
  
  real(r8), parameter :: big_cost        = 5.e+6_r8! An arbitrarily large cost
 
  !  array index when plant is fixing
  integer, parameter :: plants_are_fixing = 1
  integer, parameter :: plants_not_fixing = 2
 
  !  array index for ECM step versus AM step
  integer, parameter :: ecm_step          = 1
  integer, parameter :: am_step           = 2

  integer,  private, parameter :: nmyc            = 2             ! Number of calculation part          orginal version in FUN: [nstp]
  integer,  private, parameter :: npaths          = 5             ! Number of  N transport pathways     orginal version in FUN: [ncost6] retranslocation is handled differently
  !
  !  populated in readParamsMod
  type, public :: cnfunmimicsplus_type
  ! declare arrays and scalars here with pointers
  real(r8)          :: rootc_dens_step                                       ! root C for each PFT in each soil layer(gC/m2)
  real(r8)          :: dn                                                    ! Increment of N                        (gN/m2)  
  real(r8), pointer :: rootc_dens(:,:)                                       ! the root carbon density               (gC/m2)
  real(r8), pointer :: rootC(:)                                              ! root biomass                          (gC/m2)
  real(r8), pointer :: n_uptake_myc_frac(:)                                  ! the arrary for the ECM and AM ratio   (-)       orginal version in FUN: [permyc]
  real(r8), pointer :: plant_ndemand_pool(:)                                 ! The N demand pool (gN/m2)
  real(r8), pointer :: plant_ndemand_pool_step(:)                            ! the N demand pool (gN/m2)
  real(r8), pointer :: litterfall_n_step(:)                                  ! N loss based on the leafc to litter (gN/m2)
  real(r8), pointer :: litterfall_c_step(:)                                  ! N loss based on the leafc to litter (gN/m2)
  real(r8), pointer :: npp_remaining(:)                                      ! A temporary variable for npp_remaining(gC/m2) 
  real(r8), pointer :: n_retrans_acc(:)                                      ! N acquired by retranslocation (gN/m2)
  real(r8), pointer :: free_nretrans_acc(:)                                  ! N acquired by retranslocation (gN/m2)
  real(r8), pointer :: npp_retrans_acc(:)                                    ! NPP used for the extraction (gC/m2)
  real(r8), pointer :: nt_uptake(:)                                          ! N uptake from retrans, active, and fix(gN/m2)
  real(r8), pointer :: npp_uptake(:)                                         ! NPP used by the uptakes (gC/m2)
  real(r8), pointer :: sminfrc(:,:)                                          ! fraction of N to handle NO3 / NH4 input 
  real(r8), pointer :: sminn_to_plant(:,:)                                   ! Nitrogen to plant (to handle NO3 / NH4 input)

  !----------NITRIF_DENITRIF-------------!
  !REPLACE NO# with N
  real(r8)            :: sminn_diff                                          ! A temporary limit for N uptake                  (gN/m2)
  real(r8),  pointer  :: costs_paths(:,:,:)                                  ! Costs for all paths (gC/gN) [patch,nlev,ipath]
  real(r8),  pointer  :: npp_frac_paths(:,:,:)                               ! NPP fraction for all paths () [patch,nlev,ipath]
  real(r8),  pointer  :: npp_to_paths(:,:,:)                                 ! NPP spent all paths (gC/m2) [patch,nlev,ipath]
  real(r8),  pointer  :: n_from_paths(:,:,:)                                 ! NPP spent all paths (gN/m2) [patch,nlev,ipath]
  real(r8),  pointer  :: n_paths_acc (:,:)                                   ! accumulated N spent all paths (gN/m2) [patch,nlev,ipath]
  real(r8),  pointer  :: npp_paths_acc (:,:)                                 ! accumulated NPP spent all paths (gC/m2) [patch,nlev,ipath]
  real(r8),  pointer  :: sminn_layer(:,:)                                      ! Available N in each soil layer (gN/m2)
  real(r8),  pointer  :: sminn_layer_step(:,:)                               ! A temporary variable for soil N (gN/m2) 
  real(r8),  pointer  :: n_active_vr(:,:)                                    ! Layer mycorrhizal N uptake (gN/m2)
  real(r8),  pointer  :: n_nonmyc_vr(:,:)                                    ! Layer non-myc     N uptake (gN/m2)
  
  ! Uptake fluxes for COST_METHOD=2
  ! Update Fluxes
  real(r8),  pointer ::  no3_myc_to_plant_col(:,:)       ! No3 flux from min soil to plant  by mycorrhiza
  real(r8),  pointer ::  nh4_myc_to_plant_col(:,:)       ! NH4 flux from min soil to plant by mycorrhiza

contains
  procedure , public  :: Init
  procedure , private :: InitAllocate
  procedure , private :: SetZeros
   
end type cnfunmimicsplus_type

type, private :: params_type
   real(r8) :: ndays_on        ! number of days to complete leaf onset
   real(r8) :: ndays_off       ! number of days to complete leaf offset
end type params_type  

   type(params_type), private :: params_inst  ! params_inst is

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

   contains

  subroutine Init (this, bounds, soilbiogeochem_carbonflux_inst, soilbiogeochem_nitrogenflux_inst)
   !
   ! !DESCRIPTION:
   ! Initialize xxx
   !
   ! !USES:
   !
   ! !ARGUMENTS:
     class(cnfunmimicsplus_type)       :: this
     type(bounds_type) , intent(in)    :: bounds
     type(soilbiogeochem_carbonflux_type),   intent(inout) :: soilbiogeochem_carbonflux_inst
     type(soilbiogeochem_nitrogenflux_type), intent(inout) :: soilbiogeochem_nitrogenflux_inst

     call this%InitAllocate(bounds)
     call this%SetZeros(bounds,soilbiogeochem_carbonflux_inst,soilbiogeochem_nitrogenflux_inst)

  end subroutine Init


  subroutine InitAllocate (this, bounds)
   !
   ! !DESCRIPTION:
   ! Initialize module data structure
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   class(cnfunmimicsplus_type)       :: this
   type(bounds_type) , intent(in)    :: bounds  

   integer :: begc, endc
   integer :: begg, endg
   integer :: begp, endp

   begg = bounds%begg; endg = bounds%endg
   begc = bounds%begc; endc = bounds%endc
   begp = bounds%begp; endp = bounds%endp

     allocate(this%rootc_dens(bounds%begp:bounds%endp,1:nlevdecomp));        this%rootc_dens(:,:)= spval

     allocate(this%rootC(bounds%begp:bounds%endp));                          this%rootC(:)= spval

     allocate(this%n_uptake_myc_frac(bounds%begp:bounds%endp));              this%n_uptake_myc_frac(:) = spval

     allocate(this%plant_ndemand_pool(bounds%begp:bounds%endp));             this%plant_ndemand_pool(:) = spval

     allocate(this%plant_ndemand_pool_step(bounds%begp:bounds%endp));        this%plant_ndemand_pool_step(:) = spval
     
     allocate(this%litterfall_n_step(bounds%begp:bounds%endp));              this%litterfall_n_step(:) = spval
     
     allocate(this%litterfall_c_step(bounds%begp:bounds%endp));              this%litterfall_c_step(:) = spval
     
     allocate(this%npp_remaining(bounds%begp:bounds%endp));                  this%npp_remaining(:) = spval
     
     allocate(this%n_retrans_acc(bounds%begp:bounds%endp));                  this%n_retrans_acc(:) = spval

     allocate(this%free_nretrans_acc(bounds%begp:bounds%endp));              this%free_nretrans_acc(:) = spval
     
     allocate(this%npp_retrans_acc(bounds%begp:bounds%endp));                this%npp_retrans_acc(:) = spval
     
     allocate(this%nt_uptake(bounds%begp:bounds%endp));                      this%nt_uptake(:) = spval
     
     allocate(this%npp_uptake(bounds%begp:bounds%endp));                     this%npp_uptake(:) = spval

     allocate(this%costs_paths(bounds%begp:bounds%endp,1:nlevdecomp,1:npaths));  this%costs_paths(:,:,:) = spval
     
     allocate(this%sminn_layer(bounds%begc:bounds%endc,1:nlevdecomp));       this%sminn_layer(:,:) = spval

     allocate(this%sminn_layer_step(bounds%begp:bounds%endp,1:nlevdecomp));  this%sminn_layer_step(:,:) = spval

     allocate(this%n_active_vr(bounds%begp:bounds%endp, 1:nlevdecomp));      this%n_active_vr(:,:) = spval

     allocate(this%n_nonmyc_vr(bounds%begp:bounds%endp, 1:nlevdecomp));      this%n_nonmyc_vr(:,:) = spval

     allocate(this%sminfrc(bounds%begc:bounds%endc,1:nlevdecomp));           this%sminfrc(:,:) = spval

     allocate(this%sminn_to_plant(bounds%begc:bounds%endc,1:nlevdecomp));    this%sminn_to_plant(:,:) = spval

    ! Uptake fluxes for COST_METHOD=2
    ! actual npp to each layer for each uptake process
     allocate(this%npp_to_paths(bounds%begp:bounds%endp,1:nlevdecomp,1:npaths));     this%npp_to_paths(:,:,:) = spval

     allocate(this%npp_frac_paths(bounds%begp:bounds%endp,1:nlevdecomp,1:npaths));   this%npp_frac_paths(:,:,:) = spval

     allocate(this%n_from_paths(bounds%begp:bounds%endp,1:nlevdecomp,1:npaths));     this%n_from_paths(:,:,:) = spval

     allocate(this%n_paths_acc(bounds%begp:bounds%endp,1:npaths));                   this%n_paths_acc(:,:) = spval

     allocate(this%npp_paths_acc(bounds%begp:bounds%endp,1:npaths));                 this%npp_paths_acc(:,:) = spval

     allocate(this%no3_myc_to_plant_col(bounds%begc:bounds%endc, 1:nlevdecomp));     this%no3_myc_to_plant_col(:,:) = spval

     allocate(this%nh4_myc_to_plant_col(bounds%begc:bounds%endc, 1:nlevdecomp));     this%nh4_myc_to_plant_col(:,:) = spval

  end subroutine InitAllocate

  subroutine SetZeros (this, bounds, soilbiogeochem_carbonflux_inst,soilbiogeochem_nitrogenflux_inst)
   !
   ! !DESCRIPTION:
   ! Sets module data structure to zero
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   class(cnfunmimicsplus_type)                           :: this
   type(bounds_type) ,                     intent(in)    :: bounds 
   type(soilbiogeochem_carbonflux_type),   intent(inout) :: soilbiogeochem_carbonflux_inst
   type(soilbiogeochem_nitrogenflux_type), intent(inout) :: soilbiogeochem_nitrogenflux_inst
   integer :: begc, endc
   integer :: begg, endg
   integer :: begp, endp

   begg = bounds%begg; endg = bounds%endg
   begc = bounds%begc; endc = bounds%endc
   begp = bounds%begp; endp = bounds%endp

   associate( cf   => soilbiogeochem_carbonflux_inst , &
              nf   =>  soilbiogeochem_nitrogenflux_inst )
 
     this%rootc_dens_step                                               = 0._r8
     this%sminn_diff                                                    = 0._r8
     this%rootc_dens(begp:endp,1:nlevdecomp)                            = 0._r8
     this%rootC(begp:endp)                                              = 0._r8
     this%n_uptake_myc_frac(begp:endp)                                  = 0._r8
     this%plant_ndemand_pool(begp:endp)                                 = 0._r8
     this%plant_ndemand_pool_step(begp:endp)                            = 0._r8
     this%litterfall_n_step(begp:endp)                                  = 0._r8
     this%litterfall_c_step(begp:endp)                                  = 0._r8
     this%npp_remaining(begp:endp)                                      = 0._r8
     this%costs_paths(begp:endp,1:nlevdecomp,1:npaths)                  = big_cost
     this%npp_to_paths(begp:endp,1:nlevdecomp,1:npaths)                 = 0._r8
     this%npp_frac_paths(begp:endp,1:nlevdecomp,1:npaths)               = 0._r8
     this%n_from_paths(begp:endp,1:nlevdecomp,1:npaths)                 = 0._r8
     this%n_paths_acc(bounds%begp:bounds%endp,1:npaths)                 = 0._r8
     this%npp_paths_acc(bounds%begp:bounds%endp,1:npaths)               = 0._r8
     this%n_retrans_acc(bounds%begp:bounds%endp)                        = 0._r8
     this%free_nretrans_acc(bounds%begp:bounds%endp)                    = 0._r8
     this%npp_retrans_acc(bounds%begp:bounds%endp)                      = 0._r8
     this%nt_uptake(bounds%begp:bounds%endp)                            = 0._r8
     this%npp_uptake(bounds%begp:bounds%endp)                           = 0._r8
     this%sminn_layer(bounds%begc:bounds%endc,1:nlevdecomp)             = 0._r8
     this%sminn_layer_step(bounds%begc:bounds%endc,1:nlevdecomp)        = 0._r8
     this%n_active_vr(bounds%begp:bounds%endp, 1:nlevdecomp)            = 0._r8
     this%n_nonmyc_vr(bounds%begp:bounds%endp, 1:nlevdecomp)            = 0._r8
     this%sminfrc(bounds%begc:bounds%endc,1:nlevdecomp)                 = 0._r8
     this%sminn_to_plant(bounds%begc:bounds%endc,1:nlevdecomp)          = 0._r8

     cf%c_am_resp_vr_col(begc:endc,1:nlevdecomp)                        = 0._r8
     cf%c_ecm_resp_vr_col(begc:endc,1:nlevdecomp)                       = 0._r8
     cf%c_am_growth_vr_col(begc:endc,1:nlevdecomp)                      = 0._r8
     cf%c_ecm_growth_vr_col(begc:endc,1:nlevdecomp)                     = 0._r8
     nf%n_am_growth_vr_col(begc:endc,1:nlevdecomp)                      = 0._r8
     nf%n_ecm_growth_vr_col(begc:endc,1:nlevdecomp)                     = 0._r8
     cf%c_ecm_enz_vr_col(begc:endc,1:nlevdecomp)                        = 0._r8
     nf%n_somc2ecm_vr_col(begc:endc,1:nlevdecomp)                       = 0._r8
     nf%n_somp2ecm_vr_col(begc:endc,1:nlevdecomp)                       = 0._r8
     cf%c_somc2soma_vr_col(begc:endc,1:nlevdecomp)                      = 0._r8
     cf%c_somp2soma_vr_col(begc:endc,1:nlevdecomp)                      = 0._r8
     nf%sminno3_to_ecm_vr_col(begc:endc,1:nlevdecomp)                   = 0._r8
     nf%sminno3_to_am_vr_col(begc:endc,1:nlevdecomp)                    = 0._r8
     nf%sminnh4_to_ecm_vr_col(begc:endc,1:nlevdecomp)                   = 0._r8
     nf%sminnh4_to_am_vr_col(begc:endc,1:nlevdecomp)                    = 0._r8
     nf%sminno3_nonmyc_to_plant_col(begc:endc,1:nlevdecomp)             = 0._r8
     nf%sminnh4_nonmyc_to_plant_col(begc:endc,1:nlevdecomp)             = 0._r8
     this%no3_myc_to_plant_col(begc:endc,1:nlevdecomp)                  = 0._r8
     this%nh4_myc_to_plant_col(begc:endc,1:nlevdecomp)                  = 0._r8

   end associate
  end subroutine SetZeros

 subroutine readParams ( ncid )
  ! 
  ! !DESCRIPTION:
  ! Read in parameters
  !
  ! !USES:
  use ncdio_pio , only : file_desc_t,ncd_io
  !
  ! !ARGUMENTS:
  implicit none
  type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
  !
  ! !LOCAL VARIABLES:
  character(len=32)  :: subname = 'CNFUNParamsType'
  character(len=100) :: errCode = '-Error reading in parameters file:'
  logical            :: readv ! has variable been read in or not
  real(r8)           :: tempr ! temporary to read in parameter
  character(len=100) :: tString ! temp. var for reading

  tString='ndays_on'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
  params_inst%ndays_on=tempr

  tString='ndays_off'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
  params_inst%ndays_off=tempr
    
 end subroutine readParams


 subroutine CNFUNMIMICSplusInit (bounds,cnveg_state_inst,cnveg_carbonstate_inst,cnveg_nitrogenstate_inst)
  !
  ! !DESCRIPTION:
  ! Initialization of FUN for MIMICSplus
  !
  ! !USES:
  use clm_varcon      , only: secspday, fun_period
  use clm_time_manager, only: get_step_size_real,get_nstep,get_curr_date,get_curr_days_per_year
  use SoilBiogeochemDecompCascadeMIMICSplusMod, only: calc_myc_roi, cost_FUN, fun_fluxes_myc_update1
  !
  ! !ARGUMENTS:
  type(bounds_type)             , intent(in)    :: bounds
  type(cnveg_state_type)        , intent(inout) :: cnveg_state_inst
  type(cnveg_carbonstate_type)  , intent(inout) :: cnveg_carbonstate_inst
  type(cnveg_nitrogenstate_type), intent(inout) :: cnveg_nitrogenstate_inst
  !
  ! !LOCAL VARIABLES:         
  real(r8)          :: ndays_on        ! number of days to complete leaf onset
  real(r8)          :: ndays_off       ! number of days to complete leaf offset
  real(r8)          :: dt              ! timestep size (seconds)
  real(r8)          :: dayspyr         ! days per year (days)
  real(r8)          :: timestep_fun    ! Timestep length for FUN (s)
  real(r8)          :: numofyear       ! number of days per year
  integer           :: nstep           ! time step number
  integer           :: nstep_fun       ! Number of atmospheric timesteps between calls to FUN

  character(len=32) :: subname = 'CNFUNMIMICSplusInit'

  ! Set local pointers
  associate(ivt                 => patch%itype                                          , & ! Input:  [integer  (:)   ]  p
         leafcn                 => pftcon%leafcn                                        , & ! Input:  leaf C:N (gC/gN)
         leafcn_offset          => cnveg_state_inst%leafcn_offset_patch                 , & ! Output:
         !  [real(r8) (:)   ]  Leaf C:N used by FUN  
         leafc_storage_xfer_acc => cnveg_carbonstate_inst%leafc_storage_xfer_acc_patch  , & ! Output: [real(r8) (:)
         !   ]  Accmulated leaf C transfer (gC/m2)
         storage_cdemand        => cnveg_carbonstate_inst%storage_cdemand_patch         , & ! Output: [real(r8) (:)
         !   ]  C use from the C storage pool
         leafn_storage_xfer_acc => cnveg_nitrogenstate_inst%leafn_storage_xfer_acc_patch, & ! Output: [real(r8) (:)
         !  ]  Accmulated leaf N transfer (gC/m2)                 
         storage_ndemand        => cnveg_nitrogenstate_inst%storage_ndemand_patch         & ! Output: [real(r8) (:)
         !  ]  N demand during the offset period 
         )
 
  ! Calculate some timestep-related values.
  
  ! set time steps
  dt           = get_step_size_real()
  dayspyr      = get_curr_days_per_year()
  nstep        = get_nstep()
  timestep_fun = real(secspday * fun_period)
  nstep_fun    = int(secspday * dayspyr / dt) 

  ndays_on     = params_inst%ndays_on
  ndays_off    = params_inst%ndays_off

  ! Decide if FUN will be called on this timestep.
  
  numofyear = nstep/nstep_fun
  if (mod(nstep,nstep_fun) == 0) then
     leafcn_offset(bounds%begp:bounds%endp)          = leafcn(ivt(bounds%begp:bounds%endp))
     storage_cdemand(bounds%begp:bounds%endp)        = 0._r8
     storage_ndemand(bounds%begp:bounds%endp)        = 0._r8
     leafn_storage_xfer_acc(bounds%begp:bounds%endp) = 0._r8
     leafc_storage_xfer_acc(bounds%begp:bounds%endp) = 0._r8
  end if   

  end associate
  end subroutine CNFUNMIMICSplusInit 


subroutine CNFUNMIMICSplus (bounds, num_soilc, filter_soilc, num_soilp ,filter_soilp, &
   waterstatebulk_inst, waterfluxbulk_inst, temperature_inst,soilstate_inst, &
   cnveg_state_inst,cnveg_carbonstate_inst, cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,cnveg_nitrogenflux_inst, &
   soilbiogeochem_nitrogenflux_inst,soilbiogeochem_carbonflux_inst,canopystate_inst, soilbiogeochem_nitrogenstate_inst, &
   soilbiogeochem_carbonstate_inst, cnfunmimicsplus_inst, soilbiogeochem_state_inst)
   
! !USES:
   use clm_time_manager, only : get_step_size_real, get_curr_date
   use clm_varpar      , only : nlevdecomp
   use clm_varcon      , only : secspday, smallValue, fun_period, tfrz, dzsoi_decomp, spval
   use clm_varctl      , only : use_nitrif_denitrif
   use PatchType       , only : patch
   use subgridAveMod   , only : p2c
   use pftconMod       , only : npcropmin
   use SoilBiogeochemDecompCascadeMIMICSplusMod , only : calc_myc_roi, decomp_rates_mimicsplus, calc_myc_mortality, &
                                                         calc_myc_mining_rates, cost_FUN, fun_fluxes_myc_update1, myc_n_extraction, myc_cn_fluxes
!
! !ARGUMENTS: 
   type(bounds_type)                       , intent(in)    :: bounds
   integer                                 , intent(in)    :: num_soilc             ! number of soil columns in filter
   integer                                 , intent(in)    :: filter_soilc(:)       ! filter for soil columns
   integer                                 , intent(in)    :: num_soilp             ! number of soil patches in filter
   integer                                 , intent(in)    :: filter_soilp(:)       ! filter for soil patches
   type(waterstatebulk_type)               , intent(in)    :: waterstatebulk_inst
   type(waterfluxbulk_type)                , intent(in)    :: waterfluxbulk_inst
   type(temperature_type)                  , intent(in)    :: temperature_inst
   type(soilstate_type)                    , intent(in)    :: soilstate_inst
   type(cnveg_state_type)                  , intent(inout) :: cnveg_state_inst
   type(cnveg_carbonstate_type)            , intent(inout) :: cnveg_carbonstate_inst
   type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenstate_type)          , intent(inout) :: cnveg_nitrogenstate_inst
   type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
   type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
   type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst 
   type(canopystate_type)                  , intent(inout) :: canopystate_inst  
   type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
   type(soilbiogeochem_carbonstate_type)   , intent(inout) :: soilbiogeochem_carbonstate_inst
   type(cnfunmimicsplus_type)              , intent(inout) :: cnfunmimicsplus_inst
   type(soilbiogeochem_state_type)         , intent(in)    :: soilbiogeochem_state_inst

   ! LOCAL VARIABLES:
   integer   :: fn                    ! number of values in pft filter
   integer   :: fp                    ! lake filter pft index 
   integer   :: fc                    ! lake filter column index
   integer   :: p, c                  ! pft index
   integer   :: g, l                  ! indices
   integer   :: j, i, k               ! soil/snow level index
   integer   :: imyc                  ! Loop counters/work                                    orginal version in FUN: [istp]
   integer   :: ipath                 ! a local index                                         orginal version in FUN: [icost]
   integer   :: fixer                 ! 0 = non-fixer, 1 =fixer 
   logical   :: unmetDemand           ! True while there is still demand for N
   integer   :: FIX 
   real(r8)  :: dt                    ! timestep size (s)
   real(r8)  :: dnpp                  ! Soulution to main FUN equation, C that actually goes into pathways
   real(r8)  :: stepspday             ! number of timesteps per day (real)
   real(r8)  :: ndays_on              ! number of days to complete leaf onset
   real(r8)  :: ndays_off             ! number of days to complete leaf offset
   real(r8)  :: frac_ideal_C_use      ! How much less C do we use for 'buying' N than that needed to get to the ideal ratio?  fraction. 
   real(r8)  :: N_acquired
   real(r8)  :: N_before_corr(1:npaths)
   real(r8)  :: C_before_corr(1:npaths)
   real(r8)  :: C_spent
   real(r8)  :: leaf_narea            ! leaf n per unit leaf area in gN/m2 (averaged across canopy, which is OK for the cost calculation)
   real(r8)  :: sum_n_acquired        ! Sum N aquired from one unit of C (unitless)  
   real(r8)  :: burned_off_carbon     ! carbon wasted by poor allocation algorithm. If this is too big, we need a better iteration. 
   real(r8)  :: temp_n_flux  
   real(r8)  :: delta_cn              ! difference between 'ideal' leaf CN ration and actual leaf C:N ratio. C/N
   real(r8)  :: excess_carbon         ! how much carbon goes into the leaf C pool on account of the flexibleCN modifications.   
   real(r8)  :: excess_carbon_acc     ! excess accumulated over layers.
   real(r8)  :: roi_ecm               ! RoI nitrogen per carbon invested from vegetation to mycorrhiza [gN/gC]
   real(r8)  :: roi_am                ! RoI nitrogen per carbon invested from vegetation to mycorrhiza [gN/gC]
   real(r8)  :: fixer_mult            ! multiplier for fixer


   
   !  WITHOUT GROWTH RESP
   real(r8) ::                    fixerfrac            ! what fraction of plants can fix?
   real(r8) ::                    npp_to_spend         ! how much carbon do we need to get rid of?
   real(r8) ::                    npp_to_spend_init    ! npp that is available for nuptake
   real(r8) ::                    npp_to_spend_fix     ! temporary variable to check if npp from fixation is sammaler than all pathways
   real(r8) ::                    npp_spent            ! temporary
   real(r8) ::                    soil_n_extraction    ! calculates total N pulled from soil
   real(r8) ::                    total_N_conductance  ! inverse of C to of N for whole soil-leaf pathway
   real(r8) ::                    total_N_resistance   ! C to of N for whole soil -leaf pathway
   real(r8) ::                    paid_for_n_retrans
   real(r8) ::                    free_n_retrans
   real(r8) ::                    total_c_spent_retrans
   real(r8) ::                    total_c_accounted_retrans

   real(r8) :: frac_alloc_ecm                     ! fraction from ROI allocated to EcM

   ! soil mineral nitrogen availible
   real(r8) :: sminno3_to_paths(bounds%begp:bounds%endp, 1:nlevdecomp)   ! miniral no3 available
   real(r8) :: sminnh4_to_paths(bounds%begp:bounds%endp, 1:nlevdecomp)   ! miniral nh4 available
   real(r8) :: sminno3_extracted                                         ! 
   real(r8) :: sminnh4_extracted                                         !
   real(r8) :: sminno3_overlimit                                         !
   real(r8) :: sminnh4_overlimit                                         !

   real(r8) :: c_am_resp_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)          ! carbon respiration flux for AM mycorrhiza
   real(r8) :: c_ecm_resp_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)         ! carbon respiration flux for ECM mycorrhiza
   real(r8) :: c_am_growth_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)        ! carbon growth flux for AM mycorrhiza
   real(r8) :: c_ecm_growth_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)       ! carbon growth flux for ECM mycorrhiza
   real(r8) :: n_am_growth_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)        ! nitrogen growth flux for AM mycorrhiza
   real(r8) :: n_ecm_growth_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)       ! nitrogen growth flux for ECM mycorrhiza
   real(r8) :: c_ecm_enz_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)          ! carbon enzyme production flux for ECM mycorrhiza
   real(r8) :: n_somc2ecm_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)         ! nitrogen mining from ECM mycorrhiza
   real(r8) :: n_somp2ecm_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)         ! nitrogen mining from ECM mycorrhiza
   real(r8) :: c_somc2soma_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)        ! carbon release from mining from somc pool
   real(r8) :: c_somp2soma_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)        ! carbon release from mining from somp pool
   
   real(r8) :: n_ecm(bounds%begp:bounds%endp, 1:nlevdecomp)             ! Layer mycorrhizal no3 uptake (gN/m2)
   real(r8) :: n_nonmyc_no3_vr(bounds%begp:bounds%endp, 1:nlevdecomp)             ! Layer non-myc no3 uptake (gN/m2)
   real(r8) :: n_am(bounds%begp:bounds%endp, 1:nlevdecomp)             ! Layer mycorrhizal nh4 uptake (gN/m2)
   real(r8) :: n_nonmyc_nh4_vr(bounds%begp:bounds%endp, 1:nlevdecomp)             ! Layer non-myc nh4 uptake (gN/m2)
   real(r8) :: sminno3_to_ecm_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)     ! No3 flux from soil to ECM
   real(r8) :: sminno3_to_am_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)      ! No3 flux from soil to AM
   real(r8) :: sminnh4_to_ecm_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)     ! NH4 flux from soil to ECM
   real(r8) :: sminnh4_to_am_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp)      ! NH4 flux from soil to AM

   associate(ivt             => patch%itype                                          , & ! Input:   [integer  (:) ]  p
      leafcn                 => pftcon%leafcn                                        , & ! Input:   leaf C:N (gC/gN)
      season_decid           => pftcon%season_decid                                  , & ! Input:   binary flag for seasonal deciduous leaf habit (0 or 1)
      stress_decid           => pftcon%stress_decid                                  , & ! Input:   binary flag for stress deciduous leaf habit (0 or 1)
      a_fix                  => pftcon%a_fix                                         , & ! Input:   A BNF parameter
      b_fix                  => pftcon%b_fix                                         , & ! Input:   A BNF parameter
      c_fix                  => pftcon%c_fix                                         , & ! Input:   A BNF parameter
      s_fix                  => pftcon%s_fix                                         , & ! Input:   A BNF parameter
      kc_nonmyc              => pftcon%kc_nonmyc                                     , & ! Input:   A non-mycorrhizal uptake parameter
      kn_nonmyc              => pftcon%kn_nonmyc                                     , & ! Input:   A non-mycorrhizal uptake parameter
      perecm                 => pftcon%perecm                                        , & ! Input:   The fraction of ECM associated PFT 
      grperc                 => pftcon%grperc                                        , & ! Input:   Growth respiration per unit Carbon gained. (fraction)
      fun_cn_flex_a          => pftcon%fun_cn_flex_a                                 , & ! Parameter a of FUN-flexcn link code (def 5)
      fun_cn_flex_b          => pftcon%fun_cn_flex_b                                 , & ! Parameter b of FUN-flexcn link code (def 200)
      fun_cn_flex_c          => pftcon%fun_cn_flex_c                                 , & ! Parameter b of FUN-flexcn link code (def 80)         
      FUN_fracfixers         => pftcon%FUN_fracfixers                                , & ! Fraction of C that can be used for fixation.    
      leafcn_offset          => cnveg_state_inst%leafcn_offset_patch                 , & ! Output:  [real(r8)  (:)]  Leaf C:N used by FUN
      plantCN                => cnveg_state_inst%plantCN_patch                       , & ! Output:  [real(r8)  (:)]  Plant C:N used by FUN
      onset_flag             => cnveg_state_inst%onset_flag_patch                    , & ! Output:  [real(r8)  (:)]  onset flag
      offset_flag            => cnveg_state_inst%offset_flag_patch                   , & ! Output:  [real(r8)  (:)]  offset flag
      availc                 => cnveg_carbonflux_inst%availc_patch                   , & ! Iutput:  [real(r8)  (:)]  C flux available for allocation (gC/m2/s)
      leafc                  => cnveg_carbonstate_inst%leafc_patch                   , & ! Input:   [real(r8)  (:)]  (gC/m2) leaf C
      leafc_storage          => cnveg_carbonstate_inst%leafc_storage_patch           , & ! Input:   [real(r8) (:)]  (gC/m2) leaf C storage
      frootc                 => cnveg_carbonstate_inst%frootc_patch                  , & ! Input:   [real(r8) (:)]  (gC/m2) fine root C
      livestemc              => cnveg_carbonstate_inst%livestemc_patch               , & ! Input:   [real(r8) (:)]  (gC/m2) live stem C
      livecrootc             => cnveg_carbonstate_inst%livecrootc_patch              , & ! Input:   [real(r8) (:)]  (gC/m2) live coarse root C
      leafc_storage_xfer_acc => cnveg_carbonstate_inst%leafc_storage_xfer_acc_patch  , & ! uutput:  [real(r8) (:)]  Accmulated leaf C transfer (gC/m2)
      leafn                  => cnveg_nitrogenstate_inst%leafn_patch                 , & ! Input:   [real(r8)  (:)] (gN/m2) leaf N
      retransn               => cnveg_nitrogenstate_inst%retransn_patch              , & ! Input:   [real(r8)  (:)] (gN/m2) live coarse root N
      leafn_storage_xfer_acc => cnveg_nitrogenstate_inst%leafn_storage_xfer_acc_patch, & ! Output:  [real(r8)  (:)] Accmulated leaf N transfer (gC/m2)
      storage_ndemand        => cnveg_nitrogenstate_inst%storage_ndemand_patch       , & ! Output:  [real(r8)  (:)] N demand during the offset period
      leafc_to_litter        => cnveg_carbonflux_inst%leafc_to_litter_patch          , & ! Output:  [real(r8) (:) ]  leaf C litterfall (gC/m2/s)
      leafc_to_litter_fun    => cnveg_carbonflux_inst%leafc_to_litter_fun_patch      , & ! Output:  [real(r8) (:) ]  leaf C litterfall used by FUN (gC/m2/s)
      npp_Nactive            => cnveg_carbonflux_inst%npp_Nactive_patch              , & ! Output:  [real(r8) (:) ]  Mycorrhizal N uptake used C (gC/m2/s)
      npp_Nnonmyc            => cnveg_carbonflux_inst%npp_Nnonmyc_patch              , & ! Output:  [real(r8) (:) ]  Non-mycorrhizal N uptake use C (gC/m2/s)
      npp_Nam                => cnveg_carbonflux_inst%npp_Nam_patch                  , & ! Output:  [real(r8) (:) ]  AM uptake use C (gC/m2/s)
      npp_Necm               => cnveg_carbonflux_inst%npp_Necm_patch                 , & ! Output:  [real(r8) (:) ]  ECM uptake use C (gC/m2/s)
      npp_Nactive_no3        => cnveg_carbonflux_inst%npp_Nactive_no3_patch          , & ! Output:  [real(r8) (:) ]  Mycorrhizal N uptake used C (gC/m2/s)
      npp_Nnonmyc_no3        => cnveg_carbonflux_inst%npp_Nnonmyc_no3_patch          , & ! Output:  [real(r8) (:) ]  Non-myco uptake use C (gC/m2/s) rrhizal N uptake (gN/m2/s)         
      npp_Nam_no3            => cnveg_carbonflux_inst%npp_Nam_no3_patch              , & ! Output:  [real(r8) (:) ]  AM uptake use C (gC/m2/s)
      npp_Necm_no3           => cnveg_carbonflux_inst%npp_Necm_no3_patch             , & ! Output:  [real(r8) (:) ]  ECM uptake use C (gC/m2/s)
      npp_Nactive_nh4        => cnveg_carbonflux_inst%npp_Nactive_nh4_patch          , & ! Output:  [real(r8) (:) ]  Mycorrhizal N uptake used C (gC/m2/s)
      npp_Nnonmyc_nh4        => cnveg_carbonflux_inst%npp_Nnonmyc_nh4_patch          , & ! Output:  [real(r8) (:) ]  Non-mycorrhizal N uptake used C (gC/m2/s)
      npp_Nam_nh4            => cnveg_carbonflux_inst%npp_Nam_nh4_patch              , & ! Output:  [real(r8) (:) ]  AM uptake used C(gC/m2/s)
      npp_Necm_nh4           => cnveg_carbonflux_inst%npp_Necm_nh4_patch             , & ! Output:  [real(r8) (:) ]  ECM uptake used C (gC/m2/s)
      npp_Nfix               => cnveg_carbonflux_inst%npp_Nfix_patch                 , & ! Output:  [real(r8) (:) ]  Symbiotic BNF used C (gC/m2/s)
      npp_Nretrans           => cnveg_carbonflux_inst%npp_Nretrans_patch             , & ! Output:  [real(r8) (:) ]  Retranslocation N uptake used C (gC/m2/s)
      npp_Nuptake            => cnveg_carbonflux_inst%npp_Nuptake_patch              , & ! Output:  [real(r8) (:) ]  Total N uptake of FUN used C (gC/m2/s)
      npp_growth             => cnveg_carbonflux_inst%npp_growth_patch               , & ! Output:  [real(r8) (:) ]  Total C used for growth in FUN  (gC/m2/s)  
      npp_burnedoff          => cnveg_carbonflux_inst%npp_burnedoff_patch            , & ! Output:  [real(r8) (:) ]  C  that cannot be used for N uptake(gC/m2/s)   
      leafn_storage_to_xfer  => cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch  , & ! Output:  [real(r8) (:) ]
      plant_ndemand_col     => soilbiogeochem_state_inst%plant_ndemand_col                   , & ! Input:  [real(r8) (:)   ]  column-level plant N demand
      plant_ndemand          => cnveg_nitrogenflux_inst%plant_ndemand_patch          , & ! Iutput:  [real(r8) (:) ]  N flux required to support initial GPP (gN/m2/s)
      plant_ndemand_retrans  => cnveg_nitrogenflux_inst%plant_ndemand_retrans_patch  , & ! Output:  [real(r8) (:) ]  N demand generated for FUN (gN/m2/s)
      Nactive                => cnveg_nitrogenflux_inst%Nactive_patch                , & ! Output:  [real(r8) (:) ]  Mycorrhizal N uptake (gN/m2/s)
      Nnonmyc                => cnveg_nitrogenflux_inst%Nnonmyc_patch                , & ! Output:  [real(r8) (:) ]  Non-mycorrhizal N uptake (gN/m2/s)
      Nam                    => cnveg_nitrogenflux_inst%Nam_patch                    , & ! Output:  [real(r8) (:) ]  AM  uptake (gN/m2/s)
      Necm                   => cnveg_nitrogenflux_inst%Necm_patch                   , & ! Output:  [real(r8) (:) ]  ECM uptake (gN/m2/s)
      Nactive_no3            => cnveg_nitrogenflux_inst%Nactive_no3_patch            , & ! Output:  [real(r8) (:) ]  Mycorrhizal N uptake (gN/m2/s)
      Nnonmyc_no3            => cnveg_nitrogenflux_inst%Nnonmyc_no3_patch            , & ! Output:  [real(r8) (:) ]  Non-mycorrhizal N uptake (gN/m2/s)
      Nam_no3                => cnveg_nitrogenflux_inst%Nam_no3_patch                , & ! Output:  [real(r8) (:) ]  AM uptake (gN/m2/s)
      Necm_no3               => cnveg_nitrogenflux_inst%Necm_no3_patch               , & ! Output:  [real(r8) (:) ]  ECM uptake (gN/m2/s)
      Nactive_nh4            => cnveg_nitrogenflux_inst%Nactive_nh4_patch            , & ! Output:  [real(r8) (:) ]  Mycorrhizal N uptake (gN/m2/s)
      Nnonmyc_nh4            => cnveg_nitrogenflux_inst%Nnonmyc_nh4_patch            , & ! Output:  [real(r8) (:) ]  Non-mycorrhizal N uptake (gN/m2/s)
      Nam_nh4                => cnveg_nitrogenflux_inst%Nam_nh4_patch                , & ! Output:  [real(r8) (:) ]  AM uptake (gN/m2/s)
      Necm_nh4               => cnveg_nitrogenflux_inst%Necm_nh4_patch               , & ! Output:  [real(r8) (:) ]  ECM uptake (gN/m2/s)
      Nfix                   => cnveg_nitrogenflux_inst%Nfix_patch                   , & ! Output:  [real(r8) (:) ] Symbiotic BNF (gN/m2/s)
      cost_nfix              => cnveg_nitrogenflux_inst%cost_Nfix_patch              , & ! Output:  [real(r8) (:) ]  Cost of fixation gC:gN
      cost_nactive           => cnveg_nitrogenflux_inst%cost_Nactive_patch           , & ! Output:  [real(r8) (:) ] Cost of active uptake gC:gN         
      cost_nretrans          => cnveg_nitrogenflux_inst%cost_Nretrans_patch          , & ! Output:  [real(r8) (:) ] Cost of retranslocation gC:gN         
      nuptake_npp_fraction_patch => cnveg_nitrogenflux_inst%nuptake_npp_fraction_patch    , & ! Output:  [real(r8) (:) ]  frac of NPP in NUPTAKE 
      c_allometry            => cnveg_state_inst%c_allometry_patch                   , & ! Output: [real(r8) (:)   ]  C allocation index (DIM)                
      n_allometry            => cnveg_state_inst%n_allometry_patch                   , & ! Output: [real(r8) (:)   ]  N allocation index (DIM)                
      leafn_storage          => cnveg_nitrogenstate_inst%leafn_storage_patch         , & ! Input:  [real(r8) (:) ]  (gN/m2) leaf N store
      nfix_to_sminn          => soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col   , & ! Output:  [real(r8) (:)] symbiotic/asymbiotic N fixation to soil mineral N (gN/m2/s)
      Nretrans               => cnveg_nitrogenflux_inst%Nretrans_patch               , & ! Output:  [real(r8) (:) ]  Retranslocation N uptake (gN/m2/s)
      Nretrans_season        => cnveg_nitrogenflux_inst%Nretrans_season_patch        , & ! Output:  [real(r8) (:) ]  Retranslocation N uptake (gN/m2/s)
      Nretrans_stress        => cnveg_nitrogenflux_inst%Nretrans_stress_patch        , & ! Output:  [real(r8) (:) ]  Retranslocation N uptake (gN/m2/s)
      Nuptake                => cnveg_nitrogenflux_inst%Nuptake_patch                , & ! Output:  [real(r8) (:) ]  Total N uptake of FUN (gN/m2/s)
      retransn_to_npool      => cnveg_nitrogenflux_inst%retransn_to_npool_patch               , & ! Output: [real(r8) (:)   ]  deployment of retranslocated N (gN/m2/s)
      free_retransn_to_npool => cnveg_nitrogenflux_inst%free_retransn_to_npool_patch          , & ! Output: [real(r8) uptake of free N from leaves (needed to allow RT during the night with no NPP
      sminn_to_plant_fun     => cnveg_nitrogenflux_inst%sminn_to_plant_fun_patch              , & ! Output: [real(r8) (:) ]  Total soil N uptake of FUN (gN/m2/s)
      sminn_to_plant_fun_no3_vr  => cnveg_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_patch   , & ! Output:  [real(r8) (:) ]  Total layer no3 uptake of FUN (gN/m2/s)
      sminn_to_plant_fun_nh4_vr  => cnveg_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_patch   , & ! Output:  [real(r8) (:) ]  Total layer nh4 uptake of FUN (gN/m2/s)
      sminn_to_plant_vr      => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_vr_col        , & ! Output:  [real(r8) (: ,:) ] [gN/m3/s]  
      smin_no3_to_plant_vr   => soilbiogeochem_nitrogenflux_inst%smin_no3_to_plant_vr_col     , & ! Input:  [real(r8) (:,:) ] [gN/m3/s]  
      smin_nh4_to_plant_vr   => soilbiogeochem_nitrogenflux_inst%smin_nh4_to_plant_vr_col     , & ! Input:  [real(r8) (:,:) ] [gN/m3/s]       
      soilc_change           => cnveg_carbonflux_inst%soilc_change_patch                      , & ! Output:  [real(r8)(:) ]  Used C from the soil (gC/m2/s)
      h2osoi_liq             => waterstatebulk_inst%h2osoi_liq_col                            , & ! Input:   [real(r8) (:,:)] liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
      t_soisno               => temperature_inst%t_soisno_col                                 , & ! Input:   [real(r8) (:,:)] soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
      crootfr                => soilstate_inst%crootfr_patch                                  , & ! Input:   [real(r8) (:,:)] fraction of roots for carbon in each soil layer  (nlevgrnd)
      costs_paths            => cnfunmimicsplus_inst%costs_paths                              , &
      npp_to_paths           => cnfunmimicsplus_inst%npp_to_paths                             , &
      npp_frac_paths         => cnfunmimicsplus_inst%npp_frac_paths                           , &
      n_from_paths           => cnfunmimicsplus_inst%n_from_paths                             , &
      n_paths_acc            => cnfunmimicsplus_inst%n_paths_acc                              , &
      npp_paths_acc          => cnfunmimicsplus_inst%npp_paths_acc                            , &
      dn                     => cnfunmimicsplus_inst%dn                                       , &
      free_nretrans_acc      => cnfunmimicsplus_inst%free_nretrans_acc                        , &
      litterfall_c_step      => cnfunmimicsplus_inst%litterfall_c_step                        , &
      litterfall_n_step      => cnfunmimicsplus_inst%litterfall_n_step                        , &
      n_active_vr            => cnfunmimicsplus_inst%n_active_vr                              , &
      n_nonmyc_vr            => cnfunmimicsplus_inst%n_nonmyc_vr                              , &
      n_retrans_acc          => cnfunmimicsplus_inst%n_retrans_acc                            , &
      n_uptake_myc_frac      => cnfunmimicsplus_inst%n_uptake_myc_frac                        , &
      npp_remaining          => cnfunmimicsplus_inst%npp_remaining                            , &
      npp_retrans_acc        => cnfunmimicsplus_inst%npp_retrans_acc                          , &
      npp_uptake             => cnfunmimicsplus_inst%npp_uptake                               , &
      nt_uptake              => cnfunmimicsplus_inst%nt_uptake                                , &
      plant_ndemand_pool     => cnfunmimicsplus_inst%plant_ndemand_pool                       , &
      plant_ndemand_pool_step => cnfunmimicsplus_inst%plant_ndemand_pool_step                 , &
      rootC                  => cnfunmimicsplus_inst%rootC                                    , &
      rootc_dens             => cnfunmimicsplus_inst%rootc_dens                               , &
      rootc_dens_step        => cnfunmimicsplus_inst%rootc_dens_step                          , &
      sminfrc                => cnfunmimicsplus_inst%sminfrc                                  , &
      sminn_layer_step       => cnfunmimicsplus_inst%sminn_layer_step                         , &
      sminn_to_plant         => cnfunmimicsplus_inst%sminn_to_plant                           , &
      
      ! C and N soil pools:
      decomp_cpools_vr      => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col           , &
      decomp_npools_vr      => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col          )
      
      dt           = get_step_size_real()
      stepspday    = secspday / dt


     ! Time step of FUNMIMICSplus: once per day                
     call cnfunmimicsplus_inst%SetZeros(bounds,soilbiogeochem_carbonflux_inst,soilbiogeochem_nitrogenflux_inst) ! set everything to zero

     do fp = 1,num_soilp          ! PFT Starts
     p = filter_soilp(fp)
   
     ! Calculate plant N based on; N lost with litterfall, root C, ...
      rootC(p)        =  frootc(p)

     ! Calculate patch plant C:N ratio
      if (n_allometry(p).gt.0._r8) then 
          plantCN(p)  = c_allometry(p)/n_allometry(p)
      else
          plantCN(p)  = 0._r8 
      end if
     end do                       ! PFT ends

      ! Calculates fraction of N uptake per mycorrhizal type per PFT
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = patch%column(p)
      ! make sure the accumulated vars are zeroed before accumulating
      n_uptake_myc_frac(p)=0.0_r8
      do j = 1, nlevdecomp
         c_am_resp_vr_patch(p,j)       =0.0_r8    ! carbon respiration flux for AM mycorrhiza
         c_ecm_resp_vr_patch(p,j)      =0.0_r8    ! carbon respiration flux for ECM mycorrhiza
         c_am_growth_vr_patch(p,j)     =0.0_r8    ! carbon growth flux for AM mycorrhiza
         c_ecm_growth_vr_patch(p,j)    =0.0_r8    ! carbon growth flux for ECM mycorrhiza
         n_am_growth_vr_patch(p,j)     =0.0_r8    ! nitrogen growth flux for AM mycorrhiza
         n_ecm_growth_vr_patch(p,j)    =0.0_r8    ! nitrogen growth flux for ECM mycorrhiza
         c_ecm_enz_vr_patch(p,j)       =0.0_r8    ! carbon enzyme production flux for ECM mycorrhiza
         n_somc2ecm_vr_patch(p,j)      =0.0_r8    ! nitrogen mining from ECM mycorrhiza
         n_somp2ecm_vr_patch(p,j)      =0.0_r8    ! nitrogen mining from ECM mycorrhiza
         c_somc2soma_vr_patch(p,j)     =0.0_r8    ! carbon release from mining from somc pool
         c_somp2soma_vr_patch(p,j)     =0.0_r8    ! carbon release from mining from somp pool
       
          n_ecm(p,j)                   =0.0_r8    ! Layer mycorrhizal no3 uptake (gN/m2)
          n_am(p,j)                    =0.0_r8    ! Layer mycorrhizal nh4 uptake (gN/m2)
          n_nonmyc_no3_vr(p,j)         =0.0_r8    ! Layer non-myc     no3 uptake (gN/m2)
          n_nonmyc_nh4_vr(p,j)         =0.0_r8    ! Layer non-myc     nh4 uptake (gN/m2)
          sminno3_to_ecm_vr_patch(p,j) =0.0_r8    ! No3 flux from soil to ECM
          sminno3_to_am_vr_patch(p,j)  =0.0_r8    ! No3 flux from soil to AM
          sminnh4_to_ecm_vr_patch(p,j) =0.0_r8    ! NH4 flux from soil to ECM
          sminnh4_to_am_vr_patch(p,j)  =0.0_r8    ! NH4 flux from soil to AM

          sminno3_to_paths(p,j)        =0.0_r8
          sminnh4_to_paths(p,j)        =0.0_r8

      end do
      n_uptake_myc_frac(p) = 1.0_r8 - n_uptake_myc_frac(p)

     if(leafc(p)>0.0_r8)then                       ! N available in leaf which fell off in this timestep. Same fraction loss as C.    
          litterfall_c_step(p)         =   dt * leafc_to_litter_fun(p) 
          litterfall_n_step(p)         =   dt * leafn(p) * leafc_to_litter_fun(p)/leafc(p)
     endif 
   end do
     
     
    do j = 1, nlevdecomp
       do fp = 1,num_soilp        ! PFT Starts
          p = filter_soilp(fp)
          c = patch%column(p)
          sminn_layer_step(p,j)=max((sminn_to_plant(c,j)) * dzsoi_decomp(j) * dt,0.0_r8) ! gN/m2
          sminno3_to_paths(p,j) = max(smin_no3_to_plant_vr(c,j) * dzsoi_decomp(j) * dt,0.0_r8) ! gN/m2
          sminnh4_to_paths(p,j) = max(smin_nh4_to_plant_vr(c,j) * dzsoi_decomp(j) * dt,0.0_r8) ! gN/m2
       end do
    end do

   ! PFT loop
   ! Calculating leaf C / N storage demand for retranslocation
pft:  do fp = 1,num_soilp        ! PFT Starts
      p = filter_soilp(fp)
      c = patch%column(p)
     
      sminn_to_plant_fun_nh4_vr(p,:) = 0._r8
      sminn_to_plant_fun_no3_vr(p,:) = 0._r8   
      burned_off_carbon              = 0.0_r8
      excess_carbon_acc              = 0.0_r8
     
      ! I have turned off this retranslocation functionality for now. To be rolled back in to a new version later on once the rest of the mode is working OK. RF

      if (season_decid(ivt(p)) == 1._r8.or.stress_decid(ivt(p)) == 1._r8) then
         if (onset_flag(p) == 1._r8) then
            leafn_storage_xfer_acc(p) = leafn_storage_xfer_acc(p) + leafn_storage_to_xfer(p) * dt
         end if
         if (offset_flag(p) == 1._r8) then
            storage_ndemand(p)        = leafn_storage_xfer_acc(p) / (ndays_off * stepspday)
            storage_ndemand(p)        = max(storage_ndemand(p),0._r8)
         else
            storage_ndemand(p)        = 0._r8   
         end if
      else
          storage_ndemand(p)          = 0._r8 
      end if   ! end for deciduous

      ! Avaliable carbon for growth or Nitrogen uptake
      !availc(p)            =  availc(p)        *  dt !!

      if (availc(p) > 0._r8) then
         do j = 1, nlevdecomp
            rootc_dens(p,j)     =  crootfr(p,j) * rootC(p)
         end do
      end if

      plant_ndemand_pool(p)     =  plant_ndemand(p) *  dt
      plant_ndemand_pool(p)     =  max(plant_ndemand_pool(p),0._r8)
      plant_ndemand_retrans(p)  =  storage_ndemand(p) 

      ! Nitrogen demand of plant & remaining carbon (NPP) per timestep
      unmetDemand              = .TRUE.
      plant_ndemand_pool_step(p)   = plant_ndemand_pool(p)
      npp_remaining(p)             = availc(p) * dt ! gC/m2 !og availc(p) *dt
         
      ! COST FIXATION PATHWAY
      ! checks which photosyntetic pathway plant has (C3 / C4) and if they can do nitrogen fixation   
      do j = 1, nlevdecomp
         if(pftcon%c3psn(patch%itype(p)).eq.1)then
           fixer=1
         else
           fixer=0
         endif
         costs_paths(p,j,ipfix)     = fun_cost_fix(fixer,a_fix(ivt(p)),b_fix(ivt(p))&
         ,c_fix(ivt(p)) ,big_cost,crootfr(p,j),s_fix(ivt(p)),t_soisno(c,j))
      end do

      ! ACTIVE UPTAKE
      ! Mycorrhizal Uptake Cost
      do j = 1,nlevdecomp
         rootc_dens_step             = rootc_dens(p,j)
         costs_paths(p,j,ipecm:ipam)=big_cost
         if (rootc_dens_step > 0.0_r8) then
         call calc_myc_roi(decomp_cpools_vr(c,j,i_ecm_myc),decomp_npools_vr(c,j,i_ecm_myc) , &
         decomp_cpools_vr(c,j,i_phys_som),decomp_cpools_vr(c,j,i_avl_som),decomp_cpools_vr(c,j,i_chem_som), &
         decomp_npools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_chem_som), &
         (smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)) * dt , ecm_step , dzsoi_decomp(j),big_cost, roi_ecm)
         call calc_myc_roi(decomp_cpools_vr(c,j,i_am_myc),decomp_npools_vr(c,j,i_am_myc) , &
         decomp_cpools_vr(c,j,i_phys_som),decomp_cpools_vr(c,j,i_avl_som),decomp_cpools_vr(c,j,i_chem_som), &
         decomp_npools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_chem_som), &
         (smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)) * dt, am_step , dzsoi_decomp(j),big_cost, roi_am)
         end if
         costs_paths(p,j,ipecm)=1./costs_paths(p,j,ipecm)
         costs_paths(p,j,ipam)=1./costs_paths(p,j,ipam) ! convert to C/N
      end do
      
      ! Non-mycorrhizal Uptake Cost
      do j = 1,nlevdecomp
         rootc_dens_step               = rootc_dens(p,j)
         costs_paths(p,j,ipnmno3)      = fun_cost_nonmyc(sminno3_to_paths(p,j) &
         ,big_cost,kc_nonmyc(ivt(p)),kn_nonmyc(ivt(p)) ,rootc_dens_step,crootfr(p,j),smallValue)
         costs_paths(p,j,ipnmnh4)      = fun_cost_nonmyc(sminnh4_to_paths(p,j) &
         ,big_cost,kc_nonmyc(ivt(p)),kn_nonmyc(ivt(p)) ,rootc_dens_step,crootfr(p,j),smallValue)
      end do
      npp_to_spend = 0.0_r8
      total_N_conductance = 1.0_r8/ (npaths * nlevdecomp * big_cost)    

      fix_loop: do FIX =plants_are_fixing, plants_not_fixing !loop around percentages of fixers and nonfixers, with differnt costs. 

          if(FIX==plants_are_fixing)then ! How much of the carbon in this PFT can in principle be used for fixation? 
            ! This is analagous to fixing the % of fixers for a given PFT - may not be realistic in the long run
            ! but prevents wholesale switching to fixer dominance during e.g. CO2 fertilization.  
            fixerfrac = FUN_fracfixers(ivt(p))
          else
            fixerfrac = 1.0_r8 - FUN_fracfixers(ivt(p))
          endif 
          npp_to_spend = npp_remaining(p)  * fixerfrac ! put parameter here.
          npp_to_spend_init = npp_to_spend
         ! has to be zeroed since depend on accumula
          n_from_paths(p,:,:) = 0.0_r8 !act and nonmyc boths nh4 and no3
          npp_frac_paths(p,:,:) = 0.0_r8
          npp_to_paths(p,:,:) = 0.0_r8
          !           Calculate Integrated Resistance OF WHOLE SOIL COLUMN
          
          sum_n_acquired      = 0.0_r8
          total_N_conductance = 1.0_r8/ (npaths * nlevdecomp * big_cost) 
          do j = 1, nlevdecomp
            ! Method changed from FUN-resistors method to a method which allocates fluxs based on conductance. rosief
            ! Sum the conductances             
           total_N_conductance  = total_N_conductance + 1._r8/ costs_paths(p,j,ipecm) + 1._r8/ costs_paths(p,j,ipam) &
                                 + 1._r8/ costs_paths(p,j,ipnmno3) + 1._r8/ costs_paths(p,j,ipnmnh4) !N/C
            if(FIX==plants_are_fixing)then
                total_N_conductance  = total_N_conductance  + 1._r8/ costs_paths(p,j,ipfix) 
            end if      
          end do 
          
          do j = 1, nlevdecomp


             ! Calculate npp allocation to pathways proportional to their exchange rate (N/C) 
             rootc_dens_step             = rootc_dens(p,j) 
             if (rootc_dens_step > 0._r8) then
               npp_frac_paths(p,j,ipecm) = (1._r8/costs_paths(p,j,ipecm)) / total_N_conductance
               npp_frac_paths(p,j,ipam) = (1._r8/costs_paths(p,j,ipam)) / total_N_conductance
               npp_frac_paths(p,j,ipnmno3) = (1._r8/costs_paths(p,j,ipnmno3)) / total_N_conductance
               npp_frac_paths(p,j,ipnmnh4) = (1._r8/costs_paths(p,j,ipnmnh4)) / total_N_conductance
               if(FIX==plants_are_fixing)then
                  npp_frac_paths(p,j,ipfix) = (1.0_r8 * 1._r8/costs_paths(p,j,ipfix)) / total_N_conductance
               else
                  npp_frac_paths(p,j,ipfix) = 0.0_r8 
               end if
             else 
               npp_frac_paths(p,j,ipecm:ipfix) = 0.0_r8
             endif
           
             ! Total N aquired from one unit of carbon  (N/C)
             sum_n_acquired = 0.0_r8
             sum_n_acquired = sum_n_acquired  + npp_frac_paths(p,j,ipecm)/costs_paths(p,j,ipecm) + npp_frac_paths(p,j,ipam)/costs_paths(p,j,ipam) &
                  + npp_frac_paths(p,j,ipnmno3)/costs_paths(p,j,ipnmno3) + npp_frac_paths(p,j,ipnmnh4)/costs_paths(p,j,ipnmnh4)
                                     
             if(FIX==plants_are_fixing)then
               sum_n_acquired= sum_n_acquired +  npp_frac_paths(p,j,ipfix)/costs_paths(p,j,ipfix)
             end if
          end do
          if (sum_n_acquired>0.0_r8) then
             total_N_resistance = 1.0_r8/sum_n_acquired
          else
             total_N_resistance = 1.0_r8/total_N_conductance !C/N
          endif 

          free_n_retrans = 0.0_r8
          !  Calculate appropriate degree of retranslocation
          if(leafc(p).gt.0.0_r8.and.litterfall_n_step(p)* fixerfrac>0.0_r8.and.ivt(p) <npcropmin)then
             call fun_retranslocation(p,dt,npp_to_spend,&
                           litterfall_c_step(p)* fixerfrac,&
                           litterfall_n_step(p)* fixerfrac,&
                           total_n_resistance, total_c_spent_retrans,total_c_accounted_retrans, &
                           free_n_retrans,paid_for_n_retrans, leafcn(ivt(p)), & 
                           grperc(ivt(p)), plantCN(p))
                            
          else
              total_c_accounted_retrans = 0.0_r8
              total_c_spent_retrans     = 0.0_r8
              total_c_accounted_retrans = 0.0_r8
              paid_for_n_retrans        = 0.0_r8
              free_n_retrans            = 0.0_r8
          endif
    
          !Add retrans fluxes in to total budgets.
          ! remove C from available pool, both directly spent and accounted for by N uptake
          npp_to_spend  = npp_to_spend - total_c_spent_retrans - total_c_accounted_retrans

          npp_retrans_acc(p) = npp_retrans_acc(p) + total_c_spent_retrans 
          ! add to to C spent pool                              
          n_retrans_acc(p)     = n_retrans_acc(p)     + paid_for_n_retrans
          free_nretrans_acc(p) = free_nretrans_acc(p) + free_n_retrans
          ! add N to the acquired from retrans pool 

          ! Spend C on extracting N
          if (plant_ndemand_pool_step(p) .gt. 0._r8) then    ! unmet demand

             if(use_flexiblecn)then   
                 if (leafn(p) == 0.0_r8) then           ! to avoid division by zero
                   delta_CN = fun_cn_flex_c(ivt(p))     ! Max CN ratio over standard
                 else
                   delta_CN = (leafc(p)+leafc_storage(p))/(leafn(p)+leafn_storage(p)) - leafcn(ivt(p)) ! leaf CN ratio                                                              
                 end if
                 ! C used for uptake is reduced if the cost of N is very high
                 frac_ideal_C_use = max(0.0_r8,1.0_r8 - (total_N_resistance-fun_cn_flex_a(ivt(p)))/fun_cn_flex_b(ivt(p)) )
                 ! then, if the plant is very much in need of N, the C used for uptake is increased accordingly.
                 if(delta_CN.lt.0.0_r8)then
                   frac_ideal_C_use = frac_ideal_C_use + (1.0_r8-frac_ideal_C_use)*min(1.0_r8, delta_CN/fun_cn_flex_c(ivt(p)))
                 end if
                 ! If we have too much N (e.g. from free N retranslocation) then make frac_ideal_c_use even lower.
                 ! For a CN delta of fun_cn_flex_c, then we reduce C expendiure to the minimum of 0.5.
                 ! This seems a little intense?
                 if(delta_CN .gt.0._r8 .and. frac_ideal_C_use.lt.1.0_r8)then
                   frac_ideal_C_use = frac_ideal_C_use + 0.5_r8*(1.0_r8*delta_CN/fun_cn_flex_c(ivt(p)))
                 end if
                 ! don't let this go above 1 or below an arbitrary minimum (to prevent zero N uptake).
                 frac_ideal_C_use = max(min(1.0_r8,frac_ideal_C_use),0.5_r8) 
             else
                 frac_ideal_C_use= 1.0_r8
          end if
           
          excess_carbon                 = npp_to_spend * (1.0_r8-frac_ideal_c_use)
          if(excess_carbon*(1.0_r8+grperc(ivt(p))).gt.npp_to_spend)then !prevent negative dnpp
               excess_carbon =  npp_to_spend/(1.0_r8+grperc(ivt(p)))
          endif
          excess_carbon_acc             = excess_carbon_acc + excess_carbon

          ! spend less C than you have to to meet the target, thus allowing C:N ratios to rise. 
          npp_to_spend         = npp_to_spend - excess_carbon*(1.0_r8+grperc(ivt(p)))
     
          ! This is the main equation of FUN, which figures out how much C to spend on uptake to remain at the target CN ratio. 
          ! nb. This term assumes that cost of N is constant through the timestep, because we don't have a concept of Michealis Menten kinetics. 
          ! 
          ! Calculate the hypothetical amount of NPP that we should use to extract N over whole profile
          ! This calculation is based on the simulataneous solution of the uptake and extrction N balance. 
          ! It satisfies the criteria (spentC+growthC=availC AND spentC/cost=growthC/plantCN
          ! Had to add growth respiration here to balance carbon pool. 
          if (total_N_resistance <= 1.0e-10_r8) then
                 write(iulog,*) 'ERROR: total_N_resistance=', total_N_resistance
                 do ipath = 1, npaths
                 do j = 1,nlevdecomp
                 write(iulog,*) 'n_from_paths:', n_from_paths(p,j,ipath), npp_to_paths(p,j,ipath)
                 enddo
                 end do
                 call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, &
                             msg= errMsg(sourcefile,  __LINE__))
              end if
          dnpp  = npp_to_spend / ( (1.0_r8+grperc(ivt(p)))*(plantCN(p) / total_N_resistance) + 1._r8)  
          dnpp  = dnpp * frac_ideal_C_use 


          !hypothetical amount of N acquired. 
          dn    = dnpp / total_N_resistance !
          do j = 1,nlevdecomp
                 
           ! RF How much of this NPP carbon do we allocate to the different pathways? fraction x gC/m2/s?
           ! Could this code now be put in a matrix? 

           npp_to_paths(p,j,ipecm:ipnmnh4) = npp_frac_paths(p,j,ipecm:ipnmnh4) * dnpp
                                
           if(FIX==plants_are_fixing)then
            npp_to_paths(p,j,ipfix) = npp_frac_paths(p,j,ipfix) * dnpp
           else
            npp_to_paths(p,j,ipfix) = 0.0_r8
           end if    

           n_from_paths(p,j,ipecm:ipnmnh4) = npp_to_paths(p,j,ipecm:ipnmnh4) / costs_paths(p,j,ipecm:ipnmnh4)
                           
           if(FIX==plants_are_fixing)then
            n_from_paths(p,j,ipfix) = npp_to_paths(p,j,ipfix) / costs_paths(p,j,ipfix)
           else
            n_from_paths(p,j,ipfix) = 0.0_r8
           end if                                                   
          end do

         ! Check if LIMITS of pools were exceeded:
          layer_loop: do j = 1,nlevdecomp
             N_before_corr(1:npaths) = n_from_paths(p,j,1:npaths)
             C_before_corr(1:npaths) = npp_to_paths(p,j,1:npaths)
             sminnh4_extracted = 0.0_r8
             sminno3_extracted = 0.0_r8
             ! get mining fluxes:
             ! ecm no3 and mining
             call myc_n_extraction(dzsoi_decomp(j),(smin_no3_to_plant_vr(c,j)) * dt, &
                  decomp_cpools_vr(c,j,i_ecm_myc),decomp_cpools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_phys_som), &
                  decomp_cpools_vr(c,j,i_chem_som),decomp_npools_vr(c,j,i_chem_som), sminno3_to_ecm_vr_patch(p,j), & 
                  c_somp2soma_vr_patch(p,j),n_somp2ecm_vr_patch(p,j), &
                  c_somc2soma_vr_patch(p,j),n_somc2ecm_vr_patch(p,j))
             ! ecm nh4
             call myc_n_extraction(dzsoi_decomp(j),(smin_nh4_to_plant_vr(c,j)) * dt, &
                  decomp_cpools_vr(c,j,i_ecm_myc),decomp_cpools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_phys_som), &
                  decomp_cpools_vr(c,j,i_chem_som),decomp_npools_vr(c,j,i_chem_som), sminnh4_to_ecm_vr_patch(p,j))
             ! am no3
             call myc_n_extraction(dzsoi_decomp(j),(smin_no3_to_plant_vr(c,j)) * dt, &
                  decomp_cpools_vr(c,j,i_am_myc),decomp_cpools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_phys_som), &
                  decomp_cpools_vr(c,j,i_chem_som),decomp_npools_vr(c,j,i_chem_som), sminno3_to_am_vr_patch(p,j))
             ! am nh4
             call myc_n_extraction(dzsoi_decomp(j),(smin_nh4_to_plant_vr(c,j)) * dt, &
                  decomp_cpools_vr(c,j,i_am_myc),decomp_cpools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_phys_som), &
                  decomp_cpools_vr(c,j,i_chem_som),decomp_npools_vr(c,j,i_chem_som), sminnh4_to_am_vr_patch(p,j))


             ! nonmyc gives to plant all the mineral nitrogen it extracts.
             sminnh4_extracted = sminnh4_to_ecm_vr_patch(p,j) + sminnh4_to_am_vr_patch(p,j) + n_from_paths(p,j,ipnmnh4)
             sminno3_extracted = sminno3_to_ecm_vr_patch(p,j) + sminno3_to_am_vr_patch(p,j) + n_from_paths(p,j,ipnmno3)
             ! Limiting nh4 extraction.
             if (sminnh4_to_paths(p,j) > 0.0_r8 .and. sminnh4_extracted > 0.0_r8) then
                ! extracted more than  there is possibly available
                if (sminnh4_to_paths(p,j) < sminnh4_extracted) then
                   sminnh4_overlimit = sminnh4_extracted - sminnh4_to_paths(p,j)
                   sminnh4_to_ecm_vr_patch(p,j) = sminnh4_to_ecm_vr_patch(p,j) * (1.0_r8 - sminnh4_overlimit / sminnh4_extracted)
                   sminnh4_to_am_vr_patch(p,j)  = sminnh4_to_am_vr_patch(p,j) * (1.0_r8 - sminnh4_overlimit / sminnh4_extracted)
                   n_from_paths(p,j,ipnmnh4)    = n_from_paths(p,j,ipnmnh4) * (1.0_r8 - sminnh4_overlimit / sminnh4_extracted)
                endif

             else
                sminnh4_overlimit = 0.0_r8
                sminnh4_to_ecm_vr_patch(p,j) = 0.0_r8
                sminnh4_to_am_vr_patch(p,j)  = 0.0_r8
                n_from_paths(p,j,ipnmnh4)    = 0.0_r8
             endif
             ! Limiting no3 extraction.
             if (sminno3_to_paths(p,j) > 0.0 .and. sminno3_extracted > 0.0) then
                ! extracted more than  there is possibly available
                if (sminno3_to_paths(p,j) < sminno3_extracted) then
                   sminno3_overlimit = sminno3_extracted - sminno3_to_paths(p,j)
                   sminno3_to_ecm_vr_patch(p,j) = sminno3_to_ecm_vr_patch(p,j) * (1.0_r8 - sminno3_overlimit / sminno3_extracted)
                   sminno3_to_am_vr_patch(p,j)  = sminno3_to_am_vr_patch(p,j) * (1.0_r8 - sminno3_overlimit / sminno3_extracted)
                   n_from_paths(p,j,ipnmno3)    = n_from_paths(p,j,ipnmno3) * (1.0_r8 - sminno3_overlimit / sminno3_extracted)
                endif
             else
                sminno3_overlimit            = 0.0_r8
                sminno3_to_ecm_vr_patch(p,j) = 0.0_r8
                sminno3_to_am_vr_patch(p,j)  = 0.0_r8
                n_from_paths(p,j,ipnmno3)    = 0.0_r8
             endif


             ! Calculate actual myc fluxes now:
             ! ecm
             call myc_cn_fluxes(dzsoi_decomp(j), npp_to_paths(p,j,ipecm), sminno3_to_ecm_vr_patch(p,j) + & 
                  n_somc2ecm_vr_patch(p,j) + n_somc2ecm_vr_patch(p,j), &
                  n_from_paths(p,j,ipecm), n_ecm_growth_vr_patch(p,j), &
                  c_ecm_growth_vr_patch(p,j), c_ecm_resp_vr_patch(p,j), ecm_step, c_ecm_enz_vr_patch(p,j))

             call myc_cn_fluxes(dzsoi_decomp(j), npp_to_paths(p,j,ipam), sminno3_to_am_vr_patch(p,j), &
                  n_from_paths(p,j,ipam), n_am_growth_vr_patch(p,j), &
                  c_am_growth_vr_patch(p,j), c_am_resp_vr_patch(p,j), am_step)

          ! MVD Here should the carbon and nitrogen fluxes to the plant be summed up       
            
          end do layer_loop
            !if (carbong_spent(p)>availc(p)) then

             ! MVD here we check for availible carbon from plant and plant n demand

            !if (n_toplant>ndemand(p)) then
          do j = 1,nlevdecomp
              N_acquired     =  sum(n_from_paths(p,j,ipecm:ipnmnh4))        ! How much N did we end up with
                                
              C_spent        =   sum(npp_to_paths(p,j,ipecm:ipnmnh4))       ! How much did it actually cost? 


                         
            if(FIX==plants_are_fixing)then
               N_acquired  = N_acquired + n_from_paths(p,j,ipfix)
               C_spent     = C_spent + npp_to_paths(p,j,ipfix) 
            end if
            
            ! How much C did we allocate or spend in this layer? 
            npp_to_spend    = npp_to_spend   - C_spent - (N_acquired  * plantCN(p)*(1.0_r8+ grperc(ivt(p))))
            
            ! Accumulate those fluxes
            nt_uptake(p)             = nt_uptake(p)       + N_acquired
            npp_uptake(p)            = npp_uptake(p)      + C_spent


            !-------------------- N flux accumulation------------!
            do ipath = ipecm,ipnmnh4
               n_paths_acc(p,ipath) = n_paths_acc(p,ipath) + n_from_paths(p,j,ipath)
            end do
            !-------------------- C flux accumulation------------!
            do ipath = ipecm,ipnmnh4
               npp_paths_acc(p,ipath) = npp_paths_acc(p,ipath) + npp_to_paths(p,j,ipath)
            end do
               
            if(FIX == plants_are_fixing)then
               n_paths_acc(p,ipfix)     = n_paths_acc(p,ipfix)     + n_from_paths(p,j,ipfix)
               npp_paths_acc(p,ipfix)   = npp_paths_acc(p,ipfix)   + npp_to_paths(p,j,ipfix) 
            end if
       
          end do

         end if !unmet demand`


         npp_to_spend_fix = sum(npp_paths_acc(p,1:ipfix)) +  total_c_spent_retrans + total_c_accounted_retrans 
         npp_spent = sum(npp_paths_acc(p,1:ipfix)) +  total_c_spent_retrans + total_c_accounted_retrans
         
           if(npp_spent - npp_to_spend_init > 1.0e-10_r8) then
            write(iulog,*) 'ERROR: TO MUCH CARBON HAS BEEN SPENT ON N UPTAKE'
            write(iulog,*) 'npp_to_spend_fix, npp before', npp_spent, npp_to_spend_init
            write(iulog,*) 'npp retrans:', total_c_spent_retrans, total_c_accounted_retrans
            do ipath = ipecm,ipfix
                write(iulog,*) 'npp to path:', ipath, npp_paths_acc(p,ipath)
            end do
        end if

      end do fix_loop ! FIXER. 
             ! Turn step level quantities back into fluxes per second. 
             Nfix(p)                   = (n_paths_acc(p,ipfix)) / dt                   
             retransn_to_npool(p)      = (n_retrans_acc(p)) / dt 
             ! Without matrix solution
             if(.not. use_matrixcn)then
                free_retransn_to_npool(p) = (free_nretrans_acc(p)) / dt
             ! With matrix solution (when it comes in)
             end if
             ! this is the N that comes off leaves. 
             Nretrans(p)               = retransn_to_npool(p) + free_retransn_to_npool(p)
    
             
             !Extract active uptake N from soil pools. 
             do j = 1, nlevdecomp
             !RF change. The N fixed doesn't actually come out of the soil mineral pools, it is 'new'... 
             ! we have  these variable instead of smin*_to_plant_fun_vr, just need give them proper units
             n_ecm(p,j)  =  n_ecm(p,j)/(dzsoi_decomp(j)*dt)
             n_am(p,j)  =  n_am(p,j)/(dzsoi_decomp(j)*dt)
             n_nonmyc_no3_vr(p,j)  =  n_nonmyc_no3_vr(p,j)/(dzsoi_decomp(j)*dt)
             n_nonmyc_nh4_vr(p,j)  =  n_nonmyc_nh4_vr(p,j)/(dzsoi_decomp(j)*dt)
             end do
    
             !SPLIT TO NO3 and NH4 like in original fun
             !---------------------------N fluxes--------------------! ! total means sum(var(p,1:nmyc,ipath)),ecm/am var(p,ecm_step,ipath)
             Nactive_no3(p) = 0.0_r8
             Nactive_nh4(p) = 0.0_r8
             Nnonmyc_no3(p) = (n_paths_acc(p,ipnmno3)) / dt
             Nnonmyc_nh4(p) = (n_paths_acc(p,ipnmnh4)) / dt
             Necm(p) = (n_paths_acc(p,ipecm)) / dt
             Nam(p) = (n_paths_acc(p,ipam)) / dt
             Nnonmyc(p) = Nnonmyc_no3(p) + Nnonmyc_nh4(p)
            
             plant_ndemand_retrans(p)  = plant_ndemand_retrans(p)/dt

             Nuptake(p) = Necm(p) + Nam(p) + Nnonmyc_no3(p) + Nnonmyc_nh4(p) + Nfix(p)+ &
                                  retransn_to_npool(p)+free_retransn_to_npool(p) 

            ! Nactive(p) = Nactive_no3(p)  + Nactive_nh4(p) + Nnonmyc_no3(p) + Nnonmyc_nh4(p)
             Nactive(p) = Necm(p)  +   Nam(p) + Nnonmyc_no3(p) + Nnonmyc_nh4(p)
             if (Nuptake(p) > 10000._r8) then
              !write(iulog,*) 'ERROR: Nactive_no3 negative: ', Nactive_no3(p)
              !write(iulog,*) 'ERROR: Nactive_nh4 negative: ', Nactive_nh4(p)
              write(iulog,*) 'ERROR: Nnonmyc_no3 negative: ', Nnonmyc_no3(p)
              write(iulog,*) 'ERROR: Nnonmyc_nh4 negative: ', Nnonmyc_nh4(p)
              write(iulog,*) 'ERROR: Nfix negative: ', Nfix(p)
              write(iulog,*) 'ERROR: retransn_to_npool negative: ', retransn_to_npool(p)
              write(iulog,*) 'ERROR: free_retransn_to_npool negative: ', free_retransn_to_npool(p)

                 call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, &
                             msg= errMsg(sourcefile,  __LINE__))
             endif                                          
             ! free N goes straight to the npool, not throught Nuptake...
             sminn_to_plant_fun(p)     = Nactive_no3(p) + Nactive_nh4(p) + Nnonmyc_no3(p) + Nnonmyc_nh4(p) + Nfix(p) 
    
             !---------------------------C fluxes--------------------!

             npp_Nactive_no3(p) = 0.0_r8
             npp_Nactive_nh4(p) = 0.0_r8
             npp_Necm(p) = (npp_paths_acc(p,ipecm)) / dt
             npp_Nam(p) = (npp_paths_acc(p,ipam)) / dt
             npp_Nnonmyc_no3(p) = (npp_paths_acc(p,ipnmno3)) / dt
             npp_Nnonmyc_nh4(p) = (npp_paths_acc(p,ipnmnh4)) / dt
             npp_Nnonmyc(p) = npp_Nnonmyc_no3(p) + npp_Nnonmyc_nh4(p)
             npp_Nfix(p) = (npp_paths_acc(p,ipfix)) /dt
             npp_Nretrans(p) = (npp_retrans_acc(p))/dt
             npp_Nactive(p) = npp_Necm(p) + npp_Nam(p) + npp_Nnonmyc_no3(p) + npp_Nnonmyc_nh4(p)

             !---------------------------Extra Respiration Fluxes--------------------!      
             soilc_change(p)           = npp_Nactive(p) + npp_Nfix(p) + npp_Nretrans(p)
            ! soilc_change(p)           = npp_Nactive(p) + npp_Nfix(p) + npp_Nnonmyc(p) + npp_Nretrans(p)
             soilc_change(p)           = npp_Necm(p) + npp_Nam(p) + npp_Nnonmyc_no3(p) + npp_Nnonmyc_nh4(p) + npp_Nfix(p) + npp_Nretrans(p)
             soilc_change(p)           = soilc_change(p) + burned_off_carbon / dt                 
             npp_burnedoff(p)          = burned_off_carbon/dt          
             npp_Nuptake(p)            = soilc_change(p)
             ! how much carbon goes to growth of tissues?  
             npp_growth(p)             = (Nuptake(p)- free_retransn_to_npool(p))*plantCN(p)+(excess_carbon_acc/dt) !does not include gresp, since this is calculated from growth 
             if (availc(p) <= 0.0_r8 .and. soilc_change(p) > 0.0_r8) then
              write(iulog,*) 'ERROR: availc(p): ', availc(p)
              write(iulog,*) 'ERROR: soilc_change(p): ', soilc_change(p)
              write(iulog,*) 'ERROR: free_retransn_to_npool(p) is negative: ', free_retransn_to_npool(p)
              write(iulog,*) 'ERROR:burned_off_carbon / dt: ',burned_off_carbon / dt
              write(iulog,*) 'ERROR: excess_carbon_acc: ', (npp_uptake(p))
              write(iulog,*) 'ERROR: npp_Nretrans(p): ', npp_Nretrans(p)
              write(iulog,*) 'npp_growth(p): ',npp_growth(p)
             endif

             if (npp_growth(p) < -1.0e-7_r8 .or. npp_growth(p) > 10000._r8) then
              write(iulog,*) 'ERROR: Nuptake(p) is negative: ', Nuptake(p)
              write(iulog,*) 'ERROR: npp_Nuptake(p) is negative: ', npp_Nuptake(p)
              write(iulog,*) 'ERROR: free_retransn_to_npool(p) is negative: ', free_retransn_to_npool(p)
              write(iulog,*) 'ERROR: plantCN(p) is negative: ', plantCN(p)
              write(iulog,*) 'ERROR: excess_carbon_acc is negative: ', excess_carbon_acc
              write(iulog,*) 'ERROR: npp_Nactive(p) is negative: ', npp_Nactive(p)
              write(iulog,*) 'ERROR: npp_Nnonmyc(p) is npp_Nnonmyc: ', npp_Nnonmyc(p)
                 call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, &
                             msg= errMsg(sourcefile,  __LINE__))
             endif
             if (availc(p) - npp_Nuptake(p) - npp_growth(p) < -1.e-8_r8) then
              write(iulog,*) 'ERROR: balance Cfun is negative: '
              write(iulog,*) 'Acailc, npp_Nuptake/growth:',availc(p), npp_Nuptake(p),npp_growth(p)
              write(iulog,*)  'soilchange, burned off c', soilc_change(p), burned_off_carbon/dt
              write(iulog,*), 'Excess carbon, npp_Nretrans,freeretrans',excess_carbon_acc/dt,npp_Nretrans(p),free_retransn_to_npool(p)
             endif
             !-----------------------Diagnostic Fluxes------------------------------!
             if(availc(p).gt.0.0_r8)then !what happens in the night? 
                nuptake_npp_fraction_patch(p) = npp_Nuptake(p)/availc(p)
             else
                nuptake_npp_fraction_patch(p) = spval
             endif 
             if(npp_Nfix(p).gt.0.0_r8)then
                cost_nfix(p) = Nfix(p)/npp_Nfix(p)
             else
                cost_nfix(p) = spval
             endif 
             if((npp_Nactive(p)+npp_Nnonmyc(p)).gt.0.0_r8)then
                cost_nactive(p) = (Nactive(p)+Nnonmyc(p))/(npp_Nactive(p)+npp_Nnonmyc(p))
             else
                cost_nactive(p) = spval
             endif 
             if(npp_Nretrans(p).gt.0.0_r8)then
                cost_nretrans(p) = Nretrans(p)/npp_Nretrans(p)
             else
                cost_nretrans(p) = spval
             endif 
            
end do pft 
! go through all variables avove and see if they are used in a nothe rmodule from 1231
! patch variables have to be p2c into col variables
! compilable 
! is everything in CN fun,  0red, allocated
end associate
call t_startf( 'updateCNFUNMIMICSplus' )
call updateCNFUNMIMICSplus (bounds, num_soilc, filter_soilc, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
                           soilbiogeochem_nitrogenflux_inst, soilbiogeochem_carbonflux_inst, cnfunmimicsplus_inst, &
                           c_am_resp_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), c_ecm_resp_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), &
                           c_am_growth_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), c_ecm_growth_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), &
                           n_am_growth_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), n_ecm_growth_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), &
                           c_ecm_enz_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), n_somc2ecm_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), &
                           n_somp2ecm_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), &
                           c_somc2soma_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), c_somp2soma_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp),&
                           sminno3_to_ecm_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), sminno3_to_am_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), &
                           sminnh4_to_ecm_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), sminnh4_to_am_vr_patch(bounds%begp:bounds%endp, 1:nlevdecomp), &
                           n_ecm(bounds%begp:bounds%endp, 1:nlevdecomp), n_nonmyc_no3_vr(bounds%begp:bounds%endp, 1:nlevdecomp), &
                           n_am(bounds%begp:bounds%endp, 1:nlevdecomp), n_nonmyc_nh4_vr(bounds%begp:bounds%endp, 1:nlevdecomp) &
                           )

call t_stopf( 'updateCNFUNMIMICSplus' )
end subroutine CNFUNMIMICSplus


subroutine updateCNFUNMIMICSplus (bounds, num_soilc, filter_soilc, &
                                 cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
                                 soilbiogeochem_nitrogenflux_inst, soilbiogeochem_carbonflux_inst, cnfunmimicsplus_inst, &
                                 c_am_resp_vr_patch, c_ecm_resp_vr_patch, c_am_growth_vr_patch, c_ecm_growth_vr_patch, &
                                 n_am_growth_vr_patch, n_ecm_growth_vr_patch, c_ecm_enz_vr_patch, &
                                 n_somc2ecm_vr_patch, n_somp2ecm_vr_patch, c_somc2soma_vr_patch, c_somp2soma_vr_patch,&
                                 sminno3_to_ecm_vr_patch, sminno3_to_am_vr_patch, sminnh4_to_ecm_vr_patch, sminnh4_to_am_vr_patch, &
                                 n_ecm, n_nonmyc_no3_vr, n_am, n_nonmyc_nh4_vr &
                                 )
    !
    ! !DESCRIPTION:

   ! !USES
   use subgridAveMod   , only : p2c
   use clm_varctl      , only : use_nitrif_denitrif

   ! !ARGUMENTS
   type(bounds_type)                       , intent(in)    :: bounds
   integer                                 , intent(in)    :: num_soilc             ! number of soil columns in filter
   integer                                 , intent(in)    :: filter_soilc(:)       ! filter for soil columns
   type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
   type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
   type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
   type(cnfunmimicsplus_type)              , intent(inout) :: cnfunmimicsplus_inst
   ! mycorrhiza patch fluxes
   real(r8), intent(in) :: c_am_resp_vr_patch(:,:)          ! carbon respiration flux for AM mycorrhiza
   real(r8), intent(in) :: c_ecm_resp_vr_patch(:,:)         ! carbon respiration flux for ECM mycorrhiza
   real(r8), intent(in) :: c_am_growth_vr_patch(:,:)        ! carbon growth flux for AM mycorrhiza
   real(r8), intent(in) :: c_ecm_growth_vr_patch(:,:)       ! carbon growth flux for ECM mycorrhiza
   real(r8), intent(in) :: n_am_growth_vr_patch(:,:)        ! nitrogen growth flux for AM mycorrhiza
   real(r8), intent(in) :: n_ecm_growth_vr_patch(:,:)       ! nitrogen growth flux for ECM mycorrhiza
   real(r8), intent(in) :: c_ecm_enz_vr_patch(:,:)          ! carbon enzyme production flux for ECM mycorrhiza
   real(r8), intent(in) :: n_somc2ecm_vr_patch(:,:)         ! nitrogen mining from ECM mycorrhiza
   real(r8), intent(in) :: n_somp2ecm_vr_patch(:,:)         ! nitrogen mining from ECM mycorrhiza
   real(r8), intent(in) :: c_somc2soma_vr_patch(:,:)        ! carbon release from mining from somc pool
   real(r8), intent(in) :: c_somp2soma_vr_patch(:,:)        ! carbon release from mining from somp pool
   real(r8), intent(in) :: sminno3_to_ecm_vr_patch(:,:)     ! No3 flux from soil to ECM
   real(r8), intent(in) :: sminno3_to_am_vr_patch(:,:)      ! No3 flux from soil to AM
   real(r8), intent(in) :: sminnh4_to_ecm_vr_patch(:,:)     ! NH4 flux from soil to ECM
   real(r8), intent(in) :: sminnh4_to_am_vr_patch(:,:)      ! NH4 flux from soil to AM
   real(r8), intent(in) :: n_ecm(:,:)                       ! Layer EcM mycorrhizal N uptake (gN/m2)
   real(r8), intent(in) :: n_nonmyc_no3_vr(:,:)             ! Layer non-myc no3 uptake (gN/m2)
   real(r8), intent(in) :: n_am(:,:)                        ! Layer AM mycorrhizal N uptake (gN/m2)
   real(r8), intent(in) :: n_nonmyc_nh4_vr(:,:)             ! Layer non-myc nh4 uptake (gN/m2)

   associate( cnfun => cnfunmimicsplus_inst              , &
              scf   => soilbiogeochem_carbonflux_inst    , &
              snf   => soilbiogeochem_nitrogenflux_inst  , &
              cf    => cnveg_carbonflux_inst             , &
              nf    => cnveg_nitrogenflux_inst             )
       

   !!! soilc_change_col is not used anywhere
   call p2c(bounds, num_soilc, filter_soilc, cf%soilc_change_patch(bounds%begp:bounds%endp), &
             scf%soilc_change_col(bounds%begc:bounds%endc))
             
   call p2c(bounds, num_soilc, filter_soilc, nf%Nfix_patch(bounds%begp:bounds%endp), &
             snf%nfix_to_sminn_col(bounds%begc:bounds%endc))


   if (use_nitrif_denitrif) then
      ! plant fluxes
      call p2c(bounds,nlevdecomp, &
               n_ecm(bounds%begp:bounds%endp,1:nlevdecomp),&
               cnfun%no3_myc_to_plant_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')
      call p2c(bounds,nlevdecomp, &
               n_nonmyc_no3_vr(bounds%begp:bounds%endp,1:nlevdecomp),&
               snf%sminno3_nonmyc_to_plant_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')
      call p2c(bounds,nlevdecomp, &
               n_am(bounds%begp:bounds%endp,1:nlevdecomp),&
               cnfun%nh4_myc_to_plant_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')
      call p2c(bounds,nlevdecomp, &
               n_nonmyc_nh4_vr(bounds%begp:bounds%endp,1:nlevdecomp),&
               snf%sminnh4_nonmyc_to_plant_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')

      ! mycorrhyza fluxes
      call p2c(bounds,nlevdecomp, &
               sminno3_to_ecm_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
               snf%sminno3_to_ecm_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')
      call p2c(bounds,nlevdecomp, &
               sminnh4_to_ecm_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
               snf%sminnh4_to_ecm_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')
      call p2c(bounds,nlevdecomp, &
               sminno3_to_am_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
               snf%sminno3_to_am_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')
      call p2c(bounds,nlevdecomp, &
               sminnh4_to_am_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
               snf%sminnh4_to_am_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')

   else

   endif

   ! fluxes that will be used in decomposition:
   call p2c(bounds,nlevdecomp, &
   c_am_resp_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   scf%c_am_resp_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   call p2c(bounds,nlevdecomp, &
   c_ecm_resp_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   scf%c_ecm_resp_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   call p2c(bounds,nlevdecomp, &
   c_am_growth_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   scf%c_am_growth_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   call p2c(bounds,nlevdecomp, &
   c_ecm_growth_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   scf%c_ecm_growth_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   ! these nitrogen fluxes might not bee needed since C:N ratio for growth is constant
   call p2c(bounds,nlevdecomp, &
   n_am_growth_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   snf%n_am_growth_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   call p2c(bounds,nlevdecomp, &
   n_ecm_growth_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   snf%n_ecm_growth_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   call p2c(bounds,nlevdecomp, &
   c_ecm_enz_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   scf%c_ecm_enz_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   ! mining fluxes
   call p2c(bounds,nlevdecomp, &
   n_somc2ecm_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   snf%n_somc2ecm_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')
   call p2c(bounds,nlevdecomp, &
   n_somp2ecm_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   snf%n_somp2ecm_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   call p2c(bounds,nlevdecomp, &
   c_somc2soma_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   scf%c_somc2soma_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   call p2c(bounds,nlevdecomp, &
   c_somp2soma_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
   scf%c_somp2soma_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
   'unity')

   end associate
end subroutine updateCNFUNMIMICSplus

!=========================================================================================
real(r8) function fun_cost_fix(fixer,a_fix,b_fix,c_fix,big_cost,crootfr,s_fix, tc_soisno)

! Description:
!   Calculate the cost of fixing N by nodules.
! Code Description:
!   This code is written to CLM4CN by Mingjie Shi on 06/27/2013
use clm_varcon      , only : tfrz

  implicit none
!--------------------------------------------------------------------------
! Function result.
!--------------------------------------------------------------------------
! real(r8) , intent(out) :: cost_of_n   !!! cost of fixing N (kgC/kgN)
!--------------------------------------------------------------------------
! Scalar arguments with intent(in).
!--------------------------------------------------------------------------
  integer,  intent(in) :: fixer     ! flag indicating if plant is a fixer
                                    ! 1=yes, otherwise no.
  real(r8), intent(in) :: a_fix     ! As in Houlton et al. (Nature) 2008
  real(r8), intent(in) :: b_fix     ! As in Houlton et al. (Nature) 2008
  real(r8), intent(in) :: c_fix     ! As in Houlton et al. (Nature) 2008
  real(r8), intent(in) :: big_cost  ! an arbitrary large cost (gC/gN)
  real(r8), intent(in) :: crootfr   ! fraction of roots for carbon that are in this layer
  real(r8), intent(in) :: s_fix     ! Inverts Houlton et al. 2008 and constrains between 7.5 and 12.5
  real(r8), intent(in) :: tc_soisno ! soil temperature (K)
  ! LOCAL VARIABIABLES:
  real(r8) :: t_soil_c  ! temperragture in celcious


  t_soil_c = tc_soisno - tfrz
  if (fixer == 1 .and. crootfr > 1.e-6_r8) then
     ! New term to directly account for Ben Houlton's temperature response function.
     ! Assumes s_fix is -6.  (RF, Jan 2015)  
     ! 1.25 converts from the Houlton temp response function to a 0-1 limitation factor. 
     ! The cost of N should probably be 6 gC/gN (or 9, including maintenance costs of nodules) 
     ! for 'optimal' temperatures. This cost should increase in a way that mirrors 
     ! Houlton et al's observations of temperautre limitations on the mirboial fixation rates. 
     ! We don't actually simulate the rate of fixation (and assume that N uptake is instantaneous) 
     ! here, so instead the limitation term is here rolled into the cost function.  
     
     ! Here we invert the 'cost' to give the optimal N:C ratio (1/6 gN/gC)  The amount of N 
     ! you get for a given C goes down as it gets colder, so this can be multiplied by 
     ! the temperature function to give a temperature-limited N:C of  f/6. This number 
     ! can then be inverted to give a temperature limited C:N, as 1/(f/6). Which is the 
     ! same as 6/f, given here" 
     fun_cost_fix  = (-1*s_fix) * 1.0_r8 / (1.25_r8* (exp(a_fix + b_fix * t_soil_c * (1._r8 - 0.5_r8 * t_soil_c / c_fix)) ))
  else
     fun_cost_fix = big_cost
  end if    ! ends up with the fixer or non-fixer decision
  
  end function fun_cost_fix
!=========================================================================================

!=========================================================================================
  real(r8) function fun_cost_nonmyc(sminn_layer,big_cost,kc_nonmyc,kn_nonmyc,rootc_dens,crootfr,smallValue)         

! Description:
!    Calculate the cost of nonmyc uptake of N frm the soil.
! Code Description:
!   This code is written to CLM4 by Mingjie Shi.

  implicit none
!--------------------------------------------------------------------------
! Function result.
!--------------------------------------------------------------------------
  real(r8), intent(in) :: sminn_layer   !  Amount of N (as NH4 or NO3) in the soil that is available to plants (gN/m2).
  real(r8), intent(in) :: big_cost      !  An arbitrary large cost (gC/gN).
  real(r8), intent(in) :: kc_nonmyc     !  Constant for cost of nonmyc uptake (gC/m2).
  real(r8), intent(in) :: kn_nonmyc     !  Constant for cost of nonmyc uptake (gC/m2).
  real(r8), intent(in) :: rootc_dens   !  Root carbon density in layer (gC/m3).
  real(r8), intent(in) :: crootfr        !  Fraction of roots that are in this layer.
  real(r8), intent(in) :: smallValue    !  A small number.

  if (rootc_dens > 1.e-6_r8.and.sminn_layer > smallValue) then
    fun_cost_nonmyc =  kn_nonmyc / sminn_layer + kc_nonmyc / rootc_dens 
  else
!   There are very few roots in this layer. Set a high cost.
    fun_cost_nonmyc = big_cost
  end if

  end function fun_cost_nonmyc

!==========================================================================

 subroutine fun_retranslocation(p,dt,npp_to_spend,total_falling_leaf_c,         &
               total_falling_leaf_n, total_n_resistance, total_c_spent_retrans, &
               total_c_accounted_retrans, free_n_retrans, paid_for_n_retrans,   &
               target_leafcn, grperc, plantCN)
!
! Description:
! This subroutine (should it be a function?) calculates the amount of N absorbed and C spent 
! during retranslocation. 
! Rosie Fisher. April 2016. 
! !USES:
  implicit none 

! !ARGUMENTS:
  real(r8), intent(IN) :: total_falling_leaf_c  ! INPUT  gC/m2/timestep
  real(r8), intent(IN) :: total_falling_leaf_n  ! INPUT  gC/m2/timestep
  real(r8), intent(IN) :: total_n_resistance    ! INPUT  gC/gN
  real(r8), intent(IN) :: npp_to_spend          ! INPUT  gN/m2/timestep
  real(r8), intent(IN) :: target_leafcn         ! INPUT  gC/gN
  real(r8), intent(IN) :: dt                    ! INPUT  seconds
  real(r8), intent(IN) :: grperc                ! INPUT growth respiration fraction
  real(r8), intent(IN) :: plantCN               ! INPUT plant CN ratio
  integer, intent(IN)  :: p                     ! INPUT  patch index

  real(r8), intent(OUT) :: total_c_spent_retrans     ! OUTPUT gC/m2/timestep
  real(r8), intent(OUT) :: total_c_accounted_retrans ! OUTPUT gC/m2/timestep
  real(r8), intent(OUT) :: paid_for_n_retrans        ! OUTPUT gN/m2/timestep
  real(r8), intent(OUT) :: free_n_retrans            ! OUTPUT gN/m2/timestep

  !
  ! !LOCAL VARIABLES:
  real(r8) :: kresorb               ! INTERNAL used factor
  real(r8) :: falling_leaf_c        ! INTERNAL gC/m2/timestep
  real(r8) :: falling_leaf_n        ! INTERNAL gN/m2/timestep
  real(r8) :: falling_leaf_cn       ! INTERNAL gC/gN
  real(r8) :: cost_retrans_temp     ! INTERNAL gC/gN
  real(r8) :: leaf_n_ext            ! INTERNAL gN/m2/timestep
  real(r8) :: c_spent_retrans       ! INTERNAL gC/m2/timestep
  real(r8) :: c_accounted_retrans   ! INTERNAL gC/m2/timestep
  real(r8) :: npp_to_spend_temp     ! INTERNAL gC/m2/timestep
  real(r8) :: max_falling_leaf_cn   ! INTERNAL gC/gN
  real(r8) :: min_falling_leaf_cn   ! INTERNAL gC/gN
  real(r8) :: cost_escalation       ! INTERNAL cost function parameter
  integer  :: iter                  ! INTERNAL
  integer  :: exitloop              ! INTERNAL
  ! ------------------------------------------------------------------------------- 


   ! ------------------ Initialize total fluxes. ------------------!
   total_c_spent_retrans = 0.0_r8
   total_c_accounted_retrans = 0.0_r8
   c_accounted_retrans   = 0.0_r8
   paid_for_n_retrans    = 0.0_r8
   npp_to_spend_temp     = npp_to_spend

   ! ------------------ Initial C and N pools in falling leaves. ------------------!
   falling_leaf_c       =  total_falling_leaf_c      
   falling_leaf_n       =  total_falling_leaf_n 

   !  ------------------ PARAMETERS ------------------ 
   max_falling_leaf_cn = target_leafcn * 3.0_r8 
   min_falling_leaf_cn = target_leafcn * 1.5_r8
   cost_escalation     = 1.3_r8

   !  ------------------ Free uptake ------------------ 
   free_n_retrans  = max(falling_leaf_n -  (falling_leaf_c/min_falling_leaf_cn),0.0_r8)
   falling_leaf_n = falling_leaf_n -  free_n_retrans  

   ! ------------------ Initial CN ratio and costs ------------------!  
   falling_leaf_cn      = falling_leaf_c/falling_leaf_n 
   kresorb =  (1.0_r8/target_leafcn)
   cost_retrans_temp    = kresorb / ((1.0_r8/falling_leaf_cn )**1.3_r8)

   ! ------------------ Iteration loops to figure out extraction limit ------------!
   iter = 0
   exitloop = 0
   do while(exitloop==0.and.cost_retrans_temp .lt. total_n_resistance.and. &
            falling_leaf_n.ge.0.0_r8.and.npp_to_spend.gt.0.0_r8)
      ! ------------------ Spend some C on removing N ------------!
      ! spend enough C to increase leaf C/N by 1 unit. 
      c_spent_retrans   = cost_retrans_temp * (falling_leaf_n - falling_leaf_c / &
                          (falling_leaf_cn + 1.0_r8))
      ! don't spend more C than you have  
      c_spent_retrans   = min(npp_to_spend_temp, c_spent_retrans) 
      ! N extracted, per this amount of C expenditure
      leaf_n_ext        = c_spent_retrans / cost_retrans_temp     
      ! Do not empty N pool 
      leaf_n_ext        = min(falling_leaf_n, leaf_n_ext)    
      !How much C do you need to account for the N that got taken up? 
      c_accounted_retrans = leaf_n_ext * plantCN * (1.0_r8 + grperc)      

      ! ------------------ Update leafCN, recalculate costs ------------!
      falling_leaf_n    = falling_leaf_n - leaf_n_ext          ! remove N from falling leaves pool 
      if(falling_leaf_n.gt.0.0_r8)then
         falling_leaf_cn   = falling_leaf_c/falling_leaf_n     ! C/N ratio
         cost_retrans_temp = kresorb /((1.0_r8/falling_leaf_cn)**1.3_r8) ! cost function. PARAMETER
      else
         exitloop=1
      endif 
 
      ! ------------------ Accumulate total fluxes ------------!
      total_c_spent_retrans     = total_c_spent_retrans + c_spent_retrans 
      total_c_accounted_retrans = total_c_accounted_retrans + c_accounted_retrans 
      paid_for_n_retrans    = paid_for_n_retrans    + leaf_n_ext
      npp_to_spend_temp     = npp_to_spend_temp     - c_spent_retrans  - c_accounted_retrans
      iter = iter+1
   
      ! run out of C or N
      if(npp_to_spend_temp.le.0.0_r8)then
         exitloop=1
         ! if we made a solving error on this (expenditure and n uptake should 
         ! really be solved simultaneously)
         ! then remove the error from the expenditure. This changes the notional cost, 
         ! but only by a bit and prevents cpool errors. 

         total_c_spent_retrans  = total_c_spent_retrans + npp_to_spend_temp 
      endif 
      ! leaf CN is too high
      if(falling_leaf_cn.ge.max_falling_leaf_cn)then
         exitloop=1
      endif
      ! safety check to prevent hanging code
      if(iter.ge.150)then
          exitloop=1
      endif 
   end do

 end subroutine fun_retranslocation

!==========================================================================


end module CNFUNMIMICSplusMod
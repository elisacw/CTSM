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
  use SoilBiogeochemNitrogenStateType  , only : soilbiogeochem_nitrogenstate_type
     
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemDecompCascadeConType, only : mimicsplus_decomp, decomp_method
  use WaterStateBulkType                  , only : waterstatebulk_type
  use WaterFluxBulkType                   , only : waterfluxbulk_type
  use TemperatureType                 , only : temperature_type
  use SoilStateType                   , only : soilstate_type
  use CanopyStateType                 , only : canopystate_type
  use perf_mod                        , only : t_startf, t_stopf
  use clm_varpar                      , only : nlevdecomp
  use clm_varcon                     , only : secspday
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
  private :: fun_cost_active
  private :: fun_cost_nonmyc
  private :: fun_cost_fix
  private :: fun_retranslocation

  ! DATA STRUCTURES


  integer,  parameter :: ipano3    = 1             ! Process number for mycorrhizal uptake of NO3.
  integer,  parameter :: ipanh4    = 2             ! Process number for mycorrhizal uptake of NH4
  integer,  parameter :: ipnmno3    = 3             ! Process number for nonmyc uptake of NO3.
  integer,  parameter :: ipnmnh4    = 4             ! Process number for nonmyc uptake of NH4.
  integer,  parameter :: ipfix      = 5             ! Process number for fixing.                          orginal version in FUN: [icost...]
  integer,  parameter :: ipret      = 6             ! Process number for retranslocation.
  
  real(r8), parameter :: big_cost        = 1000000000._r8! An arbitrarily large cost
 
  !  array index when plant is fixing
  integer, parameter :: plants_are_fixing = 1
  integer, parameter :: plants_not_fixing = 2
 
  !  array index for ECM step versus AM step
  integer, parameter :: ecm_step          = 1
  integer, parameter :: am_step           = 2
  !  arbitrary large cost (gC/gN).

  integer,  private, parameter :: nmyc            = 2             ! Number of calculation part          orginal version in FUN: [nstp]
  integer,  private, parameter :: npath6          = 6             ! Number of  N transport pathways     orginal version in FUN: [ncost6]
  
  !
  !  populated in readParamsMod
  
  type, public :: cnfunmimicsplus_type
  ! declare arrays and scalars here with pointers

  real(r8)           :: rootc_dens_step                                       ! root C for each PFT in each soil layer(gC/m2)
  real(r8)           :: dn                                                    ! Increment of N                        (gN/m2)  
  real(r8), pointer  :: rootc_dens(:,:)                                       ! the root carbon density               (gC/m2)
  real(r8), pointer  :: rootC(:)                                              ! root biomass                          (gC/m2)
  real(r8), pointer  :: n_uptake_myc_frac(:,:)                                ! the arrary for the ECM and AM ratio   (-)       orginal version in FUN: [permyc]
  real(r8), pointer  :: kc_active(:,:)                                        ! the kc_active parameter               (gC/m2)
  real(r8), pointer  :: kn_active(:,:)                                        ! the kn_active parameter               (gC/m2)
  real(r8), pointer :: plant_ndemand_pool(:)                                  ! The N demand pool (gN/m2)
  real(r8), pointer :: plant_ndemand_pool_step(:,:)                           ! the N demand pool (gN/m2)
  real(r8), pointer :: litterfall_n_step(:,:)                                 ! N loss based on the leafc to litter (gN/m2)
  real(r8), pointer :: litterfall_c_step(:,:)                                 ! N loss based on the leafc to litter (gN/m2)
  real(r8), pointer :: tc_soisno(:,:)                                         ! Soil temperature (degrees Celsius)
  real(r8), pointer :: npp_remaining(:,:)                                     ! A temporary variable for npp_remaining(gC/m2) 
  real(r8), pointer :: n_passive_step(:,:)                                    ! N taken up by transpiration at substep(gN/m2)
  real(r8), pointer :: n_passive_acc(:)                                       ! N acquired by passive uptake (gN/m2)
  real(r8), pointer :: cost_fix(:,:)                                          ! cost of fixation (gC/gN)
  real(r8), pointer :: n_fix_acc(:,:)                                         ! N acquired by fixation (gN/m2)
  real(r8), pointer :: n_fix_acc_total(:)                                     ! N acquired by fixation (gN/m2)
  real(r8), pointer :: npp_fix_acc(:,:)                                       ! Amount of NPP used by fixation (gC/m2)
  real(r8), pointer :: npp_fix_acc_total(:)                                   ! Amount of NPP used by fixation (gC/m2)
  real(r8), pointer :: n_retrans_acc(:,:)                                     ! N acquired by retranslocation (gN/m2)
  real(r8), pointer :: n_retrans_acc_total(:)                                 ! N acquired by retranslocation (gN/m2)
  real(r8), pointer :: free_nretrans_acc(:,:)                                 ! N acquired by retranslocation (gN/m2)
  real(r8), pointer :: npp_retrans_acc(:,:)                                   ! NPP used for the extraction (gC/m2)
  real(r8), pointer :: npp_retrans_acc_total(:)                               ! NPP used for the extraction (gC/m2)
  real(r8), pointer :: nt_uptake(:,:)                                         ! N uptake from retrans, active, and fix(gN/m2)
  real(r8), pointer :: npp_uptake(:,:)                                        ! NPP used by the uptakes (gC/m2)
  real(r8), pointer :: sminfrc(:,:)                                           ! fraction of N to handle NO3 / NH4 input 
  real(r8), pointer :: sminn_to_plant(:,:)                                    ! Nitrogen to plant (to handle NO3 / NH4 input)

  !----------NITRIF_DENITRIF-------------!
  !REPLACE NO# with N
  real(r8)            :: sminn_diff                                        ! A temporary limit for N uptake                  (gN/m2)
  real(r8)            :: active_limit1                                     ! A temporary limit for N uptake                  (gN/m2)
  real(r8),  pointer  :: costs_paths(:,:,:)                                 ! Costs for all paths (gC/gN) [patch,nlev,ipath]
  real(r8),  pointer  :: npp_frac_paths(:,:)                                 ! NPP fraction for all paths () [patch,nlev,ipath]
  real(r8),  pointer  :: npp_to_paths(:,:,:)                                 ! NPP spent all paths (gC/gN) [patch,nlev,ipath]
  real(r8),  pointer  :: n_from_paths(:,:,:)                                 ! NPP spent all paths (gC/gN) [patch,nlev,ipath]

  real(r8),  pointer  :: cost_active(:,:) ! cost of mycorrhizal                             (gC/gN)
  real(r8),  pointer  :: cost_nonmyc(:,:) ! cost of nonmyc                                  (gC/gN)

  real(r8),  pointer  :: sminn_layer(:,:)            ! Available N in each soil layer (gN/m2)
  real(r8),  pointer  :: sminn_layer_step(:,:,:)     ! A temporary variable for soil N (gN/m2) 
  real(r8),  pointer  :: n_am_acc(:)                            ! AM N uptake (gN/m2)
  real(r8),  pointer  :: n_ecm_acc(:)                           ! ECM N uptake (gN/m2)
  real(r8),  pointer  :: n_active_acc(:,:)                 ! Mycorrhizal N uptake (gN/m2)
  real(r8),  pointer  :: n_nonmyc_acc(:,:)                 ! Non-myc     N uptake (gN/m2)
  real(r8),  pointer  :: n_active_acc_total(:)                  ! Mycorrhizal N uptake (gN/m2)
  real(r8),  pointer  :: n_nonmyc_acc_total(:)                  ! Non-myc     N uptake (gN/m2)
  real(r8),  pointer  :: npp_active_acc(:,:)               ! Mycorrhizal N uptake used C (gC/m2)
  real(r8),  pointer  :: npp_nonmyc_acc(:,:)               ! Non-myc     N uptake used C (gC/m2)
  real(r8),  pointer  :: npp_active_acc_total(:)                ! Mycorrhizal N uptake used C (gC/m2)
  real(r8),  pointer  :: npp_nonmyc_acc_total(:)                ! Non-myc     N uptake used C (gC/m2)
  real(r8),  pointer  :: n_am_retrans(:)                        ! AM N uptake for offset (gN/m2)
  real(r8),  pointer  :: n_ecm_retrans(:)                       ! ECM N uptake for offset (gN/m2)
  real(r8),  pointer  :: n_active_retrans_total(:)              ! Mycorrhizal N for offset (gN/m2)
  real(r8),  pointer  :: n_nonmyc_retrans_total(:)              ! Non-myc     N for offset (gN/m2)
  real(r8),  pointer  :: n_passive_vr(:,:)           ! Layer passive N uptake (gN/m2)
  real(r8),  pointer  :: n_active_vr(:,:)            ! Layer mycorrhizal N uptake (gN/m2)
  real(r8),  pointer  :: n_nonmyc_vr(:,:)            ! Layer non-myc     N uptake (gN/m2)
  real(r8),  pointer  :: npp_active_retrans_total(:)            ! Mycorrhizal N uptake used C for offset (gN/m2)
  real(r8),  pointer  :: npp_nonmyc_retrans_total(:)             ! Non-myc N uptake used C for offset (gN/m2)
  


  ! Uptake fluxes for COST_METHOD=2
  ! actual npp to each layer for each uptake process
  real(r8),  pointer  ::                   npp_to_fixation(:) 
 ! real(r8),  pointer  ::                   npp_to_retrans(:)
  real(r8),  pointer  ::                   npp_to_active(:)

  ! fraction of carbon to each uptake process 
  real(r8),  pointer  ::                   npp_frac_to_fixation(:) 
  !real(r8),  pointer  ::                   npp_frac_to_retrans(:)
  real(r8),  pointer  ::                   npp_frac_to_nonmyc (:)  
   
  ! hypothetical fluxes on N in each layer 
  real(r8),  pointer  ::                  n_exch_fixation(:)        ! N aquired from one unit of C for fixation (unitless)
  real(r8),  pointer  ::                  n_exch_retrans(:)         ! N aquired from one unit of C for retrans (unitless)
  real(r8),  pointer  ::                  n_exch_active(:)          ! N aquired from one unit of C for act no3 (unitless)
  real(r8),  pointer  ::                  n_exch_nonmyc(:)          ! N aquired from one unit of C for nonmyc no3 (unitless) 
  
   !actual fluxes of N in each layer
  real(r8),  pointer  ::                  n_from_retrans(:)         ! N aquired in each layer of C for retrans (gN m-2 s-)
  real(r8),  pointer  ::                  free_Nretrans(:)          ! the total amount of NO3 and NH4                 (gN/m3/s) 

contains
  procedure , public  :: Init
  procedure , private  :: InitAllocate
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


  subroutine Init (this, bounds)
   !
   ! !DESCRIPTION:
   ! Initialize xxx
   !
   ! !USES:
   !
   ! !ARGUMENTS:
     class(cnfunmimicsplus_type)       :: this
     type(bounds_type) , intent(in)    :: bounds  

     call this%InitAllocate(bounds)
     call this%SetZeros(bounds)

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

     allocate(this%rootc_dens(bounds%begp:bounds%endp,1:nlevdecomp));        this%rootc_dens(:,:)= nan

     allocate(this%rootC(bounds%begp:bounds%endp));                          this%rootC(:)= nan

     allocate(this%n_uptake_myc_frac(bounds%begp:bounds%endp,1:nmyc));       this%n_uptake_myc_frac(:,:) = nan

     allocate(this%kc_active(bounds%begp:bounds%endp,1:nmyc));               this%kc_active(:,:) = nan

     allocate(this%kn_active(bounds%begp:bounds%endp,1:nmyc));               this%kn_active(:,:) = nan

     allocate(this%plant_ndemand_pool(bounds%begp:bounds%endp));             this%plant_ndemand_pool(:) = nan

     allocate(this%plant_ndemand_pool_step(bounds%begp:bounds%endp,1:nmyc)); this%plant_ndemand_pool_step(:,:) = nan
     
     allocate(this%litterfall_n_step(bounds%begp:bounds%endp,1:nmyc));       this%litterfall_n_step(:,:) = nan
     
     allocate(this%litterfall_c_step(bounds%begp:bounds%endp,1:nmyc));       this%litterfall_c_step(:,:) = nan
     
     allocate(this%tc_soisno(bounds%begc:bounds%endc,1:nlevdecomp));         this%tc_soisno(:,:) = nan
     
     allocate(this%npp_remaining(bounds%begp:bounds%endp,1:nmyc));           this%npp_remaining(:,:) = nan
     
     allocate(this%n_passive_step(bounds%begp:bounds%endp,1:nmyc));          this%n_passive_step(:,:) = nan
     
     allocate(this%n_passive_acc(bounds%begp:bounds%endp));                  this%n_passive_acc(:) = nan
     
     allocate(this%cost_fix(bounds%begp:bounds%endp,1:nlevdecomp));          this%cost_fix(:,:) = nan
     
     allocate(this%n_fix_acc(bounds%begp:bounds%endp,1:nmyc));               this%n_fix_acc(:,:) = nan
     
     allocate(this%n_fix_acc_total(bounds%begp:bounds%endp));                this%n_fix_acc_total(:) = nan
     
     allocate(this%npp_fix_acc(bounds%begp:bounds%endp,1:nmyc));             this%npp_fix_acc(:,:) = nan
     
     allocate(this%npp_fix_acc_total(bounds%begp:bounds%endp));              this%npp_fix_acc_total(:) = nan
     
     allocate(this%n_retrans_acc(bounds%begp:bounds%endp,1:nmyc));           this%n_retrans_acc(:,:) = nan
     
     allocate(this%n_retrans_acc_total(bounds%begp:bounds%endp));            this%n_retrans_acc_total(:) = nan
     
     allocate(this%free_nretrans_acc(bounds%begp:bounds%endp,1:nmyc));       this%free_nretrans_acc(:,:) = nan
     
     allocate(this%npp_retrans_acc(bounds%begp:bounds%endp,1:nmyc));         this%npp_retrans_acc(:,:) = nan
     
     allocate(this%npp_retrans_acc_total(bounds%begp:bounds%endp));          this%npp_retrans_acc_total(:) = nan
     
     allocate(this%nt_uptake(bounds%begp:bounds%endp,1:nmyc));               this%nt_uptake(:,:) = nan
     
     allocate(this%npp_uptake(bounds%begp:bounds%endp,1:nmyc));              this%npp_uptake(:,:) = nan

     allocate(this%costs_paths(bounds%begp:bounds%endp,1:nlevdecomp,1:npath6));  this%costs_paths(:,:,:) = nan
     
     allocate(this%cost_active(bounds%begp:bounds%endp,1:nlevdecomp));       this%cost_active(:,:) = nan

     allocate(this%cost_nonmyc(bounds%begp:bounds%endp,1:nlevdecomp));       this%cost_nonmyc(:,:) = nan

     allocate(this%sminn_layer(bounds%begc:bounds%endc,1:nlevdecomp));             this%sminn_layer(:,:) = nan

     allocate(this%sminn_layer_step(bounds%begp:bounds%endp,1:nlevdecomp,1:nmyc)); this%sminn_layer_step(:,:,:) = nan

     allocate(this%n_am_acc(bounds%begp:bounds%endp));                       this%n_am_acc = nan

     allocate(this%n_ecm_acc(bounds%begp:bounds%endp));                      this%n_ecm_acc = nan

     allocate(this%n_active_acc(bounds%begp:bounds%endp, 1:nmyc));           this%n_active_acc = nan

     allocate(this%n_nonmyc_acc(bounds%begp:bounds%endp, 1:nmyc));           this%n_nonmyc_acc = nan

     allocate(this%n_active_acc_total(bounds%begp:bounds%endp));             this%n_active_acc_total = nan

     allocate(this%n_nonmyc_acc_total(bounds%begp:bounds%endp));             this%n_nonmyc_acc_total(:) = nan

     allocate(this%npp_active_acc(bounds%begp:bounds%endp, 1:nmyc));         this%npp_active_acc(:,:) = nan

     allocate(this%npp_nonmyc_acc(bounds%begp:bounds%endp, 1:nmyc));         this%npp_nonmyc_acc(:,:) = nan

     allocate(this%npp_active_acc_total(bounds%begp:bounds%endp));           this%npp_active_acc_total(:) = nan

     allocate(this%npp_nonmyc_acc_total(bounds%begp:bounds%endp));           this%npp_nonmyc_acc_total(:) = nan

     allocate(this%n_am_retrans(bounds%begp:bounds%endp));                   this%n_am_retrans(:) = nan

     allocate(this%n_ecm_retrans(bounds%begp:bounds%endp));                  this%n_ecm_retrans(:) = nan

     allocate(this%n_active_retrans_total(bounds%begp:bounds%endp));         this%n_active_retrans_total(:) = nan

     allocate(this%n_nonmyc_retrans_total(bounds%begp:bounds%endp));         this%n_nonmyc_retrans_total(:) = nan

     allocate(this%n_passive_vr(bounds%begp:bounds%endp, 1:nlevdecomp));     this%n_passive_vr(:,:) = nan

     allocate(this%n_active_vr(bounds%begp:bounds%endp, 1:nlevdecomp));      this%n_active_vr(:,:) = nan

     allocate(this%n_nonmyc_vr(bounds%begp:bounds%endp, 1:nlevdecomp));      this%n_nonmyc_vr(:,:) = nan

     allocate(this%npp_active_retrans_total(bounds%begp:bounds%endp));       this%npp_active_retrans_total(:) = nan

     allocate(this%npp_nonmyc_retrans_total(bounds%begp:bounds%endp));       this%npp_nonmyc_retrans_total(:) = nan

     allocate(this%sminfrc(bounds%begc:bounds%endc,1:nlevdecomp));           this%sminfrc(:,:) = nan

     allocate(this%sminn_to_plant(bounds%begc:bounds%endc,1:nlevdecomp));    this%sminn_to_plant(:,:) = nan

    ! Uptake fluxes for COST_METHOD=2
    ! actual npp to each layer for each uptake process

     allocate(this%npp_to_paths(bounds%begp:bounds%endp,1:nlevdecomp,1:npath6));  this%npp_to_paths(:,:,:) = nan

     allocate(this%npp_frac_paths(bounds%begp:bounds%endp,1:nlevdecomp,1:npath6));  this%npp_frac_paths(:,:,:) = nan

     allocate(this%n_from_paths(bounds%begp:bounds%endp,1:nlevdecomp,1:npath6));  this%n_from_paths(:,:,:) = nan

     allocate(this%npp_to_fixation(1:nlevdecomp));                           this%npp_to_fixation(:) = nan

    ! allocate(this%npp_to_retrans(1:nlevdecomp));                            this%npp_to_retrans(:) = nan

     allocate(this%npp_to_active(1:nlevdecomp));                             this%npp_to_active(:) = nan

     ! fraction of carbon to each uptake process 
     allocate(this%npp_frac_to_fixation(1:nlevdecomp));                      this%npp_frac_to_fixation(:) = nan 

     !allocate(this%npp_frac_to_retrans(1:nlevdecomp));                       this%npp_frac_to_retrans(:) = nan

     allocate(this%npp_frac_to_nonmyc (1:nlevdecomp));                       this%npp_frac_to_nonmyc(:) = nan  
      
     ! hypothetical fluxes on N in each layer 
     allocate(this%n_exch_fixation(1:nlevdecomp));                           this%n_exch_fixation(:) = nan

     allocate(this%n_exch_retrans(1:nlevdecomp));                            this%n_exch_retrans(:) = nan

     allocate(this%n_exch_active(1:nlevdecomp));                             this%n_exch_active(:) = nan

     allocate(this%n_exch_nonmyc(1:nlevdecomp));                             this%n_exch_nonmyc(:) = nan
     
      !actual fluxes of N in each layer

     allocate(this%n_from_retrans(1:nlevdecomp));                            this%n_from_retrans(:) = nan

     allocate(this% free_Nretrans(bounds%begp:bounds%endp));                 this%free_Nretrans(:) = nan 
    
  
   end subroutine InitAllocate

  

  subroutine SetZeros (this, bounds)
   !
   ! !DESCRIPTION:
   ! Sets module data structure to zero
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
   

     
     this%rootc_dens_step                 = 0._r8
     
     !this%excess_carbon_acc               = 0._r8
     !this%burned_off_carbon               = 0._r8
     this%sminn_diff                      = 0._r8
     this%active_limit1                   = 0._r8
    ! this%total_N_conductance             = 0._r8

     this%rootc_dens(:,:)                = 0._r8
     this%rootC(:)                       = 0._r8
     this%n_uptake_myc_frac(:,:)         = 0._r8
     this%kc_active(:,:)                 = 0._r8
     this%kn_active(:,:)                 = 0._r8
     this%plant_ndemand_pool(:)          = 0._r8
     this%plant_ndemand_pool_step(:,:)   = 0._r8
     this%litterfall_n_step(:,:)         = 0._r8
     this%litterfall_c_step(:,:)         = 0._r8
     this%tc_soisno(:,:)                 = 0._r8 
     this%npp_remaining(:,:)             = 0._r8
     this%n_passive_step(:,:)            = 0._r8
     this%n_passive_acc(:)               = 0._r8


     this%costs_paths(:,:,:)              = big_cost
     this%npp_to_paths(:,:,:)             = 0.0_r8
     this%npp_frac_paths(:,:,:)           = 0.0_r8
     this%n_from_paths(:,:,:)             = 0.0_r8
     
     
     
     this%cost_fix(:,:)                  = 0._r8
     this%n_fix_acc(:,:)                 = 0._r8
     this%n_fix_acc_total(:)             = 0._r8
     this%npp_fix_acc(:,:)               = 0._r8
     this%npp_fix_acc_total(:)           = 0._r8
     this%n_retrans_acc(:,:)             = 0._r8
     this%n_retrans_acc_total(:)         = 0._r8
     this%free_nretrans_acc(:,:)         = 0._r8
     this%npp_retrans_acc(:,:)           = 0._r8
     this%npp_retrans_acc_total(:)       = 0._r8
     this%nt_uptake(:,:)                 = 0._r8
     this%npp_uptake(:,:)                = 0._r8
     this%cost_active(:,:)               = 0._r8
     this%cost_nonmyc(:,:)               = 0._r8
     this%sminn_layer(:,:)               = 0._r8
     this%sminn_layer_step(:,:,:)        = 0._r8
     this%n_am_acc                       = 0._r8
     this%n_ecm_acc                      = 0._r8
     this%n_active_acc                   = 0._r8
     this%n_nonmyc_acc                   = 0._r8
     this%n_active_acc_total             = 0._r8
     this%n_nonmyc_acc_total(:)          = 0._r8
     this%npp_active_acc(:,:)            = 0._r8
     this%npp_nonmyc_acc(:,:)            = 0._r8
     this%npp_active_acc_total(:)        = 0._r8
     this%npp_nonmyc_acc_total(:)        = 0._r8
     this%n_am_retrans(:)                = 0._r8
     this%n_ecm_retrans(:)               = 0._r8
     this%n_active_retrans_total(:)      = 0._r8
     this%n_nonmyc_retrans_total(:)      = 0._r8
     this%n_passive_vr(:,:)              = 0._r8
     this%n_active_vr(:,:)               = 0._r8
     this%n_nonmyc_vr(:,:)               = 0._r8
     this%npp_active_retrans_total(:)    = 0._r8
     this%npp_nonmyc_retrans_total(:)    = 0._r8
     this%sminfrc(:,:)                   = 0._r8
     this%sminn_to_plant(:,:)            = 0._r8

     this%npp_to_fixation(:)             = 0._r8
     this%npp_to_active(:)               = 0._r8
     
     this%npp_frac_to_fixation(:)        = 0._r8 
     this%npp_frac_to_nonmyc(:)          = 0._r8  
     
     this%n_exch_fixation(:)             = 0._r8
     this%n_exch_retrans(:)              = 0._r8
     this%n_exch_active(:)               = 0._r8
     this%n_exch_nonmyc(:)               = 0._r8
     
     this%n_from_retrans(:)              = 0._r8
     this%free_Nretrans(:)               = 0._r8 
   

  end subroutine SetZeros

!ECW move things from down

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
!--------------------------------------------------------------------
  

  ! read in parameters

  tString='ndays_on'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
  params_inst%ndays_on=tempr

  tString='ndays_off'
  call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
  if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
  params_inst%ndays_off=tempr
    
 end subroutine readParams

 !--------------------------------------------------------------------
  !---

 subroutine CNFUNMIMICSplusInit (bounds,cnveg_state_inst,cnveg_carbonstate_inst,cnveg_nitrogenstate_inst)
  !
  ! !DESCRIPTION:
  ! Initialization of FUN for MIMICSplus
  !
  ! !USES:
  use clm_varcon      , only: secspday, fun_period
  use clm_time_manager, only: get_step_size_real,get_nstep,get_curr_date,get_curr_days_per_year
  use SoilBiogeochemDecompCascadeMIMICSplusMod, only: calc_myc_roi
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
  real(r8)          :: dt                       ! timestep size (seconds)
  real(r8)          :: dayspyr                  ! days per year (days)
  real(r8)          :: timestep_fun             ! Timestep length for
  !  FUN (s)
  real(r8)          :: numofyear                ! number of days per
  !  year
  integer           :: nstep                    ! time step number
  integer           :: nstep_fun                ! Number of

  !  atmospheric timesteps between calls to FUN
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
  !--------------------------------------------------------------------
  !---
  ! Calculate some timestep-related values.
  !--------------------------------------------------------------------
  !---
  ! set time steps
  dt           = get_step_size_real()
  dayspyr      = get_curr_days_per_year()
  nstep        = get_nstep()
  timestep_fun = real(secspday * fun_period)
  nstep_fun    = int(secspday * dayspyr / dt) 

  ndays_on     = params_inst%ndays_on
  ndays_off    = params_inst%ndays_off

  !--------------------------------------------------------------------
  !---
  ! Decide if FUN will be called on this timestep.
  !--------------------------------------------------------------------
  !---
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
   soilbiogeochem_carbonstate_inst, cnfunmimicsplus_inst)
   
! !USES:
   use clm_time_manager, only : get_step_size_real, get_curr_date
   use clm_varpar      , only : nlevdecomp
   use clm_varcon      , only : secspday, smallValue, fun_period, tfrz, dzsoi_decomp, spval
   use clm_varctl      , only : use_nitrif_denitrif
   use PatchType       , only : patch
   use subgridAveMod   , only : p2c
   use pftconMod       , only : npcropmin
   use SoilBiogeochemDecompCascadeMIMICSplusMod , only : calc_myc_roi, decomp_rates_mimicsplus, calc_myc_mortality, calc_myc_mining_rates
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


   ! LOCAL VARIABLES:
   integer   :: fn                                ! number of values in pft filter
   integer   :: fp                                ! lake filter pft index 
   integer   :: fc                                ! lake filter column index
   integer   :: p, c                              ! pft index
   integer   :: g, l                              ! indices
   integer   :: j, i, k                           ! soil/snow level index
   integer   :: imyc                              ! Loop counters/work                                    orginal version in FUN: [istp]
   integer   :: ipath                             ! a local index                                         orginal version in FUN: [icost]
   integer   :: fixer                             ! 0 = non-fixer, 1 =fixer 
   logical   :: unmetDemand                       ! True while there is still demand for N
   integer   :: FIX 
   real(r8)  :: dt                               ! timestep size (s)
   real(r8)  :: dnpp                             ! Soulution to main FUN equation, C that actually goes into pathways
   real(r8)  :: stepspday                        ! number of timesteps per day (real)
   real(r8)  :: ndays_on        ! number of days to complete leaf onset
   real(r8)  :: ndays_off       ! number of days to complete leaf offset
   real(r8)  :: frac_ideal_C_use      ! How much less C do we use for 'buying' N than that needed to get to the ideal ratio?  fraction. 
  
   real(r8)  ::                  N_acquired
   real(r8)  ::                  C_spent
   real(r8)  ::                  leaf_narea            ! leaf n per unit leaf area in gN/m2 (averaged across canopy, which is OK for the cost calculation)
   
   real(r8)  ::                  sum_n_acquired        ! Sum N aquired from one unit of C (unitless)  
   real(r8)  ::                  burned_off_carbon     ! carbon wasted by poor allocation algorithm. If this is too big, we need a better iteration. 
   real(r8)  ::                  temp_n_flux  
   real(r8)  ::                  delta_cn              ! difference between 'ideal' leaf CN ration and actual leaf C:N ratio. C/N
   real(r8)  ::                  excess_carbon        ! how much carbon goes into the leaf C pool on account of the flexibleCN modifications.   
   real(r8)  ::                  excess_carbon_acc    ! excess accumulated over layers.
   
   !  WITHOUT GROWTH RESP
   real(r8) ::                    fixerfrac            ! what fraction of plants can fix?
   real(r8) ::                    npp_to_spend         ! how much carbon do we need to get rid of? 
   real(r8) ::                    soil_n_extraction    ! calculates total N pullled from soil
   real(r8) ::                    total_N_conductance  !inverse of C to of N for whole soil-leaf pathway
   real(r8) ::                    total_N_resistance   ! C to of N for whole soil -leaf pathway
   
   real(r8) ::                    paid_for_n_retrans
   real(r8) ::                    free_n_retrans
   real(r8) ::                    total_c_spent_retrans
   real(r8) ::                    total_c_accounted_retrans
 

   associate(ivt             => patch%itype                                          , & ! Input:   [integer  (:) ]  p
      leafcn                 => pftcon%leafcn                                        , & ! Input:   leaf C:N (gC/gN)
      season_decid           => pftcon%season_decid                                  , & ! Input:   binary flag for seasonal deciduous leaf habit (0 or 1)
      stress_decid           => pftcon%stress_decid                                  , & ! Input:   binary flag for stress deciduous leaf habit (0 or 1)
      a_fix                  => pftcon%a_fix                                         , & ! Input:   A BNF parameter
      b_fix                  => pftcon%b_fix                                         , & ! Input:   A BNF parameter
      c_fix                  => pftcon%c_fix                                         , & ! Input:   A BNF parameter
      s_fix                  => pftcon%s_fix                                         , & ! Input:   A BNF parameter
      akc_active             => pftcon%akc_active                                    , & ! Input:   A mycorrhizal uptake parameter
      akn_active             => pftcon%akn_active                                    , & ! Input:   A mycorrhizal uptake parameter
      ekc_active             => pftcon%ekc_active                                    , & ! Input:   A mycorrhizal uptake parameter
      ekn_active             => pftcon%ekn_active                                    , & ! Input:   A mycorrhizal upatke parameter
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
      frootc_storage         => cnveg_carbonstate_inst%frootc_storage_patch          , & ! Input:   [real(r8) (:)]  (gC/m2) fine root C storage
      livestemc              => cnveg_carbonstate_inst%livestemc_patch               , & ! Input:   [real(r8) (:)]  (gC/m2) live stem C
      livecrootc             => cnveg_carbonstate_inst%livecrootc_patch              , & ! Input:   [real(r8) (:)]  (gC/m2) live coarse root C
      leafc_storage_xfer_acc => cnveg_carbonstate_inst%leafc_storage_xfer_acc_patch  , & ! uutput:  [real(r8) (:)]  Accmulated leaf C transfer (gC/m2)
      storage_cdemand        => cnveg_carbonstate_inst%storage_cdemand_patch         , & ! Output:  [real(r8) (:)]  C use f rom the C storage pool
      tlai                   => canopystate_inst%tlai_patch                          , & ! Input:   [real(r8) (:)   ] one sided leaf area index
      leafn                  => cnveg_nitrogenstate_inst%leafn_patch                 , & ! Input:   [real(r8)  (:)] (gN/m2) leaf N
      frootn                 => cnveg_nitrogenstate_inst%frootn_patch                , & ! Input:   [real(r8)  (:)] (gN/m2) fine root N
      livestemn              => cnveg_nitrogenstate_inst%livestemn_patch             , & ! Input:   [real(r8)  (:)] (gN/m2) live stem N
      livecrootn             => cnveg_nitrogenstate_inst%livecrootn_patch            , & ! Input:   [real(r8)  (:)] (gN/m2) retranslocation N
      retransn               => cnveg_nitrogenstate_inst%retransn_patch              , & ! Input:   [real(r8)  (:)] (gN/m2) live coarse root N
      leafn_storage_xfer_acc => cnveg_nitrogenstate_inst%leafn_storage_xfer_acc_patch, & ! Output:  [real(r8)  (:)] Accmulated leaf N transfer (gC/m2)
      storage_ndemand        => cnveg_nitrogenstate_inst%storage_ndemand_patch       , & ! Output:  [real(r8)  (:)] N demand during the offset period
      leafc_to_litter        => cnveg_carbonflux_inst%leafc_to_litter_patch          , & ! Output:  [real(r8) (:) ]  leaf C litterfall (gC/m2/s)
      leafc_to_litter_fun    => cnveg_carbonflux_inst%leafc_to_litter_fun_patch      , & ! Output:  [real(r8) (:) ]  leaf C litterfall used by FUN (gC/m2/s)
      prev_leafc_to_litter   => cnveg_carbonflux_inst%prev_leafc_to_litter_patch     , & ! Output:  [real(r8) (:)] previous timestep leaf C litterfall flux (gC/m2/s)
      leafc_storage_to_xfer  => cnveg_carbonflux_inst%leafc_storage_to_xfer_patch    , & ! Output:  [real(r8) (:) ] 
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
      npp_growth             => cnveg_carbonflux_inst%npp_growth_patch               , & ! Output:  [real(r8) (:) ]  Total N uptake of FUN used C (gC/m2/s) 
      npp_burnedoff       => cnveg_carbonflux_inst%npp_burnedoff_patch            , & ! Output:  [real(r8) (:) ]  C  that cannot be used for N uptake(gC/m2/s)   
      leafc_change           => cnveg_carbonflux_inst%leafc_change_patch             , & ! Output:  [real(r8) (:) ]  Used C from the leaf (gC/m2/s)
      leafn_storage_to_xfer  => cnveg_nitrogenflux_inst%leafn_storage_to_xfer_patch  , & ! Output:  [real(r8) (:) ]
      plant_ndemand          => cnveg_nitrogenflux_inst%plant_ndemand_patch          , & ! Iutput:  [real(r8) (:) ]  N flux required to support initial GPP (gN/m2/s)
      plant_ndemand_retrans  => cnveg_nitrogenflux_inst%plant_ndemand_retrans_patch  , & ! Output:  [real(r8) (:) ]  N demand generated for FUN (gN/m2/s)
      plant_ndemand_season   => cnveg_nitrogenflux_inst%plant_ndemand_season_patch   , & ! Output:  [real(r8) (:) ]  N demand for seasonal deciduous forest (gN/m2/s)
      plant_ndemand_stress   => cnveg_nitrogenflux_inst%plant_ndemand_stress_patch   , & ! Output:  [real(r8) (:) ]  N demand for stress deciduous forest   (gN/m2/s)
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
      Npassive               => cnveg_nitrogenflux_inst%Npassive_patch               , & ! Output:  [real(r8) (:) ]  Passive N uptake (gN/m2/s)
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
      sminn_to_plant_fun_vr  => cnveg_nitrogenflux_inst%sminn_to_plant_fun_vr_patch           , & ! Output: [real(r8) (:) ]  Total layer soil N uptake of FUN (gN/m2/s) 
      sminn_to_plant_fun_no3_vr  => cnveg_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_patch   , & ! Output:  [real(r8) (:) ]  Total layer no3 uptake of FUN (gN/m2/s)
      sminn_to_plant_fun_nh4_vr  => cnveg_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_patch   , & ! Output:  [real(r8) (:) ]  Total layer nh4 uptake of FUN (gN/m2/s)
      sminn_to_plant_vr      => soilbiogeochem_nitrogenflux_inst%sminn_to_plant_vr_col        , & ! Output:  [real(r8) (: ,:) ]
      smin_no3_to_plant_vr   => soilbiogeochem_nitrogenflux_inst%smin_no3_to_plant_vr_col     , & ! Input:  [real(r8) (:,:) ]
      smin_nh4_to_plant_vr   => soilbiogeochem_nitrogenflux_inst%smin_nh4_to_plant_vr_col     , & ! Input:  [real(r8) (:,:) ]        
      soilc_change           => cnveg_carbonflux_inst%soilc_change_patch                      , & ! Output:  [real(r8)(:) ]  Used C from the soil (gC/m2/s)
      h2osoi_liq             => waterstatebulk_inst%h2osoi_liq_col                            , & ! Input:   [real(r8) (:,:)] liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)
      t_soisno               => temperature_inst%t_soisno_col                                 , & ! Input:   [real(r8) (:,:)] soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
      crootfr                => soilstate_inst%crootfr_patch                                  , & ! Input:   [real(r8) (:,:)] fraction of roots for carbon in each soil layer  (nlevgrnd)

      active_limit1          => cnfunmimicsplus_inst%active_limit1                            , &
      !burned_off_carbon      => cnfunmimicsplus_inst%burned_off_carbon                        , &





      costs_paths            => cnfunmimicsplus_inst%costs_paths                              , &
      npp_to_paths           => cnfunmimicsplus_inst%npp_to_paths                             , &
      npp_frac_paths         => cnfunmimicsplus_inst%npp_frac_paths                           , &
      n_from_paths            => cnfunmimicsplus_inst%n_from_paths                             , &
      
      
      cost_active            => cnfunmimicsplus_inst%cost_active                              , &
      cost_fix               => cnfunmimicsplus_inst%cost_fix                                 , &
      cost_nonmyc            => cnfunmimicsplus_inst%cost_nonmyc                              , &
      dn                     => cnfunmimicsplus_inst%dn                                       , &
      !excess_carbon_acc      => cnfunmimicsplus_inst%excess_carbon_acc                        , &
      free_Nretrans          => cnfunmimicsplus_inst%free_Nretrans                            , &
      free_nretrans_acc      => cnfunmimicsplus_inst%free_nretrans_acc                        , &
      kc_active              => cnfunmimicsplus_inst%kc_active                                , &
      kn_active              => cnfunmimicsplus_inst%kn_active                                , &
      litterfall_c_step      => cnfunmimicsplus_inst%litterfall_c_step                        , &
      litterfall_n_step      => cnfunmimicsplus_inst%litterfall_n_step                        , &
      n_active_acc           => cnfunmimicsplus_inst%n_active_acc                             , &
      n_active_acc_total     => cnfunmimicsplus_inst%n_active_acc_total                       , &
      n_active_retrans_total => cnfunmimicsplus_inst%n_active_retrans_total                   , &
      n_active_vr            => cnfunmimicsplus_inst%n_active_vr                              , &
      n_am_acc               => cnfunmimicsplus_inst%n_am_acc                                 , &
      n_am_retrans           => cnfunmimicsplus_inst%n_am_retrans                             , &
      n_ecm_acc              => cnfunmimicsplus_inst%n_ecm_acc                                , &
      n_ecm_retrans          => cnfunmimicsplus_inst%n_ecm_retrans                            , &
      n_exch_active          => cnfunmimicsplus_inst%n_exch_active                            , &
      n_exch_fixation        => cnfunmimicsplus_inst%n_exch_fixation                          , &
      n_exch_nonmyc          => cnfunmimicsplus_inst%n_exch_nonmyc                            , &
      n_exch_retrans         => cnfunmimicsplus_inst%n_exch_retrans                           , &
      n_fix_acc              => cnfunmimicsplus_inst%n_fix_acc                                , &
      n_fix_acc_total        => cnfunmimicsplus_inst%n_fix_acc_total                          , &
      n_from_retrans         => cnfunmimicsplus_inst%n_from_retrans                           , &
      n_nonmyc_acc           => cnfunmimicsplus_inst%n_nonmyc_acc                             , &
      n_nonmyc_acc_total     => cnfunmimicsplus_inst%n_nonmyc_acc_total                       , &
      n_nonmyc_retrans_total => cnfunmimicsplus_inst%n_nonmyc_retrans_total                   , &
      n_nonmyc_vr            => cnfunmimicsplus_inst%n_nonmyc_vr                              , &
      n_passive_acc          => cnfunmimicsplus_inst%n_passive_acc                            , &
      n_passive_step         => cnfunmimicsplus_inst%n_passive_step                           , &
      n_passive_vr           => cnfunmimicsplus_inst%n_passive_vr                             , &
      n_retrans_acc          => cnfunmimicsplus_inst%n_retrans_acc                            , &
      n_retrans_acc_total    => cnfunmimicsplus_inst%n_retrans_acc_total                      , &
      n_uptake_myc_frac      => cnfunmimicsplus_inst%n_uptake_myc_frac                        , &
      npp_active_acc         => cnfunmimicsplus_inst%npp_active_acc                           , &
      npp_active_acc_total   => cnfunmimicsplus_inst%npp_active_acc_total                     , &
      npp_active_retrans_total => cnfunmimicsplus_inst%npp_active_retrans_total               , &
      npp_fix_acc            => cnfunmimicsplus_inst%npp_fix_acc                              , &
      npp_fix_acc_total      => cnfunmimicsplus_inst%npp_fix_acc_total                        , &
      npp_frac_to_fixation   => cnfunmimicsplus_inst%npp_frac_to_fixation                     , &
      npp_frac_to_nonmyc     => cnfunmimicsplus_inst%npp_frac_to_nonmyc                       , &
      npp_nonmyc_acc         => cnfunmimicsplus_inst%npp_nonmyc_acc                           , &
      npp_nonmyc_acc_total   => cnfunmimicsplus_inst%npp_nonmyc_acc_total                     , &
      npp_nonmyc_retrans_total => cnfunmimicsplus_inst%npp_nonmyc_retrans_total               , &
      npp_remaining          => cnfunmimicsplus_inst%npp_remaining                            , &
      npp_retrans_acc        => cnfunmimicsplus_inst%npp_retrans_acc                          , &
      npp_retrans_acc_total  => cnfunmimicsplus_inst%npp_retrans_acc_total                    , &
      npp_to_active          => cnfunmimicsplus_inst%npp_to_active                            , &
      npp_to_fixation        => cnfunmimicsplus_inst%npp_to_fixation                          , &
      npp_uptake             => cnfunmimicsplus_inst%npp_uptake                               , &
      nt_uptake              => cnfunmimicsplus_inst%nt_uptake                                , &
      plant_ndemand_pool     => cnfunmimicsplus_inst%plant_ndemand_pool                       , &
      plant_ndemand_pool_step => cnfunmimicsplus_inst%plant_ndemand_pool_step                 , &
      rootC                  => cnfunmimicsplus_inst%rootC                                    , &
      rootc_dens             => cnfunmimicsplus_inst%rootc_dens                               , &
      rootc_dens_step        => cnfunmimicsplus_inst%rootc_dens_step                          , &
      sminfrc                => cnfunmimicsplus_inst%sminfrc                                  , &
      sminn_diff             => cnfunmimicsplus_inst%sminn_diff                               , &
      sminn_layer_step       => cnfunmimicsplus_inst%sminn_layer_step                         , &
      sminn_to_plant         => cnfunmimicsplus_inst%sminn_to_plant                           , &
      !total_N_conductance    => cnfunmimicsplus_inst%total_N_conductance                           , &
      ! C and N soil pools:
      decomp_cpools_vr      => soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col          , &
      decomp_npools_vr      => soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col          , &
      tc_soisno              => cnfunmimicsplus_inst%tc_soisno                                )
      

      dt           = get_step_size_real()
      stepspday    = secspday / dt


     ! Time step of FUNMIMICSplus: once per day                
     call cnfunmimicsplus_inst%SetZeros(bounds)                                        ! set everything to zero

     do fp = 1,num_soilp          ! PFT Starts
     p = filter_soilp(fp)
   
     ! Calculate plant N based on; N lost with litterfall, root C, ...
     ! litterfall_n(p) =  (leafc_to_litter_fun(p) / leafcn_offset(p))  * dt ! NOT USED
      rootC(p)        =  frootc(p)
    
     ! plantN(p)       =  leafn(p) + frootn(p) + livestemn(p) + livecrootn(p)    ! NOT USED
     ! Calculate patch plant C:N ratio
      if (n_allometry(p).gt.0._r8) then 
          plantCN(p)  = c_allometry(p)/n_allometry(p)
      else
          plantCN(p)  = 0._r8 
      end if
     end do                       ! PFT ends

      ! Calculates fraction of N uptake per mycorrhizal type per PFT
  do imyc = 1, nmyc
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      if (imyc.eq.ecm_step) then
         n_uptake_myc_frac(p,imyc)      = perecm(ivt(p))
         kc_active(p,imyc)              = ekc_active(ivt(p))
         kn_active(p,imyc)              = ekn_active(ivt(p))
      else
         n_uptake_myc_frac(p,imyc)      = 1._r8 - perecm(ivt(p))
         kc_active(p,imyc)              = akc_active(ivt(p))
         kn_active(p,imyc)              = akn_active(ivt(p))
      end if

      ! 824 N avaliabale throug litterfall
     if(leafc(p)>0.0_r8)then                       ! N available in leaf which fell off in this timestep. Same fraction loss as C.    
      litterfall_c_step(p,imyc)         =   dt * n_uptake_myc_frac(p,imyc) * leafc_to_litter_fun(p) 
      litterfall_n_step(p,imyc)         =   dt * n_uptake_myc_frac(p,imyc) * leafn(p) * leafc_to_litter_fun(p)/leafc(p) 
     endif 

   end do
  end do     
    
  ! Vertically-resolved plant uptake of soil NO3 / NH4 is an INPUT variable from SoilBiogeochemNitrogenFluxType,
  ! therefore, needs to be kept separate 
  do j = 1, nlevdecomp                  ! Loop through soil column 
   do fp = 1,num_soilp       
     p = filter_soilp(fp)
     c = patch%column(p)
     
     if ((smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)) > 0)  then           
     sminfrc(c,j) = smin_no3_to_plant_vr(c,j)/(smin_no3_to_plant_vr(c,j)+smin_nh4_to_plant_vr(c,j))      ! sminn_nh4_frac = 1._r8 - sminn_no3_frac
     else 
     sminfrc(c,j)=0._r8
     sminn_to_plant(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
     end if
   end do 
  end do    
     

  ! Calculaings available N in each soil layer (gN/m2)   ! 840 - 876
   do imyc=1,nmyc
    do j = 1, nlevdecomp
       do fp = 1,num_soilp        ! PFT Starts
     p = filter_soilp(fp)
     c = patch%column(p)
     sminn_layer_step(p,j,imyc)=max((smin_no3_to_plant_vr(c,j)+smin_nh4_to_plant_vr(c,j)) * dzsoi_decomp(j) * dt,0.0_r8) * n_uptake_myc_frac(p,imyc)
      end do
    end do
   end do


   ! PFT loop
   ! Calculating leaf C / N storage demand for retranslocation
pft:  do fp = 1,num_soilp        ! PFT Starts
      p = filter_soilp(fp)
      c = patch%column(p)
     
      sminn_to_plant_fun_nh4_vr(p,:) = 0._r8
      sminn_to_plant_fun_no3_vr(p,:) = 0._r8   
      burned_off_carbon = 0.0_r8 !ECW GET RID OF THIS AND USE npp_burndoff  from the getgo
      excess_carbon_acc = 0.0_r8
     
      ! I have turned off this retranslocation functionality for now. To be rolled back in to a new version later on once the rest of the mode is working OK. RF

      if (season_decid(ivt(p)) == 1._r8.or.stress_decid(ivt(p)) == 1._r8) then
         if (onset_flag(p) == 1._r8) then
           ! leafc_storage_xfer_acc(p) = leafc_storage_xfer_acc(p) + leafc_storage_to_xfer(p) * dt
            leafn_storage_xfer_acc(p) = leafn_storage_xfer_acc(p) + leafn_storage_to_xfer(p) * dt
         end if
         if (offset_flag(p) == 1._r8) then
            !storage_cdemand(p)        = leafc_storage(p)          / (ndays_off * steppday)
            storage_ndemand(p)        = leafn_storage_xfer_acc(p) / (ndays_off * stepspday)
            storage_ndemand(p)        = max(storage_ndemand(p),0._r8)
         else
            !storage_cdemand(p)        = 0._r8    
            storage_ndemand(p)        = 0._r8   
         end if
      else
          !storage_cdemand(p)          = 0._r8
          storage_ndemand(p)          = 0._r8 
      end if   ! end for deciduous

      
      ! Avaliable carbon for growth or Nitrogen uptake
      !availc_pool(p)            =  availc(p)        *  dt

      if (availc(p) > 0._r8) then
         do j = 1, nlevdecomp
            rootc_dens(p,j)     =  crootfr(p,j) * rootC(p)
         end do
      end if

      plant_ndemand_pool(p)     =  plant_ndemand(p) *  dt
      plant_ndemand_pool(p)     =  max(plant_ndemand_pool(p),0._r8)
      plant_ndemand_retrans(p)  =  storage_ndemand(p) 

      ! Nitrogen demand of plant & remaining carbon (NPP) per timestep      1048 -1050
stp:  do imyc = ecm_step, am_step        ! TWO STEPS

      unmetDemand              = .TRUE.
      plant_ndemand_pool_step(p,imyc)   = plant_ndemand_pool(p)    * n_uptake_myc_frac(p,imyc) 
      npp_remaining(p,imyc)             = availc(p)*dt             * n_uptake_myc_frac(p,imyc)
         

      ! COST FIXATION PATHWAY             1056 - 1066
      ! checks which photosyntetic pathway plant has (C3 / C4) and if they can do nitrogen fixation   
      do j = 1, nlevdecomp
         tc_soisno(c,j)          = t_soisno(c,j)  -   tfrz     ! Soil temperature
         if(pftcon%c3psn(patch%itype(p)).eq.1)then
           fixer=1
         else
           fixer=0
         endif
         costs_paths(p,j,ipfix)     = fun_cost_fix(fixer,a_fix(ivt(p)),b_fix(ivt(p))&
         ,c_fix(ivt(p)) ,big_cost,crootfr(p,j),s_fix(ivt(p)),tc_soisno(c,j))
      end do

      ! ACTIVE UPTAKE                  1081 - 1102
      ! Mycorrhizal Uptake Cost
      do j = 1,nlevdecomp
         rootc_dens_step             = rootc_dens(p,j) *  n_uptake_myc_frac(p,imyc)
         if (imyc == ecm_step) then
            call calc_myc_roi(decomp_cpools_vr(c,j,i_ecm_myc),decomp_npools_vr(c,j,i_ecm_myc) , &
                              decomp_cpools_vr(c,j,i_phys_som),decomp_cpools_vr(c,j,i_avl_som),decomp_cpools_vr(c,j,i_chem_som), &
                              decomp_npools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_chem_som), &
                              sminn_layer_step(p,j,imyc), ecm_step, dzsoi_decomp(j), costs_paths(p,j,ipano3))
            call calc_myc_roi(decomp_cpools_vr(c,j,i_ecm_myc),decomp_npools_vr(c,j,i_ecm_myc) , &
                              decomp_cpools_vr(c,j,i_phys_som),decomp_cpools_vr(c,j,i_avl_som),decomp_cpools_vr(c,j,i_chem_som), &
                              decomp_npools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_chem_som), &
                              sminn_layer_step(p,j,imyc), ecm_step, dzsoi_decomp(j), costs_paths(p,j,ipanh4))
         elseif (imyc == am_step) then
            call calc_myc_roi(decomp_cpools_vr(c,j,i_am_myc),decomp_npools_vr(c,j,i_am_myc) , &  ! sminn_frac no3 no4
                              decomp_cpools_vr(c,j,i_phys_som),decomp_cpools_vr(c,j,i_avl_som),decomp_cpools_vr(c,j,i_chem_som), &
                              decomp_npools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_chem_som), &
                              sminn_layer_step(p,j,imyc), am_step, dzsoi_decomp(j), costs_paths(p,j,ipano3))
            call calc_myc_roi(decomp_cpools_vr(c,j,i_am_myc),decomp_npools_vr(c,j,i_am_myc) , &
                              decomp_cpools_vr(c,j,i_phys_som),decomp_cpools_vr(c,j,i_avl_som),decomp_cpools_vr(c,j,i_chem_som), &
                              decomp_npools_vr(c,j,i_phys_som),decomp_npools_vr(c,j,i_chem_som), &
                              sminn_layer_step(p,j,imyc), am_step, dzsoi_decomp(j), costs_paths(p,j,ipanh4))
         else 
            call endrun(msg='ERROR: number of mycorrhiza in FUN is more than number of mycorrhizal pools in MimicsPlus')
      end do
      
      ! Non-mycorrhizal Uptake Cost
      do j = 1,nlevdecomp
         rootc_dens_step             = rootc_dens(p,j)  *  n_uptake_myc_frac(p,imyc)
         costs_paths(p,j,ipnmno3)      = fun_cost_nonmyc(sminn_layer_step(p,j,imyc) &
         ,big_cost,kc_nonmyc(ivt(p)),kn_nonmyc(ivt(p)) ,rootc_dens_step,crootfr(p,j),smallValue)
         costs_paths(p,j,ipnmnh4)      = fun_cost_nonmyc(sminn_layer_step(p,j,imyc) &
         ,big_cost,kc_nonmyc(ivt(p)),kn_nonmyc(ivt(p)) ,rootc_dens_step,crootfr(p,j),smallValue)
      end do


      ! Remove C required to pair with N from passive uptake from the available pool. 
      npp_remaining(p,imyc)  =   npp_remaining(p,imyc) - n_passive_step(p,imyc)*plantCN(p)

      total_N_conductance = 0.0_r8 
      fix_loop: do FIX =plants_are_fixing, plants_not_fixing !loop around percentages of fixers and nonfixers, with differnt costs. 

          if(FIX==plants_are_fixing)then ! How much of the carbon in this PFT can in principle be used for fixation? 
            ! This is analagous to fixing the % of fixers for a given PFT - may not be realistic in the long run
            ! but prevents wholesale switching to fixer dominance during e.g. CO2 fertilization.  
            fixerfrac = FUN_fracfixers(ivt(p))
          else
            fixerfrac = 1.0_r8 - FUN_fracfixers(ivt(p))
          endif 
          npp_to_spend = npp_remaining(p,imyc)  * fixerfrac !put parameter here.
          ! has to be zeroed since depend on accumula
          n_from_active(1:nlevdecomp) = 0._r8
          n_from_nonmyc(1:nlevdecomp) = 0._r8
          n_from_paths(p,:,ipano3:ipnmnh4) !act and nonmyc boths nh4 and no3

          !--------------------------------------------------------------------
          !-----------
          !           Calculate Integrated Resistance OF WHOLE SOIL COLUMN
          !--------------------------------------------------------------------
          !----------- 

          sum_n_acquired      = 0.0_r8
          total_N_conductance = 0.0_r8
          do j = 1, nlevdecomp
            ! Method changed from FUN-resistors method to a method which allocates fluxs based on conductance. rosief
            ! Sum the conductances             
           total_N_conductance  = total_N_conductance + 1._r8/ costs_paths(p,j,ipano3)  &
                                 + 1._r8/ costs_paths(p,j,ipanh4) + 1._r8/ costs_paths(p,j,ipnmno3)  &
                                 + 1._r8/ costs_paths(p,j,ipnmnh4) 
            if(FIX==plants_are_fixing)then
                total_N_conductance  = total_N_conductance  + 1._r8/+ 1._r8/ costs_paths(p,j,ipfix) 
            end if      
          end do 
          
          do j = 1, nlevdecomp     
             ! Calculate npp allocation to pathways proportional to their exchange rate (N/C) 
             npp_frac_paths(p,j,ipanh4) = (1._r8/costs_paths(p,j,ipanh4)) / total_N_conductance
             npp_frac_paths(p,j,ipnmnh4) = (1._r8/costs_paths(p,j,ipnmnh4)) / total_N_conductance
             npp_frac_paths(p,j,ipano3) = (1._r8/costs_paths(p,j,ipano3)) / total_N_conductance
             npp_frac_paths(p,j,ipnmno3) = (1._r8/costs_paths(p,j,ipnmno3)) / total_N_conductance
             
             if(FIX==plants_are_fixing)then
               npp_frac_paths(p,j,ipfix) = (1.0_r8 * 1._r8/costs_paths(p,j,ipfix)) / total_N_conductance
             else
               npp_frac_paths(p,j,ipfix) = 0.0_r8 
             end if           
             ! Total N aquired from one unit of carbon  (N/C)
             sum_n_acquired = sum_n_acquired  + npp_frac_paths(p,j,ipanh4)/costs_paths(p,j,ipanh4) + &
                           npp_frac_paths(p,j,ipano3)/costs_paths(p,j,ipano3) + &
                        npp_frac_paths(p,j,ipnmnh4)/costs_paths(p,j,ipnmnh4) + npp_frac_paths(p,j,ipnmno3)/costs_paths(p,j,ipnmno3)

                                     
             if(FIX==plants_are_fixing)then
               sum_n_acquired= sum_n_acquired +  npp_frac_paths(p,j,ipfix)/costs_paths(p,j,ipfix)
             end if                                                          
          end do
 
          total_N_resistance = 1.0_r8/sum_n_acquired !gC/gN
    
          !  Calculate appropriate degree of retranslocation
          if(leafc(p).gt.0.0_r8.and.litterfall_n_step(p,imyc)* fixerfrac>0.0_r8.and.ivt(p) <npcropmin)then
             call fun_retranslocation(p,dt,npp_to_spend,&
                           litterfall_c_step(p,imyc)* fixerfrac,&
                           litterfall_n_step(p,imyc)* fixerfrac,&
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
          npp_retrans_acc(p,imyc) = npp_retrans_acc(p,imyc) + total_c_spent_retrans 
          ! add to to C spent pool                              
          n_retrans_acc(p,imyc)     = n_retrans_acc(p,imyc)     + paid_for_n_retrans
          free_nretrans_acc(p,imyc) = free_nretrans_acc(p,imyc) + free_n_retrans
          ! add N to the acquired from retrans pool 
    

          ! Spend C on extracting N
          if (plant_ndemand_pool_step(p,imyc) .gt. 0._r8) then    ! unmet demand

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


          dnpp  = npp_to_spend / ( (1.0_r8+grperc(ivt(p)))*(plantCN(p) / total_N_resistance) + 1._r8)  
          dnpp  = dnpp * frac_ideal_C_use
      
          !hypothetical amount of N acquired. 
          dn    = dnpp / total_N_resistance
          do j = 1,nlevdecomp
                 
           ! RF How much of this NPP carbon do we allocate to the different pathways? fraction x gC/m2/s?
           ! Could this code now be put in a matrix? 

           npp_to_paths(p,j,ipano3:ipnmnh4) = npp_frac_paths(p,j,ipano3:ipnmnh4) * dNPP
                                
           if(FIX==plants_are_fixing)then
            npp_to_paths(p,j,ipfix) = npp_frac_paths(p,j,ipfix) * dNPP
           else
            npp_to_paths(p,j,ipfix) = 0.0_r8
           end if    

           n_from_paths(p,j,ipano3:ipnmnh4) = npp_to_paths(p,j,ipano3:ipnmnh4) / costs_paths(p,j,ipano3:ipnmnh4)
                           
           if(FIX==plants_are_fixing)then
            n_from_paths(p,j,ipfix) = npp_to_paths(p,j,ipfix) / costs_paths(p,j,ipfix)
           else
            n_from_paths(p,j,ipfix) = 0.0_r8
           end if                                                   
          end do

                 
                ! Check if LIMITS of pools were exceeded:
          !!!!MVD THIS IS WHERE WE ENDED LAST TIME IN THE PROCESS OF REPLACING VARIABLES!!!
          do j = 1,nlevdecomp    
            ! ACTIVE UPTAKE LIMIT 
            active_limit1          = sminn_layer_step(p,j,imyc) * fixerfrac 
             ! trying to remove too much nh4 from soil. 
             if (n_from_active(j) + n_from_nonmyc(j).gt.active_limit1) then 
                sminn_diff          = n_from_active(j) + n_from_nonmyc(j) - active_limit1
                temp_n_flux = n_from_active(j)
                ! divide discrepancy between sources
                n_from_active(j)     = n_from_active(j) - sminn_diff &
                                           * (n_from_active(j) /(n_from_active(j) + n_from_nonmyc(j)))
                n_from_nonmyc(j)     = n_from_nonmyc(j) - sminn_diff &
                                      * (n_from_nonmyc(j) /(temp_n_flux + n_from_nonmyc(j)))
                npp_to_active(j)     = n_from_active(j) * cost_active(p,j) 
                npp_to_nonmyc(j)     = n_from_nonmyc(j) * cost_nonmyc(p,j)
                           
            end if

                  
            N_acquired     =  n_from_active(j)+n_from_nonmyc(j)          ! How much N did we end up with
            C_spent        =   npp_to_active(j)+npp_to_nonmyc(j)         ! How much did it actually cost? 
                         
            if(FIX==plants_are_fixing)then
               N_acquired  = N_acquired + n_from_fixation(j)
               C_spent     = C_spent + npp_to_fixation(j) 
            end if
            
            ! How much C did we allocate or spend in this layer? 
            npp_to_spend    = npp_to_spend   - C_spent - (N_acquired  * plantCN(p)*(1.0_r8+ grperc(ivt(p))))
            
            ! Accumulate those fluxes
            nt_uptake(p,imyc)             = nt_uptake(p,imyc)       + N_acquired
            npp_uptake(p,imyc)            = npp_uptake(p,imyc)      + C_spent


            !-------------------- N flux accumulation------------!
            n_active_acc(p,imyc)      = n_active_acc(p,imyc) + n_from_active(j)
            n_nonmyc_acc(p,imyc)      = n_nonmyc_acc(p,imyc) + n_from_nonmyc(j)
            !-------------------- C flux accumulation------------!
            npp_active_acc(p,imyc) = npp_active_acc(p,imyc)  + npp_to_active(j)
            npp_nonmyc_acc(p,imyc) = npp_nonmyc_acc(p,imyc)  + npp_to_nonmyc(j)          
       
            if(FIX == plants_are_fixing)then
              n_fix_acc(p,imyc)          = n_fix_acc(p,imyc)         + n_from_fixation(j)
              npp_fix_acc(p,imyc)        = npp_fix_acc(p,imyc)       + npp_to_fixation(j)
            end if
       
            end do

            ! CHECK Nitrogen amount:
            if (npp_to_spend .ge. 0.0000000000001_r8)then
               burned_off_carbon =  burned_off_carbon + npp_to_spend
            end if

             ! add vertical fluxes to patch arrays 
            do j = 1,nlevdecomp
               n_active_vr(p,j)      =  n_active_vr(p,j)      + n_from_active(j)
               n_nonmyc_vr(p,j)      =  n_nonmyc_vr(p,j)      + n_from_nonmyc(j)
            
            end do
          end if !unmet demand`
         
      end do fix_loop ! FIXER. 
                      
             if (imyc.eq.ecm_step) then
                n_ecm_acc(p)          =  n_active_acc(p,imyc)
             else
                n_am_acc(p)           =  n_active_acc(p,imyc)
             end if
       
       
                 ! Accumulate total column N fluxes over imyc
             n_active_acc_total(p)        =  n_active_acc_total(p)    + n_active_acc(p,imyc)
             n_nonmyc_acc_total(p)        =  n_nonmyc_acc_total(p)    + n_nonmyc_acc(p,imyc)
             n_fix_acc_total(p)           =  n_fix_acc_total(p)           + n_fix_acc(p,imyc)
             n_retrans_acc_total(p)       =  n_retrans_acc_total(p)       + n_retrans_acc(p,imyc)
             free_nretrans(p)             =  free_nretrans(p)            + free_nretrans_acc(p,imyc)
       
             ! Accumulate total column C fluxes over imyc
             npp_active_acc_total(p)      =  npp_active_acc_total(p)  + npp_active_acc(p,imyc)
             npp_nonmyc_acc_total(p)      =  npp_nonmyc_acc_total(p)  + npp_nonmyc_acc(p,imyc)
             npp_fix_acc_total(p)         =  npp_fix_acc_total(p)         + npp_fix_acc(p,imyc)
             npp_retrans_acc_total(p)     =  npp_retrans_acc_total(p)     + npp_retrans_acc(p,imyc) 
                    
             !end if  ! plant_ndemand_pool_step > 0._r8
end do stp  
             ! Turn step level quantities back into fluxes per second. 
             Npassive(p)               = n_passive_acc(p)/dt
             Nfix(p)                   = n_fix_acc_total(p)/dt                   
             retransn_to_npool(p)      = n_retrans_acc_total(p)/dt
             ! Without matrix solution
             if(.not. use_matrixcn)then
                free_retransn_to_npool(p) = free_nretrans(p)/dt
             ! With matrix solution (when it comes in)
             else
             end if
             ! this is the N that comes off leaves. 
             Nretrans(p)               = retransn_to_npool(p) + free_retransn_to_npool(p)
    
             
             !Extract active uptake N from soil pools. 
             do j = 1, nlevdecomp
             !RF change. The N fixed doesn't actually come out of the soil mineral pools, it is 'new'... 
             sminn_to_plant_fun_vr(p,j)    = (n_passive_vr(p,j)  + n_active_vr(p,j) + n_nonmyc_vr(p,j))/(dzsoi_decomp(j)*dt)
             end do
    
    
             !---------------------------N fluxes--------------------!
             Nactive(p)            = n_active_acc_total(p)/dt   + n_active_retrans_total(p)/dt
             Necm(p)               = n_ecm_acc(p)/dt            + n_ecm_retrans(p)/dt
             Nam(p)                = n_am_acc(p)/dt             + n_am_retrans(p)/dt
             Nnonmyc(p)            = n_nonmyc_acc_total(p)/dt   + n_nonmyc_retrans_total(p)/dt
             plant_ndemand_retrans(p)  = plant_ndemand_retrans(p)/dt
             Nuptake(p)                = Nactive(p) + Nnonmyc(p) + Nfix(p) + Npassive(p) + &
                                         retransn_to_npool(p)+free_retransn_to_npool(p) 
                                          
             ! free N goes straight to the npool, not throught Nuptake...
             sminn_to_plant_fun(p)     = Nactive(p) + Nnonmyc(p) + Nfix(p) + Npassive(p) 
    
             soil_n_extraction = ( sum(n_active_vr(p,1: nlevdecomp))+sum(n_nonmyc_vr(p,1: nlevdecomp)))
    
             !---------------------------C fluxes--------------------!
             npp_Nactive(p)        = npp_active_acc_total(p)/dt + npp_active_retrans_total(p)/dt
             npp_Nnonmyc(p)        = npp_nonmyc_acc_total(p)/dt + npp_nonmyc_retrans_total(p)/dt               
    
             npp_Nfix(p)               = npp_fix_acc_total(p)/dt      
             npp_Nretrans(p)           = npp_retrans_acc_total(p)/dt  
            
             !---------------------------Extra Respiration Fluxes--------------------!      
             soilc_change(p)           = (npp_active_acc_total(p) + npp_nonmyc_acc_total(p) + npp_fix_acc_total(p))/dt + npp_Nretrans(p)
             soilc_change(p)           = soilc_change(p) + burned_off_carbon / dt                 
             npp_burnedoff(p)       = burned_off_carbon/dt          
             npp_Nuptake(p)            = soilc_change(p)
             ! how much carbon goes to growth of tissues?  
             npp_growth(p)             = (Nuptake(p)- free_retransn_to_npool(p))*plantCN(p)+(excess_carbon_acc/dt) !does not include gresp, since this is calculated from growth
       
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

             ! Splitting Nitrogen fluxes again into NO3 and NH4

             npp_NO3active(p)        = npp_Nactive(p) * sminfrc(c,j)
             npp_NH4active(p)        = npp_Nactive(p) * (1.0_r8 - sminfrc(c,j))

             npp_NO3nonmyc(p)        = npp_Nnonmyc(p) * sminfrc(c,j)
             npp_NH4nonmyc(p)        = npp_Nnonmyc(p) * (1.0_r8 - sminfrc(c,j))       
            
             Nactive_no3(p)          = Nactive(p) * sminfrc(c,j)         
             Nactive_nh4(p)          = Nactive(p) * * (1.0_r8 - sminfrc(c,j))                  
                               
             Nnonmyc_no3(p)          = Nnonmyc(p) * sminfrc(c,j)
             Nnonmyc_nh4(p)          = Nnonmyc(p) * * (1.0_r8 - sminfrc(c,j))

            
end do pft 
        


end associate
  
end subroutine CNFUNMIMICSplus


subroutine updateCNFUNMIMICSplus (bounds, num_soilc, filter_soilc, cnveg_carbonflux_inst, cnveg_nitrogenflux_inst, &
                                  soilbiogeochem_nitrogenflux_inst, soilbiogeochem_carbonflux_inst)

   use subgridAveMod   , only : p2c
   use clm_varctl      , only : use_nitrif_denitrif

   ! ARGUMENTS
   type(bounds_type)                       , intent(in)    :: bounds
   integer                                 , intent(in)    :: num_soilc             ! number of soil columns in filter
   integer                                 , intent(in)    :: filter_soilc(:)       ! filter for soil columns
   type(cnveg_carbonflux_type)             , intent(inout) :: cnveg_carbonflux_inst
   type(cnveg_nitrogenflux_type)           , intent(inout) :: cnveg_nitrogenflux_inst
   type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst 
   type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst 
       
   call p2c(bounds, num_soilc, filter_soilc,                               &
             cnveg_carbonflux_inst%soilc_change_patch(bounds%begp:bounds%endp), &
             soilbiogeochem_carbonflux_inst%soilc_change_col(bounds%begc:bounds%endc))
             
   call p2c(bounds, num_soilc, filter_soilc,                               &
             cnveg_nitrogenflux_inst%Nfix_patch(bounds%begp:bounds%endp), &
             soilbiogeochem_nitrogenflux_inst%nfix_to_sminn_col(bounds%begc:bounds%endc))


   if (use_nitrif_denitrif) then
      call p2c(bounds,nlevdecomp, &
               cnveg_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
               soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')

      call p2c(bounds,nlevdecomp, &
               cnveg_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
               soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_col(bounds%begc:bounds%endc,1:nlevdecomp),&
               'unity')
   else

      call p2c(bounds, nlevdecomp, &
               cnveg_nitrogenflux_inst%sminn_to_plant_fun_vr_patch(bounds%begp:bounds%endp,1:nlevdecomp),&
               soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_vr_col(bounds%begc:bounds%endc,1:nlevdecomp), &
               'unity')
   endif

   end subroutine updateCNFUNMIMICSplus




!=========================================================================================
real(r8) function fun_cost_fix(fixer,a_fix,b_fix,c_fix,big_cost,crootfr,s_fix, tc_soisno)

! Description:
!   Calculate the cost of fixing N by nodules.
! Code Description:
!   This code is written to CLM4CN by Mingjie Shi on 06/27/2013

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
  real(r8), intent(in) :: tc_soisno ! soil temperature (degrees Celsius)

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
     fun_cost_fix  = (-1*s_fix) * 1.0_r8 / (1.25_r8* (exp(a_fix + b_fix * tc_soisno * (1._r8 - 0.5_r8 * tc_soisno / c_fix)) ))
  else
     fun_cost_fix = big_cost
  end if    ! ends up with the fixer or non-fixer decision
  
  end function fun_cost_fix
!=========================================================================================
  real(r8) function fun_cost_active(sminn_layer,big_cost,kc_active,kn_active,rootc_dens,crootfr,smallValue)         

! Description:
!    Calculate the cost of active uptake of N frm the soil.
! Code Description:
!   This code is written to CLM4 by Mingjie Shi.

  implicit none
!--------------------------------------------------------------------------
! Function result.
!--------------------------------------------------------------------------
  real(r8), intent(in) :: sminn_layer   !  Amount of N (as NH4 or NO3) in the soil that is available to plants (gN/m2).
  real(r8), intent(in) :: big_cost      !  An arbitrary large cost (gC/gN).
  real(r8), intent(in) :: kc_active     !  Constant for cost of active uptake (gC/m2).
  real(r8), intent(in) :: kn_active     !  Constant for cost of active uptake (gC/m2).
  real(r8), intent(in) :: rootc_dens    !  Root carbon density in layer (gC/m3).
  real(r8), intent(in) :: crootfr        !  Fraction of roots that are in this layer.
  real(r8), intent(in) :: smallValue    !  A small number.

  if (rootc_dens > 1.e-6_r8.and.sminn_layer > smallValue) then
     fun_cost_active =  kn_active/sminn_layer + kc_active/rootc_dens 
  else
!    There are very few roots in this layer. Set a high cost.
     fun_cost_active =  big_cost
  end if
 
  end function fun_cost_active
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
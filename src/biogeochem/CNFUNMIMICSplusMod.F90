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
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  !
  ! TYPES


  ! !PUBLIC MEMBER FUNCTIONS:
  public:: readParams            ! Read in parameters needed for FUNMimicplus
  public:: CNFUNMIMICSplusInit   ! FUNMIMICSplus calculation initialization
  public:: CNFUNMIMICSplus      ! FUNMIMICSplus calculation itself

  ! !LOCAL VARIABLES
  integer,  parameter :: icostFix        = 1             ! Process number for fixing.
  integer,  parameter :: icostRetrans    = 2             ! Process number for retranslocation.
  integer,  parameter :: icostActiveNO3  = 3             ! Process number for mycorrhizal uptake of NO3.
  integer,  parameter :: icostActiveNH4  = 4             ! Process number for mycorrhizal uptake of NH4
  integer,  parameter :: icostnonmyc_no3 = 5             ! Process number for nonmyc uptake of NO3.
  integer,  parameter :: icostnonmyc_nh4 = 6             ! Process number for nonmyc uptake of NH4.
  
  real(r8), parameter :: big_cost        = 1000000000._r8! An arbitrarily large cost
 
  !  array index when plant is fixing
  integer, parameter :: plants_are_fixing = 1
  integer, parameter :: plants_not_fixing = 2
 
  !  array index for ECM step versus AM step
  integer, parameter :: ecm_step          = 1
  integer, parameter :: am_step           = 2
  !  arbitrary large cost (gC/gN).

  integer,  private, parameter :: nstp            = 2             ! Number of calculation part
  integer,  private, parameter :: ncost6          = 6             ! Number of  N transport pathways



  ! DATA STRUCTURES

  type, private :: params_type
   real(r8) :: ndays_on        ! number of days to complete leaf onset
   real(r8) :: ndays_off       ! number of days to complete leaf offset
  end type params_type   

   
  !
  type(params_type), private :: params_inst  ! params_inst is
  !  populated in readParamsMod
  
  type, public :: cnfunmimicsplus_type
  ! declare arrays and scalars here with pointers

  real(r8)           :: excess                                                ! excess N taken up by transpiration    (gN/m2) 
  real(r8)           :: steppday                                              ! model time steps in each day          (-)
  real(r8)           :: rootc_dens_step                                       ! root C for each PFT in each soil layer(gC/m2)
  real(r8)           :: retrans_limit1                                        ! a temporary variable for leafn        (gN/m2)
  real(r8)           :: qflx_tran_veg_layer                                   ! transpiration in each soil layer      (mm H2O/S)
  real(r8)           :: dn                                                    ! Increment of N                        (gN/m2)  
  real(r8)           :: dn_retrans                                            ! Increment of N                        (gN/m2)
  real(r8)           :: dnpp                                                  ! Increment of NPP                      (gC/m2)
  real(r8)           :: dnpp_retrans                                          ! Increment of NPP                      (gC/m2)
  real(r8), pointer  :: rootc_dens(:,:)                                       ! the root carbon density               (gC/m2)
  real(r8), pointer  :: rootC(:)                                              ! root biomass                          (gC/m2)
  real(r8), pointer  :: permyc(:,:)                                           ! the arrary for the ECM and AM ratio   (-) 
  real(r8), pointer  :: kc_active(:,:)                                        ! the kc_active parameter               (gC/m2)
  real(r8), pointer  :: kn_active(:,:)                                        ! the kn_active parameter               (gC/m2)
  real(r8), pointer  :: availc_pool(:)                                        ! The avaible C pool for allocation     (gC/m2)
  real(r8), pointer :: plantN(:)                                              ! Plant N (gN/m2)
  real(r8), pointer :: plant_ndemand_pool(:)                                  ! The N demand pool (gN/m2)
  real(r8), pointer :: plant_ndemand_pool_step(:,:)                           ! the N demand pool (gN/m2)
  real(r8), pointer :: leafn_step(:,:)                                        ! N loss based for deciduous trees (gN/m2)
  real(r8), pointer :: leafn_retrans_step(:,:)                                ! N loss based for deciduous trees (gN/m2) 
  real(r8), pointer :: litterfall_n(:)                                        ! N loss based on the leafc to litter (gN/m2) 
  real(r8), pointer :: litterfall_n_step(:,:)                                 ! N loss based on the leafc to litter (gN/m2)
  real(r8), pointer :: litterfall_c_step(:,:)                                 ! N loss based on the leafc to litter (gN/m2)
  real(r8), pointer :: tc_soisno(:,:)                                         ! Soil temperature (degrees Celsius)
  real(r8), pointer :: npp_remaining(:,:)                                     ! A temporary variable for npp_remaining(gC/m2) 
  real(r8), pointer :: n_passive_step(:,:)                                    ! N taken up by transpiration at substep(gN/m2)
  real(r8), pointer :: n_passive_acc(:)                                       ! N acquired by passive uptake (gN/m2)
  real(r8), pointer :: cost_retran(:,:)                                       ! cost of retran (gC/gN)
  real(r8), pointer :: cost_fix(:,:)                                          ! cost of fixation (gC/gN)
  real(r8), pointer :: cost_resis(:,:)                                        ! cost of resis (gC/gN)
  real(r8), pointer :: cost_res_resis(:,:)                                    ! The cost of resis (gN/gC)
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

  !----------NITRIF_DENITRIF-------------!
  !REPLACE NO# with N
  real(r8)            :: sminn_diff                                        ! A temporary limit for N uptake                  (gN/m2)
  real(r8)            :: active_limit1                                     ! A temporary limit for N uptake                  (gN/m2)
  real(r8),  pointer  :: cost_active(:,:) ! cost of mycorrhizal                             (gC/gN)
  real(r8),  pointer  :: cost_nonmyc(:,:) ! cost of nonmyc                                  (gC/gN)

  real(r8),  pointer  :: sminn_conc(:,:)             ! Concentration of N in soil water (gN/gH2O)
  real(r8),  pointer  :: sminn_conc_step(:,:,:)      ! A temporary variable for soil mineral N (gN/gH2O)
  real(r8),  pointer  :: sminn_layer(:,:)            ! Available N in each soil layer (gN/m2)
  real(r8),  pointer  :: sminn_layer_step(:,:,:)     ! A temporary variable for soil N (gN/m2) 
  real(r8),  pointer  :: sminn_uptake(:,:,:)         ! A temporary variable for soil mineral N (gN/m2/s)
  
  real(r8),  pointer  :: active_uptake1(:,:)         ! N mycorrhizal uptake (gN/m2)
  real(r8),  pointer  :: nonmyc_uptake1(:,:)         ! N non-mycorrhizal uptake (gN/m2) 
  real(r8),  pointer  :: active_uptake2(:,:)         ! N mycorrhizal uptake (gN/m2) 
  real(r8),  pointer  :: nonmyc_uptake2(:,:)         ! N non-mycorrhizal uptake (gN/m2) 
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
  real(r8),  pointer  :: n_active_retrans(:,:)             ! Mycorrhizal N for offset (gN/m2)
  real(r8),  pointer  :: n_nonmyc_retrans(:,:)             ! Non-myc     N for offset (gN/m2)
  real(r8),  pointer  :: n_active_retrans_total(:)              ! Mycorrhizal N for offset (gN/m2)
  real(r8),  pointer  :: n_nonmyc_retrans_total(:)              ! Non-myc     N for offset (gN/m2)
  real(r8),  pointer  :: n_passive_vr(:,:)           ! Layer passive N uptake (gN/m2)
  real(r8),  pointer  :: n_fix_vr(:,:)               ! Layer fixation N uptake (gN/m2)
  real(r8),  pointer  :: n_active_vr(:,:)            ! Layer mycorrhizal N uptake (gN/m2)
  real(r8),  pointer  :: n_nonmyc_vr(:,:)            ! Layer non-myc     N uptake (gN/m2)
  real(r8),  pointer  :: npp_active_retrans(:,:)           ! Mycorrhizal N uptake used C for offset (gN/m2)
  real(r8),  pointer  :: npp_nonmyc_retrans(:,:)           ! Non-myc N uptake used C for offset (gN/m2)
  real(r8),  pointer  :: npp_active_retrans_total(:)            ! Mycorrhizal N uptake used C for offset (gN/m2)
  real(r8),  pointer  :: npp_nonmyc_retrans_total(:)             ! Non-myc N uptake used C for offset (gN/m2)
  

  real(r8),  pointer  :: costNit(:,:)                          ! Cost of N via each process                      (gC/gN)

  ! Uptake fluxes for COST_METHOD=2
  ! actual npp to each layer for each uptake process
  real(r8),  pointer  ::                   npp_to_fixation(:) 
  real(r8),  pointer  ::                   npp_to_retrans(:)
  real(r8),  pointer  ::                   npp_to_active(:)
  real(r8),  pointer  ::                   npp_to_nonmyc(:)  

  ! fraction of carbon to each uptake process 
  real(r8),  pointer  ::                   npp_frac_to_fixation(:) 
  real(r8),  pointer  ::                   npp_frac_to_retrans(:)
  real(r8),  pointer  ::                   npp_frac_to_active(:)
  real(r8),  pointer  ::                   npp_frac_to_nonmyc (:)  
   
  ! hypothetical fluxes on N in each layer 
  real(r8),  pointer  ::                  n_exch_fixation(:)        ! N aquired from one unit of C for fixation (unitless)
  real(r8),  pointer  ::                  n_exch_retrans(:)         ! N aquired from one unit of C for retrans (unitless)
  real(r8),  pointer  ::                  n_exch_active(:)      ! N aquired from one unit of C for act no3 (unitless)
  real(r8),  pointer  ::                  n_exch_nonmyc(:)      ! N aquired from one unit of C for nonmyc no3 (unitless) 
  
   !actual fluxes of N in each layer
  real(r8),  pointer  ::                  n_from_fixation(:)        ! N aquired in each layer for fixation       (gN m-2 s-)
  real(r8),  pointer  ::                  n_from_retrans(:)         ! N aquired in each layer of C for retrans (gN m-2 s-)
  real(r8),  pointer  ::                  n_from_active(:)      ! N aquired in each layer of C for act no3 (gN m-2 s-)
  real(r8),  pointer  ::                  n_from_nonmyc(:)      ! N aquired in each layer of C for nonmyc no3 (gN m-2 s-) 

  real(r8),  pointer  ::                  free_Nretrans(:)          ! the total amount of NO3 and NH4                 (gN/m3/s)


  ! Uptake fluxes for COST_METHOD=2
  ! actual fluxes of N in each layer
  real(r8)  ::                  frac_ideal_C_use      ! How much less C do we use for 'buying' N than that needed to get to the ideal ratio?  fraction. 
  
  real(r8)  ::                  N_acquired
  real(r8)  ::                  C_spent
  real(r8)  ::                  leaf_narea            ! leaf n per unit leaf area in gN/m2 (averaged across canopy, which is OK for the cost calculation)
  
  real(r8)  ::                  sum_n_acquired        ! Sum N aquired from one unit of C (unitless)  
  real(r8)  ::                  burned_off_carbon     ! carbon wasted by poor allocation algorithm. If this is too big, we need a better iteration. 
  real(r8)   ::                 temp_n_flux  
  real(r8)  ::                  delta_cn              ! difference between 'ideal' leaf CN ration and actual leaf C:N ratio. C/N
  real(r8) :: excess_carbon        ! how much carbon goes into the leaf C pool on account of the flexibleCN modifications.   
  real(r8) :: excess_carbon_acc    ! excess accumulated over layers.
  
  !  WITHOUT GROWTH RESP
  real(r8) :: fixerfrac            ! what fraction of plants can fix?
  real(r8) :: npp_to_spend         ! how much carbon do we need to get rid of? 
  real(r8) :: soil_n_extraction    ! calculates total N pullled from soil
  real(r8) :: total_N_conductance  !inverse of C to of N for whole soil-leaf pathway
  real(r8) :: total_N_resistance   ! C to of N for whole soil -leaf pathway
  real(r8) :: free_RT_frac=0.0_r8  !fraction of N retranslocation which is automatic/free.
  !  SHould be made into a PFT parameter. 
  
  real(r8) :: paid_for_n_retrans
  real(r8) :: free_n_retrans
  real(r8) :: total_c_spent_retrans
  real(r8) :: total_c_accounted_retrans
  
  
  !------end of not_use_nitrif_denitrif------!
  !--------------------------------------------------------------------
  !------------

  


contains
  procedure , public  :: Init
  procedure , private  :: InitAllocate
  procedure , private :: SetZeros
   
end type cnfunmimicsplus_type


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

     allocate(this%kc_active(bounds%begp:bounds%endp,1:nstp));               this%kc_active(:,:) = nan

     allocate(this%kn_active(bounds%begp:bounds%endp,1:nstp));               this%kn_active(:,:) = nan

     allocate(this%availc_pool(bounds%begp:bounds%endp));                    this%availc_pool(:) = nan

     allocate(this%plantN(bounds%begp:bounds%endp));                         this%plantN(:) = nan

     allocate(this%plant_ndemand_pool(bounds%begp:bounds%endp));             this%plant_ndemand_pool(:) = nan

     allocate(this%plant_ndemand_pool_step(bounds%begp:bounds%endp,1:nstp)); this%plant_ndemand_pool_step(:,:) = nan
     
     allocate(this%leafn_step(bounds%begp:bounds%endp,1:nstp));              this%leafn_step(:,:) = nan
     
     allocate(this%leafn_retrans_step(bounds%begp:bounds%endp,1:nstp));      this%leafn_retrans_step(:,:) = nan
     
     allocate(this%litterfall_n(bounds%begp:bounds%endp));                   this%litterfall_n(:) = nan
     
     allocate(this%litterfall_n_step(bounds%begp:bounds%endp,1:nstp));       this%litterfall_n_step(:,:) = nan
     
     allocate(this%litterfall_c_step(bounds%begp:bounds%endp,1:nstp));       this%litterfall_c_step(:,:) = nan
     
     allocate(this%tc_soisno(bounds%begc:bounds%endc,1:nlevdecomp));         this%tc_soisno(:,:) = nan
     
     allocate(this%npp_remaining(bounds%begp:bounds%endp,1:nstp));           this%npp_remaining(:,:) = nan
     
     allocate(this%n_passive_step(bounds%begp:bounds%endp,1:nstp));          this%n_passive_step(:,:) = nan
     
     allocate(this%n_passive_acc(bounds%begp:bounds%endp));                  this%n_passive_acc(:) = nan
     
     allocate(this%cost_retran(bounds%begp:bounds%endp,1:nlevdecomp));       this%cost_retran(:,:) = nan
     
     allocate(this%cost_fix(bounds%begp:bounds%endp,1:nlevdecomp));          this%cost_fix(:,:) = nan
     
     allocate(this%cost_resis(bounds%begp:bounds%endp,1:nlevdecomp));        this%cost_resis(:,:) = nan
     
     allocate(this%cost_res_resis(bounds%begp:bounds%endp,1:nlevdecomp));    this%cost_res_resis(:,:) = nan
     
     allocate(this%n_fix_acc(bounds%begp:bounds%endp,1:nstp));               this%n_fix_acc(:,:) = nan
     
     allocate(this%n_fix_acc_total(bounds%begp:bounds%endp));                this%n_fix_acc_total(:) = nan
     
     allocate(this%npp_fix_acc(bounds%begp:bounds%endp,1:nstp));             this%npp_fix_acc(:,:) = nan
     
     allocate(this%npp_fix_acc_total(bounds%begp:bounds%endp));              this%npp_fix_acc_total(:) = nan
     
     allocate(this%n_retrans_acc(bounds%begp:bounds%endp,1:nstp));           this%n_retrans_acc(:,:) = nan
     
     allocate(this%n_retrans_acc_total(bounds%begp:bounds%endp));            this%n_retrans_acc_total(:) = nan
     
     allocate(this%free_nretrans_acc(bounds%begp:bounds%endp,1:nstp));       this%free_nretrans_acc(:,:) = nan
     
     allocate(this%npp_retrans_acc(bounds%begp:bounds%endp,1:nstp));         this%npp_retrans_acc(:,:) = nan
     
     allocate(this%npp_retrans_acc_total(bounds%begp:bounds%endp));          this%npp_retrans_acc_total(:) = nan
     
     allocate(this%nt_uptake(bounds%begp:bounds%endp,1:nstp));               this%nt_uptake(:,:) = nan
     
     allocate(this%npp_uptake(bounds%begp:bounds%endp,1:nstp));              this%npp_uptake(:,:) = nan
     
     allocate(this%cost_active(bounds%begp:bounds%endp,1:nlevdecomp));       this%cost_active(:,:) = nan

     allocate(this%cost_nonmyc(bounds%begp:bounds%endp,1:nlevdecomp));       this%cost_nonmyc(:,:) = nan

     allocate(this%sminn_conc(bounds%begc:bounds%endc,1:nlevdecomp));        this%sminn_conc(:,:) = nan

     allocate(this%sminn_conc_step(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp));  this%sminn_conc_step(:,:,:) = nan

     allocate(this%sminn_layer(bounds%begc:bounds%endc,1:nlevdecomp));             this%sminn_layer(:,:) = nan

     allocate(this%sminn_layer_step(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp)); this%sminn_layer_step(:,:,:) = nan

     allocate(this%sminn_uptake(bounds%begp:bounds%endp,1:nlevdecomp,1:nstp));     this%sminn_uptake(:,:,:) = nan

     allocate(this%active_uptake1(bounds%begp:bounds%endp, 1:nlevdecomp));   this%active_uptake1 = nan

     allocate(this%nonmyc_uptake1(bounds%begp:bounds%endp, 1:nlevdecomp));   this%nonmyc_uptake1 = nan

     allocate(this%active_uptake2(bounds%begp:bounds%endp, 1:nlevdecomp));   this%active_uptake2 = nan

     allocate(this%nonmyc_uptake2(bounds%begp:bounds%endp, 1:nlevdecomp));   this%nonmyc_uptake2 = nan

     allocate(this%n_am_acc(bounds%begp:bounds%endp));                       this%n_am_acc = nan

     allocate(this%n_ecm_acc(bounds%begp:bounds%endp));                      this%n_ecm_acc = nan

     allocate(this%n_active_acc(bounds%begp:bounds%endp, 1:nstp));           this%n_active_acc = nan

     allocate(this%n_nonmyc_acc(bounds%begp:bounds%endp, 1:nstp));           this%n_nonmyc_acc = nan

     allocate(this%n_active_acc_total(bounds%begp:bounds%endp));             this%n_active_acc_total = nan

     allocate(this%n_nonmyc_acc_total(bounds%begp:bounds%endp));             this%n_nonmyc_acc_total(:) = nan

     allocate(this%npp_active_acc(bounds%begp:bounds%endp, 1:nstp));         this%npp_active_acc(:,:) = nan

     allocate(this%npp_nonmyc_acc(bounds%begp:bounds%endp, 1:nstp));         this%npp_nonmyc_acc(:,:) = nan

     allocate(this%npp_active_acc_total(bounds%begp:bounds%endp));           this%npp_active_acc_total(:) = nan

     allocate(this%npp_nonmyc_acc_total(bounds%begp:bounds%endp));           this%npp_nonmyc_acc_total(:) = nan

     allocate(this%n_am_retrans(bounds%begp:bounds%endp));                   this%n_am_retrans(:) = nan

     allocate(this%n_ecm_retrans(bounds%begp:bounds%endp));                  this%n_ecm_retrans(:) = nan

     allocate(this%n_active_retrans(bounds%begp:bounds%endp, 1:nstp));       this%n_active_retrans(:,:) = nan

     allocate(this%n_nonmyc_retrans(bounds%begp:bounds%endp, 1:nstp));       this%n_nonmyc_retrans(:,:) = nan

     allocate(this%n_active_retrans_total(bounds%begp:bounds%endp));         this%n_active_retrans_total(:) = nan

     allocate(this%n_nonmyc_retrans_total(bounds%begp:bounds%endp));         this%n_nonmyc_retrans_total(:) = nan

     allocate(this%n_passive_vr(bounds%begp:bounds%endp, 1:nlevdecomp));     this%n_passive_vr(:,:) = nan

     allocate(this%n_fix_vr(bounds%begp:bounds%endp, 1:nlevdecomp));         this%n_fix_vr(:,:) = nan

     allocate(this%n_active_vr(bounds%begp:bounds%endp, 1:nlevdecomp));      this%n_active_vr(:,:) = nan

     allocate(this%n_nonmyc_vr(bounds%begp:bounds%endp, 1:nlevdecomp));      this%n_nonmyc_vr(:,:) = nan

     allocate(this%npp_active_retrans(bounds%begp:bounds%endp, 1:nstp));     this%npp_active_retrans(:,:) = nan

     allocate(this%npp_nonmyc_retrans(bounds%begp:bounds%endp, 1:nstp));     this%npp_nonmyc_retrans(:,:) = nan

     allocate(this%npp_active_retrans_total(bounds%begp:bounds%endp));       this%npp_active_retrans_total(:) = nan

     allocate(this%npp_nonmyc_retrans_total(bounds%begp:bounds%endp));       this%npp_nonmyc_retrans_total(:) = nan

     allocate(this%costNit(1:nlevdecomp,ncost6));                            this%costNit(:,:) = nan

    ! Uptake fluxes for COST_METHOD=2
  ! actual npp to each layer for each uptake process
     allocate(this%npp_to_fixation(1:nlevdecomp));                           this%npp_to_fixation(:) = nan

     allocate(this%npp_to_retrans(1:nlevdecomp));                            this%npp_to_retrans(:) = nan

     allocate(this%npp_to_active(1:nlevdecomp));                             this%npp_to_active(:) = nan

     allocate(this%npp_to_nonmyc (1:nlevdecomp));                            this%npp_to_nonmyc(:) = nan  
    
     ! fraction of carbon to each uptake process 
     allocate(this%npp_frac_to_fixation(1:nlevdecomp));                      this%npp_frac_to_fixation(:) = nan 

     allocate(this%npp_frac_to_retrans(1:nlevdecomp));                       this%npp_frac_to_retrans(:) = nan

     allocate(this%npp_frac_to_active(1:nlevdecomp));                        this%npp_frac_to_active(:) = nan

     allocate(this%npp_frac_to_nonmyc (1:nlevdecomp));                       this%npp_frac_to_nonmyc(:) = nan  
      
     ! hypothetical fluxes on N in each layer 
     allocate(this%n_exch_fixation(1:nlevdecomp));                           this%n_exch_fixation(:) = nan

     allocate(this%n_exch_retrans(1:nlevdecomp));                            this%n_exch_retrans(:) = nan

     allocate(this%n_exch_active(1:nlevdecomp));                             this%n_exch_active(:) = nan

     allocate(this%n_exch_nonmyc(1:nlevdecomp));                             this%n_exch_nonmyc(:) = nan
     
      !actual fluxes of N in each layer
     allocate(this%n_from_fixation(1:nlevdecomp));                           this%n_from_fixation(:) = nan

     allocate(this%n_from_retrans(1:nlevdecomp));                            this%n_from_retrans(:) = nan

     allocate(this%n_from_active(1:nlevdecomp));                             this%n_from_active(:) = nan

     allocate(this%n_from_nonmyc(1:nlevdecomp));                             this%n_from_nonmyc(:) = nan
    
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
   
     !local_use_flexibleCN                = use_flexibleCN
     !steppday                            = 48._r8
     !big_cost
     !qflx_tran_veg_layer                 = 0._r8
     !rootc_dens_step                     = 0._r8
     !plant_ndemand_pool                  = 0._r8

     this%rootc_dens(:,:)                = 0._r8
     this%rootC(:)                       = 0._r8
     this%kc_active(:,:)                 = 0._r8
     this%kn_active(:,:)                 = 0._r8
     this%availc_pool(:)                 = 0._r8
     this%plantN(:)                      = 0._r8
     this%plant_ndemand_pool(:)          = 0._r8
     this%plant_ndemand_pool_step(:,:)   = 0._r8
     this%leafn_step(:,:)                = 0._r8
     this%leafn_retrans_step(:,:)        = 0._r8
     this%litterfall_n(:)                = 0._r8
     this%litterfall_n_step(:,:)         = 0._r8
     this%litterfall_c_step(:,:)         = 0._r8
     this%tc_soisno(:,:)                 = 0._r8 
     this%npp_remaining(:,:)             = 0._r8
     this%n_passive_step(:,:)            = 0._r8
     this%n_passive_acc(:)               = 0._r8
     this%cost_retran(:,:)               = 0._r8
     this%cost_fix(:,:)                  = 0._r8
     this%cost_resis(:,:)                = 0._r8
     this%cost_res_resis(:,:)            = 0._r8
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
     this%sminn_conc(:,:)                = 0._r8
     this%sminn_conc_step(:,:,:)         = 0._r8
     this%sminn_layer(:,:)               = 0._r8
     this%sminn_layer_step(:,:,:)        = 0._r8
     this%sminn_uptake(:,:,:)            = 0._r8
     this%active_uptake1                 = 0._r8
     this%nonmyc_uptake1                 = 0._r8
     this%active_uptake2                 = 0._r8
     this%nonmyc_uptake2                 = 0._r8
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
     this%n_active_retrans(:,:)          = 0._r8
     this%n_nonmyc_retrans(:,:)          = 0._r8
     this%n_active_retrans_total(:)      = 0._r8
     this%n_nonmyc_retrans_total(:)      = 0._r8
     this%n_passive_vr(:,:)              = 0._r8
     this%n_fix_vr(:,:)                  = 0._r8
     this%n_active_vr(:,:)               = 0._r8
     this%n_nonmyc_vr(:,:)               = 0._r8
     this%npp_active_retrans(:,:)        = 0._r8
     this%npp_nonmyc_retrans(:,:)        = 0._r8
     this%npp_active_retrans_total(:)    = 0._r8
     this%npp_nonmyc_retrans_total(:)    = 0._r8

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
  !
  ! !ARGUMENTS:
  type(bounds_type)             , intent(in)    :: bounds
  type(cnveg_state_type)        , intent(inout) :: cnveg_state_inst
  type(cnveg_carbonstate_type)  , intent(inout) :: cnveg_carbonstate_inst
  type(cnveg_nitrogenstate_type), intent(inout) :: cnveg_nitrogenstate_inst
  !
  ! !LOCAL VARIABLES:         
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
!--------------------------------------------------------------------
  !---

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
   soilbiogeochem_carbonstate_inst)
   
! !USES:
   use clm_time_manager, only : get_step_size_real, get_curr_date
   use clm_varpar      , only : nlevdecomp
   use clm_varcon      , only : secspday, smallValue, fun_period, tfrz, dzsoi_decomp, spval
   use clm_varctl      , only : use_nitrif_denitrif
   use PatchType       , only : patch
   use subgridAveMod   , only : p2c
   use pftconMod       , only : npcropmin
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
  ! 

end subroutine CNFUNMIMICSplus


end module CNFUNMIMICSplusMod
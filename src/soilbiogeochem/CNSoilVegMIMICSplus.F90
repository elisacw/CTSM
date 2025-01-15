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
  public :: readParams                       ! Read in parameters from params file

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: init_mimicsplus_veg_myc         !
  private :: CNallocation                    !
  private :: calc_roi                        !

  type, public :: params_type
     real(r8) :: mimicsplus_v_nh4         !Maximum NH4+ immobilization rate [year-1]
     real(r8) :: mimicsplus_v_no3         !Maximum NO3- immobilization rate [year-1]
     real(r8) :: mimicsplus_vmax_denit    !Maximum denitrification decomposition rate at reference temperature [year-1]
     real(r8) :: mimicsplus_fden          !Maximum denitrification decomposition rate at reference temperature [unitless]
     real(r8) :: mimicsplus_kdenit        !Half-saturation constant for nitrate concentration in denitrification [kg NO3-N kg NO3-N demand-1 year-1]
     real(r8) :: mimicsplus_root_no3      !Maximum root active nitrate uptake rate [kg N m-3 year-1]
     real(r8) :: mimicsplus_root_nh4      !Maximum root active ammonium uptake rate [kg N m-3 year-1]
     real(r8) :: mimicsplus_km_no3        !Half-saturation nitrate concentration for root active uptake [kg N m-3]
     real(r8) :: mimicsplus_km_nh4        !Half-saturation nitrate concentration for root active uptake [kg N m-3]
     real(r8) :: mimicsplus_r_rhiz        !Radius of the rhizosphere [m]
     real(r8) :: mimicsplus_v_scav        !Maximum N uptake rate by scavenging mycorrhizae [kg N m-3 year-1]
     real(r8) :: mimicsplus_k_scav_inorg  !Half-saturation inorganic N concentration for mycorrhizal uptake [kg N m-3]
     real(r8) :: mimicsplus_k_scav        !Half-saturation mycorrhizal biomass concentration for scavenging [kg C m-3]
     real(r8) :: mimicsplus_km_mine       !Half-saturation mycorrhizal biomass concentration for scavenging [kg C m-3]
     real(r8) :: mimicsplus_cue_mine      !Carbon use efficiency of mycorrhizal mining [fraction]
     real(r8) :: mimicsplus_nue_mine      !Nitrogen use efficiency of mycorrhizal mining [fraction]
     real(r8) :: mimicsplus_vmax_ref_mine !Maximum decomposition rate at reference temperature for mycorrhizal mining [year-1]
     real(r8) :: mimicsplus_rfix          !N fixation rate per unit symbiotic biomass [kg N kg biomass C-1 year-1]
     real(r8) :: mimicsplus_kgrowth       !Half-saturation of intermediate C pool for symbiotic growth [kg C m -2]
     real(r8) :: mimicsplus_rgrowth       !Maximum symbiont growth rate [kg C m-2 year-1]
     real(r8) :: mimicsplus_tau_sym       !Fraction of symbiotic biomass turnover not used for maintenance respiration [fraction]
     real(r8) :: mimicsplus_growth_fix    !N fixer growth efficiency [unitless]
     real(r8) :: mimicsplus_tau_fix       !N fixer turnover time [year-1]
     real(r8) :: mimicsplus_cn_fix        !N fixer C:N [unitless]
     real(r8) :: mimicsplus_tau_int       !Turnover time of intermediate C pool [year-1]
     real(r8) :: mimicsplus_rup_veg       !Vegetation N uptake rate from intermediate N pool [year-1]
     real(r8) :: mimicsplus_fnalloc       !Fraction of NPP allocated to N uptake per unit N stress [fraction]
   contains
     procedure, private :: allocParams    ! Allocate the parameters
     procedure, private :: cleanParams    ! Deallocate parameters from member
  end type params_type
  !
  type(params_type), public, protected :: params_inst  ! params_inst is populated in readParamsMod 

 
  type, public :: symbiont_type

  ! real(r8), pointer, private :: ac_phs_patch      (:,:,:) ! patch Rubisco-limited gross photosynthesis (umol CO2/m**2/s)
  
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

    class(photosyns_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate (bounds)
    call this%InitHistory  (bounds)
    call this%InitCold     (bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
   !
   ! !ARGUMENTS:
   class(photosyns_type) :: this
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
      class(photosyns_type) :: this
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
     class(photosyns_type) :: this
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

   !-----------------------------------------------------------------------
   subroutine allocParams ( this )
     !
     implicit none
  
     ! !ARGUMENTS:
     class(symbiont_type) :: this
     !
     ! !LOCAL VARIABLES:
     character(len=32)  :: subname = 'allocParams'
     !-----------------------------------------------------------------------
  
     ! allocate parameters
  
     !allocate( this%krmax       (0:mxpft) )          ; this%krmax(:)        = nan
     
   end subroutine allocParams

   !-----------------------------------------------------------------------
   subroutine cleanParams ( this )
     !
     implicit none
  
     ! !ARGUMENTS:
     class(symbiont_type) :: this
     !
     ! !LOCAL VARIABLES:
     character(len=32)  :: subname = 'cleanParams'
     !-----------------------------------------------------------------------
  
     ! deallocate parameters
  
     !deallocate( this%krmax       )
     
  end subroutine cleanParams



    !-----------------------------------------------------------------------
    subroutine readParams ( this, ncid )
      !
      ! !USES:
      use ncdio_pio ,   only : file_desc_t,ncd_io
      use paramUtilMod, only: readNcdioScalar
      implicit none
  
      ! !ARGUMENTS:
      class(symbiont_type) :: this
      type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
      !
      ! !LOCAL VARIABLES:
      character(len=32)  :: subname = 'readParams'
      character(len=100) :: errCode = '-Error reading in parameters file:'
      logical            :: readv ! has variable been read in or not
      real(r8)           :: temp1d(0:mxpft) ! temporary to read in parameter
      real(r8)           :: temp2d(0:mxpft,nvegwcs) ! temporary to read in parameter
      character(len=100) :: tString ! temp. var for reading
      !-----------------------------------------------------------------------

      ! read in parameters
  
      call params_inst%allocParams()
     
      tString='mimicsplus_v_nh4'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_v_nh4=tempr
      
      tString='mimicsplus_v_no3'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_v_no3=tempr
      
      tString='mimicsplus_vmax_denit'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_vmax_denit=tempr
      
      tString='mimicsplus_fden'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_fden=tempr
      
      tString='mimicsplus_kdenit'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_kdenit=tempr
      
      tString='mimicsplus_root_no3'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_root_no3=tempr
      
      tString='mimicsplus_root_nh4'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_root_nh4=tempr
      
      tString='mimicsplus_km_no3'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_km_no3=tempr
      
      tString='mimicsplus_km_nh4'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_km_nh4=tempr
      
      tString='mimicsplus_r_rhiz'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_r_rhiz=tempr
      
      tString='mimicsplus_v_scav'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_v_scav=tempr
      
      tString='mimicsplus_k_scav_inorg'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_k_scav_inorg=tempr
      
      tString='mimicsplus_k_scav'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_k_scav=tempr
      
      tString='mimicsplus_km_mine'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_km_mine=tempr
      
      tString='mimicsplus_cue_mine'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_cue_mine=tempr
      
      tString='mimicsplus_nue_mine'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_nue_mine=tempr
      
      tString='mimicsplus_vmax_ref_mine'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_vmax_ref_mine=tempr
      
      tString='mimicsplus_rfix'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_rfix=tempr
      
      tString='mimicsplus_kgrowth'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_kgrowth=tempr
      
      tString='mimicsplus_rgrowth'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_rgrowth=tempr
      
      tString='mimicsplus_tau_sym'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_tau_sym=tempr
      
      tString='mimicsplus_growth_fix'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_growth_fix=tempr
      
      tString='mimicsplus_tau_fix'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_tau_fix=tempr
      
      tString='mimicsplus_cn_fix'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_cn_fix=tempr
      
      tString='mimicsplus_tau_int'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_tau_int=tempr
      
      tString='mimicsplus_rup_veg'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_rup_veg=tempr
      
      tString='mimicsplus_fnalloc'
      call ncd_io(trim(tString), tempr, 'read', ncid, readvar=readv)
      if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
      params_inst%mimicsplus_fnalloc=tempr
     
   end subroutine readParams 

   !------------------------------------------------------------------------
   subroutine Restart(this, bounds, ncid, flag)
     !
     ! !USES:
     use ncdio_pio  , only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
     use restUtilMod
     !
     ! !ARGUMENTS:
     class(photosyns_type) :: this
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

      contains
   

   ! !DESCRIPTION:

   ! !USES:
    use clm_varcon, only: pct_to_frac
   !
   ! !ARGUMENTS:
   type(bounds_type)               , intent(in)    :: bounds 
   !
   
   ! !LOCAL VARIABLES

   mimicsplus_vmax_denit_fast = params_inst%mimicsplus_vmax_denit(1)
   mimicsplus_vmax_denit_slow = params_inst%mimicsplus_vmax_denit(2)
   mimicsplus_vmax_denit_necr = params_inst%mimicsplus_vmax_denit(3)

   mimicsplus_cue_mine_fast = params_inst%mimicsplus_cue_mine(1)
   mimicsplus_cue_mine_slow = params_inst%mimicsplus_cue_mine(2)
   mimicsplus_cue_mine_necr = params_inst%mimicsplus_cue_mine(3)

   mimicsplus_nue_mine_fast = params_inst%mimicsplus_nue_mine(1)
   mimicsplus_nue_mine_slow = params_inst%mimicsplus_nue_mine(2)
   mimicsplus_nue_mine_necr = params_inst%mimicsplus_nue_mine(3)

   mimicsplus_vmax_ref_mine_fast = params_inst%mimicsplus_vmax_ref_mine(1)
   mimicsplus_vmax_ref_mine_slow = params_inst%mimicsplus_vmax_ref_mine(2)
   mimicsplus_vmax_ref_mine_necr = params_inst%mimicsplus_vmax_ref_mine(3)

   !ECW get parameters from the mimics module
      associate(                                                                                     &
      mimicsplus_k_myc_som             => params_inst%mimicsplus_k_myc_som      , & ! 
      mimicsplus_k_mo                  => params_inst%mimicsplus_k_mo             & !          
      )
      end associate


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

   subroutine CNallocation (cnveg_nitrogenstate_inst, leaf_prof_patch, froot_prof_patch, croot_prof_patch)
   
   ! !DESCRIPTION:
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   type(cnveg_nitrogenstate_type)       , intent(in)    :: cnveg_nitrogenstate_inst
   real(r8) :: N_stress                                     ! N demand of plant, based on current N amount in plant []

   real(r8)                             , intent(in)    :: leaf_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: froot_prof_patch(bounds%begp:,1:)
   real(r8)                             , intent(in)    :: croot_prof_patch(bounds%begp:,1:) 
   
   associate(  
   leafn                               => cnveg_nitrogenstate_inst%leafn_patch                              , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N                                    
   leafn_storage                       => cnveg_nitrogenstate_inst%leafn_storage_patch                      , & ! Input:  [real(r8) (:)     ]  (gN/m2) leaf N storage                            
   frootn                              => cnveg_nitrogenstate_inst%frootn_patch                             , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N                               
   frootn_storage                      => cnveg_nitrogenstate_inst%frootn_storage_patch                     , & ! Input:  [real(r8) (:)     ]  (gN/m2) fine root N storage                       
   livecrootn                          => cnveg_nitrogenstate_inst%livecrootn_patch                         , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N                        
   livecrootn_storage                  => cnveg_nitrogenstate_inst%livecrootn_storage_patch                 , & ! Input:  [real(r8) (:)     ]  (gN/m2) live coarse root N storage                
   )

   !N_stress = (2(Nleaf + Nroot)-Nstorage) / (Nleaf + Nroot)
   N_stress = (2*(leafn + frootn + livecrootn) - (leafn_storage + frootn_storage + livecrootn_storage)) / &
               (leafn + frootn + livecrootn)

   


   ! Calculate the amount of C transferent to Naquisation
   ! Ctransfer = max(NPP,0)fNallocNstress
   real(r8) :: C_transfer                                   ! Carbon allocated to N acquisition, higher N_stress leads to higher C_transfer

   ! Dynamic allocation of a fraction of NPP to root exudation
   C_transfer = max(NPP,0) * params_inst%mimicsplus_fnalloc * N_stress
      

   end associate
   end subroutine

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

   ! Growth and Turnover of Symbiotic Biomass



end module CNSoilVegMIMICSplus
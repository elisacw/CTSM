module SoilBiogeochemDecompMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used in litter and soil decomposition model
  !
  ! !USES:
  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use decompMod                          , only : bounds_type
  use clm_varpar                         , only : nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools
  use clm_varctl                         , only : use_nitrif_denitrif, use_lch4, iulog
  use clm_varcon                         , only : dzsoi_decomp
  use clm_varpar                         , only : i_phys_som, i_chem_som, i_avl_som, i_ecm_myc, i_am_myc
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con, mimics_decomp, mimicsplus_decomp, decomp_method, use_soil_matrixcn
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemNitrogenStateType    , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenFluxType     , only : soilbiogeochem_nitrogenflux_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams
  public :: SoilBiogeochemDecomp
  !
  type, private :: params_type
     real(r8) :: dnp         !denitrification proportion
  end type params_type
  !
  type(params_type), private ::  params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read parameters
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io
    use abortutils   , only: endrun
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    tString='dnp'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%dnp=tempr 

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemDecomp (bounds, num_bgc_soilc, filter_bgc_soilc,                                &
       soilbiogeochem_state_inst, soilbiogeochem_carbonstate_inst, soilbiogeochem_carbonflux_inst, &
       soilbiogeochem_nitrogenstate_inst, soilbiogeochem_nitrogenflux_inst, &
       cn_decomp_pools, p_decomp_cpool_loss, pmnf_decomp_cascade, &
       p_decomp_npool_to_din, cnfunmimicsplus_inst)
    !
    ! !USES:
    use SoilBiogeochemDecompCascadeConType, only : i_atm
    use CNFUNMIMICSplusMod,                 only : cnfunmimicsplus_type
    use clm_time_manager,                   only : get_step_size_real
    !
    ! !ARGUMENT:
    type(bounds_type)                       , intent(in)    :: bounds   
    integer                                 , intent(in)    :: num_bgc_soilc          ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_bgc_soilc(:)    ! filter for soil columns
    type(soilbiogeochem_state_type)         , intent(inout) :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonstate_type)   , intent(in)    :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    real(r8)                                , intent(inout) :: cn_decomp_pools(bounds%begc:,1:,1:)     ! c:n ratios of applicable pools
    real(r8)                                , intent(inout) :: p_decomp_cpool_loss(bounds%begc:,1:,1:) ! potential C loss from one pool to another
    real(r8)                                , intent(inout) :: pmnf_decomp_cascade(bounds%begc:,1:,1:) ! potential mineral N flux from one pool to another
    real(r8)                                , intent(in)    :: p_decomp_npool_to_din(bounds%begc:,1:,1:) ! potential flux to dissolved inorganic N
    type(cnfunmimicsplus_type)              , intent(inout) :: cnfunmimicsplus_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l,m                                    ! indices
    integer  :: fc                                           ! lake filter column index
    integer  :: begc,endc                                    ! bounds 
    real(r8) :: dt                                           ! timestep

    !  For methane code
    real(r8):: hrsum(bounds%begc:bounds%endc,1:nlevdecomp)  ! sum of HR (gC/m2/s) 
    !-----------------------------------------------------------------------
   
    begc = bounds%begc; endc = bounds%endc
    
    SHR_ASSERT_ALL_FL((ubound(cn_decomp_pools)     == (/endc,nlevdecomp,ndecomp_pools/))               , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(p_decomp_cpool_loss) == (/endc,nlevdecomp,ndecomp_cascade_transitions/)) , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(pmnf_decomp_cascade) == (/endc,nlevdecomp,ndecomp_cascade_transitions/)) , sourcefile, __LINE__)
    SHR_ASSERT_ALL_FL((ubound(p_decomp_npool_to_din) == (/endc,nlevdecomp,ndecomp_cascade_transitions/)) , sourcefile, __LINE__)

    associate(                                                                                                          & 
         cascade_donor_pool               =>    decomp_cascade_con%cascade_donor_pool                                 , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool            =>    decomp_cascade_con%cascade_receiver_pool                              , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step
         floating_cn_ratio_decomp_pools   =>    decomp_cascade_con%floating_cn_ratio_decomp_pools                     , & ! Input:  [logical  (:)     ]  TRUE => pool has fixed C:N ratio                   
         initial_cn_ratio                 =>    decomp_cascade_con%initial_cn_ratio                                   , & ! Input:  [real(r8) (:)     ]  c:n ratio for initialization of pools             

         fpi_vr                           =>    soilbiogeochem_state_inst%fpi_vr_col                                  , & ! Input:  [real(r8) (:,:)   ]  fraction of potential immobilization (no units) 
         rf_decomp_cascade                =>    soilbiogeochem_carbonflux_inst%rf_decomp_cascade_col                  , & ! Input:  [real(r8) (:,:,:) ]  respired fraction in decomposition step (frac)

         decomp_npools_vr                 =>    soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col                , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
         decomp_cpools_vr                 =>    soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col                  , & ! Input:  [real(r8) (:,:,:) ]  (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) c pools

         decomp_cascade_ntransfer_vr      =>    soilbiogeochem_nitrogenflux_inst%decomp_cascade_ntransfer_vr_col      , & ! Output: [real(r8) (:,:,:) ]  vert-res transfer of N from donor to receiver pool along decomp. cascade (gN/m3/s)
         decomp_cascade_sminn_flux_vr     =>    soilbiogeochem_nitrogenflux_inst%decomp_cascade_sminn_flux_vr_col     , & ! Output: [real(r8) (:,:,:) ]  vert-res mineral N flux for transition along decomposition cascade (gN/m3/s)
         potential_immob_vr               =>    soilbiogeochem_nitrogenflux_inst%potential_immob_vr_col               , & ! Output: [real(r8) (:,:)   ]                                                  
         sminn_to_denit_decomp_cascade_vr =>    soilbiogeochem_nitrogenflux_inst%sminn_to_denit_decomp_cascade_vr_col , & ! Output: [real(r8) (:,:,:) ] 
         gross_nmin_vr                    =>    soilbiogeochem_nitrogenflux_inst%gross_nmin_vr_col                    , & ! Input: [real(r8) (:,:)   ]
         net_nmin_vr                      =>    soilbiogeochem_nitrogenflux_inst%net_nmin_vr_col                      , & ! Output: [real(r8) (:,:)   ]                                                  
         gross_nmin                       =>    soilbiogeochem_nitrogenflux_inst%gross_nmin_col                       , & ! Output: [real(r8) (:)     ]  gross rate of N mineralization (gN/m2/s)          
         net_nmin                         =>    soilbiogeochem_nitrogenflux_inst%net_nmin_col                         , & ! Output: [real(r8) (:)     ]  net rate of N mineralization (gN/m2/s)            
         
         w_scalar                         =>    soilbiogeochem_carbonflux_inst%w_scalar_col                           , & ! Input:  [real(r8) (:,:)   ]  fraction by which decomposition is limited by moisture availability
         c_overflow_vr                    =>    soilbiogeochem_carbonflux_inst%c_overflow_vr                          , & ! Input:  [real(r8) (:,:,:) ]  vertically-resolved C rejected by microbes that cannot process it (gC/m3/s)
         decomp_cascade_hr_vr             =>    soilbiogeochem_carbonflux_inst%decomp_cascade_hr_vr_col               , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         decomp_cascade_ctransfer_vr      =>    soilbiogeochem_carbonflux_inst%decomp_cascade_ctransfer_vr_col        , & ! Output: [real(r8) (:,:,:) ]  vertically-resolved het. resp. from decomposing C pools (gC/m3/s)
         phr_vr                           =>    soilbiogeochem_carbonflux_inst%phr_vr_col                             , & ! Input:  [real(r8) (:,:)   ]  potential HR (gC/m3/s)                           
         fphr                             =>    soilbiogeochem_carbonflux_inst%fphr_col                               , & ! Output: [real(r8) (:,:)   ]  fraction of potential SOM + LITTER heterotrophic

         c_am_growth_vr                   => soilbiogeochem_carbonflux_inst%c_am_growth_vr_col                        , & ! Input: [real(r8) (:,:)    ]  vertically resolved C growth flux for AM mycorrhiza (gC/m3/s)
         c_ecm_growth_vr                  => soilbiogeochem_carbonflux_inst%c_ecm_growth_vr_col                       , & ! Input: [real(r8) (:,:)    ]  vertically resolved C growth flux for ECM mycorrhiza (gC/m3/s)
         n_am_growth_vr                   => soilbiogeochem_nitrogenflux_inst%n_am_growth_vr_col                      , & ! Input: [real(r8) (:,:)    ]  vertically resolved N growth flux for AM mycorrhiza (gC/m3/s)
         n_ecm_growth_vr                  => soilbiogeochem_nitrogenflux_inst%n_ecm_growth_vr_col                     , & ! Input: [real(r8) (:,:)    ]  vertically resolved N growth flux for ECM mycorrhiza (gC/m3/s)
         c_ecm_enz_vr                     => soilbiogeochem_carbonflux_inst%c_ecm_enz_vr_col                          , & ! Input: [real(r8) (:,:)    ]  vertically resolved C enzyme flux for ECM mycorrhiza (goes from plant) (gC/m3/s)
         n_somc2ecm_vr                    => soilbiogeochem_nitrogenflux_inst%n_somc2ecm_vr_col(:,:)                  , & ! Input: [real(r8) (:,:)    ]nitrogen mining from ECM mycorrhiza
         n_somp2ecm_vr                    => soilbiogeochem_nitrogenflux_inst%n_somp2ecm_vr_col(:,:)                  , & ! Input: [real(r8) (:,:)    ]nitrogen mining from ECM mycorrhiza
         c_somc2soma_vr                   => soilbiogeochem_carbonflux_inst%c_somc2soma_vr_col(:,:)                   , & ! Input: [real(r8) (:,:) carbon release from mining from somс pool
         c_somp2soma_vr                   => soilbiogeochem_carbonflux_inst%c_somp2soma_vr_col(:,:)                   , & ! Input: [real(r8) (:,:) carbon release from mining from somp pool
         sminno3_to_ecm_vr                => soilbiogeochem_nitrogenflux_inst%sminno3_to_ecm_vr_col(:,:)              , & ! Input: [real(r8) (:,:)No3 flux from soil NO3 to ECM
         sminno3_to_am_vr                 => soilbiogeochem_nitrogenflux_inst%sminno3_to_am_vr_col(:,:)               , & ! Input: [real(r8) (:,:)No3 flux from soil NO3 to AM
         sminnh4_to_ecm_vr                => soilbiogeochem_nitrogenflux_inst%sminnh4_to_ecm_vr_col(:,:)              , & ! Input: [real(r8) (:,:)No3 flux from soil NO3 to ECM
         sminnh4_to_am_vr                 => soilbiogeochem_nitrogenflux_inst%sminnh4_to_am_vr_col(:,:)                 & ! Input: [real(r8) (:,:)No3 flux from soil NO3 to A
          )

      dt = get_step_size_real()
      ! column loop to calculate actual immobilization and decomp rates, following
      ! resolution of plant/heterotroph  competition for mineral N
      dt = get_step_size_real()
      ! calculate c:n ratios of applicable pools !ECW CN ratio calculated based on pool sizes
      do l = 1, ndecomp_pools
         if ( floating_cn_ratio_decomp_pools(l) ) then
            do j = 1,nlevdecomp
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  if ( decomp_npools_vr(c,j,l) > 0._r8 ) then
                     cn_decomp_pools(c,j,l) = decomp_cpools_vr(c,j,l) / decomp_npools_vr(c,j,l)
                  end if
               end do
            end do
         else
            do j = 1,nlevdecomp
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  cn_decomp_pools(c,j,l) = initial_cn_ratio(l)
               end do
            end do
         end if
      end do

      ! column loop to calculate actual immobilization and decomp rates, following
      ! resolution of plant/heterotroph  competition for mineral N
      
      ! upon return from SoilBiogeochemCompetition, the fraction of potential immobilization
      ! has been set (soilbiogeochem_state_inst%fpi_vr_col). now finish the decomp calculations.
      ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
      ! Also calculate denitrification losses as a simple proportion
      ! of mineralization flux.

      do k = 1, ndecomp_cascade_transitions
         do j = 1,nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)

               if (decomp_cpools_vr(c,j,cascade_donor_pool(k)) > 0._r8) then
                  if ( pmnf_decomp_cascade(c,j,k) > 0._r8 ) then
                     p_decomp_cpool_loss(c,j,k) = p_decomp_cpool_loss(c,j,k) * fpi_vr(c,j)
                     pmnf_decomp_cascade(c,j,k) = pmnf_decomp_cascade(c,j,k) * fpi_vr(c,j)
                     if (use_soil_matrixcn)then ! correct only when one transfer from each litter pool
                     end if
                     if (.not. use_nitrif_denitrif) then
                        sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
                     end if
                  else
                     if (.not. use_nitrif_denitrif) then
                        sminn_to_denit_decomp_cascade_vr(c,j,k) = -params_inst%dnp * pmnf_decomp_cascade(c,j,k)
                     end if
                  end if
                  decomp_cascade_hr_vr(c,j,k) = rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
                  decomp_cascade_ctransfer_vr(c,j,k) = (1._r8 - rf_decomp_cascade(c,j,k)) * p_decomp_cpool_loss(c,j,k)
                  if (decomp_method == mimics_decomp) then
                     decomp_cascade_hr_vr(c,j,k) = min( &
                        p_decomp_cpool_loss(c,j,k), &
                        decomp_cascade_hr_vr(c,j,k) + c_overflow_vr(c,j,k))
                     decomp_cascade_ctransfer_vr(c,j,k) = max(0.0_r8, p_decomp_cpool_loss(c,j,k) - decomp_cascade_hr_vr(c,j,k))
                  else if (decomp_method == mimicsplus_decomp) then !ECW maybe .or.
                     decomp_cascade_hr_vr(c,j,k) = min( &
                        p_decomp_cpool_loss(c,j,k), &
                        decomp_cascade_hr_vr(c,j,k) + c_overflow_vr(c,j,k))
                     decomp_cascade_ctransfer_vr(c,j,k) = max(0.0_r8, p_decomp_cpool_loss(c,j,k) - decomp_cascade_hr_vr(c,j,k))
                  end if
                  if (decomp_npools_vr(c,j,cascade_donor_pool(k)) > 0._r8 .and. cascade_receiver_pool(k) /= i_atm) then
                     decomp_cascade_ntransfer_vr(c,j,k) = p_decomp_cpool_loss(c,j,k) / cn_decomp_pools(c,j,cascade_donor_pool(k))
                  else
                     decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
                  endif
                  if ( cascade_receiver_pool(k) /= 0 ) then
                     decomp_cascade_sminn_flux_vr(c,j,k) = pmnf_decomp_cascade(c,j,k)
                  else  ! keep sign convention negative for terminal pools
                     decomp_cascade_sminn_flux_vr(c,j,k) = - pmnf_decomp_cascade(c,j,k)
                  endif
                  net_nmin_vr(c,j) = net_nmin_vr(c,j) - pmnf_decomp_cascade(c,j,k)
                  if (decomp_method == mimics_decomp) then
                     decomp_cascade_sminn_flux_vr(c,j,k) = decomp_cascade_sminn_flux_vr(c,j,k) - p_decomp_npool_to_din(c,j,k)
                     net_nmin_vr(c,j) = net_nmin_vr(c,j) + p_decomp_npool_to_din(c,j,k)
                  else if (decomp_method == mimicsplus_decomp) then !ECW maybe .or.
                     decomp_cascade_sminn_flux_vr(c,j,k) = decomp_cascade_sminn_flux_vr(c,j,k) - p_decomp_npool_to_din(c,j,k)
                     net_nmin_vr(c,j) = net_nmin_vr(c,j) + p_decomp_npool_to_din(c,j,k)
                  end if
               else
                  decomp_cascade_ntransfer_vr(c,j,k) = 0._r8
                  if (.not. use_nitrif_denitrif) then
                     sminn_to_denit_decomp_cascade_vr(c,j,k) = 0._r8
                  end if
                  decomp_cascade_sminn_flux_vr(c,j,k) = 0._r8
               end if

            end do
         end do
      end do



      if (use_lch4) then
         ! Calculate total fraction of potential HR, for methane code
         do j = 1,nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               hrsum(c,j) = 0._r8
            end do
         end do
         do k = 1, ndecomp_cascade_transitions
            do j = 1,nlevdecomp
               do fc = 1,num_bgc_soilc
                  c = filter_bgc_soilc(fc)
                  hrsum(c,j) = hrsum(c,j) + rf_decomp_cascade(c,j,k) * p_decomp_cpool_loss(c,j,k)
               end do
            end do
         end do


        ! Nitrogen limitation / (low)-moisture limitation                                                                    
         do j = 1,nlevdecomp
            do fc = 1,num_bgc_soilc
               c = filter_bgc_soilc(fc)
               if (phr_vr(c,j) > 0._r8) then
                  fphr(c,j) = hrsum(c,j) / phr_vr(c,j) * w_scalar(c,j)
                  fphr(c,j) = max(fphr(c,j), 0.01_r8) ! Prevent overflow errors for 0 respiration                          
               else
                  fphr(c,j) = 1._r8
               end if
            end do
         end do
      end if

      
      ! vertically integrate net and gross mineralization fluxes for diagnostic output                                     

     do fc = 1,num_bgc_soilc
       c = filter_bgc_soilc(fc)
         do j = 1,nlevdecomp
              net_nmin(c) = net_nmin(c) + net_nmin_vr(c,j) * dzsoi_decomp(j)
              gross_nmin(c) = gross_nmin(c) + gross_nmin_vr(c,j) * dzsoi_decomp(j)
         end do
      end do

      ! update fluxes from mimicsplus after FUN


      if (decomp_method ==  mimicsplus_decomp) then

         do fc = 1,num_bgc_soilc
            c = filter_bgc_soilc(fc)
            do j = 1,nlevdecomp
               do k = 1,ndecomp_cascade_transitions

                  if (cascade_donor_pool(k) == i_ecm_myc) then
                     ! mining fluxes
                     if (cascade_receiver_pool(k) == i_chem_som) then
                     ! mortality fluxes have been calculated in the mimicsplus_decomp_rates()
                        decomp_cascade_ntransfer_vr(c,j,k) = decomp_cascade_ntransfer_vr(c,j,k) - n_somc2ecm_vr(c,j)
                     else if (cascade_receiver_pool(k) == i_phys_som) then
                        decomp_cascade_ntransfer_vr(c,j,k) = decomp_cascade_ntransfer_vr(c,j,k) - n_somp2ecm_vr(c,j)
                     end if
                  endif
                  if  (cascade_receiver_pool(k) == i_avl_som) then
                     ! carbon release associated with mining
                     if (cascade_donor_pool(k) == i_chem_som) then
                        decomp_cascade_ctransfer_vr(c,j,k) = decomp_cascade_ctransfer_vr(c,j,k) + c_somc2soma_vr(c,j)
                     else if (cascade_donor_pool(k) == i_phys_som) then
                        decomp_cascade_ctransfer_vr(c,j,k) = decomp_cascade_ctransfer_vr(c,j,k) + c_somp2soma_vr(c,j)
                     end if
                  end if

                end do ! transitions
                ! fluxes that are not part of the cascade.
                decomp_cpools_vr(c,j,i_ecm_myc) = decomp_cpools_vr(c,j,i_ecm_myc) + c_ecm_growth_vr(c,j) * dt
                decomp_npools_vr(c,j,i_ecm_myc) = decomp_npools_vr(c,j,i_ecm_myc) + (n_ecm_growth_vr(c,j) + &
                                                  sminno3_to_ecm_vr(c,j) + sminnh4_to_ecm_vr(c,j)) * dt
                decomp_cpools_vr(c,j,i_am_myc) = decomp_cpools_vr(c,j,i_am_myc) + c_am_growth_vr(c,j) * dt
                decomp_npools_vr(c,j,i_am_myc) = decomp_npools_vr(c,j,i_am_myc) + (n_ecm_growth_vr(c,j) + &
                                                   sminno3_to_am_vr(c,j) + sminnh4_to_am_vr(c,j)) * dt
                decomp_cpools_vr(c,j,i_avl_som) = decomp_cpools_vr(c,j,i_avl_som) + c_ecm_enz_vr(c,j) * dt
            end do ! layer
         enddo !column




      end if

    end associate

  end subroutine SoilBiogeochemDecomp
 
end module SoilBiogeochemDecompMod

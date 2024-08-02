module CNFUNMIMICSplusMod
   !--------------------------------------------------------------------
     !---
   ! ! DESCRIPTION
   ! ! Updated version of theThe FUN model developed by Fisher et al. 2010
   ! ! for coupling with MIMICSplus (Aas et al., 2023) decomposition model

   
   ! !USES: 
     use shr_kind_mod                    , only : r8 => shr_kind_r8
     use shr_log_mod                     , only : errMsg => shr_log_errMsg
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




     implicit none
     private
     !
     ! TYPES


     ! !PUBLIC MEMBER FUNCTIONS:

     



end module CNFUNMIMICSplusMod
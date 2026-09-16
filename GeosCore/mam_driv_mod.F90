! FAB think about ifdef
!#ifdef  MAM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: apm_driv_mod
!
! !DESCRIPTION: Module APM\_DRIV\_MOD contains variables and routines to drive
!  the Advanced Particle Microphysics (APM) model.  It serves as the
!  interface between APM module and the 3D model.
!\\
!\\
! !INTERFACE:
!
MODULE MAM_DRIV_MOD

USE precision_mod, only :r8 => f8, fp, f8  
USE physics_buffer, only: physics_buffer_desc
USE physics_types, only : physics_state, physics_ptend
USE mam_utils, only : masterproc, pcols, pver, l_h2so4g, l_soag, l_hno3g, l_hclg, l_nh3g
USE constituents, only : pcnst
use radconstants, only : nswbands, nlwbands
!
IMPLICIT NONE

PRIVATE

!PUBLIC MEMBER FUNCTIONS:

PUBLIC :: MAM_DRIV, MAM_INIT, MAM_APPLY_RAINOUT_EFF, MAM_OPT_to_RRTMG, MAM_ACTIVATE_EVAP, &
          MAM_OPT_to_PHOTOL, MAM_to_HETRATES

! !REMARKS:
!  The MAM model was designed and developed for implementation into GEOS-Chem
!
!
!EOP
!------------------------------------------------------------------------------
!BOC
!PUBLIC DATA


! sulf production rate calculates in   
! perhaps use AeroMass state variables 
REAL(fp), pointer, public :: PSO4AQ_RATE(:,:,:) ! Cld chem sulfate prod rate [kg s-1]
! FAB (MAM-decouple-std, Step 2b): SO4 from SO2+O3 in alkaline sea-salt water
! (KPP K_MT(1) fine, K_MT(4) coarse), set in fullchem_mod.F90 after KPP.
! [kg SO4 box-1 per chemistry step]; 4th dim: 1 = fine, 2 = coarse
REAL(fp), pointer, public :: PSO4SS_RATE(:,:,:,:)
REAL(fp), pointer, public :: H2SO4_RATE(:,:,:) ! H2SO4 prod rate [kg s-1]
! 
REAL(fp), pointer, public :: PSO4_SO2MAM(:,:,:)
!

TYPE(physics_buffer_desc), pointer :: pbuf(:)
TYPE(physics_state) :: physta
TYPE(physics_ptend) :: ptend

! define a specific type to handle MAM/GC prognostic species information.
! it is a bit similar to State_Chm%SpcData(N)%Info
! Convenient for communication between MAM and GC worlds for instance 
! Potentially it could live as a specific subtype of chem_state ? just an idea..
! Perhaps this type should also be declared in Headers for consistency with GC ? 

TYPE, public ::  mamspec ! 
  INTEGER :: gcind  ! sp index relative to Spc ( GC chemstate species )
  INTEGER :: mamind ! sp index relative to both q and qqcw MAM states 
  INTEGER :: modID  ! MAM mode index to which this sp belongs (also used to point to chemstate%GCMAM(mode)%xx)
  LOGICAL :: isnum  ! True if number concentration vs mass concentration
  LOGICAL :: iscb
  CHARACTER* 12  :: name !GC species name for MAM tracers ( ABSOLUTLY must be consistent with speciesdat.yml)  
  CHARACTER* 12  :: namecb !GC species name for MAM cloudborne sp ( ABSOLUTLY must be consistent with speciesdat.yml)
END TYPE mamspec

TYPE(mamspec), pointer :: mamgc(:)

INTEGER nmamgc ! number of GC advected MAM tracers 


    INTEGER:: loffset, lchnk

    REAL(r8) :: deltat
    INTEGER, save  :: mamstep ! number of elapsed mam call since first call
    LOGICAL, save  :: lfirstcall
    LOGICAL, save  :: is_cbsim

    ! Interstitial-to-CB partner index: cb_partner(s) = index in mamgc of the
    ! CB species with same (mamind, modId) as interstitial species s; -1 if none.
    INTEGER, ALLOCATABLE, SAVE :: cb_partner(:)
    ! Cloud fraction from previous MAM_DRIV call, for activation/evaporation tendency.
    REAL(fp), ALLOCATABLE, SAVE :: CLDF_prev(:,:,:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-/Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: mam_driv
!
! !DESCRIPTION: Subroutine MAM\_DRIV is the interface between MAM and
!  the GEOS-Chem model.
!\\
!\\
! !INTERFACE:
!

SUBROUTINE MAM_DRIV( Input_Opt,  State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : SpcConc
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Diag_Mod, ONLY : DgnState
    USE UnitConv_Mod,   ONLY : Check_Units 

    USE mam_utils, only: begchunk, endrun,  mdo_coldstart, & 
                         mdo_gaschem, mdo_cloudchem,  mdo_gasaerexch,mdo_rename,          &
                         mdo_newnuc,  mdo_coag           
    USE chem_mods, only: adv_mass, gas_pcnst, imozart
    USE physconst, only: mwdry, rga
    USE modal_aero_data, only: numptr_amode, lptr_so4_a_amode, &
                               lptr_bc_a_amode, lptr_nacl_a_amode,&
                               lptr_pom_a_amode, lptr_soa_a_amode,&
                               lptr_dust_a_amode, lptr_so4_cw_amode,&
                               lptr_nh4_a_amode,lptr_no3_a_amode,&
                               lptr_ca_a_amode,lptr_cl_a_amode,&
                               lptr_co3_a_amode,lptr_mom_a_amode,   &
                               modeptr_accum, modeptr_aitken,alnsg_amode, voltonumb_amode, &
                               ntot_amode, modeptr_coarse
    USE modal_aero_initialize_data, only: MAM_cold_start
    USE modal_aero_calcsize, only: modal_aero_calcsize_sub
    USE modal_aero_wateruptake, only: modal_aero_wateruptake_dr, &
                                      load_pbuf, unload_pbuf
    USE modal_aero_amicphys, only: modal_aero_amicphys_intr
    use mam_opt , only : mam_aero_sw,mam_aero_lw, mamoptdiag
    !
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
!
! !LOCAL VARIABLES:
!
    INTEGER :: I,J,L,S,l2,K,N,M,SIZENUM,MDAY

    ! Make a pointer to the tracer array
    TYPE(SpcConc), POINTER :: Spc(:)

! !LOCAL VARIABLES:
!
      REAL(r8) :: vmr(pcols,pver,gas_pcnst)     ! gas & aerosol volume mixing ratios
      REAL(r8) :: vmrcw(pcols,pver,gas_pcnst)   ! gas & aerosol cloud-borne volume mixing ratios
      REAL(r8) :: vmr_svaa(pcols,pver,gas_pcnst) ! temp save before gas chem
      REAL(r8) :: vmr_svbb(pcols,pver,gas_pcnst) ! temp save before cloud chem
      REAL(r8) :: vmrcw_svbb(pcols,pver,gas_pcnst)!temp save before cloud chem 
      REAL(r8) :: aircon(pcols,pver) !  air concentration (kmol/m3)

      !-----FAB  this mode avg local variables could be perhaps suppressed since mamoptdiag carries the information--- 
      real(r8)  :: tauxar(pcols,pver,nswbands)  ! aerosol extinction optical depth
      real(r8)  :: wa(pcols,pver,nswbands)      ! aerosol single scattering albedo * tau
      real(r8)  :: ga(pcols,pver,nswbands)      ! aerosol asymmetry parameter * wa
      real(r8)  :: fa(pcols,pver,nswbands)      ! aerosol forward scattered fraction * ga
      real(r8) :: taux_lw(pcols,pver,nlwbands)
      real(r8) :: hplus_aer_out(pcols,pver,ntot_amode)
      real(r8) :: relhum_loc(pcols,pver)   ! clear-sky RH for wateruptake
      real(r8) :: pso4ss(pcols,pver,2)     ! FAB Step 2b: sea-salt aq. SO4 [kg/kg per step]
      !-------------
      INTEGER :: latndx(pcols),lonndx(pcols)  ! required by amicphys interface
!--------------------------------------------------------------------------


    IF ( lfirstcall ) then
           mamstep = 1
    ELSE 
           mamstep =mamstep+1
    END IF        
    ! Point to Spc
    Spc => State_Chm%Species

    latndx(:) = 0
    pso4ss(:,:,:) = 0.0_r8   ! FAB Step 2b
    lonndx(:) = 0

!     if (masterproc) then

    lchnk = begchunk
    loffset = imozart -1
    ! load the mam met state   
    physta%lchnk = lchnk
    physta%ncol  = pcols
    
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX 
      n = J + (I-1)*State_Grid%NY 
      physta%t(n,L) = State_Met%T(I,J,L)
      physta%pdel(n,l) = State_Met%DELP(I,J,L) * 100.0e+0_f8 ! hPa to Pa 
      physta%pdeldry(n,l) = State_Met%DELP_DRY(I,J,L) * 100.0e+0_f8 !
      physta%pmid(n,l) = State_Met%PMID(I,J,L) * 100.0e+0_f8
      physta%cld(n,l)= State_Met%CLDF(I,J,L)
      physta%relhum(n,l)= State_Met%RH(I,J,L)/100.0e+0_f8 ! used actually for Optics for now , could replace mam calc..
      physta%qv(n,l) = State_Met%SPHU(I,J,L) * 1.0e-3_r8 ! in kg/kgair  Caution here make sure
      physta%zm(n,l) =  sum(State_Met%BXHEIGHT(I,J,1:l))- 0.5 * State_Met%BXHEIGHT(I,J,l)

      IF (L==1) physta%pblh(n) = State_Met%PBLH(I,J)

      ! first element water vapor mr (used in water )  
      physta%q(n,l,1) = physta%qv(n,l) / (1.0e+0_f8 - physta%qv(n,l))

      ! load sulf production rate and convert from Kg.s-1  to kg.kg-1.s-1 
      ! needs fullchem activated 
      physta%ph2so4(n,l) = Spc(Ind_('PH2SO4'))%Conc(I,J,L) / State_Met%AD(I,J,L)/ deltat              

      physta%paqso4(n,l) =  Spc(Ind_('PSO4AQ'))%Conc(I,J,L)  / State_Met%AD(I,J,L) /deltat

      ! FAB Step 2b: sea-salt aqueous SO4 from previous KPP call, kg -> kg/kg
      pso4ss(n,l,1) = PSO4SS_RATE(I,J,L,1) / State_Met%AD(I,J,L)
      pso4ss(n,l,2) = PSO4SS_RATE(I,J,L,2) / State_Met%AD(I,J,L)
      ! load q gas ...
      ! the gas phase species SO2,DMS,H2O2 in q are not used/modified 
      ! if we use GC production rate for SO4 instead of MAM simple chem 
      ! ,maybe get rid of them later or consider haveing some oxidaton routines
      ! which could run dindependant of fullchem e.g. from prescribed oxidants..

    
      physta%q(n,l,l_h2so4g) = Spc(Ind_('H2SO4'))%Conc(I,J,L) / State_Met%AD(I,J,L)
      physta%q(n,l,l_soag) = Spc(Ind_('SOAP'))%Conc(I,J,L) / State_Met%AD(I,J,L) 
      ! Rq in GC standard lumped SOAP is not treated as semi_volatil but in MAM yes
      ! is this reqsonqble ? also develop options with advanced SOA scheme  
      !
      IF(l_nh3g > 0) physta%q(n,l,l_nh3g) = Spc(Ind_('NH3'))%Conc(I,J,L) / State_Met%AD(I,J,L)
      IF(l_hno3g > 0) physta%q(n,l,l_hno3g) = Spc(Ind_('HNO3'))%Conc(I,J,L) / State_Met%AD(I,J,L)
      IF(l_hclg > 0) physta%q(n,l,l_hclg) = Spc(Ind_('HCL'))%Conc(I,J,L) / State_Met%AD(I,J,L)
     END DO
     END DO
     END DO


        ! This cold start init is temporary until handling a proper GC restart file
    IF (lfirstcall ) then
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX ! 
      n = J + (I-1)*State_Grid%NY
      physta%aircon(n,l) =  State_Met%AIRDEN(I,J,L) * 0.0345e+0_f8 !Kmol/m3 
    END DO
    END DO
    END DO
    ! initialise q aerosol state (cold start) / for testing phase  

     CALL MAM_cold_start (physta)
 
     IF ( mdo_coldstart == 1) then
     ! mdo_coldstart initialies in Mam_cold_start
     DO L = 1, State_Grid%NZ
     DO J = 1, State_Grid%NY
     DO I = 1, State_Grid%NX ! 
       n = J + (I-1)*State_Grid%NY
! Cold-start overrides of MAM accumulation-mode SO4/NO3/NH4 mass + number from
! the equivalent standard GC species, so MAM and standard chem start from
! the same aerosol burden.
       physta%q(n,l,lptr_so4_a_amode(1)) = Spc(IND_('SO4'))%Conc(I,J,L)/State_Met%AD(I,J,L)
       physta%q(n,l,numptr_amode(1)) =    physta%q(n,l,lptr_so4_a_amode(1)) /1700. * voltonumb_amode(1)

       physta%q(n,l,lptr_no3_a_amode(1)) = Spc(IND_('NIT'))%Conc(I,J,L)/State_Met%AD(I,J,L)
       physta%q(n,l,numptr_amode(1)) =   physta%q(n,l,numptr_amode(1))+  physta%q(n,l,lptr_no3_a_amode(1)) /1700. * voltonumb_amode(1)

       physta%q(n,l,lptr_nh4_a_amode(1)) = Spc(IND_('NH4'))%Conc(I,J,L)/State_Met%AD(I,J,L)
       physta%q(n,l,numptr_amode(1)) =   physta%q(n,l,numptr_amode(1))+  physta%q(n,l,lptr_nh4_a_amode(1)) /1700. * voltonumb_amode(1)

!      Spc(IND_('MAMDEV'))%Conc(I,J,L) = Spc(IND_('SO4'))%Conc(I,J,L)
    END DO
    END DO
    END DO
    ENDIF
    ENDIF


    IF (.not. lfirstcall .or. mdo_coldstart < 1) then 
     ! load q and qqcw mam state aerosol variables 
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX ! GEOS-Chem vertical grid is bottom-up
      n = J + (I-1)*State_Grid%NY         
      !number concentrations #/gridbox and convert to #/kg 
      !mass concentrations Kg/gridbox to Kg/Kg (mixing ratios) 

       IF(.not. is_cbsim) then !cloud borne state is not considered    
         DO s = 1 , nmamgc
          physta%q(n,l,mamgc(s)%mamind) = Spc(mamgc(s)%gcind)%Conc(I,J,L)/State_Met%AD(I,J,L) 
         END DO
       ELSE !cloud born and interstitial states are considered       
         DO s = 1 , nmamgc
           IF (.not. mamgc(s)%iscb ) then
             physta%q(n,l,mamgc(s)%mamind) = Spc(mamgc(s)%gcind)%Conc(I,J,L)/State_Met%AD(I,J,L) !
           ELSE
             physta%qqcw(n,l,mamgc(s)%mamind) = Spc(mamgc(s)%gcind)%Conc(I,J,L)/State_Met%AD(I,J,L) !
           END IF
         END DO  
       END IF  
    END DO
    END DO
    END DO
    END IF ! .

!FAB TEST 
  physta%qaerwat(:,:,:) = 0.
! CALCSIZE INTERFACE     
     
CALL load_pbuf( pbuf, lchnk, pcols, &
        physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens,physta%hygro )

! call calcsize
    ptend%lq = .false.
    ptend%q = 0.0e+0_f8
    CALL modal_aero_calcsize_sub( physta, ptend, deltat, pbuf, &
         do_adjust_in=.true., do_aitacc_transfer_in=.true. )

! unload pbuf
      CALL unload_pbuf( pbuf, lchnk, pcols, &
         physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens,physta%hygro )
! note also that the cloudborne state is supposed to be directly updated in calcsize
! (perhaps because unlinke q , qqcw is not an advected state in cesm and the tendencie does not need to be 
!  passed up ...
!  
      DO l = 1, pcnst
         IF ( .not. ptend%lq(l) ) cycle
         DO k = 1, pver
         DO i = 1, pcols 
            physta%q(i,k,l) = physta%q(i,k,l) + ptend%q(i,k,l)*deltat  
            physta%q(i,k,l) = max( physta%q(i,k,l), 0.0e+0_f8 )
         END DO
         END DO
      END DO
      
      physta%lchnk = lchnk      
      physta%ncol = pcols      




! WATER UPTAKE
     CALL load_pbuf( pbuf, lchnk, pcols, &
        physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens, physta%hygro )
!
     ! Clear-sky RH for water uptake: RH_clear = (RH_mean - fcld) / (1 - fcld)
     do l = 1, pver
       do n = 1, pcols
         if (physta%cld(n,l) < 1.0_r8) then
           relhum_loc(n,l) = max( 0.0_r8, &
              (physta%relhum(n,l) - physta%cld(n,l)) / (1.0_r8 - physta%cld(n,l)) )
         else
           relhum_loc(n,l) = 0.0_r8
         end if
       end do
     end do
     CALL modal_aero_wateruptake_dr( physta, pbuf, deltat, mamstep, clear_rh_in = relhum_loc)
     
     CALL unload_pbuf( pbuf, lchnk, pcols, &
         physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens,physta%hygro )

!--------------------------------------------------------------------------------
! switch from q & qqcw mass mixing ratios to volume mixing ratios  vmr and vmrcw
! only adress the gas/aerosol variables in q, qqcw 
! Rq : not all gases are modified by MAM routines so q size could be reduced  
! Rq2: the flow of unit changes could be perhaps simpler
      vmr = 0.0e+0_f8
      vmrcw = 0.0e+0_f8
      DO l = imozart, pcnst
         l2 = l - loffset
         vmr(  1:pcols,1:pver,l2) =physta%q(  1:pcols,1:pver,l)*mwdry/adv_mass(l2)
         vmrcw(1:pcols,1:pver,l2) =physta%qqcw(1:pcols,1:pver,l)*mwdry/adv_mass(l2)
      END DO
!-------------------------------
! GASCHEM interface 

!-------------------------------
! FAB (MAM-decouple-std, Step 2b): SO4 produced by SO2 + O3 in alkaline
! sea-salt aerosol water (KPP K_MT(1)/K_MT(4), diagnosed in fullchem_mod).
! Added to INTERSTITIAL SO4 (the reaction takes place in aerosol water, not
! cloud water): fine -> accumulation mode (Aitken sea-salt mass is small),
! coarse -> coarse mode. Added BEFORE amicphys so that MOSAIC re-equilibrates
! the new sulfate with the sea-salt Cl-/NO3- (acid displacement -> HCl).
! Added before vmr_svaa so it is treated as state, not as a gas/cloud-chem
! tendency. Number is unchanged (condensed mass; calcsize adjusts size).
      l2 = lptr_so4_a_amode(modeptr_accum) - loffset
      vmr(1:pcols,1:pver,l2) = vmr(1:pcols,1:pver,l2) + pso4ss(1:pcols,1:pver,1)*mwdry/adv_mass(l2)
      l2 = lptr_so4_a_amode(modeptr_coarse) - loffset
      vmr(1:pcols,1:pver,l2) = vmr(1:pcols,1:pver,l2) + pso4ss(1:pcols,1:pver,2)*mwdry/adv_mass(l2)

      vmr_svaa   = vmr  !save before gas chem , this is how the mam code proceed

      IF (mdo_gaschem > 0) then
        !
        l2 = l_h2so4g-loffset
        vmr(1:pcols,1:pver,l2) = vmr(1:pcols,1:pver,l2) + physta%ph2so4(1:pcols,1:pver)*mwdry/adv_mass(l2)*deltat
      END IF 
      vmr_svbb = vmr    ! save before cloud chem
      vmrcw_svbb = vmrcw! save before cloud chem

!CLOUDCHEM
!rq vmrcw / qcw are not advected in MAM /CESM

      IF (mdo_cloudchem > 0 .and. is_cbsim) then
        !updating sulfate from aq.chem mass in the cb accumulation mode
        
        l2 = lptr_so4_cw_amode(modeptr_accum) - loffset
        vmrcw(1:pcols,1:pver,l2) = vmrcw(1:pcols,1:pver,l2) +  &
                                   physta%paqso4(1:pcols,1:pver)*mwdry/adv_mass(l2)*deltat
      END IF
      
!------------------------------
  IF(.true.) then ! .and. masterproc ) then
     CALL modal_aero_amicphys_intr(               &
         mdo_gasaerexch,     mdo_rename,          &
         mdo_newnuc,         mdo_coag,            &
         lchnk,    pcols,     mamstep   ,           &
         loffset,  deltat,                        &
         latndx,   lonndx,                        &
         physta%t,   physta%pmid, physta%pdel,    &
         physta%zm,  physta%pblh,                 &
         physta%qv,  physta%cld ,                 &
         vmr,                vmrcw,               &   ! after  cloud chem
         vmr_svaa,                                &   ! before gas chem
         vmr_svbb,           vmrcw_svbb,          &   ! before cloud chem!
!         nqtendbb,           nqqcwtendbb,         &  ! ifdef cambox not enabled for now
!         dvmrdt_bb,          dvmrcwdt_bb,         &  ! in the interface maybe conssider for diag
         physta%dgncur_a,     physta%dgncur_awet,  &
         physta%wetdens,      physta%qaerwat,             &
         hplus_aer_out = hplus_aer_out,           &
         relhum_in     = physta%relhum             )
    END IF
! vmr and vmrcw have been updated in modal_aero_amicphys_intr


! if not cb_sim, apply aqueous SO4 production after gas-aerosol exchange so it does not
! feed into condensation/renaming/nucleation (it was produced in cloud droplets,
! not in the interstitial phase)
      IF (mdo_cloudchem > 0 .and. .not. is_cbsim) then
        l2 = lptr_so4_a_amode(modeptr_accum) - loffset
        vmr(1:pcols,1:pver,l2) = vmr(1:pcols,1:pver,l2) +  &
                                 physta%paqso4(1:pcols,1:pver)*mwdry/adv_mass(l2)*deltat
      END IF

! switch back from vmr & vmrcw to q & qqcw state
!
      DO l = imozart, pcnst
         l2 = l - loffset
         physta%q(    1:pcols,1:pver,l)  = vmr(  1:pcols,1:pver,l2) * adv_mass(l2)/mwdry 
         physta%qqcw( 1:pcols,1:pver,l)  = vmrcw(1:pcols,1:pver,l2) * adv_mass(l2)/mwdry
      END DO


      
      
 ! 
! calculate optical properties 
      IF(.true.) then 
       call  mam_aero_sw(physta, tauxar, wa, ga, fa, mamoptdiag)
       call  mam_aero_lw(physta, taux_lw, mamoptdiag)      
      end if 
! Update GC/MAM species 
     DO L = 1, State_Grid%NZ
     DO J = 1, State_Grid%NY
     DO I = 1, State_Grid%NX
        n = J + (I-1)*State_Grid%NY
       
        IF (.not.is_cbsim) then ! there is no cloud borne state 
          DO s =1, nmamgc  
            Spc(mamgc(s)%gcind)%Conc(I,J,L) = physta%q(n,l,mamgc(s)%mamind) &
                                            * State_Met%AD(I,J,L)
          END DO
        ELSE  ! there is a cloud-borne state
            DO s =1, nmamgc                                  
              IF ( .not. mamgc(s)%iscb) then
                Spc(mamgc(s)%gcind)%Conc(I,J,L) = physta%q(n,l,mamgc(s)%mamind) &
                                              * State_Met%AD(I,J,L)
               ELSE ! treat cloudborne    
                 Spc(mamgc(s)%gcind)%Conc(I,J,L) = physta%qqcw(n,l,mamgc(s)%mamind) &
                      * State_Met%AD(I,J,L)
               END IF
           END DO
         END IF   
        !gas species affected by mam 
        Spc(Ind_('SOAP'))%Conc(I,J,L) = physta%q(n,l,l_soag) &
                                       * State_Met%AD(I,J,L)  

        Spc(Ind_('H2SO4'))%Conc(I,J,L) = physta%q(n,l,l_h2so4g) &
                                       * State_Met%AD(I,J,L)

        IF(l_nh3g > 0) Spc(Ind_('NH3'))%Conc(I,J,L) = physta%q(n,l,l_nh3g) &
                                       * State_Met%AD(I,J,L)

        IF(l_hno3g > 0) Spc(Ind_('HNO3'))%Conc(I,J,L) = physta%q(n,l,l_hno3g) &
                                       * State_Met%AD(I,J,L)

        IF(l_hclg > 0) Spc(Ind_('HCL'))%Conc(I,J,L) = physta%q(n,l,l_hclg) &
                                       * State_Met%AD(I,J,L)
      ENDDO
      ENDDO
      ENDDO

! fill state GCMAM chem state variables,  used in e.g. drydep  nd diags
! harmonize mamgc and GCMAM 
! try to optimize the if statements within loops

IF(1==1) THEN
      DO L = 1, State_Grid%NZ
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX
        n = J + (I-1)*State_Grid%NY

      DO m= 1, size(State_Chm%GCMAM) ! loop on modes     

        State_Chm%GCMAM(m)%nudryrad(I,J,L) = 0.5e+0_f8*physta%dgncur_a(n,L,m)   

        State_Chm%GCMAM(m)%nuwetrad(I,J,L) = 0.5e+0_f8*physta%dgncur_awet(n,L,m)   

        ! volume mean geo radius
        State_Chm%GCMAM(m)%dryrad(I,J,L) = 0.5e+0_f8*physta%dgncur_a(n,L,m)     &
                                            *exp(3.0e+0_f8*(alnsg_amode(m)**2))
        State_Chm%GCMAM(m)%wetrad(I,J,L) = 0.5e+0_f8*physta%dgncur_awet(n,L,m)  &
                                            *exp(3.0e+0_f8*(alnsg_amode(m)**2))
        ! wet aer density
        State_Chm%GCMAM(m)%aerdens(I,J,L) =  physta%wetdens(n,L,m)         
       
        ! volume mean hygroscopicity 
        State_Chm%GCMAM(m)%hygro(I,J,L) =  physta%hygro(n,L,m)  

        ! aerosol water concentration
        State_Chm%GCMAM(m)%aerwat(I,J,L) =  physta%qaerwat(n,L,m)

        !modal number concentrations in #.m-3
        State_Chm%GCMAM(m)%Nu(I,J,L) =              &
                                physta%q(n,L,numptr_amode(m))*State_Met%AIRDEN(I,J,L)

        ! number-based wet surface area [cm2/cm3]
        ! SA = pi * N[#/m3] * Dgn_wet[m]^2 * exp(2*ln(sg)^2) * 1e-2 (m^-1 -> cm2/cm3)
        State_Chm%GCMAM(m)%saer(I,J,L) = acos(-1.0_f8)                       &
                                        * State_Chm%GCMAM(m)%Nu(I,J,L)        &
                                        * physta%dgncur_awet(n,L,m)**2         &
                                        * exp(2.0_f8*alnsg_amode(m)**2)        &
                                        * 1.0e-2_f8
        ! modal mass concentrations in Kg.m-3  
        IF(lptr_so4_a_amode(m) > 0 ) State_Chm%GCMAM(m)%so4(I,J,L) =              & 
                                physta%q(n,L,lptr_so4_a_amode(m))*State_Met%AIRDEN(I,J,L) 
        !  + physta%qqcw(n,L,lptr_so4_a_amode(m))*State_Met%AIRDEN(I,J,L) 
                        !  
        !FAB TEMP add the cloud borne sulf to chm state for diag // change that once 
        ! transfer from qqcw to q is properly trated !!
        ! physta%qqcw(n,L,lptr_so4_a_amode(m))*State_Met%AIRDEN(I,J,L)    

         
        IF(lptr_bc_a_amode(m) > 0 ) State_Chm%GCMAM(m)%bc(I,J,L) =              &
                                physta%q(n,L,lptr_bc_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_pom_a_amode(m) > 0 ) State_Chm%GCMAM(m)%pom(I,J,L) =              &
                                physta%q(n,L,lptr_pom_a_amode(m))*State_Met%AIRDEN(I,J,L)
 
        IF(lptr_soa_a_amode(m) > 0 ) State_Chm%GCMAM(m)%soa(I,J,L) =              &
                                physta%q(n,L,lptr_soa_a_amode(m))*State_Met%AIRDEN(I,J,L)
        
        IF(lptr_nacl_a_amode(m) > 0 ) State_Chm%GCMAM(m)%sslt(I,J,L) =              &
                                physta%q(n,L,lptr_nacl_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_dust_a_amode(m) > 0 ) State_Chm%GCMAM(m)%dust(I,J,L) =              &
                                physta%q(n,L,lptr_dust_a_amode(m))*State_Met%AIRDEN(I,J,L)
 
        !FAB: lptr_xxx_a_amode(m)>0 is structural (mechanism-level), but
        ! GCMAM(m)%xxx is only allocated if the species is actually advected
        ! (see GCMAM(m)%lxxx, mam_container_mod.F90) -- guard on both, else
        ! dropping a species (e.g. MAMMOM) segfaults on a null-pointer write.
        IF(lptr_nh4_a_amode(m) > 0 .AND. State_Chm%GCMAM(m)%lnh4 ) State_Chm%GCMAM(m)%nh4(I,J,L) =  &
                               physta%q(n,L,lptr_nh4_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_no3_a_amode(m) > 0 .AND. State_Chm%GCMAM(m)%lno3 ) State_Chm%GCMAM(m)%no3(I,J,L) =  &
                               physta%q(n,L,lptr_no3_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_ca_a_amode(m) > 0 .AND. State_Chm%GCMAM(m)%lca ) State_Chm%GCMAM(m)%ca(I,J,L) =     &
                               physta%q(n,L,lptr_ca_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_co3_a_amode(m) > 0 .AND. State_Chm%GCMAM(m)%lco3 ) State_Chm%GCMAM(m)%co3(I,J,L) =  &
                               physta%q(n,L,lptr_co3_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_cl_a_amode(m) > 0 .AND. State_Chm%GCMAM(m)%lcl ) State_Chm%GCMAM(m)%cl(I,J,L) =     &
                               physta%q(n,L,lptr_cl_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_mom_a_amode(m) > 0 .AND. State_Chm%GCMAM(m)%lmom ) State_Chm%GCMAM(m)%mom(I,J,L) =  &
                               physta%q(n,L,lptr_mom_a_amode(m))*State_Met%AIRDEN(I,J,L)

! optics
        State_Chm%GCMAM(m)%tauxar(I,J,L,:) = mamoptdiag(m)%tauxar(n,L,:)
        State_Chm%GCMAM(m)%ssa(I,J,L,:) = mamoptdiag(m)%ssa(n,L,:)
        State_Chm%GCMAM(m)%g(I,J,L,:) = mamoptdiag(m)%g(n,L,:)
        State_Chm%GCMAM(m)%hplus(I,J,L) = hplus_aer_out(n,L,m)
      END DO
     ENDDO
     ENDDO
     ENDDO
END IF



     ! call MAM aerosol update for gravitational settling (this call could be somewhere else )        
     CALL  MAM_SETTL( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
   
                
!Fill out MAM diags (cf Headers/state_diag)    

    CALL Set_MAM_Diagnostic( Input_Opt, State_Chm, State_Diag, &
                                     State_Grid, State_Met, RC )

!

    ! Activation / evaporation: transfer between interstitial and CB based on
    ! CLDF tendency. On first call just snapshot CLDF, no transfer.
    IF (is_cbsim) THEN
      CALL MAM_ACTIVATE_EVAP( State_Grid, State_Met, State_Chm, Spc, lfirstcall )
    END IF

    Spc => NULL()
    IF (lfirstcall) lfirstcall = .false.


  END SUBROUTINE MAM_DRIV

!------------------------------------------------------------------------------
! !IROUTINE: MAM_ACTIVATE_EVAP
!
! !DESCRIPTION:
!  Transfer mass between interstitial (q) and cloud-borne (qqcw / MAMCB*) MAM
!  tracers based on the cloud-fraction tendency over one chemistry timestep.
!
!  ACTIVATION  (dCLDF > 0):  the newly clouded volume fraction dCLDF of the
!    grid cell is assumed to contain interstitial aerosol at the grid-mean
!    concentration.  The fraction f_act activates into cloud droplets, where
!    f_act = MIN(1, κ · 27 · Dd³ · Ss² / (4·A³)) with Ss = 0.3 % (stratiform).
!    Both κ and Dd are mode- and grid-point-specific, updated by calcsize /
!    wateruptake earlier in the same MAM_DRIV call.
!
!  EVAPORATION (dCLDF < 0):  the fraction |dCLDF| / CLDF_prev of the cloud
!    volume evaporates, releasing that fraction of the CB mass back to the
!    interstitial pool.
!
!  On the first call (is_firstcall = .TRUE.) only the CLDF snapshot is stored;
!  no transfer is performed to avoid a cold-start spike.
!------------------------------------------------------------------------------
SUBROUTINE MAM_ACTIVATE_EVAP( State_Grid, State_Met, State_Chm, Spc, is_firstcall )

    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE Species_Mod,    ONLY : SpcConc

    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(ChmState), INTENT(IN)    :: State_Chm
    TYPE(SpcConc),  INTENT(INOUT) :: Spc(:)
    LOGICAL,        INTENT(IN)    :: is_firstcall

    ! Köhler constants: A = 2.1e-9 m → 4·A³ = 3.7044e-26 m³; Ss_strat = 0.003
    REAL(fp), PARAMETER :: A3_x4 = 3.7044e-26_fp   ! 4·(2.1e-9)³  [m³]
    REAL(fp), PARAMETER :: Ss2_s = 9.0e-6_fp        ! (0.003)²

    INTEGER  :: I, J, L, s, s_cb, m
    REAL(fp) :: CF, CF_prev, dCF, f_act, delta

    ! On first call: snapshot CLDF and return without any transfer.
    IF (is_firstcall) THEN
      CLDF_prev(:,:,:) = State_Met%CLDF(:,:,:)
      RETURN
    END IF

    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

      CF      = State_Met%CLDF(I,J,L)
      CF_prev = CLDF_prev(I,J,L)
      dCF     = CF - CF_prev

      IF (ABS(dCF) < 1.0e-10_fp) CYCLE

      DO s = 1, nmamgc
        IF (mamgc(s)%iscb) CYCLE
        s_cb = cb_partner(s)
        IF (s_cb < 0) CYCLE

        m = mamgc(s)%modId

        IF (dCF > 0.0_fp) THEN
          ! Activation: f_act fraction of interstitial in the newly clouded
          ! volume dCF transfers to CB.
          f_act = MIN(1.0_fp,                                               &
                      State_Chm%GCMAM(m)%hygro(I,J,L)                     &
                      * 27.0_fp                                             &
                      * (2.0_fp * State_Chm%GCMAM(m)%nudryrad(I,J,L))**3  &
                      * Ss2_s / A3_x4 )
          delta = f_act * MAX(Spc(mamgc(s)%gcind)%Conc(I,J,L), 0.0_fp) * dCF
          Spc(mamgc(s)%gcind)%Conc(I,J,L)    = Spc(mamgc(s)%gcind)%Conc(I,J,L)    - delta
          Spc(mamgc(s_cb)%gcind)%Conc(I,J,L) = Spc(mamgc(s_cb)%gcind)%Conc(I,J,L) + delta

        ELSE
          ! Evaporation: fraction |dCF|/CF_prev of CB releases back to interstitial.
          delta = MAX(Spc(mamgc(s_cb)%gcind)%Conc(I,J,L), 0.0_fp)        &
                  * ABS(dCF) / MAX(CF_prev, 1.0e-6_fp)
          Spc(mamgc(s_cb)%gcind)%Conc(I,J,L) = Spc(mamgc(s_cb)%gcind)%Conc(I,J,L) - delta
          Spc(mamgc(s)%gcind)%Conc(I,J,L)    = Spc(mamgc(s)%gcind)%Conc(I,J,L)    + delta
        END IF

      END DO

      CLDF_prev(I,J,L) = CF

    END DO
    END DO
    END DO

END SUBROUTINE MAM_ACTIVATE_EVAP

!-----------------------------------------------------------------

SUBROUTINE MAM_INIT( Input_Opt, State_Chm,  State_Diag, State_Grid, RC )

    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : SpcConc
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Diag_Mod, ONLY : DgnState

    USE TIME_MOD,     ONLY : GET_TS_CHEM

    USE mam_utils, ONLY :endrun, plev, ncol_for_outfld  
    USE physics_buffer, only: physics_buffer_desc 
    USE physics_types, only : physics_state
    USE modal_aero_data, only: numptr_amode, lptr_so4_a_amode, &
                               lptr_bc_a_amode, lptr_nacl_a_amode,&
                               lptr_pom_a_amode, lptr_soa_a_amode,&
                               lptr_dust_a_amode, lptr_nh4_a_amode,&
                               lptr_no3_a_amode,lptr_ca_a_amode,&
                               lptr_cl_a_amode,lptr_co3_a_amode,&
                               lptr_mom_a_amode 

    USE modal_aero_initialize_data, only: MAM_init_basics, MAM_ALLOCATE
    USE mam_opt, only: mam_init_opt
    ! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object

!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!!

    INTEGER :: species_class(pcnst) = -1
    INTEGER  :: m , i ,j, s 
    CHARACTER * 12 :: tmp

    
    
!-------initialise namelist parameters

!   mdo_gaschem=1
!   mdo_cloudchem =1

!   mdo_gasaerexch=1
!   mdo_rename=1
!   mdo_newnuc=1
!   mdo_coag=1

   is_cbsim = .true.


   masterproc = Input_Opt%amIRoot

!FAB test restart functionning      
   lfirstcall = .true.

! Chemistry timestep [s]
   deltat  = GET_TS_CHEM()


pcols = State_Grid%NX *  State_Grid%Ny
pver  = State_Grid%NZ
plev = pver

!this control diag file printing in mam ...
!set to 1 for limiting printing
!anyway you should revisit these printing in GC context.
ncol_for_outfld = 1

CALL MAM_init_basics(pbuf)

!allocate MAM state object
CALL MAM_ALLOCATE (physta,ptend )

! allocate and init OPTICS

CALL MAM_INIT_OPT()


!allocate specific GC diqg usefull for mam  
ALLOCATE( PSO4AQ_RATE(State_Grid%NX,State_Grid%NY,State_Grid%NZ) )
ALLOCATE( H2SO4_RATE(State_Grid%NX,State_Grid%NY,State_Grid%NZ) )
ALLOCATE( PSO4_SO2MAM(State_Grid%NX,State_Grid%NY,State_Grid%NZ) )
!FAB Step 2b: sea-salt aqueous SO4 production (1 = fine, 2 = coarse)
ALLOCATE( PSO4SS_RATE(State_Grid%NX,State_Grid%NY,State_Grid%NZ,2) )
PSO4SS_RATE = 0.0_fp



!Initialize MAM4/GC transported tracer info  

nmamgc = State_Chm%nMam
allocate(mamgc(State_Chm%nMam)) 

DO i = 1, State_Chm%nMam
 m = State_Chm%map_Mam(i) 
 mamgc(i)%name = State_Chm%SpcData(m )%Info%name 
 mamgc(i)%gcind = m 
 mamgc(i)%modId =State_Chm%SpcData(m )%Info%MamModId
 mamgc(i)%isnum = State_Chm%SpcData(m )%Info%MP_SizeResNum
 mamgc(i)%iscb = State_Chm%SpcData(m )%Info%Is_CloudBorne
! add maping info for mam q and qqcw states 
! test name and  for cloud borne simulation - important species in q and qqcw have the same indexing
 s=4 ! work only for MAMXXX naming convention
 IF (is_cbsim .and. mamgc(i)%iscb) s=6
 mamgc(i)%mamind =-1
 !mam modal indices should match existing species and be > 0 
 !should also be consistent with the allocation state of 
 !the chem state%GCMAM(mode)%XXX(:,:,:) in /headers
 !t.b.d make a security test with explicit error message
 IF (mamgc(i)%name(s:s+1) == 'Nu') mamgc(i)%mamind = numptr_amode(mamgc(i)%modId)
 IF (mamgc(i)%name(s:s+2) == 'SO4') mamgc(i)%mamind = lptr_so4_a_amode(mamgc(i)%modId)
 IF (mamgc(i)%name(s:s+1) == 'BC') mamgc(i)%mamind =  lptr_bc_a_amode(mamgc(i)%modId)
 IF (mamgc(i)%name(s:s+2) == 'POM') mamgc(i)%mamind = lptr_pom_a_amode(mamgc(i)%modId)
 IF (mamgc(i)%name(s:s+2) == 'SOA') mamgc(i)%mamind = lptr_soa_a_amode(mamgc(i)%modId)
 IF (mamgc(i)%name(s:s+3) == 'SSLT') mamgc(i)%mamind = lptr_nacl_a_amode(mamgc(i)%modId)
 IF (mamgc(i)%name(s:s+3) == 'DUST') mamgc(i)%mamind = lptr_dust_a_amode(mamgc(i)%modId)
! MOSAIC species / consider NH4 for default 
#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
IF (mamgc(i)%name(s:s+2) == 'NH4') mamgc(i)%mamind = lptr_nh4_a_amode(mamgc(i)%modId)
IF (mamgc(i)%name(s:s+2) == 'NO3') mamgc(i)%mamind = lptr_no3_a_amode(mamgc(i)%modId)
IF (mamgc(i)%name(s:s+2) == 'MOM') mamgc(i)%mamind = lptr_mom_a_amode(mamgc(i)%modId)
IF (mamgc(i)%name(s:s+1) == 'CA') mamgc(i)%mamind = lptr_ca_a_amode(mamgc(i)%modId)
IF (mamgc(i)%name(s:s+1) == 'CL') mamgc(i)%mamind = lptr_cl_a_amode(mamgc(i)%modId)
IF (mamgc(i)%name(s:s+2) == 'CO3') mamgc(i)%mamind = lptr_co3_a_amode(mamgc(i)%modId)
#endif
 
IF  (masterproc .and. mamgc(i)%mamind < 0.0e+0_f8 ) then 
        CALL endrun('Stoping in MAM_INIT, GC species not consistent with MAM species, check species_database.yml') 
END IF          
! to be updated when adding species to MAM
END DO

 IF (is_cbsim .and. .not.any(mamgc(:)%iscb)) then
     stop
 END IF

! Build interstitial → CB partner mapping (matched by mamind + modId).
 ALLOCATE(cb_partner(nmamgc))
 cb_partner = -1
 IF (is_cbsim) THEN
   DO i = 1, nmamgc
     IF (mamgc(i)%iscb) CYCLE
     DO j = 1, nmamgc
       IF (mamgc(j)%iscb .and. &
           mamgc(j)%mamind == mamgc(i)%mamind .and. &
           mamgc(j)%modId  == mamgc(i)%modId ) THEN
         cb_partner(i) = j
         EXIT
       END IF
     END DO
   END DO
 END IF

! Allocate cloud-fraction history for activation/evaporation tendency.
 ALLOCATE(CLDF_prev(State_Grid%NX, State_Grid%NY, State_Grid%NZ))
 CLDF_prev = 0.0_fp

END SUBROUTINE MAM_INIT



!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: aero_drydep
!
! !DESCRIPTION: Subroutine AERO\_DRYDEP removes size-resolved aerosol number
!  and mass by dry deposition.  The deposition velocities are calcualted from
!  drydep_mod.f and only aerosol number NK01-NK30 are really treated as dry
!  depositing species while each of the mass species are depositing accordingly
!  with number.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MAM_SETTL( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : g0
    USE PhysConstants,      ONLY : AVO
    USE PRECISION_MOD
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object

! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  22 Jul 2007 - Win T. - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:


! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL,  SAVE     :: DOSETTLING = .True.
    LOGICAL,  SAVE     :: FIRST      = .TRUE.

    ! Scalars
    INTEGER            :: I,      J,        L,    n      
    REAL(fp)           :: DTCHEM, AREA_CM2, FLUX,  X,    Y
    REAL(fp)           :: DEN,   DP,   PDP
    REAL(fp)           :: TEMP,   P,        CONST, SLIP, VISC
    REAL(fp)           :: DELZ,   DELZ1,    TOT1,  TOT2

    ! Strings
    CHARACTER(LEN=255) :: LOC, MSG

    ! Arrays
    REAL(fp)           :: TC(State_Grid%NZ)
    REAL(fp)           :: TC0(State_Grid%NZ)
    REAL(fp)           :: VTS(State_Grid%NZ) ! Settling V [m/s]

    ! Pointers
    TYPE(SpcConc), POINTER  :: Spc     (:      )
    REAL(fp),      POINTER  :: BXHEIGHT(:,:,:  )
    REAL(fp),      POINTER  :: T       (:,:,:  )

    !=================================================================
    ! MAM_SETTL begins here!
    !=================================================================

    ! DTCHEM is the chemistry timestep in seconds
    DTCHEM    = GET_TS_CHEM()

    ! Initialize pointers
    Spc      => State_Chm%Species
    BXHEIGHT => State_Met%BXHEIGHT
    T        => State_Met%T

    !---------- GRAVITATIONAL SETTLING -------------
    !/ma
    ! First calculate vertical movement and removal by
    ! gravitational settling
    !
    ! Clarify units:
    !
    !      v_settling = rho   * Dp**2  *  g    *  C
    !                  -----------------------------
    !                   18    *  visc
    ! [units]
    !         m/s    = kg/m^3 *  m^2   * m/s^2  * -
    !                  -----------------------------
    !                    -    * kg/m/s
    !
    ! NOTES:
    ! (1 ) Pa s = kg/m/s
    ! (2 ) Slip correction factor is unitless, however, the
    !      equation from Hinds' Aerosol Technology that is
    !      a function of P and Dp needs the correct units
    !      P [=] kPa and Dp [=] um


    IF ( DOSETTLING ) THEN

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( N, I, J, DP, DEN, CONST, L, P, TEMP )   &
       !$OMP PRIVATE( PDP, SLIP, VISC, VTS, JC, ID, TC0, TC )   &
       !$OMP PRIVATE( DELZ, DELZ1, AREA_CM2, TOT1, TOT2, FLUX ) &
       !$OMP SCHEDULE( DYNAMIC )
       DO I = 1, State_Grid%NX
       DO J = 1, State_Grid%NY
       DO n = 1 , size(mamgc)

          ! CB species are not gravitationally settled here: they have no
          ! surface sink in this routine (Is_DryDep=false by design), and
          ! their removal is handled exclusively by stratiform rainout in
          ! wetscav_mod.F90. Settling them anyway drains mass into a
          ! sink-less L=1, producing unbounded fake accumulation there.
          IF (mamgc(n)%iscb) CYCLE

          DO L = 1, State_Grid%NZ
             
               IF(mamgc(n)%isnum) THEN
                   DP =State_Chm%GCMAM(mamgc(n)%ModId)  &
                                          %nuwetrad(I,J,L)*2.D6 ![=] um
               ELSE
                   DP  = State_Chm%GCMAM(mamgc(n)%ModId)  &
                                          %wetrad(I,J,L)*2.D6 ![=] um
               END IF
                   DEN   = State_Chm%GCMAM(mamgc(n)%ModId)  &
                                           %aerdens(I,J,L)

               CONST = DEN *  (DP*1.d-6)**2.d0 * g0 / 18.d0
             ! Get P [kPa], T [K], and P*DP
             ! Use moist pressure for mean free path (ewl, 3/2/2015)
               P    = State_Met%PMID(I,J,L) * 0.1d0  ![=] kPa
               TEMP = T(I,J,L)          ![=] K
               PDP  = P * DP

             !=====================================================
             ! # air molecule number density
             ! num = P * 1d3 * 6.023d23 / (8.314 * Temp)
             !
             ! # gas mean free path
             ! lamda = 1.d6 /
             !     &   ( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
             !
             ! # Slip correction
             ! Slip = 1. + 2. * lamda * (1.257 + 0.4 *
             !      &  exp( -1.1 * Dp / (2. * lamda))) / Dp
             !=====================================================
             ! NOTE, Slip correction factor calculations following
             !       Seinfeld, pp464 which is thought to be more
             !       accurate but more computation required.
             !=====================================================

             ! Slip correction factor as function of (P*dp)
             SLIP = 1d0 + ( 15.60d0 + 7.0d0 * EXP(-0.059d0*PDP) ) / PDP

             !=====================================================
             ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
             ! which produce slip correction factor with small
             ! error compared to the above with less computation.
             !=====================================================

             ! Viscosity [Pa s] of air as a function of temp (K)
             ! Sutherland eqn. (ref. pp 25 in Hinds (Aerosol Technology)
             VISC = 1.458d-6 * (TEMP)**(1.5d0) / ( TEMP + 110.4d0 )

             ! Settling velocity [m/s]
             VTS(L) = CONST * SLIP / VISC

             ! Method is to solve bidiagonal matrix
             ! which is implicit and first order accurate in Z

             TC0(L) = Spc(mamgc(n)%gcind)%Conc(I,J,L)
             TC(L)  = TC0(L)
          ENDDO  !L-loop

          ! We know the boundary condition at L = model top
  !FAB TEST        L     = State_Grid%MaxChemLev
          L     = State_Grid%NZ
          DELZ  = BXHEIGHT(I,J,L)           ![=] meter, model top
          TC(L) = TC(L) / ( 1.d0 + DTCHEM * VTS(L) / DELZ )

          DO L = State_Grid%MaxChemLev-1, 2, -1
                DELZ  = BXHEIGHT(I,J,L)
                DELZ1 = BXHEIGHT(I,J,L+1)
                TC(L) = 1.d0 / &
                      ( 1.d0   + DTCHEM * VTS(L)   / DELZ ) * &
                      ( TC(L)  + DTCHEM * VTS(L+1) / DELZ1  *  TC(L+1) )
          ENDDO

          ! Surface layer: only receive the settling flux from L=2, with
          ! no additional removal term here -- actual surface deposition
          ! is handled separately by the standard dry deposition scheme
          ! (drydep_mod.F90 / mixing_mod.F90), which already includes a
          ! gravitational settling contribution to Vd. Removing mass here
          ! too would double-count that surface loss.
          DELZ1   = BXHEIGHT(I,J,2)
          TC(1)   = TC(1) + DTCHEM * VTS(2) / DELZ1 * TC(2)

          DO L = 1, State_Grid%NZ
                Spc(mamgc(n)%gcind)%Conc(I,J,L) = TC(L)
          ENDDO

       ENDDO  ! MAMGC species (transported MAM species) 
       ENDDO  ! I-loop
       ENDDO  ! J-loop
       !$OMP END PARALLEL DO
  
    ENDIF  ! DOSETTLING

END SUBROUTINE MAM_SETTL

!---------------------------------------------------------------------------------------------

SUBROUTINE Set_MAM_Diagnostic( Input_Opt,  State_Chm, State_Diag, &
                                     State_Grid, State_Met, RC )
USE ErrCode_Mod
USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : Species, SpcConc
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag  ! Diagnostic State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    ! Scalars
    INTEGER                  :: I, J, L, M, iSlot
    ! Strings
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=512)       :: ErrMsg
    ! Convert [kg/m3] to [ug/m3]
    REAL(fp),      PARAMETER :: kgm3_to_ugm3 = 1.0e+9_fp
! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Set_AerMass_Diagnostic (in module GeosCore/aerosol_mod.F90)'
    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED   ) &
    !$OMP PRIVATE( I, J, L, M, iSlot )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
! Start with state_chm%GCMAM diag     
    
! modal
      IF ( State_Diag%Archive_MamNu ) THEN
          DO m = 1, size(State_Chm%GCMAM)
            iSlot = State_Diag%Map_MamNu%id2slot(m)
            IF (iSlot > 0)  &
            State_Diag%MamNu(I,J,L,iSlot) = State_Chm%GCMAM(m)%Nu(I,J,L) * 1.0e-6_fp !#m-3 to #cm-3
           END DO
      ENDIF
!now everything is set up to have modal species concentration diag as well
!total concentrations
      IF ( State_Diag%Archive_MamSO4Mass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamSO4Mass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lso4) &
                   State_Diag%MamSO4Mass(I,J,L,iSlot) = State_Chm%GCMAM(m)%so4(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamBCMass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamBCMass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lbc) &
                   State_Diag%MamBCMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%bc(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamPOMMass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamPOMMass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lpom) &
                   State_Diag%MamPOMMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%pom(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamSOAMass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamSOAMass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lsoa) &
                   State_Diag%MamSOAMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%soa(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamSSLTMass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamSSLTMass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lsslt) &
                   State_Diag%MamSSLTMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%sslt(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamDUSTMass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamDUSTMass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%ldust) &
                   State_Diag%MamDUSTMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%dust(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamNH4Mass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamNH4Mass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lnh4) &
                   State_Diag%MamNH4Mass(I,J,L,iSlot) = State_Chm%GCMAM(m)%nh4(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamNO3Mass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamNO3Mass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lno3) &
                   State_Diag%MamNO3Mass(I,J,L,iSlot) = State_Chm%GCMAM(m)%no3(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamCAMass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamCAMass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lca) &
                   State_Diag%MamCAMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%ca(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamCO3Mass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamCO3Mass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lco3) &
                   State_Diag%MamCO3Mass(I,J,L,iSlot) = State_Chm%GCMAM(m)%co3(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamCLMass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamCLMass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lcl) &
                   State_Diag%MamCLMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%cl(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamMOMMass ) THEN
           DO m = 1, size(State_Chm%GCMAM)
               iSlot = State_Diag%Map_MamMOMMass%id2slot(m)
               IF (iSlot > 0 .AND. State_Chm%GCMAM(m)%lmom) &
                   State_Diag%MamMOMMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%mom(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF

  IF ( State_Diag%Archive_MamWATMass ) THEN
       DO m = 1, size(State_Chm%GCMAM)
           iSlot = State_Diag%Map_MamWATMass%id2slot(m)
           IF (iSlot > 0) &
               State_Diag%MamWATMass(I,J,L,iSlot) = State_Chm%GCMAM(m)%aerwat(I,J,L) * kgm3_to_ugm3
       END DO
  ENDIF

     IF ( State_Diag%Archive_Mamwetrad ) THEN
         DO m = 1, size(State_Chm%GCMAM)
             iSlot = State_Diag%Map_Mamwetrad%id2slot(m)
             IF (iSlot > 0) &
                 State_Diag%Mamwetrad(I,J,L,iSlot) = State_Chm%GCMAM(m)%wetrad(I,J,L)
         END DO
     ENDIF

      IF ( State_Diag%Archive_Mamdryrad ) THEN
          DO m = 1, size(State_Chm%GCMAM)
              iSlot = State_Diag%Map_Mamdryrad%id2slot(m)
              IF (iSlot > 0) &
                  State_Diag%Mamdryrad(I,J,L,iSlot) = State_Chm%GCMAM(m)%dryrad(I,J,L)
          END DO
      ENDIF

      IF ( State_Diag%Archive_Mamhygro ) THEN
          DO m = 1, size(State_Chm%GCMAM)
              iSlot = State_Diag%Map_Mamhygro%id2slot(m)
              IF (iSlot > 0) &
                  State_Diag%Mamhygro(I,J,L,iSlot) = State_Chm%GCMAM(m)%hygro(I,J,L)
          END DO
      ENDIF

      IF ( State_Diag%Archive_MamPH ) THEN
          DO m = 1, size(State_Chm%GCMAM)
              iSlot = State_Diag%Map_MamPH%id2slot(m)
              IF (iSlot > 0) &
                  State_Diag%MamPH(I,J,L,iSlot) = State_Chm%GCMAM(m)%hplus(I,J,L)
          END DO
      ENDIF
! aerosol optical properties

        ! modal
      IF ( State_Diag%Archive_MamTauxarv ) THEN
          DO m = 1, size(State_Chm%GCMAM)
            iSlot = State_Diag%Map_MamTauxarv%id2slot(m)
            IF (iSlot > 0)  &
            State_Diag%MamTauxarv(I,J,L,iSlot) = State_Chm%GCMAM(m)%tauxar(I,J,L,10) !visible band
           END DO
      ENDIF
      IF ( State_Diag%Archive_Mamssav ) THEN
          DO m = 1, size(State_Chm%GCMAM)
            iSlot = State_Diag%Map_Mamssav%id2slot(m)
            IF (iSlot > 0)  &
            State_Diag%Mamssav(I,J,L,iSlot) = State_Chm%GCMAM(m)%ssa(I,J,L,10) !visible band
           END DO
      ENDIF
      IF ( State_Diag%Archive_Mamgv ) THEN
          DO m = 1, size(State_Chm%GCMAM)
            iSlot = State_Diag%Map_Mamgv%id2slot(m)
            IF (iSlot > 0)  &
            State_Diag%Mamgv(I,J,L,iSlot) = State_Chm%GCMAM(m)%g(I,J,L,10) !visible band
           END DO
      ENDIF

   ENDDO
   ENDDO
   ENDDO
END SUBROUTINE Set_MAM_Diagnostic
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: apply_rainout_eff
!
! !DESCRIPTION: Subroutine APPLY\_RAINOUT\_EFF multiplies the rainout fraction
!  computed by RAINOUT with the rainout efficiency for one of 3 temperature
!  ranges: (1) T < 237 K; (2) 237 K <= T < 258 K; (3) T > 258 K. The rainout
!  efficiencies for each aerosol species are defined in the species database
!  object (i.e. State\_Chm%SpcData(:)%Info).
!\\
!\\
!  This allows us to apply the impaction scavenging of certain aerosol species
!  (BC, dust, HNO3) as implemented by Qiaoqiao Wang, while also suppressing
!  rainout for other aerosol species.  The prior code achieved this by using
!  a large and confusing IF statement, whose logic was hard to understand.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MAM_APPLY_RAINOUT_EFF( hygro, nudryrad, TK, SpcInfo, RainFrac )
!
! !USES:
!
    USE Species_Mod, ONLY : Species
!
! !INPUT PARAMETERS:
!
    REAL(fp),      INTENT(IN)    :: TK         ! Temperature [K]
    REAL(fp),      INTENT(IN)    :: hygro      ! mode hygroscopicity
    REAL(fp),      INTENT(IN)    :: nudryrad   ! number-mode dry radius [m]
    TYPE(Species), INTENT(IN)    :: SpcInfo    ! Species Database object
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),      INTENT(INOUT) :: RainFrac   ! Rainout fraction

! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: RainoutEff, dgn_dry, kappa_ref
    ! Kohler theory constants: A=(2.1e-9 m)^3, S_s=0.15% -> Ss2=(1.5e-3)^2
    REAL(fp), PARAMETER :: A3       = 9.261e-27_fp
    REAL(fp), PARAMETER :: Ss2      = 2.25e-6_fp
    REAL(fp), PARAMETER :: luo_norm = 153.5_fp

    Rainouteff = 0.0_fp

    IF ( TK < 237.0_fp ) THEN
       ! Pure ice regime: ice nucleation controls, kappa irrelevant.
       ! Retain empirical fixed values (Luo et al., 2020).
       SELECT CASE ( SpcInfo%MamModId )
         CASE (1)  ! Accumulation
           IF ( hygro > 0.5_fp ) THEN
             Rainouteff = 0.4_fp
           ELSE
             Rainouteff = 0.6_fp
           END IF
         CASE (2)  ! Aitken
           Rainouteff = 0.08_fp
         CASE (3)  ! Coarse (dust-dominated)
           IF ( hygro > 0.4_fp ) THEN
             Rainouteff = 0.4_fp
           ELSE
             Rainouteff = 1.0_fp
           END IF
         CASE (4)  ! Primary carbon (hydrophobic)
           Rainouteff = 0.5_fp
       END SELECT

    ELSE IF ( TK < 258.0_fp ) THEN
       ! Mixed-phase regime (237-258 K).
       ! Coarse (dust) and primary carbon: ice nucleation fraction from
       ! DeMott et al. (2015) via Luo et al. (2020). Primary carbon uses
       ! 50% of dust value (Luo et al., 2020).
       ! Accum and Aitken: CCN activation into supercooled liquid droplets
       ! still dominates -> kappa-based mapping, no Luo correction.
       SELECT CASE ( SpcInfo%MamModId )
         CASE (3)  ! Coarse: dust ice nucleation efficiency x Luo correction
           IF ( hygro > 0.4_fp ) THEN
             Rainouteff = 0.4_fp
           ELSE
             Rainouteff = 1.0_fp
           END IF
           Rainouteff = Rainouteff * &
              ( EXP( 0.46_fp * ( 273.16_fp - TK ) - 11.6_fp ) / luo_norm )
         CASE (4)  ! Primary carbon: 50% of dust ice nucleation fraction
           Rainouteff = 0.5_fp * &
              ( EXP( 0.46_fp * ( 273.16_fp - TK ) - 11.6_fp ) / luo_norm )
         CASE DEFAULT  ! Accum, Aitken: kappa-based liquid CCN activation
           dgn_dry   = nudryrad * 2.0_fp
           kappa_ref = 4.0_fp * A3 / ( 27.0_fp * dgn_dry**3 * Ss2 )
           Rainouteff = MIN( 1.0_fp, hygro / kappa_ref )
       END SELECT

    ELSE
       ! Warm liquid clouds (T >= 258 K): CCN activation controls.
       ! kappa_ref from Kohler theory at mode dry diameter and S_s=0.15%.
       ! Coarse mode: kappa_ref << kappa -> always activates.
       SELECT CASE ( SpcInfo%MamModId )
         CASE (3)
           Rainouteff = 1.0_fp
         CASE DEFAULT
           dgn_dry   = nudryrad * 2.0_fp
           kappa_ref = 4.0_fp * A3 / ( 27.0_fp * dgn_dry**3 * Ss2 )
           Rainouteff = MIN( 1.0_fp, hygro / kappa_ref )
       END SELECT

    END IF

    RainFrac = RainFrac * Rainouteff

  END SUBROUTINE MAM_APPLY_RAINOUT_EFF
SUBROUTINE MAM_OPT_to_RRTMG( Input_Opt,  State_Chm,  State_Diag, &
                                State_Grid)

    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE Input_Opt_Mod,  ONLY : OptInput
    USE mam_opt,        ONLY : mamoptdiag
    USE radconstants,   ONLY : nswbands, nlwbands
    USE physconst,      ONLY : rga
    USE precision_mod,  ONLY : f8

    ! =========================================================================
    ! Arguments
    ! =========================================================================
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(DgnState), INTENT(INOUT) :: State_Diag

    ! =========================================================================
    ! FAB : Here is the tricky part : MAM Species index mapping in RRTMG RTODAER,SSA,ASY  
    ! MAM species must be consistent with GC RRTMG standard aerosol species ordering and aggregation 
    ! see SPECMASK in rrtmg_rad_transfer. Here we calculate and pass directly the MAM modal mean 
    ! of each relevant species to the RT tables. To avoid multiple count in rrtmg_rad_trnsfer 
    ! when standard species are aggregated ( e.g 7 bins of dust,2 bins of SSALT, etc), we just 
    ! fill one of the index with MAM values and zero out the others.  
    ! Again see  SPECMASK definition  in rrtmg_rad_transfer.  
    ! =========================================================================
    !RTXX, MASK, IS index     MAM speices ( averaged over modes) 
    !  = 1                    ! Sulfate
    !  = 2                    ! Nitrate    (MOSAIC_SPECIES only)
    !  = 3                    ! Ammonium   (MOSAIC_SPECIES only)
    !  = 4                    ! Black carbon
    !  = 5                    ! Organic aerosol (here POM + SOA [+ MOM], revisit this equivalence)
    !  = 6-7                  ! Sea salt : only 6 filled
    !  = 10-16                ! Dust: only 10 filled
    !indices 8,9 corresponds to stratospheric aerosol : need more thoughts on that 
    !for now zeroed. Rq: MAM is active in stratosphere.  
    ! =========================================================================
    ! Local variables
    ! =========================================================================
    REAL(f8), POINTER :: RTODAER   (:,:,:,:,:)
    REAL(f8), POINTER :: RTSSAER   (:,:,:,:,:)
    REAL(f8), POINTER :: RTASYMAER (:,:,:,:,:)

    INTEGER  :: NBNDS, IB, IBX, IB_SW  
    INTEGER  :: I, J, L, n, m, nmodes

    ! Accumulators for bulk optical properties across modes (tau, tau*ssa, tau*ssa*g)
    ! Using explicit variables rather than filling RTODAER directly so we can
    ! normalise SSA and g correctly before storing.
    REAL(f8) :: tau_s, wa_s, ga_s     ! extinction, single-scatter, asymmetry accumulators
    REAL(f8) :: EXT_MIN               ! small number to avoid divide-by-zero

    PARAMETER ( EXT_MIN = 1.0e-40_f8 )

    ! =========================================================================
    ! Set up pointers and dimensions
    ! =========================================================================
    NBNDS  = nswbands + nlwbands   ! must match RRTMG configuration

    RTODAER    => State_Chm%Phot%RTODAER
    RTSSAER    => State_Chm%Phot%RTSSAER
    RTASYMAER  => State_Chm%Phot%RTASYMAER

    nmodes = SIZE( mamoptdiag )

    ! =========================================================================
    ! Loop over all RRTMG bands (LW first, then SW — skip LW here)
    ! =========================================================================
    DO IB = 1, NBNDS

       ! RRTMG waveband slots start after NWVAA0 standard wavelengths in GC arrays
       IBX = IB + State_Chm%Phot%NWVAA0

       IF ( IB > nlwbands ) THEN

          ! SW band index into mamoptdiag arrays (1-based)
          IB_SW = IB - nlwbands

          !$OMP PARALLEL DO                 &
          !$OMP DEFAULT( SHARED )           &
          !$OMP PRIVATE( I, J, L, n, m,    &
          !$OMP          tau_s, wa_s, ga_s )&
          !$OMP SCHEDULE( DYNAMIC )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             ! zero out all aerosol species in RT tables 
             RTODAER  (I,J,L,IBX,1:State_Chm%Phot%NASPECRAD) =  0.0_f8 
             RTSSAER  (I,J,L,IBX,1:State_Chm%Phot%NASPECRAD) =  0.0_f8
             RTASYMAER(I,J,L,IBX,1:State_Chm%Phot%NASPECRAD) =  0.0_f8

             ! Linear column index used by mamoptdiag (pcols ordering)
             n = J + ( I - 1 ) * State_Grid%NY
             ! In order to use the existing GC/RRTMG interface, 
             ! Rad. Species are expected in the exact same order as set in
             ! Set_SpecMask from the GeosCore/rrtm_rad_transfer_mode.F90  
             ! ==============================================================
             ! IS = 1 : Sulfate
             ! ==============================================================
             tau_s = 0.0_f8 ;  wa_s = 0.0_f8 ;  ga_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_sulfate(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_sulfate(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_sulfate(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_sulfate(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_sulfate(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_sulfate(n,L,IB_SW)
             END DO
             ! Convert specific extinction (m2/kg air) -> AOD for the layer
             RTODAER  (I,J,L,IBX,1) = tau_s * physta%pdeldry(n,L) * rga
             RTSSAER  (I,J,L,IBX,1) = wa_s  / MAX( tau_s, EXT_MIN )
             RTASYMAER(I,J,L,IBX,1) = ga_s  / MAX( wa_s,  EXT_MIN )

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
             ! =============================================================
             ! IS = 2 : Nitrate  (MOSAIC_SPECIES only)
             ! =============================================================
             tau_s = 0.0_f8 ;  wa_s = 0.0_f8 ;  ga_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_no3(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_no3(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_no3(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_no3(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_no3(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_no3(n,L,IB_SW)
             END DO
             RTODAER  (I,J,L,IBX,2) = tau_s * physta%pdeldry(n,L) * rga
             RTSSAER  (I,J,L,IBX,2) = wa_s  / MAX( tau_s, EXT_MIN )
             RTASYMAER(I,J,L,IBX,2) = ga_s  / MAX( wa_s,  EXT_MIN )

             ! =============================================================
             ! IS = 3 : Ammonium  (MOSAIC_SPECIES only)
             ! =============================================================
             tau_s = 0.0_f8 ;  wa_s = 0.0_f8 ;  ga_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_nh4(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_nh4(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_nh4(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_nh4(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_nh4(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_nh4(n,L,IB_SW)
             END DO
             RTODAER  (I,J,L,IBX,3) = tau_s * physta%pdeldry(n,L) * rga
             RTSSAER  (I,J,L,IBX,3) = wa_s  / MAX( tau_s, EXT_MIN )
             RTASYMAER(I,J,L,IBX,3) = ga_s  / MAX( wa_s,  EXT_MIN )
#else
             ! Without MOSAIC_SPECIES, zero out NO3 and NH4 slots
             RTODAER  (I,J,L,IBX,3) = 0.0_f8
             RTSSAER  (I,J,L,IBX,3) = 0.0_f8
             RTASYMAER(I,J,L,IBX,3) = 0.0_f8
             RTODAER  (I,J,L,IBX,2) = 0.0_f8
             RTSSAER  (I,J,L,IBX,2) = 0.0_f8
             RTASYMAER(I,J,L,IBX,2) = 0.0_f8
#endif

             ! =============================================================
             ! IS = 4 : Black Carbon
             ! =============================================================
             tau_s = 0.0_f8 ;  wa_s = 0.0_f8 ;  ga_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_bc(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_bc(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_bc(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_bc(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_bc(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_bc(n,L,IB_SW)
             END DO
             RTODAER  (I,J,L,IBX,4) = tau_s * physta%pdeldry(n,L) * rga
             RTSSAER  (I,J,L,IBX,4) = wa_s  / MAX( tau_s, EXT_MIN )
             RTASYMAER(I,J,L,IBX,4) = ga_s  / MAX( wa_s,  EXT_MIN )

             ! =============================================================
             ! IS = 5 : Organic Aerosol  (POM + SOA + MOM if available)
             ! All sub-types are aggregated before normalization so that the
             ! bulk SSA and g are extinction-weighted correctly.
             ! =============================================================
             tau_s = 0.0_f8 ;  wa_s = 0.0_f8 ;  ga_s = 0.0_f8
             DO m = 1, nmodes
                ! POM
                tau_s = tau_s + mamoptdiag(m)%vext_pom(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_pom(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_pom(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_pom(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_pom(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_pom(n,L,IB_SW)
                ! SOA
                tau_s = tau_s + mamoptdiag(m)%vext_soa(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_soa(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_soa(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_soa(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_soa(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_soa(n,L,IB_SW)
                ! MOM (marine organic matter) 
                tau_s = tau_s + mamoptdiag(m)%vext_mom(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_mom(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_mom(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_mom(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_mom(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_mom(n,L,IB_SW)
             END DO
             RTODAER  (I,J,L,IBX,5) = tau_s * physta%pdeldry(n,L) * rga
             RTSSAER  (I,J,L,IBX,5) = wa_s  / MAX( tau_s, EXT_MIN )
             RTASYMAER(I,J,L,IBX,5) = ga_s  / MAX( wa_s,  EXT_MIN )

             ! =============================================================
             ! IS = 6-7 : Sea Salt- only 6 is filled
             ! =============================================================
             tau_s = 0.0_f8 ;  wa_s = 0.0_f8 ;  ga_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_seasalt(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_seasalt(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_seasalt(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_seasalt(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_seasalt(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_seasalt(n,L,IB_SW)
             END DO
             RTODAER  (I,J,L,IBX,6) = tau_s * physta%pdeldry(n,L) * rga
             RTSSAER  (I,J,L,IBX,6) = wa_s  / MAX( tau_s, EXT_MIN )
             RTASYMAER(I,J,L,IBX,6) = ga_s  / MAX( wa_s,  EXT_MIN )

             ! =============================================================
             ! IS = 10-16 : Dust - only 10 is filled
             ! =============================================================
             tau_s = 0.0_f8 ;  wa_s = 0.0_f8 ;  ga_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_dust(n,L,IB_SW)
                wa_s  = wa_s  + mamoptdiag(m)%vext_dust(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_dust(n,L,IB_SW)
                ga_s  = ga_s  + mamoptdiag(m)%vext_dust(n,L,IB_SW) &
                                * mamoptdiag(m)%vssa_dust(n,L,IB_SW) &
                                * mamoptdiag(m)%vasm_dust(n,L,IB_SW)
             END DO
             RTODAER  (I,J,L,IBX,10) = tau_s * physta%pdeldry(n,L) * rga
             RTSSAER  (I,J,L,IBX,10) = wa_s  / MAX( tau_s, EXT_MIN )
             RTASYMAER(I,J,L,IBX,10) = ga_s  / MAX( wa_s,  EXT_MIN )

          END DO  ! I
          END DO  ! J
          END DO  ! L
          !$OMP END PARALLEL DO  !

       ELSE
          ! ==================================================================
          ! LW bands  ( IB = 1 .. nlwbands )
          ! In the LW, RRTMG only requires absorption optical depth.
          ! RTSSAER and RTASYMAER are zeroed; only RTODAER is filled.
          ! vext_lw_* [m2/kg_air] = abs_lw_interp * specmmr, accumulated in
          ! mam_aero_lw over modes.  The same pdeldry*rga conversion as SW
          ! gives layer AOD.  Species IS-index mapping is identical to SW.
          ! ==================================================================

          !$OMP PARALLEL DO                 &
          !$OMP DEFAULT( SHARED )           &
          !$OMP PRIVATE( I, J, L, n, m,    &
          !$OMP          tau_s )            &
          !$OMP SCHEDULE( DYNAMIC )
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Zero all RT tables for this LW band (SSA/g not used in LW)
             RTODAER  (I,J,L,IBX,1:State_Chm%Phot%NASPECRAD) = 0.0_f8
             RTSSAER  (I,J,L,IBX,1:State_Chm%Phot%NASPECRAD) = 0.0_f8
             RTASYMAER(I,J,L,IBX,1:State_Chm%Phot%NASPECRAD) = 0.0_f8

             ! Linear column index used by mamoptdiag (pcols ordering)
             n = J + ( I - 1 ) * State_Grid%NY

             ! =============================================================
             ! IS = 1 : Sulfate
             ! =============================================================
             tau_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_lw_sulfate(n,L,IB)
             END DO
             RTODAER(I,J,L,IBX,1) = tau_s * physta%pdeldry(n,L) * rga

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined MOSAIC_SPECIES ) )
             ! =============================================================
             ! IS = 2 : Nitrate  (MOSAIC_SPECIES only)
             ! =============================================================
             tau_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_lw_no3(n,L,IB)
             END DO
             RTODAER(I,J,L,IBX,2) = tau_s * physta%pdeldry(n,L) * rga

             ! =============================================================
             ! IS = 3 : Ammonium  (MOSAIC_SPECIES only)
             ! =============================================================
             tau_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_lw_nh4(n,L,IB)
             END DO
             RTODAER(I,J,L,IBX,3) = tau_s * physta%pdeldry(n,L) * rga
#endif
             ! IS = 2 and 3 remain 0.0 in the non-MOSAIC case (set above)

             ! =============================================================
             ! IS = 4 : Black Carbon
             ! =============================================================
             tau_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_lw_bc(n,L,IB)
             END DO
             RTODAER(I,J,L,IBX,4) = tau_s * physta%pdeldry(n,L) * rga

             ! =============================================================
             ! IS = 5 : Organic Aerosol  (POM + SOA + MOM if available)
             ! =============================================================
             tau_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_lw_pom(n,L,IB) &
                              + mamoptdiag(m)%vext_lw_soa(n,L,IB)
                tau_s = tau_s + mamoptdiag(m)%vext_lw_mom(n,L,IB)
             END DO
             RTODAER(I,J,L,IBX,5) = tau_s * physta%pdeldry(n,L) * rga

             ! =============================================================
             ! IS = 6-7 : Sea Salt — only 6 is filled
             ! =============================================================
             tau_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_lw_seasalt(n,L,IB)
             END DO
             RTODAER(I,J,L,IBX,6) = tau_s * physta%pdeldry(n,L) * rga

             ! =============================================================
             ! IS = 10-16 : Dust — only 10 is filled
             ! =============================================================
             tau_s = 0.0_f8
             DO m = 1, nmodes
                tau_s = tau_s + mamoptdiag(m)%vext_lw_dust(n,L,IB)
             END DO
             RTODAER(I,J,L,IBX,10) = tau_s * physta%pdeldry(n,L) * rga

          END DO  ! I
          END DO  ! J
          END DO  ! L
          !$OMP END PARALLEL DO
       END IF  ! SW / LW bands

    END DO  ! IB

    RTODAER    => NULL()
    RTSSAER    => NULL()
    RTASYMAER  => NULL()

  END SUBROUTINE MAM_OPT_to_RRTMG


  FUNCTION MamMassFrac( r_min, r_max, r_g, sigma_g ) RESULT( f )
  USE, INTRINSIC :: ISO_C_BINDING
  REAL(fp), INTENT(IN) :: r_min, r_max  ! bin edges [same units as r_g]
  REAL(fp), INTENT(IN) :: r_g           ! geometric mean radius
  REAL(fp), INTENT(IN) :: sigma_g       ! geometric standard deviation
  REAL(fp)             :: f
  REAL(fp)             :: ln_sg, cdf_max, cdf_min

  ln_sg   = LOG( sigma_g )
  cdf_min = 0.5_fp * ( 1.0_fp + ERF( LOG(r_min/r_g) / (SQRT(2.0_fp) * ln_sg) ) )
  cdf_max = 0.5_fp * ( 1.0_fp + ERF( LOG(r_max/r_g) / (SQRT(2.0_fp) * ln_sg) ) )
  f       = cdf_max - cdf_min

  END FUNCTION






!------------------------------------------------------------------------------
SUBROUTINE MAM_to_HETRATES( Input_Opt, State_Chm, State_Grid, State_Met )
!
! Replace legacy aerosol fields in State_Chm that feed KPP heterogeneous
! chemistry with MAM-derived equivalents.  Called from chemistry_mod.F90
! after RDAER and RDust_Online, under
! #if ( defined MODAL_AERO_4MODE || defined MODAL_AERO_4MODE_MOM )
! (FAB, MAM-decouple-std Step 1: was MODAL_AERO_4MODE_MOM only).
!
! FAB (MAM-decouple-std, Step 2): composition-resolved surface areas.
!   Previously all MAM wet SA went into the SUL slot and every other slot
!   (dust, BC, OC, SSA, SSC) was zeroed, which switched off halogen het chem
!   on sea salt, N2O5/HO2 uptake on dust, NO3 hydrolysis, etc.
!   Now, for each INTERSTITIAL mode m (cloud-borne mass has no particle
!   surface: it sits in droplets), the wet SA S_m (GCMAM%saer) is split among
!   the GEOS-Chem aerosol-type slots (fullchem_RateLawFuncs.F90:30-42) by the
!   component dry-volume fraction f_c,m = V_c,m / sum_c V_c,m,
!   V = mass/density (MAM specdens_amode):
!     DU1..DU7 (1-7) <- dust, Ca, CO3      (bin with R_eff closest to mode r_eff)
!     SUL      (8)   <- SO4 (non-sea-salt part), NH4, NO3
!     BKC      (9)   <- BC
!     ORC      (10)  <- POM, SOA, MOM
!     SSA      (11)  <- Na (seasalt), Cl, sea-salt SO4  in accum + Aitken
!     SSC      (12)  <- Na (seasalt), Cl, sea-salt SO4  in coarse mode
!   NOTE: in the MOSAIC build the MAM "seasalt" type (GCMAM%sslt, MAMSSLT*
!   tracers, lptr_nacl_a_amode) holds Na+ only (MAM4 naming heritage; to be
!   renamed one day).
!   Sea-salt SO4 (FAB, Step 2 amendment): MAM sea salt is emitted with its own
!   primary SO4 (HEMCO SeaSalt ext., 'MAM SSA SO4 mass fraction'), which STD
!   counts inside SALA/SALC area. The primary part of each mode's SO4 is
!   estimated from non-volatile Na and moved to the SS category:
!     SO4_ss = min( SO4, SO4_to_Na_SS * Na ),  SO4_nss = SO4 - SO4_ss (-> SUL)
!   Secondary SO4 (gas/cloud/K_MT production) stays in SUL.
!   Slots 13-14 (strat. liquid aerosol, ice) are not touched.
!   Radius of each slot = SA-weighted wet EFFECTIVE radius of contributing
!   modes, r_eff = 3V/S = r_gv,wet * exp(-0.5 ln^2 sigma) (wetrad is the
!   volume-mean geometric radius), so xVol = xArea*xRadi/3 in
!   fullchem_HetStateFuncs is the true wet volume. A slot with zero area
!   keeps its legacy radius (rate is zero anyway).
!   Mode water is apportioned among components by kappa*V (ZSR additivity,
!   kappa = MAM spechygro).
!   Fine-inorganic fields keep their STD meaning (ISORROPIA/HETP "fine"
!   system = SNA + fine sea salt):
!     aClArea/aClRadi, AeroH2O(SUL), IsorropAeroH2O(1) <- SNA + SS share of
!     accumulation + Aitken modes; IsorropAeroH2O(2) <- SNA + SS share of
!     the coarse mode.
!
! Fields overwritten:
!   AeroArea / WetAeroArea (1:NDUST+NRHAER), AeroRadi / WetAeroRadi (same,
!                             where area > 0)
!   aClArea / aClRadi       -- fine inorganic SA and r_eff (N2O5_InorgOrg,
!                              NO3 hydrolysis on SALA)
!   AeroH2O(SUL,BKC,ORC,SSA,SSC) [g/m3] -- H%xH2O (N2O5 gamma)
!   IsorropAeroH2O(1:2)     [ug/m3]     -- H%AWATER (NO3 hydrolysis on SS)
!   IsorropHplus(1)         -- accumulation-mode H+ [mol/L]
!   IsorropAeropH(1)        -- pH = -log10(H+); sets H%H_conc_Sul for
!                              HOBr/HOCl/ClNO3 het reactions
!------------------------------------------------------------------------------
    USE CMN_SIZE_Mod,    ONLY : NDUST, NRHAER
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE precision_mod,   ONLY : fp
    USE modal_aero_data, ONLY : ntot_aspectype, specname_amode,            &
                                specdens_amode, spechygro, alnsg_amode,    &
                                modeptr_accum, modeptr_aitken, modeptr_coarse

    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met

    ! Het-chem slot indices (fullchem_RateLawFuncs.F90:30-42)
    INTEGER,  PARAMETER :: iSUL = NDUST+1, iBKC = NDUST+2, iORC = NDUST+3
    INTEGER,  PARAMETER :: iSSA = NDUST+4, iSSC = NDUST+5
    INTEGER,  PARAMETER :: NSLOT = NDUST+NRHAER
    ! Composition categories
    INTEGER,  PARAMETER :: cDST = 1, cSNA = 2, cBC = 3, cOM = 4, cSS = 5
    INTEGER,  PARAMETER :: NCAT = 5
    ! GEOS-Chem dust bin effective radii [um] (fullchem_RateLawFuncs.F90:29-35)
    REAL(fp), PARAMETER :: DU_Reff_um(NDUST) = (/ 0.151_fp, 0.253_fp,       &
                           0.402_fp, 0.818_fp, 1.491_fp, 2.417_fp, 3.721_fp /)
    ! MAM species types carried in GCMAM(m)%<field>
    INTEGER,  PARAMETER :: NSPC = 12
    CHARACTER(LEN=10), PARAMETER :: spcType(NSPC) = (/                      &
         'sulfate   ', 'ammonium  ', 'nitrate   ', 'black-c   ',            &
         'p-organic ', 's-organic ', 'm-organic ', 'seasalt   ',            &
         'chloride  ', 'dust      ', 'calcium   ', 'carbonate ' /)
    INTEGER,  PARAMETER :: spcCat(NSPC) = (/ cSNA, cSNA, cSNA, cBC,         &
                           cOM, cOM, cOM, cSS, cSS, cDST, cDST, cDST /)
    ! FAB: sea-salt SO4/Na emitted mass ratio = 'MAM SSA SO4 mass fraction' /
    ! 'MAM SSA Na mass fraction' in HEMCO_Config.rc (defaults 0.077/0.385,
    ! hcox_seasalt_mod.F90). Hard-coded: KEEP IN SYNC if those options change.
    REAL(fp), PARAMETER :: SO4_to_Na_SS = 0.077_fp / 0.385_fp

    REAL(fp) :: rDens(NSPC), kappa(NSPC)
    REAL(fp) :: vol(NCAT), kvol(NCAT), mass(NSPC)
    REAL(fp) :: area(NSLOT), areaR(NSLOT), wat(NSLOT)
    REAL(fp) :: vtot, kvtot, S_m, reff, wm, a, w, dmin, d
    REAL(fp) :: aclA, aclAR, watFine, watCoarse
    INTEGER  :: I, J, L, m, k, t, c, iDU, nmodes
    LOGICAL  :: isFine

    nmodes = SIZE( State_Chm%GCMAM )

    ! Inverse density [m3/kg] and kappa of each carried species type, taken
    ! from MAM (rad_cnst physprop); types absent from this build get 0.
    rDens = 0.0_fp
    kappa = 0.0_fp
    DO k = 1, NSPC
       DO t = 1, ntot_aspectype
          IF ( TRIM(specname_amode(t)) == TRIM(spcType(k)) ) THEN
             IF ( specdens_amode(t) > 0.0_fp .AND.                          &
                  specdens_amode(t) < 1.0e+20_fp ) THEN
                rDens(k) = 1.0_fp / specdens_amode(t)
                kappa(k) = MAX( spechygro(t), 0.0_fp )
             ENDIF
             EXIT
          ENDIF
       END DO
    END DO

    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED )                                                  &
    !$OMP PRIVATE( I, J, L, m, k, c, iDU, isFine, vol, kvol, mass, area )    &
    !$OMP PRIVATE( areaR, wat, vtot, kvtot, S_m, reff, wm, a, w, dmin, d )   &
    !$OMP PRIVATE( aclA, aclAR, watFine, watCoarse )                         &
    !$OMP COLLAPSE( 3 )                                                      &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       area      = 0.0_fp
       areaR     = 0.0_fp
       wat       = 0.0_fp
       aclA      = 0.0_fp
       aclAR     = 0.0_fp
       watFine   = 0.0_fp
       watCoarse = 0.0_fp

       DO m = 1, nmodes

          S_m = State_Chm%GCMAM(m)%saer(I,J,L)          ! [cm2/cm3]
          IF ( S_m <= 0.0_fp ) CYCLE

          ! Interstitial component masses [kg/m3] (allocated only if carried)
          mass = 0.0_fp
          IF ( State_Chm%GCMAM(m)%lso4  ) mass(1)  = State_Chm%GCMAM(m)%so4 (I,J,L)
          IF ( State_Chm%GCMAM(m)%lnh4  ) mass(2)  = State_Chm%GCMAM(m)%nh4 (I,J,L)
          IF ( State_Chm%GCMAM(m)%lno3  ) mass(3)  = State_Chm%GCMAM(m)%no3 (I,J,L)
          IF ( State_Chm%GCMAM(m)%lbc   ) mass(4)  = State_Chm%GCMAM(m)%bc  (I,J,L)
          IF ( State_Chm%GCMAM(m)%lpom  ) mass(5)  = State_Chm%GCMAM(m)%pom (I,J,L)
          IF ( State_Chm%GCMAM(m)%lsoa  ) mass(6)  = State_Chm%GCMAM(m)%soa (I,J,L)
          IF ( State_Chm%GCMAM(m)%lmom  ) mass(7)  = State_Chm%GCMAM(m)%mom (I,J,L)
          IF ( State_Chm%GCMAM(m)%lsslt ) mass(8)  = State_Chm%GCMAM(m)%sslt(I,J,L)
          IF ( State_Chm%GCMAM(m)%lcl   ) mass(9)  = State_Chm%GCMAM(m)%cl  (I,J,L)
          IF ( State_Chm%GCMAM(m)%ldust ) mass(10) = State_Chm%GCMAM(m)%dust(I,J,L)
          IF ( State_Chm%GCMAM(m)%lca   ) mass(11) = State_Chm%GCMAM(m)%ca  (I,J,L)
          IF ( State_Chm%GCMAM(m)%lco3  ) mass(12) = State_Chm%GCMAM(m)%co3 (I,J,L)

          ! Dry volume and kappa*volume per category
          vol  = 0.0_fp
          kvol = 0.0_fp
          DO k = 1, NSPC
             c       = spcCat(k)
             a       = MAX( mass(k), 0.0_fp ) * rDens(k)
             vol(c)  = vol(c)  + a
             kvol(c) = kvol(c) + kappa(k) * a
          END DO

          ! FAB: move primary sea-salt SO4 (estimated from Na = mass(8)) from
          ! the SNA to the SS category (see header)
          a = MIN( MAX( mass(1), 0.0_fp ),                                   &
                   SO4_to_Na_SS * MAX( mass(8), 0.0_fp ) ) * rDens(1)
          a = MIN( a, vol(cSNA) )
          vol(cSNA)  = vol(cSNA)  - a
          vol(cSS)   = vol(cSS)   + a
          kvol(cSNA) = MAX( kvol(cSNA) - kappa(1) * a, 0.0_fp )
          kvol(cSS)  = kvol(cSS)  + kappa(1) * a

          vtot  = SUM( vol  )
          kvtot = SUM( kvol )
          IF ( vtot <= 0.0_fp ) CYCLE

          ! Wet effective radius [cm]: wetrad [m] is volume-mean geometric
          reff = State_Chm%GCMAM(m)%wetrad(I,J,L)                           &
               * EXP( -0.5_fp * alnsg_amode(m)**2 ) * 1.0e+2_fp

          ! Mode aerosol water [kg/m3]: aerwat [kg/kg] * AIRDEN [kg/m3]
          wm = MAX( State_Chm%GCMAM(m)%aerwat(I,J,L), 0.0_fp )              &
             * State_Met%AIRDEN(I,J,L)

          isFine = ( m == modeptr_accum .OR. m == modeptr_aitken )

          ! Dust bin whose effective radius is closest (log space)
          iDU  = 1
          dmin = HUGE( 1.0_fp )
          DO k = 1, NDUST
             d = ABS( LOG( MAX( reff, 1.0e-10_fp ) * 1.0e+4_fp / DU_Reff_um(k) ) )
             IF ( d < dmin ) THEN
                dmin = d
                iDU  = k
             ENDIF
          END DO

          DO c = 1, NCAT
             IF ( vol(c) <= 0.0_fp ) CYCLE
             a = S_m * vol(c) / vtot
             IF ( kvtot > 0.0_fp ) THEN
                w = wm * kvol(c) / kvtot
             ELSE
                w = wm * vol(c) / vtot
             ENDIF
             SELECT CASE ( c )
                CASE ( cDST )
                   k = iDU
                CASE ( cSNA )
                   k = iSUL
                CASE ( cBC )
                   k = iBKC
                CASE ( cOM )
                   k = iORC
                CASE ( cSS )
                   k = iSSA
                   IF ( m == modeptr_coarse ) k = iSSC
             END SELECT
             area(k)  = area(k)  + a
             areaR(k) = areaR(k) + a * reff
             wat(k)   = wat(k)   + w

             ! Fine / coarse inorganic system (STD ISORROPIA meaning)
             IF ( c == cSNA .OR. c == cSS ) THEN
                IF ( isFine ) THEN
                   aclA    = aclA    + a
                   aclAR   = aclAR   + a * reff
                   watFine = watFine + w
                ELSE IF ( m == modeptr_coarse ) THEN
                   watCoarse = watCoarse + w
                ENDIF
             ENDIF
          END DO

       END DO   ! modes

       ! Surface areas [cm2/cm3] and effective radii [cm]
       State_Chm%AeroArea   (I,J,L,1:NSLOT) = area
       State_Chm%WetAeroArea(I,J,L,1:NSLOT) = area
       DO k = 1, NSLOT
          IF ( area(k) > 0.0_fp ) THEN
             State_Chm%AeroRadi   (I,J,L,k) = areaR(k) / area(k)
             State_Chm%WetAeroRadi(I,J,L,k) = areaR(k) / area(k)
          ENDIF
       END DO

       ! Aerosol water per slot [g/m3] (dust slots are not set in STD either)
       State_Chm%AeroH2O(I,J,L,iSUL) = watFine    * 1.0e+3_fp
       State_Chm%AeroH2O(I,J,L,iBKC) = wat(iBKC)  * 1.0e+3_fp
       State_Chm%AeroH2O(I,J,L,iORC) = wat(iORC)  * 1.0e+3_fp
       State_Chm%AeroH2O(I,J,L,iSSA) = wat(iSSA)  * 1.0e+3_fp
       State_Chm%AeroH2O(I,J,L,iSSC) = wat(iSSC)  * 1.0e+3_fp

       ! Fine inorganic SA [cm2/cm3] and r_eff [cm]
       State_Chm%aClArea(I,J,L) = aclA
       IF ( aclA > 0.0_fp ) State_Chm%aClRadi(I,J,L) = aclAR / aclA

       ! Inorganic aerosol water, fine / coarse [ug/m3] (H%AWATER)
       State_Chm%IsorropAeroH2O(I,J,L,1) = watFine   * 1.0e+9_fp
       State_Chm%IsorropAeroH2O(I,J,L,2) = watCoarse * 1.0e+9_fp

    END DO
    END DO
    END DO
    !$OMP END PARALLEL DO

    ! -----------------------------------------------------------------------
    ! IsorropHplus / IsorropAeropH: aerosol H+ and pH from MAM accumulation
    ! mode, replacing the ISORROPIA values (zero when bypassed).
    ! IsorropHplus(1) [mol/L] -> H%H_plus in fullchem_HetStateFuncs.
    ! IsorropAeropH(1) = -log10(H+) -> H%pHSSA(1) -> H%H_conc_Sul used by
    ! HOBr/HOCl/ClNO3 heterogeneous reactions.
    ! -----------------------------------------------------------------------
    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( I, J, L )  &
    !$OMP SCHEDULE( DYNAMIC )
    DO I = 1, State_Grid%NX
    DO J = 1, State_Grid%NY
    DO L = 1, State_Grid%NZ
       State_Chm%IsorropHplus(I,J,L,1)  = State_Chm%GCMAM(1)%hplus(I,J,L)
       State_Chm%IsorropAeropH(I,J,L,1) =                            &
            -LOG10( MAX( State_Chm%GCMAM(1)%hplus(I,J,L), 1.0e-30_fp ) )
    END DO
    END DO
    END DO
    !$OMP END PARALLEL DO

END SUBROUTINE MAM_to_HETRATES

!------------------------------------------------------------------------------
SUBROUTINE MAM_OPT_to_PHOTOL( Input_Opt, State_Chm, State_Grid )
!
! Replace ODAER and ODMDUST at 1000 nm (Fast-JX photolysis reference
! wavelength) with MAM optical depths from mamoptdiag.
!
! RDAER fills ODAER slots 1:NRHAER from legacy species masses; RDUST_ONLINE
! fills ODMDUST(1:NDUST) from standard dust bins.  This routine zeroes both
! sets and fills ODAER slot 1 with the total MAM extinction optical depth
! (tauxar summed over all modes).  tauxar already integrates every species in
! each mode (sulfate, BC, OC, sea salt, dust, ...) so no explicit species
! mapping is needed and coarse-mode dust is automatically included, avoiding
! double-counting with ODMDUST.  RDUST_ONLINE still runs for its surface-area
! and heterogeneous-chemistry outputs; only its optical depth is suppressed.
!
! Slot 1 retains the SNA Mie entry for Fast-JX UV spectral scaling — a
! reasonable bulk aerosol proxy.  Strat aerosol slots (NRHAER+1:NAER) are
! left untouched; they are filled afterwards by the strat aerosol loop.
!
! Precision note: the legacy code looks up Mie coefficients at exactly 1000 nm
! (one of the 11 discrete wavelengths in the GC aerosol optics dat files).
! MAM tauxar is a band-average over SW band 8 (8050-12850 cm-1, 778-1242 nm),
! a broader spectral range.  This is inherent to the MAM optics data structure
! which is discretised to RRTMG-SW bands, not individual wavelengths.
!------------------------------------------------------------------------------
    USE CMN_SIZE_Mod,   ONLY : NRHAER, NDUST
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
    USE mam_opt,        ONLY : mamoptdiag
    USE precision_mod,  ONLY : f8

    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid

    ! Band 8 of the 14 SW bands in the MAM optics data files (8050-12850 cm-1,
    ! 778-1242 nm).  Fixed by the optics file format; independent of whether
    ! RRTMG radiation is active at runtime.
    INTEGER, PARAMETER :: IB_1000 = 8

    REAL(f8), POINTER :: ODAER  (:,:,:,:,:)
    REAL(f8), POINTER :: ODMDUST(:,:,:,:,:)
    INTEGER           :: I, J, L, m, n, nmodes, IWV1000

    ODAER   => State_Chm%Phot%ODAER
    ODMDUST => State_Chm%Phot%ODMDUST
    IWV1000  = State_Chm%Phot%IWV1000
    nmodes   = SIZE( mamoptdiag )

    ! Zero legacy hygroscopic and dust optical depths at the 1000-nm index
    ODAER  (:,:,:,IWV1000,1:NRHAER) = 0.0_f8
    ODMDUST(:,:,:,IWV1000,1:NDUST)  = 0.0_f8

    ! Accumulate total MAM optical depth (all modes) into slot 1
    DO m = 1, nmodes
       !$OMP PARALLEL DO                &
       !$OMP DEFAULT( SHARED )          &
       !$OMP PRIVATE( I, J, L, n )     &
       !$OMP SCHEDULE( DYNAMIC )
       DO I = 1, State_Grid%NX
       DO J = 1, State_Grid%NY
          n = J + ( I - 1 ) * State_Grid%NY
          DO L = 1, State_Grid%NZ
             ODAER(I,J,L,IWV1000,1) = ODAER(I,J,L,IWV1000,1) &
                                     + mamoptdiag(m)%tauxar(n,L,IB_1000)
          END DO
       END DO
       END DO
       !$OMP END PARALLEL DO
    END DO

    ODAER   => NULL()
    ODMDUST => NULL()

END SUBROUTINE MAM_OPT_to_PHOTOL

  END MODULE MAM_DRIV_MOD
!FAB#endif

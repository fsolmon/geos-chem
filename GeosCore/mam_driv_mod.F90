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

PUBLIC :: MAM_DRIV, MAM_INIT, MAM_APPLY_RAINOUT_EFF 

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

    USE mam_utils, only: begchunk, endrun, mdo_coldstart,                            & 
                         mdo_gaschem, mdo_cloudchem,  mdo_gasaerexch,     mdo_rename,          &
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
                               modeptr_accum, alnsg_amode, voltonumb_amode
    USE modal_aero_initialize_data, only: MAM_cold_start 
    USE modal_aero_calcsize, only: modal_aero_calcsize_sub
    USE modal_aero_wateruptake, only: modal_aero_wateruptake_dr
    USE modal_aero_amicphys, only: modal_aero_amicphys_intr
    use mam_opt , only : mam_aero_sw, mamoptdiag
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
      real(r8)  :: tauxar(pcols,pver,nswbands)  ! aerosol extinction optical depth
      real(r8)  :: wa(pcols,pver,nswbands)      ! aerosol single scattering albedo * tau
      real(r8)  :: ga(pcols,pver,nswbands)      ! aerosol asymmetry parameter * wa
      real(r8)  :: fa(pcols,pver,nswbands)      ! aerosol forward scattered fraction * ga


      INTEGER :: latndx(pcols),lonndx(pcols)                 !required by the mam interface
                                                !not used now potentiall usefull for diags 

                                          
      CHARACTER(len=8) :: spcnam
!--------------------------------------------------------------------------


    IF ( lfirstcall ) then
           mamstep = 1
    ELSE 
           mamstep =mamstep+1
    END IF        
    ! Point to Spc
    Spc => State_Chm%Species
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
      physta%relhum(n,l)= State_Met%RH(I,J,L)
      physta%qv(n,l) = State_Met%SPHU(I,J,L) * 1.0e-3_r8 ! in kg/kgair  Caution here make sure
      physta%zm(n,l) =  sum(State_Met%BXHEIGHT(I,J,1:l))- 0.5 * State_Met%BXHEIGHT(I,J,l)

      IF (L==1) physta%pblh(n) = State_Met%PBLH(I,J)

      ! first element water vapor mr (used in water )  
      physta%q(n,l,1) = physta%qv(n,l) / (1.0e+0_f8 - physta%qv(n,l))

      ! load sulf production rate and convert from Kg.s-1  to kg.kg-1.s-1 
      ! needs fullchem activated 
      physta%ph2so4(n,l) = Spc(Ind_('PH2SO4'))%Conc(I,J,L) / State_Met%AD(I,J,L)/ deltat              

      physta%paqso4(n,l) =  Spc(Ind_('PSO4AQ'))%Conc(I,J,L)  / State_Met%AD(I,J,L) /deltat
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
     ! mdo_coldstart initialies in Mam_cold_start 
     if (mdo_coldstart == 1 .and. masterproc) print*, 'MAM q from namelist in coldstart '
!    DO L = 1, State_Grid%NZ
!    DO J = 1, State_Grid%NY
!    DO I = 1, State_Grid%NX ! 
!      n = J + (I-1)*State_Grid%NY
!FAB  try something temporary 
 !     physta%q(n,l,lptr_so4_a_amode(1)) = Spc(IND_('SO4'))%Conc(I,J,L)/State_Met%AD(I,J,L)
 !     physta%q(n,l,numptr_amode(1)) =    physta%q(n,l,lptr_so4_a_amode(1)) /1700. * voltonumb_amode(1)
!      Spc(IND_('MAMDEV'))%Conc(I,J,L) = Spc(IND_('SO4'))%Conc(I,J,L)
!    END DO
!    END DO
!    END DO
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
!
! apply tendencies ! note ptend.lq is modified by calcsize
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
     CALL modal_aero_wateruptake_dr( physta, pbuf, deltat, mamstep)
     
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

      IF (mdo_cloudchem > 0) then
        !start by updating mass in the accumulation mode 
        !(consider partitioning with aitken )
        l2 = lptr_so4_cw_amode(modeptr_accum) - loffset
        IF (is_cbsim) then
          vmrcw(1:pcols,1:pver,l2) = vmrcw(1:pcols,1:pver,l2) +  & 
                                   physta%paqso4(1:pcols,1:pver)*mwdry/adv_mass(l2)*deltat
        ELSE  ! consider all aersol intersticial
          vmr(1:pcols,1:pver,l2) = vmr(1:pcols,1:pver,l2) +  &
                                   physta%paqso4(1:pcols,1:pver)*mwdry/adv_mass(l2)*deltat
        END IF        
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
         physta%wetdens,      physta%qaerwat              )
    END IF
! vmr and vmrcw have been updated in modal_aero_amicphys_intr  
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
        ! modal mass concentrations in Kg.m-3  
        IF(lptr_so4_a_amode(m) > 0 ) State_Chm%GCMAM(m)%so4(I,J,L) =              & 
                                physta%q(n,L,lptr_so4_a_amode(m))*State_Met%AIRDEN(I,J,L) + &  
        !FAB TEMP add the cloud borne sulf to chm state for diag // change that once 
        ! transfer from qqcw to q is properly trated !!
                                physta%qqcw(n,L,lptr_so4_a_amode(m))*State_Met%AIRDEN(I,J,L)    

         
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
 
        IF(lptr_nh4_a_amode(m) > 0 ) State_Chm%GCMAM(m)%nh4(I,J,L) =                &
                               physta%q(n,L,lptr_nh4_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_no3_a_amode(m) > 0 ) State_Chm%GCMAM(m)%no3(I,J,L) =                &
                               physta%q(n,L,lptr_no3_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_ca_a_amode(m) > 0 ) State_Chm%GCMAM(m)%ca(I,J,L) =                  &
                               physta%q(n,L,lptr_ca_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_co3_a_amode(m) > 0 ) State_Chm%GCMAM(m)%co3(I,J,L) =                &
                               physta%q(n,L,lptr_co3_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_cl_a_amode(m) > 0 ) State_Chm%GCMAM(m)%cl(I,J,L) =                  &
                               physta%q(n,L,lptr_cl_a_amode(m))*State_Met%AIRDEN(I,J,L)

        IF(lptr_mom_a_amode(m) > 0 ) State_Chm%GCMAM(m)%mom(I,J,L) =                &
                               physta%q(n,L,lptr_mom_a_amode(m))*State_Met%AIRDEN(I,J,L)  

! optics                  
        State_Chm%GCMAM(m)%tauxar(I,J,L,:) = mamoptdiag(m)%tauxar(n,L,:)
        State_Chm%GCMAM(m)%ssa(I,J,L,:) = mamoptdiag(m)%ssa(n,L,:)
        State_Chm%GCMAM(m)%g(I,J,L,:) = mamoptdiag(m)%g(n,L,:)        
                 
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

    Spc => NULL() 
    IF (lfirstcall) lfirstcall = .false.
    
    
  END SUBROUTINE MAM_DRIV 

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
    INTEGER  :: m , i , s 
    CHARACTER * 12 :: tmp

    
    
!-------initialise namelist parameters

!   mdo_gaschem=1
!   mdo_cloudchem =1

!   mdo_gasaerexch=1
!   mdo_rename=1
!   mdo_newnuc=1
!   mdo_coag=1

   is_cbsim = .false.


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

 IF ( masterproc) then 
  print*, 'MAM species INFO'
  print*, mamgc(:)%name
  print*, mamgc(:)%gcind
  print*, mamgc(:)%mamind
  print*, mamgc(:)%modId
  print*, mamgc(:)%isnum
  print*, mamgc(:)%iscb
 ENDIF 
 IF (is_cbsim .and. .not.any(mamgc(:)%name(1:5)=='MAMCB')) then
     print*, 'CloudBorne aerosol simulation enabled but no MAMCBxx species present in Species list !' 
     stop    
 END IF         


END SUBROUTINE MAM_INIT        


!-------------------------------------------------------------------------------
SUBROUTINE load_pbuf( pbuf, lchnk, ncol,  &
         cld, qqcw, dgncur_a, dgncur_awet, qaerwat, wetdens,hygro)


      USE mam_utils, only: pcols,pver
      USE constituents, only : pcnst
      USE chem_mods, only: adv_mass, gas_pcnst, imozart
      USE physconst, only: mwdry

      USE modal_aero_data, only:  &
         lmassptrcw_amode, nspec_amode, numptrcw_amode, &
         qqcw_get_field, ntot_amode

      USE physics_buffer, only: physics_buffer_desc, &
         pbuf_get_index, pbuf_get_field

      TYPE(physics_buffer_desc), pointer :: pbuf(:)  ! physics buffer for a chunk

      INTEGER,  intent(in   ) :: lchnk, ncol

      REAL(r8), intent(in   ) :: cld(pcols,pver)    ! stratiform cloud fraction
      REAL(r8), intent(in   ) :: qqcw(pcols,pver,pcnst)  ! Cloudborne aerosol MR array
      REAL(r8), intent(in   ) :: dgncur_a(pcols,pver,ntot_amode)
      REAL(r8), intent(in   ) :: dgncur_awet(pcols,pver,ntot_amode)
      REAL(r8), intent(in   ) :: qaerwat(pcols,pver,ntot_amode)
      REAL(r8), intent(in   ) :: wetdens(pcols,pver,ntot_amode)
      REAL(r8), intent(in   ) :: hygro(pcols,pver,ntot_amode)
      INTEGER :: idx, l, ll, n

      REAL(r8), pointer :: fldcw(:,:)
      REAL(r8), pointer :: ycld(:,:)
      REAL(r8), pointer :: ydgnum(:,:,:)
      REAL(r8), pointer :: ydgnumwet(:,:,:)
      REAL(r8), pointer :: yqaerwat(:,:,:)
      REAL(r8), pointer :: ywetdens(:,:,:)
      REAL(r8), pointer :: yhygro(:,:,:)
 

 ! FAB ncol = pcols , maybe getrif of it   
      idx = pbuf_get_index( 'CLD' )
      CALL pbuf_get_field( pbuf, idx, ycld )
      ycld(:,:) = 0.0e+0_f8
      ycld(1:ncol,:) = cld(1:ncol,:)
      
      idx = pbuf_get_index( 'DGNUM' )
      CALL pbuf_get_field( pbuf, idx, ydgnum )
      ydgnum(:,:,:) = 0.0e+0_f8
      ydgnum(1:ncol,:,:) = dgncur_a(1:ncol,:,:)
      
      idx = pbuf_get_index( 'DGNUMWET' )
      CALL pbuf_get_field( pbuf, idx, ydgnumwet )
      ydgnumwet(:,:,:) = 0.0e+0_f8
      ydgnumwet(1:ncol,:,:) = dgncur_awet(1:ncol,:,:)
      
      idx = pbuf_get_index( 'QAERWAT' )
      CALL pbuf_get_field( pbuf, idx, yqaerwat )
      yqaerwat(:,:,:) = 0.0e+0_f8
      yqaerwat(1:ncol,:,:) = qaerwat(1:ncol,:,:)
      
      idx = pbuf_get_index( 'WETDENS_AP' )
      CALL pbuf_get_field( pbuf, idx, ywetdens )
      ywetdens(:,:,:) = 0.0e+0_f8
      ywetdens(1:ncol,:,:) = wetdens(1:ncol,:,:)
      
      idx = pbuf_get_index( 'HYGROM' )
      CALL pbuf_get_field( pbuf, idx, yhygro )
      yhygro(:,:,:) = 0.0e+0_f8
      yhygro(1:ncol,:,:) = hygro(1:ncol,:,:)

      DO n = 1, ntot_amode
      DO ll = 0, nspec_amode(n)
         l = numptrcw_amode(n)
         IF (ll > 0) l = lmassptrcw_amode(ll,n)
         fldcw => qqcw_get_field( pbuf, l, lchnk )
         fldcw(:,:) = 0.0e+0_f8
         fldcw(1:ncol,:) = qqcw(1:ncol,:,l)
      END DO
      END DO


      RETURN
      END SUBROUTINE load_pbuf


!-------------------------------------------------------------------------------
      SUBROUTINE unload_pbuf( pbuf, lchnk, ncol, &
         cld, qqcw, dgncur_a, dgncur_awet, qaerwat, wetdens,hygro )

      USE mam_utils, only: pcols,pver
      USE constituents, only : pcnst
      USE chem_mods, only: adv_mass, gas_pcnst, imozart
      USE physconst, only: mwdry

      USE modal_aero_data, only:  &
         lmassptrcw_amode, nspec_amode, numptrcw_amode, &
         qqcw_get_field, ntot_amode

      USE physics_buffer, only: physics_buffer_desc, &
         pbuf_get_index, pbuf_get_field

      TYPE(physics_buffer_desc), pointer :: pbuf(:)  ! physics buffer for a chunk

      INTEGER,  intent(in   ) :: lchnk, ncol

      REAL(r8), intent(in   ) :: cld(pcols,pver)    ! stratiform cloud fraction

      REAL(r8), intent(inout) :: qqcw(pcols,pver,pcnst)  ! Cloudborne aerosol MR array
      REAL(r8), intent(inout) :: dgncur_a(pcols,pver,ntot_amode)
      REAL(r8), intent(inout) :: dgncur_awet(pcols,pver,ntot_amode)
      REAL(r8), intent(inout) :: qaerwat(pcols,pver,ntot_amode)
      REAL(r8), intent(inout) :: wetdens(pcols,pver,ntot_amode)
      REAL(r8), intent(inout) :: hygro(pcols,pver,ntot_amode)

      INTEGER :: i, idx, k, l, ll, n
      REAL(r8) :: tmpa

      REAL(r8), pointer :: fldcw(:,:)
      REAL(r8), pointer :: ycld(:,:)
      REAL(r8), pointer :: ydgnum(:,:,:)
      REAL(r8), pointer :: ydgnumwet(:,:,:)
      REAL(r8), pointer :: yqaerwat(:,:,:)
      REAL(r8), pointer :: ywetdens(:,:,:)
      REAL(r8), pointer :: yhygro(:,:,:)


      idx = pbuf_get_index( 'CLD' )
      CALL pbuf_get_field( pbuf, idx, ycld )
! cld should not have changed, so check for changes rather than unloading it
!     cld(1:ncol,:) = ycld(1:ncol,:)
      tmpa = maxval( abs( cld(1:ncol,:) - ycld(1:ncol,:) ) )
      IF (tmpa /= 0.0e+0_f8) then
         write(*,*) '*** unload_pbuf cld change error - ', tmpa
         stop
      END IF

      idx = pbuf_get_index( 'DGNUM' )
      CALL pbuf_get_field( pbuf, idx, ydgnum )
      dgncur_a(1:ncol,:,:) = ydgnum(1:ncol,:,:)

      idx = pbuf_get_index( 'DGNUMWET' )
      CALL pbuf_get_field( pbuf, idx, ydgnumwet )
      dgncur_awet(1:ncol,:,:) = ydgnumwet(1:ncol,:,:)

      idx = pbuf_get_index( 'QAERWAT' )
      CALL pbuf_get_field( pbuf, idx, yqaerwat )
      qaerwat(1:ncol,:,:) = yqaerwat(1:ncol,:,:)

      idx = pbuf_get_index( 'WETDENS_AP' )
      CALL pbuf_get_field( pbuf, idx, ywetdens )
      wetdens(1:ncol,:,:) = ywetdens(1:ncol,:,:)

      idx = pbuf_get_index( 'HYGROM' )
      CALL pbuf_get_field( pbuf, idx, yhygro )
      hygro(1:ncol,:,:) = yhygro(1:ncol,:,:)


      DO n = 1, ntot_amode
      DO ll = 0, nspec_amode(n)
         l = numptrcw_amode(n)
         IF (ll > 0) l = lmassptrcw_amode(ll,n)
         fldcw => qqcw_get_field( pbuf, l, lchnk )
         qqcw(1:ncol,:,l) = fldcw(1:ncol,:)
      END DO
      END DO


      RETURN
      END SUBROUTINE unload_pbuf

!---------------------------------------------------------------------

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

          DO L = State_Grid%MaxChemLev-1, 1, -1
                DELZ  = BXHEIGHT(I,J,L)
                DELZ1 = BXHEIGHT(I,J,L+1)
                TC(L) = 1.d0 / &
                      ( 1.d0   + DTCHEM * VTS(L)   / DELZ ) * &
                      ( TC(L)  + DTCHEM * VTS(L+1) / DELZ1  *  TC(L+1) )
          ENDDO
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
    INTEGER                  :: I, J, L, M
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
    !$OMP PRIVATE( I, J, L  )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
! Start with state_chm%GCMAM diag     
    
! modal 
      IF ( State_Diag%Archive_MamNu ) THEN
          DO m = 1, size(State_Chm%GCMAM)
            IF (State_Diag%Map_MamNu%id2slot(m) > 0)  &
            State_Diag%MamNu(I,J,L,m) = State_Chm%GCMAM(m)%Nu(I,J,L) * 1.0e-6_fp !#m-3 to #cm-3
           END DO
      ENDIF
!now everything is set up to have modal species concentration diag as well 
!total concentrations
      IF ( State_Diag%Archive_MamSO4Mass ) THEN
           State_Diag%MamSO4Mass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lso4)   State_Diag%MamSO4Mass(I,J,L) = &
                   State_Diag%MamSO4Mass(I,J,L) + State_Chm%GCMAM(m)%so4(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamBCMass ) THEN
           State_Diag%MamBCMass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lbc)    State_Diag%MamBCMass(I,J,L) = &
                   State_Diag%MamBCMass(I,J,L) + State_Chm%GCMAM(m)%bc(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamPOMMass ) THEN
           State_Diag%MamPOMMass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lpom)   State_Diag%MamPOMMass(I,J,L) = &
                   State_Diag%MamPOMMass(I,J,L) + State_Chm%GCMAM(m)%pom(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamSOAMass ) THEN
           State_Diag%MamSOAMass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lsoa)  State_Diag%MamSOAMass(I,J,L) = &
                   State_Diag%MamSOAMass(I,J,L) + State_Chm%GCMAM(m)%soa(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamSSLTMass ) THEN
           State_Diag%MamSSLTMass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lsslt) State_Diag%MamSSLTMass(I,J,L) = &
                   State_Diag%MamSSLTMass(I,J,L) + State_Chm%GCMAM(m)%sslt(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamDUSTMass ) THEN
           State_Diag%MamDUSTMass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%ldust)  State_Diag%MamDUSTMass(I,J,L) = &
                   State_Diag%MamDUSTMass(I,J,L) + State_Chm%GCMAM(m)%dust(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamNH4Mass ) THEN
           State_Diag%MamNH4Mass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lnh4)  State_Diag%MamNH4Mass(I,J,L) = &
                   State_Diag%MamNH4Mass(I,J,L) + State_Chm%GCMAM(m)%nh4(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamNO3Mass ) THEN
           State_Diag%MamNO3Mass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lno3)  State_Diag%MamNO3Mass(I,J,L) = &
                   State_Diag%MamNO3Mass(I,J,L) + State_Chm%GCMAM(m)%no3(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamCAMass ) THEN
           State_Diag%MamCAMass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lca)  State_Diag%MamCAMass(I,J,L) = &
                   State_Diag%MamCAMass(I,J,L) + State_Chm%GCMAM(m)%ca(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamCO3Mass ) THEN
           State_Diag%MamCO3Mass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lco3)  State_Diag%MamCO3Mass(I,J,L) = &
                   State_Diag%MamCO3Mass(I,J,L) + State_Chm%GCMAM(m)%co3(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamCLMass ) THEN
           State_Diag%MamCLMass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lcl)  State_Diag%MamCLMass(I,J,L) = &
                   State_Diag%MamCLMass(I,J,L) + State_Chm%GCMAM(m)%cl(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF
      IF ( State_Diag%Archive_MamMOMMass ) THEN
           State_Diag%MamMOMMass(I,J,L) = 0.
           DO m = 1, size(State_Chm%GCMAM)
               IF (State_Chm%GCMAM(m)%lmom)  State_Diag%MamMOMMass(I,J,L) = &
                   State_Diag%MamMOMMass(I,J,L) + State_Chm%GCMAM(m)%mom(I,J,L) * kgm3_to_ugm3
           END DO
      ENDIF

  IF ( State_Diag%Archive_MamWATMass ) THEN
       State_Diag%MamWATMass(I,J,L) = 0.
       DO m = 1, size(State_Chm%GCMAM)
             State_Diag%MamWATMass(I,J,L) = &
               State_Diag%MamWATMass(I,J,L) + State_Chm%GCMAM(m)%aerwat(I,J,L) * kgm3_to_ugm3
       END DO
  ENDIF

     IF ( State_Diag%Archive_Mamwetrad ) THEN
       State_Diag%Mamwetrad(I,J,L) = 0.
       DO m = 1, 1 ! FAB just output the accum mode fro now
           State_Diag%Mamwetrad(I,J,L) = &
               State_Diag%Mamwetrad(I,J,L) + State_Chm%GCMAM(m)%wetrad(I,J,L)
       END DO
     ENDIF

      IF ( State_Diag%Archive_Mamdryrad ) THEN
       State_Diag%Mamdryrad(I,J,L) = 0.
       DO m = 1, 1
           State_Diag%Mamdryrad(I,J,L) = &
               State_Diag%Mamdryrad(I,J,L) + State_Chm%GCMAM(m)%dryrad(I,J,L)
       END DO
      ENDIF

      IF ( State_Diag%Archive_Mamhygro ) THEN
       State_Diag%Mamhygro(I,J,L) = 0.
       DO m = 1, 1
           State_Diag%Mamhygro(I,J,L) = &
               State_Diag%Mamhygro(I,J,L) + State_Chm%GCMAM(m)%hygro(I,J,L)
       END DO
       ENDIF
! aerosol optical properties 
     
        ! modal 
      IF ( State_Diag%Archive_MamTauxarv ) THEN
          DO m = 1, size(State_Chm%GCMAM)
            IF (State_Diag%Map_MamTauxarv%id2slot(m) > 0)  &
            State_Diag%MamTauxarv(I,J,L,m) = State_Chm%GCMAM(m)%tauxar(I,J,L,10) !visible band  
           END DO
      ENDIF
      IF ( State_Diag%Archive_Mamssav ) THEN
          DO m = 1, size(State_Chm%GCMAM)
            IF (State_Diag%Map_Mamssav%id2slot(m) > 0)  &
            State_Diag%Mamssav(I,J,L,m) = State_Chm%GCMAM(m)%ssa(I,J,L,10) !visible band  
           END DO
      ENDIF
      IF ( State_Diag%Archive_Mamgv ) THEN
          DO m = 1, size(State_Chm%GCMAM)
            IF (State_Diag%Map_Mamgv%id2slot(m) > 0)  &
            State_Diag%Mamgv(I,J,L,m) = State_Chm%GCMAM(m)%g(I,J,L,10) !visible band  
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
  SUBROUTINE MAM_APPLY_RAINOUT_EFF( hygro, TK, SpcInfo, RainFrac )
!
! !USES:
!
    USE Species_Mod, ONLY : Species
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
!    INTEGER        INTENT(IN)    :: I,J,L      ! loop indices 
    REAL(fp),      INTENT(IN)    :: TK         ! Temperature [K]
    REAL(fp),      INTENT(IN)    :: hygro      ! mode hygroscopicty
    TYPE(Species), INTENT(IN)    :: SpcInfo    ! Species Database object
                                               ! mode index notably
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(fp),      INTENT(INOUT) :: RainFrac   ! Rainout fraction

! LOCAL 
    REAL(fp)                   :: RainoutEff    
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    Rainouteff = 0.0e+0_fp
    
    ! Apply hygro-temperature-mode dependent rainout efficiencies
 
    IF ( TK < 258.0_fp ) THEN ! cold and mixed phase clouds 
       ! Ice: T < 237 K
       ! first approach : value set following arguments developped in Luo et al., 2020
       ! perhaps find more appropriate rational ( use mode modal species composition)  
      IF (SpcInfo%MamModId == 1) then ! Accum 
         IF(hygro > 0.5e+0_fp ) then 
            Rainouteff = 0.4e+0_fp
         ELSE 
            Rainouteff = 0.6e+0_fp 
         END IF
      END IF

      IF (SpcInfo%MammodId == 2)  Rainouteff = 0.08e+0_fp  !Aitken assumed to    

      IF (SpcInfo%MamModId == 3) then !coarse
         IF(hygro > 0.4e+0_fp ) then 
             Rainouteff = 0.4e+0_fp ! scavenge like hydrophilic 
         ELSE
             Rainouteff = 1.0e+0_fp ! likely coarse dust dominated      
         END IF
       END IF    

       IF (SpcInfo%MamModId == 4)  Rainouteff = 0.5e+0_fp !MAM primary carbon : hydrophobic
         
       IF ( TK >= 237.0_fp )  then  ! mixed phase , temp. correction ( Luo et al., 2020)
             Rainouteff = Rainouteff * &
             ( EXP( 0.46e+0_fp * ( 273.16_fp - TK ) - 11.6_fp ) / 153.5_fp )
       END IF 

    ELSE ! warm clouds, Liquid rain: T > 258 K
         ! consider that rainout efficiency equals MAM hygroscopcity, limited to 1. 
           Rainouteff = min(hygro,1.0e+0_fp) 
    ENDIF

     RainFrac = RainFrac * Rainouteff

  END SUBROUTINE MAM_APPLY_RAINOUT_EFF



END MODULE MAM_DRIV_MOD

!FAB#endif

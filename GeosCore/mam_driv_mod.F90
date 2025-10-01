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

use precision_mod, only :r8 => f8, fp, f8  
use physics_buffer, only: physics_buffer_desc
use physics_types, only : physics_state, physics_ptend
use mam_utils, only : masterproc, pcols, pver, l_h2so4g, l_soag
use constituents, only : pcnst

!
IMPLICIT NONE

PRIVATE

!PUBLIC MEMBER FUNCTIONS:

PUBLIC :: MAM_DRIV, MAM_INIT 

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
!

type(physics_buffer_desc), pointer :: pbuf(:)
type(physics_state) :: physta
type(physics_ptend) :: ptend

! define a specific type to handle MAM/GC prognostic species information.
! it is a bit similar to State_Chm%SpcData(N)%Info
! Convenient for communication between MAM and GC worlds for instance 
! Potentially it could live as a specific subtype of chem_state ? just an idea..
! Perhaps this type should also be declared in Headers for consistency with GC ? 

type, public ::  mamspec ! 
  integer :: gcind  ! sp index relative to Spc ( GC chemstate species )
  integer :: mamind ! sp index relative to both q and qqcw MAM states 
  integer :: modID  ! MAM mode index to which this sp belongs (also used to point to chemstate%GCMAM(mode)%xx)
  logical :: isnum  ! True if number concentration vs mass concentration
  logical :: iscb
  character* 12  :: name !GC species name for MAM tracers ( ABSOLUTLY must be consistent with speciesdat.yml)  
  character* 12  :: namecb !GC species name for MAM cloudborne sp ( ABSOLUTLY must be consistent with speciesdat.yml)
end type mamspec

type(mamspec), pointer :: mamgc(:)

integer nmamgc ! number of GC advected MAM tracers 

!MAM control ( will go in namelist) 
    integer :: mdo_gasaerexch,     mdo_rename,          &
               mdo_newnuc,         mdo_coag
    integer :: mdo_gaschem, mdo_cloudchem

    integer:: loffset, lchnk

    real(r8) :: deltat

    logical, save  :: lfirstcall


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

    use mam_utils, only: begchunk, endrun
    use chem_mods, only: adv_mass, gas_pcnst, imozart
    use physconst, only: mwdry
    use modal_aero_data, only: numptr_amode, lptr_so4_a_amode, &
                               lptr_bc_a_amode, lptr_nacl_a_amode,&
                               lptr_pom_a_amode, lptr_soa_a_amode,&
                               lptr_dust_a_amode, lptr_so4_cw_amode,&
                               modeptr_accum, alnsg_amode
    use modal_aero_initialize_data, only: MAM_cold_start 
    use modal_aero_calcsize, only: modal_aero_calcsize_sub
    use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
    use modal_aero_amicphys, only: modal_aero_amicphys_intr

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
! !DEFINED PARAMETERS:
!
    real(kind=8), parameter :: cpdair = 1.004e3  ! Specific heat capacity
                                                 ! of dry air at constant
                                                 ! pressure at 273 K
                                                 ! (J kg-1 K-1)
!
! !LOCAL VARIABLES:
!
    INTEGER :: I,J,L,S,l2,K,N,M,SIZENUM,MDAY

    ! Make a pointer to the tracer array
    TYPE(SpcConc), POINTER :: Spc(:)

! !LOCAL VARIABLES:
!
      real(r8) :: vmr(pcols,pver,gas_pcnst)     ! gas & aerosol volume mixing ratios
      real(r8) :: vmrcw(pcols,pver,gas_pcnst)   ! gas & aerosol cloud-borne volume mixing ratios
      real(r8) :: vmr_svaa(pcols,pver,gas_pcnst) ! temp save before gas chem
      real(r8) :: vmr_svbb(pcols,pver,gas_pcnst) ! temp save before cloud chem
      real(r8) :: vmrcw_svbb(pcols,pver,gas_pcnst)!temp save before cloud chem 
      real(r8) :: aircon(pcols,pver) !  air concentration (kmol/m3)
      real(r8) :: tau_gaschem_simple(pcols,pver)
      integer :: latndx(pcols),lonndx(pcols)                 !required by the mam interface
                                                !not used now potentiall usefull for diags 
      integer :: nstep                          !same
      character(len=8) :: spcnam
!--------------------------------------------------------------------------

       
    ! Point to Spc
    Spc => State_Chm%Species
     if (masterproc) then
     print*, 'FAB comp rates',maxval(h2so4_rate), maxval(Spc(Ind_('PH2SO4'))%Conc(:,:,:))/deltat ,&
               maxval(PSO4AQ_RATE)/deltat, maxval(Spc(Ind_('PSO4AQ'))%Conc(:,:,:))/deltat  
     end if  

    lchnk = begchunk
    loffset = imozart -1
    ! load the mam met state   
    physta%lchnk = lchnk
    physta%ncol  = pcols

    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX ! GEOS-Chem vertical grid is bottom-up
      n = J + (I-1)*State_Grid%NY 
      physta%t(n,L) = State_Met%T(I,J,L)
      physta%pdel(n,l) = State_Met%DELP(I,J,L) * 100._r8 ! hPa to Pa 
      physta%pmid(n,l) = State_Met%PMID(I,J,L) * 100._r8
      physta%cld(n,l)= State_Met%CLDF(I,J,L)
      physta%relhum(n,l)= State_Met%RH(I,J,L)
      physta%qv(n,l) = State_Met%SPHU(I,J,L) * 1.E-3_r8 ! in kg/kgair  Caution here make sure
      physta%zm(n,l) =  sum(State_Met%BXHEIGHT(I,J,1:l))- 0.5 * State_Met%BXHEIGHT(I,J,l)

      if (L==1) physta%pblh(n) = State_Met%PBLH(I,J)

      ! first element water vapor mr (used in water )  
      physta%q(n,l,1) = physta%qv(n,l) / (1._r8 - physta%qv(n,l))

      ! load sulf production rate and convert from Kg.s-1  to kg.kg-1.s-1 
      ! needs fullchem activated 
       
      physta%ph2so4(n,l) = h2so4_rate(I,J,L) / State_Met%AD(I,J,L)              
      !paqso4 is the prod per time step ! consider using PSO4AQ, PH2SO4 instead!
      ! and remove the interface in fullchem ?? 
      physta%paqso4(n,l) = PSO4AQ_RATE(I,J,L)/deltat / State_Met%AD(I,J,L)

      ! load q gas ...
      ! the gas phase species SO2,DMS,H2O2 in q are not used/modified 
      ! if we use GC production rate for SO4 instead of MAM simple chem 
      ! ,maybe get rid of them later or consider haveing some oxidaton routines
      ! which could run dindependant of fullchem e.g. from prescribed oxidants..

    
       physta%q(n,l,l_h2so4g) = Spc(Ind_('H2SO4'))%Conc(I,J,L) / State_Met%AD(I,J,L)

       
      physta%q(n,l,l_soag) = Spc(Ind_('SOAP'))%Conc(I,J,L) / State_Met%AD(I,J,L) 
      ! Rq in GC standard lumped SOAP is not treated as semi_volatil but in MAM yes
      ! is this reqsonqble ? 
      ! develop options with advanced SOA scheme  
      ! other gas might be considered if MAM7 and/or MOSAIC are implemented

       ! FAB TEST 
       !Spc(Ind_('MAMDEV'))%Conc(I,J,L)= Spc(Ind_('MAMDEV'))%Conc(I,J,L) + h2so4_rate(i,j,l) *deltat + PSO4AQ_RATE(i,j,l)
 
!         Spc(Ind_('MAMDEV'))%Conc(I,J,L)= Spc(Ind_('MAMDEV'))%Conc(I,J,L) + Spc(Ind_('PH2SO4'))%Conc(I,J,L) + Spc(Ind_('PSO4AQ'))%Conc(I,J,L)  


     END DO
     END DO
     END DO

    if (.not. lfirstcall) then 
     ! load q and qqcw mam state aerosol variables 
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX ! GEOS-Chem vertical grid is bottom-up
      n = J + (I-1)*State_Grid%NY         
     do s = 1 , nmamgc
        !number concentrations #/gridbox and convert to #/kg 
        !mass concentrations Kg/gridbox to Kg/Kg (mixing ratios) 
        if (.not. mamgc(s)%iscb ) then
          physta%q(n,l,mamgc(s)%mamind) = Spc(mamgc(s)%gcind)%Conc(I,J,L)/State_Met%AD(I,J,L) !
        else
          physta%qqcw(n,l,mamgc(s)%mamind) = Spc(mamgc(s)%gcind)%Conc(I,J,L)/State_Met%AD(I,J,L) !
        end if
     end do  
    END DO
    END DO
    END DO
    end if ! lfirst call.

    ! This colde start init is temporary until handling a proper GC restart file
    if (lfirstcall) then 
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX ! 
      n = J + (I-1)*State_Grid%NY
      physta%aircon(n,l) =  State_Met%AIRDEN(I,J,L) * 0.0345_r8 !Kmol/m3 
    END DO
    END DO
    END DO
    ! initialise q aerosol state (cold start) / for testing phase  
    call MAM_cold_start (physta) 
    endif         



! CALCSIZE INTERFACE     
     
call load_pbuf( pbuf, lchnk, pcols, &
        physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens )

! call calcsize
    ptend%lq = .false.
    ptend%q = 0._r8
    call modal_aero_calcsize_sub( physta, ptend, deltat, pbuf, &
         do_adjust_in=.true., do_aitacc_transfer_in=.true. )

! unload pbuf
      call unload_pbuf( pbuf, lchnk, pcols, &
         physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens )
!
! apply tendencies ! note ptend.lq is modified by calcsize
! note also that the cloudborne state is supposed to be directly updated in calcsize
! (perhaps because unlinke q , qqcw is not an advected state in cesm and the tendencie does not need to be 
!  passed up ...
!  
      do l = 1, pcnst
         if ( .not. ptend%lq(l) ) cycle
         do k = 1, pver
         do i = 1, pcols 
            physta%q(i,k,l) = physta%q(i,k,l) + ptend%q(i,k,l)*deltat  
            physta%q(i,k,l) = max( physta%q(i,k,l), 0.0_r8 )
         end do
         end do
      end do
      
      physta%lchnk = lchnk      
      physta%ncol = pcols      

! WATER UPTAKE    
     call load_pbuf( pbuf, lchnk, pcols, &
        physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens )
!
     call modal_aero_wateruptake_dr( physta, pbuf )
     
     call unload_pbuf( pbuf, lchnk, pcols, &
         physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens )


 
 
IF( .false.) THEN
  print*, lptr_soa_a_amode(2)
  print*, 'q sox avant ',l_soag, physta%q(5,2,l_soag),physta%q(5,2,lptr_soa_a_amode(1)), physta%q(5,2,lptr_soa_a_amode(2)) 
  print*, 'q num avant', physta%q(5,2,17), physta%q(5,2,21), physta%q(5,2,25)
END IF
 
 
!-------------------------------------------------------------------------------- 
! switch from q & qqcw mass mixing ratios to volume mixing ratios  vmr and vmrcw
! only adress the gas/aerosol variables in q, qqcw 
! Rq : not all gases are modified by MAM routines so q size could be reduced  
! Rq2: the flow of unit changes could be perhaps simpler
      vmr = 0.0_r8
      vmrcw = 0.0_r8
      do l = imozart, pcnst
         l2 = l - loffset
         vmr(  1:pcols,1:pver,l2) =physta%q(  1:pcols,1:pver,l)*mwdry/adv_mass(l2)
         vmrcw(1:pcols,1:pver,l2) =physta%qqcw(1:pcols,1:pver,l)*mwdry/adv_mass(l2)
      end do
!-------------------------------
! GASCHEM interface 

      vmr_svaa   = vmr  !save before gas chem , this is how the mam code proceed

      if (mdo_gaschem > 0) then
        !
        l2 = l_h2so4g-loffset
        vmr(1:pcols,1:pver,l2) = vmr(1:pcols,1:pver,l2) + physta%ph2so4(1:pcols,1:pver)*mwdry/adv_mass(l2)*deltat          

      end if 

      vmr_svbb = vmr    ! save before cloud chem
      vmrcw_svbb = vmrcw! save before cloud chem

!CLOUDCHEM
!rq vmrcw / qcw are not advected in MAM /CESM   

      if (mdo_cloudchem > 0) then
        !start by updating mass in the accumulation mode 
        !(consider partitioning with aitken )
        l2 = lptr_so4_cw_amode(modeptr_accum) - loffset
        vmrcw(1:pcols,1:pver,l2) = vmrcw(1:pcols,1:pver,l2) +  & 
                                   physta%paqso4(1:pcols,1:pver)*mwdry/adv_mass(l2)*deltat

      end if

!------------------------------
  nstep = 1 ! get the GC integration step counter ( used only for print out in fact) 
  if(1==1) then ! .and. masterproc ) then
     call modal_aero_amicphys_intr(               &
         mdo_gasaerexch,     mdo_rename,          &
         mdo_newnuc,         mdo_coag,            &
         lchnk,    pcols,     nstep   ,           &
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
    end if
 
! vmr and vmrcw have been updated in modal_aero_amicphys_intr  
! switch back from vmr & vmrcw to q & qqcw state
!
      do l = imozart, pcnst
         l2 = l - loffset
         physta%q(    1:pcols,1:pver,l)  = vmr(  1:pcols,1:pver,l2) * adv_mass(l2)/mwdry 
         physta%qqcw( 1:pcols,1:pver,l)  = vmrcw(1:pcols,1:pver,l2) * adv_mass(l2)/mwdry
      end do

      
      
IF(.false.) THEN
      print*,'FAB DRIVER'
      print*, 'q apres sox ',l_soag, physta%q(5,2,l_soag),physta%q(5,2,lptr_soa_a_amode(1)), physta%q(5,2,lptr_soa_a_amode(2)) 
      print*, 'q apres num ', physta%q(5,2,17), physta%q(5,2,21), physta%q(5,2,25)
END IF
      
! Update GC/MAM species 
     DO L = 1, State_Grid%NZ
     DO J = 1, State_Grid%NY
     DO I = 1, State_Grid%NX
        n = J + (I-1)*State_Grid%NY
        !both for interstitial and cloudborne states
        do s =1, nmamgc
        if(.not. mamgc(s)%iscb) then   
            Spc(mamgc(s)%gcind)%Conc(I,J,L) = physta%q(n,l,mamgc(s)%mamind) &
                * State_Met%AD(I,J,L)
        else
            Spc(mamgc(s)%gcind)%Conc(I,J,L) = physta%qqcw(n,l,mamgc(s)%mamind) &
                                          * State_Met%AD(I,J,L)
        end if
        end do
        !gas species affected by mam 
        Spc(Ind_('SOAP'))%Conc(I,J,L) = physta%q(n,l,l_soag) &
                                       * State_Met%AD(I,J,L)  

        Spc(Ind_('H2SO4'))%Conc(I,J,L) = physta%q(n,l,l_h2so4g) &
                                       * State_Met%AD(I,J,L)

      ENDDO
      ENDDO
      ENDDO  
                               
! fill state GCMAM chem state variables,  used in e.g. drydep  nd diags 
! harmonize mamgc and GCMAM 
! claude AI suggest to keep m loop inside l,j,i loop as long as m is small (which is the case)
! try to optimize the if statements within loops

      DO L = 1, State_Grid%NZ
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX
        n = J + (I-1)*State_Grid%NY


      do m= 1, size(State_Chm%GCMAM) ! loop on modes     

        State_Chm%GCMAM(m)%nudryrad(I,J,L) = 0.5_r8*physta%dgncur_a(n,L,m)   

        State_Chm%GCMAM(m)%nuwetrad(I,J,L) = 0.5_r8*physta%dgncur_awet(n,L,m)   

        ! volume mean geo radius
        State_Chm%GCMAM(m)%dryrad(I,J,L) = 0.5_r8*physta%dgncur_a(n,L,m)     &
                                            *exp(3._r8*(alnsg_amode(m)**2))
        State_Chm%GCMAM(m)%wetrad(I,J,L) = 0.5_r8*physta%dgncur_awet(n,L,m)  &
                                            *exp(3._r8*(alnsg_amode(m)**2))
        ! wet aer density
        State_Chm%GCMAM(m)%aerdens(I,J,L) =  physta%wetdens(n,L,m)         

        !modal number concentrations in #.m-3
        State_Chm%GCMAM(m)%Nu(I,J,L) =              &
                                physta%q(n,L,numptr_amode(m))*State_Met%AIRDEN(I,J,L)
        ! modal mass concentrations in Kg.m-3  
        if(lptr_so4_a_amode(m) > 0 ) State_Chm%GCMAM(m)%so4(I,J,L) =              & 
                                physta%q(n,L,lptr_so4_a_amode(m))*State_Met%AIRDEN(I,J,L) + &  
        !FAB TEMP add the cloud borne sulf to chm state for diag // change that once 
        ! transfer from qqcw to q is properly trated !!
                                physta%qqcw(n,L,lptr_so4_a_amode(m))*State_Met%AIRDEN(I,J,L)    

         
        if(lptr_bc_a_amode(m) > 0 ) State_Chm%GCMAM(m)%bc(I,J,L) =              &
                                physta%q(n,L,lptr_bc_a_amode(m))*State_Met%AIRDEN(I,J,L)

        if(lptr_pom_a_amode(m) > 0 ) State_Chm%GCMAM(m)%pom(I,J,L) =              &
                                physta%q(n,L,lptr_pom_a_amode(m))*State_Met%AIRDEN(I,J,L)
 
        if(lptr_soa_a_amode(m) > 0 ) State_Chm%GCMAM(m)%soa(I,J,L) =              &
                                physta%q(n,L,lptr_soa_a_amode(m))*State_Met%AIRDEN(I,J,L)
  
        if(lptr_nacl_a_amode(m) > 0 ) State_Chm%GCMAM(m)%sslt(I,J,L) =              &
                                physta%q(n,L,lptr_nacl_a_amode(m))*State_Met%AIRDEN(I,J,L)

        if(lptr_dust_a_amode(m) > 0 ) State_Chm%GCMAM(m)%dust(I,J,L) =              &
                                physta%q(n,L,lptr_dust_a_amode(m))*State_Met%AIRDEN(I,J,L)
 

      end do 

     ENDDO
     ENDDO
     ENDDO
          
     ! call MAM aerosol update for gravitational settling (this call could be somewhere else )        
     call  MAM_SETTL( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )

!Fill out MAM diags (cf Headers/state_diag)    

    call Set_MAM_Diagnostic( Input_Opt, State_Chm, State_Diag, &
                                     State_Grid, State_Met, RC )

!

    Spc => NULL() 
    if (lfirstcall) lfirstcall = .false.

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

    use mam_utils, ONLY : plev 
    use physics_buffer, only: physics_buffer_desc 
    use physics_types, only : physics_state
    use modal_aero_data, only: numptr_amode, lptr_so4_a_amode, &
                               lptr_bc_a_amode, lptr_nacl_a_amode,&
                               lptr_pom_a_amode, lptr_soa_a_amode,&
                               lptr_dust_a_amode
    use modal_aero_initialize_data, only: MAM_init_basics, MAM_ALLOCATE

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

    integer :: species_class(pcnst) = -1
    integer  :: m , i , s 
    character * 12 :: tmp

    
    
!-------initialise namelist parameters

   mdo_gaschem=1
   mdo_cloudchem =1

   mdo_gasaerexch=1
   mdo_rename=1
   mdo_newnuc=1
   mdo_coag=1

   masterproc = Input_Opt%amIRoot
   lfirstcall = .true.
! Chemistry timestep [s]
   deltat  = GET_TS_CHEM()


pcols = State_Grid%NX *  State_Grid%Ny
pver  = State_Grid%NZ
plev = pver

!pcols pver pourraient just passer par module plutot que par argument
call MAM_init_basics(pbuf)

!allocate MAM state object
call MAM_ALLOCATE (physta,ptend )

!allocate specific GC diqg usefull for mam  
ALLOCATE( PSO4AQ_RATE(State_Grid%NX,State_Grid%NY,State_Grid%NZ) )
ALLOCATE( H2SO4_RATE(State_Grid%NX,State_Grid%NY,State_Grid%NZ) )

!Initialize MAM4/GC transported tracer info  
!Rq allocate according to MAM option , start with MAM 4 , to be refined 
nmamgc = State_Chm%nMam
allocate(mamgc(State_Chm%nMam)) 

do i = 1, State_Chm%nMam
 m = State_Chm%map_Mam(i) 
 mamgc(i)%name = State_Chm%SpcData(m )%Info%name 
 mamgc(i)%gcind = m 
 mamgc(i)%modId =State_Chm%SpcData(m )%Info%MamModId
 mamgc(i)%isnum = State_Chm%SpcData(m )%Info%MP_SizeResNum
 mamgc(i)%iscb = State_Chm%SpcData(m )%Info%Is_CloudBorne
! add maping info for mam q and qqcw states 
! test name and  for cloud borne - important species in q and qqcw have the same indexing
 s=4 ! work only for MAMXXX naming convention
 if (mamgc(i)%iscb) s=6
 mamgc(i)%mamind =-1
 !mam modal indices should match existing species and be > 0 
 !should also be consistent with the allocation state of 
 !the chem state%GCMAM(mode)%XXX(:,:,:) in /headers
 !t.b.d make a security test with explicit error message
 if (mamgc(i)%name(s:s+1) == 'Nu') mamgc(i)%mamind = numptr_amode(mamgc(i)%modId)
 if (mamgc(i)%name(s:s+2) == 'SO4') mamgc(i)%mamind = lptr_so4_a_amode(mamgc(i)%modId)
 if (mamgc(i)%name(s:s+1) == 'BC') mamgc(i)%mamind =  lptr_bc_a_amode(mamgc(i)%modId)
 if (mamgc(i)%name(s:s+2) == 'POM') mamgc(i)%mamind = lptr_pom_a_amode(mamgc(i)%modId)
 if (mamgc(i)%name(s:s+2) == 'SOA') mamgc(i)%mamind = lptr_soa_a_amode(mamgc(i)%modId)
 if (mamgc(i)%name(s:s+3) == 'SSLT') mamgc(i)%mamind = lptr_nacl_a_amode(mamgc(i)%modId)
 if (mamgc(i)%name(s:s+3) == 'DUST') mamgc(i)%mamind = lptr_dust_a_amode(mamgc(i)%modId)
! to be updated when adding species to MAM
end do

  print*, 'MAM species INFO'
  print*, mamgc(:)%name
  print*, mamgc(:)%gcind
  print*, mamgc(:)%mamind
  print*, mamgc(:)%modId
  print*, mamgc(:)%isnum
   print*, mamgc(:)%iscb
END SUBROUTINE MAM_INIT        


!-------------------------------------------------------------------------------
subroutine load_pbuf( pbuf, lchnk, ncol,  &
         cld, qqcw, dgncur_a, dgncur_awet, qaerwat, wetdens )


      use mam_utils, only: pcols,pver
      use constituents, only : pcnst
      use chem_mods, only: adv_mass, gas_pcnst, imozart
      use physconst, only: mwdry

      use modal_aero_data, only:  &
         lmassptrcw_amode, nspec_amode, numptrcw_amode, &
         qqcw_get_field, ntot_amode

      use physics_buffer, only: physics_buffer_desc, &
         pbuf_get_index, pbuf_get_field

      type(physics_buffer_desc), pointer :: pbuf(:)  ! physics buffer for a chunk

      integer,  intent(in   ) :: lchnk, ncol

      real(r8), intent(in   ) :: cld(pcols,pver)    ! stratiform cloud fraction
      real(r8), intent(in   ) :: qqcw(pcols,pver,pcnst)  ! Cloudborne aerosol MR array
      real(r8), intent(in   ) :: dgncur_a(pcols,pver,ntot_amode)
      real(r8), intent(in   ) :: dgncur_awet(pcols,pver,ntot_amode)
      real(r8), intent(in   ) :: qaerwat(pcols,pver,ntot_amode)
      real(r8), intent(in   ) :: wetdens(pcols,pver,ntot_amode)

      integer :: idx, l, ll, n

      real(r8), pointer :: fldcw(:,:)
      real(r8), pointer :: ycld(:,:)
      real(r8), pointer :: ydgnum(:,:,:)
      real(r8), pointer :: ydgnumwet(:,:,:)
      real(r8), pointer :: yqaerwat(:,:,:)
      real(r8), pointer :: ywetdens(:,:,:)

 ! FAB ncol = pcols , maybe getrif of it   
      idx = pbuf_get_index( 'CLD' )
      call pbuf_get_field( pbuf, idx, ycld )
      ycld(:,:) = 0.0_r8
      ycld(1:ncol,:) = cld(1:ncol,:)
      
      idx = pbuf_get_index( 'DGNUM' )
      call pbuf_get_field( pbuf, idx, ydgnum )
      ydgnum(:,:,:) = 0.0_r8
      ydgnum(1:ncol,:,:) = dgncur_a(1:ncol,:,:)
      
      idx = pbuf_get_index( 'DGNUMWET' )
      call pbuf_get_field( pbuf, idx, ydgnumwet )
      ydgnumwet(:,:,:) = 0.0_r8
      ydgnumwet(1:ncol,:,:) = dgncur_awet(1:ncol,:,:)
      
      idx = pbuf_get_index( 'QAERWAT' )
      call pbuf_get_field( pbuf, idx, yqaerwat )
      yqaerwat(:,:,:) = 0.0_r8
      yqaerwat(1:ncol,:,:) = qaerwat(1:ncol,:,:)
      
      idx = pbuf_get_index( 'WETDENS_AP' )
      call pbuf_get_field( pbuf, idx, ywetdens )
      ywetdens(:,:,:) = 0.0_r8
      ywetdens(1:ncol,:,:) = wetdens(1:ncol,:,:)
      
      do n = 1, ntot_amode
      do ll = 0, nspec_amode(n)
         l = numptrcw_amode(n)
         if (ll > 0) l = lmassptrcw_amode(ll,n)
         fldcw => qqcw_get_field( pbuf, l, lchnk )
         fldcw(:,:) = 0.0_r8
         fldcw(1:ncol,:) = qqcw(1:ncol,:,l)
      end do
      end do


      return
      end subroutine load_pbuf


!-------------------------------------------------------------------------------
      subroutine unload_pbuf( pbuf, lchnk, ncol, &
         cld, qqcw, dgncur_a, dgncur_awet, qaerwat, wetdens )

      use mam_utils, only: pcols,pver
      use constituents, only : pcnst
      use chem_mods, only: adv_mass, gas_pcnst, imozart
      use physconst, only: mwdry

      use modal_aero_data, only:  &
         lmassptrcw_amode, nspec_amode, numptrcw_amode, &
         qqcw_get_field, ntot_amode

      use physics_buffer, only: physics_buffer_desc, &
         pbuf_get_index, pbuf_get_field

      type(physics_buffer_desc), pointer :: pbuf(:)  ! physics buffer for a chunk

      integer,  intent(in   ) :: lchnk, ncol

      real(r8), intent(in   ) :: cld(pcols,pver)    ! stratiform cloud fraction

      real(r8), intent(inout) :: qqcw(pcols,pver,pcnst)  ! Cloudborne aerosol MR array
      real(r8), intent(inout) :: dgncur_a(pcols,pver,ntot_amode)
      real(r8), intent(inout) :: dgncur_awet(pcols,pver,ntot_amode)
      real(r8), intent(inout) :: qaerwat(pcols,pver,ntot_amode)
      real(r8), intent(inout) :: wetdens(pcols,pver,ntot_amode)

      integer :: i, idx, k, l, ll, n
      real(r8) :: tmpa

      real(r8), pointer :: fldcw(:,:)
      real(r8), pointer :: ycld(:,:)
      real(r8), pointer :: ydgnum(:,:,:)
      real(r8), pointer :: ydgnumwet(:,:,:)
      real(r8), pointer :: yqaerwat(:,:,:)
      real(r8), pointer :: ywetdens(:,:,:)


      idx = pbuf_get_index( 'CLD' )
      call pbuf_get_field( pbuf, idx, ycld )
! cld should not have changed, so check for changes rather than unloading it
!     cld(1:ncol,:) = ycld(1:ncol,:)
      tmpa = maxval( abs( cld(1:ncol,:) - ycld(1:ncol,:) ) )
      if (tmpa /= 0.0_r8) then
         write(*,*) '*** unload_pbuf cld change error - ', tmpa
         stop
      end if

      idx = pbuf_get_index( 'DGNUM' )
      call pbuf_get_field( pbuf, idx, ydgnum )
      dgncur_a(1:ncol,:,:) = ydgnum(1:ncol,:,:)

      idx = pbuf_get_index( 'DGNUMWET' )
      call pbuf_get_field( pbuf, idx, ydgnumwet )
      dgncur_awet(1:ncol,:,:) = ydgnumwet(1:ncol,:,:)

      idx = pbuf_get_index( 'QAERWAT' )
      call pbuf_get_field( pbuf, idx, yqaerwat )
      qaerwat(1:ncol,:,:) = yqaerwat(1:ncol,:,:)

      idx = pbuf_get_index( 'WETDENS_AP' )
      call pbuf_get_field( pbuf, idx, ywetdens )
      wetdens(1:ncol,:,:) = ywetdens(1:ncol,:,:)

      do n = 1, ntot_amode
      do ll = 0, nspec_amode(n)
         l = numptrcw_amode(n)
         if (ll > 0) l = lmassptrcw_amode(ll,n)
         fldcw => qqcw_get_field( pbuf, l, lchnk )
         qqcw(1:ncol,:,l) = fldcw(1:ncol,:)
      end do
      end do


      return
      end subroutine unload_pbuf

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
      if (masterproc) print*, 'END MAM SETTLING', maxval(Spc(mamgc(1)%gcind)%Conc), maxval(Spc(mamgc(3)%gcind)%Conc)

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
          L     = State_Grid%MaxChemLev
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
  
       if (masterproc) print*, 'END MAM SETTLING', maxval(Spc(mamgc(1)%gcind)%Conc), maxval(Spc(mamgc(3)%gcind)%Conc) 
    ENDIF  ! DOSETTLING

END SUBROUTINE MAM_SETTL


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


! 

    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED   ) &
    !$OMP PRIVATE( I, J, L  )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX


! Start with state_chm%GCMAM diag     
    
! modal 
      IF ( State_Diag%Archive_MamNu ) THEN
          do m = 1, size(State_Chm%GCMAM)
            if (State_Diag%Map_MamNu%id2slot(m) > 0)  &
            State_Diag%MamNu(I,J,L,m) = State_Chm%GCMAM(m)%Nu(I,J,L) * 1.E-6_fp !#m-3 to #cm-3
           end do
      ENDIF

!now everything is set up to have modal species concentration diag as well 

!total concentrations

      IF ( State_Diag%Archive_MamSO4Mass ) THEN
           State_Diag%MamSO4Mass(I,J,L) = 0.
           do m = 1, size(State_Chm%GCMAM)
               if (State_Chm%GCMAM(m)%lso4)   State_Diag%MamSO4Mass(I,J,L) = &
                   State_Diag%MamSO4Mass(I,J,L) + State_Chm%GCMAM(m)%so4(I,J,L) * kgm3_to_ugm3
           end do
      ENDIF

      IF ( State_Diag%Archive_MamBCMass ) THEN
           State_Diag%MamBCMass(I,J,L) = 0.
           do m = 1, size(State_Chm%GCMAM)
               if (State_Chm%GCMAM(m)%lbc)    State_Diag%MamBCMass(I,J,L) = &
                   State_Diag%MamBCMass(I,J,L) + State_Chm%GCMAM(m)%bc(I,J,L) * kgm3_to_ugm3
           end do
      ENDIF

      IF ( State_Diag%Archive_MamPOMMass ) THEN
           State_Diag%MamPOMMass(I,J,L) = 0.
           do m = 1, size(State_Chm%GCMAM)
               if (State_Chm%GCMAM(m)%lpom)   State_Diag%MamPOMMass(I,J,L) = &
                   State_Diag%MamPOMMass(I,J,L) + State_Chm%GCMAM(m)%pom(I,J,L) * kgm3_to_ugm3
           end do
      ENDIF

      IF ( State_Diag%Archive_MamSOAMass ) THEN
           State_Diag%MamSOAMass(I,J,L) = 0.
           do m = 1, size(State_Chm%GCMAM)
               if (State_Chm%GCMAM(m)%lsoa)  State_Diag%MamSOAMass(I,J,L) = &
                   State_Diag%MamSOAMass(I,J,L) + State_Chm%GCMAM(m)%soa(I,J,L) * kgm3_to_ugm3
           end do
      ENDIF

      IF ( State_Diag%Archive_MamSSLTMass ) THEN
           State_Diag%MamSSLTMass(I,J,L) = 0.
           do m = 1, size(State_Chm%GCMAM)
               if (State_Chm%GCMAM(m)%lsslt) State_Diag%MamSSLTMass(I,J,L) = &
                   State_Diag%MamSSLTMass(I,J,L) + State_Chm%GCMAM(m)%sslt(I,J,L) * kgm3_to_ugm3
           end do
      ENDIF

      IF ( State_Diag%Archive_MamDUSTMass ) THEN
           State_Diag%MamDUSTMass(I,J,L) = 0.
           do m = 1, size(State_Chm%GCMAM)
               if (State_Chm%GCMAM(m)%ldust)  State_Diag%MamDUSTMass(I,J,L) = &
                   State_Diag%MamDUSTMass(I,J,L) + State_Chm%GCMAM(m)%dust(I,J,L) * kgm3_to_ugm3
           end do
      ENDIF
   ENDDO
   ENDDO
   ENDDO

END SUBROUTINE Set_MAM_Diagnostic

END MODULE MAM_DRIV_MOD

!FAB#endif

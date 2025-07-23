! FAB think about ifdef
!#ifdef APM
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
use mam_utils, only : masterproc, pcols, pver
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
!
! 
!

type(physics_buffer_desc), pointer :: pbuf2d(:,:)
type(physics_buffer_desc), pointer :: pbuf(:)
type(physics_state) :: physta
type(physics_ptend) :: ptend


!MAM control ( will go in namelist) 
    integer :: mdo_gasaerexch,     mdo_rename,          &
               mdo_newnuc,         mdo_coag
    integer :: mdo_gaschem
    integer :: lchnk, loffset

!species indices pointing to constituents in the physta%q variable 
    integer :: l_nh3g, l_so2g, l_soag, l_hno3g, l_hclg, l_h2so4g 
    integer :: l_num_a1, l_num_a2, l_nh4_a1, l_nh4_a2, &
                 l_so4_a1, l_so4_a2, l_soa_a1, l_soa_a2
    integer :: l_numa, l_so4a, l_nh4a, l_soaa, l_poma, l_bcxa, l_ncla, &
                 l_dsta, l_no3a, l_clxa, l_caxa, l_co3a, l_moma


    real(r8) :: deltat


CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  
    use mam_utils, only: begchunk, endrun
    use chem_mods, only: adv_mass, gas_pcnst, imozart
    use physconst, only: mwdry
    use modal_aero_calcsize, only: modal_aero_calcsize_sub
    use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
    use modal_aero_amicphys, only: modal_aero_amicphys_intr, &
          gaexch_h2so4_uptake_optaa, newnuc_h2so4_conc_optaa, mosaic

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
    INTEGER :: I,J,L,l2,K,N,M,SIZENUM,MDAY

    ! Make a pointer to the tracer array
    TYPE(SpcConc), POINTER :: Spc(:)

! !LOCAL VARIABLES:
!
    REAL(fp)                   :: ANO3, GNO3
    REAL(f8)                   :: RHI, TEMPI

      real(r8) :: vmr(pcols,pver,gas_pcnst)     ! gas & aerosol volume mixing ratios
      real(r8) :: vmrcw(pcols,pver,gas_pcnst)   ! gas & aerosol cloud-borne volume mixing ratios
      real(r8) :: vmr_svaa(pcols,pver,gas_pcnst) ! temp save before gas chem
      real(r8) :: vmr_svbb(pcols,pver,gas_pcnst) ! temp save before cloud chem
      real(r8) :: vmrcw_svbb(pcols,pver,gas_pcnst)!temp save before cloud chem 

      real(r8) :: tau_gaschem_simple(pcols,pver)
      integer :: latndx(pcols),lonndx(pcols)                 !required by the mam interface
                                                !not used now potentiall usefull for diags 
      integer :: nstep                          !same

!--------------------------------------------------------------------------

       
    ! Point to Spc
    Spc => State_Chm%Species

    lchnk = begchunk
    ! load the mam met state   
    physta%lchnk = lchnk
    physta%ncol  = pcols

    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
      n = J + (I-1)*State_Grid%NY 
      physta%t(n,L) = State_Met%T(I,J,L)
      physta%pdel(n,l) = State_Met%DELP(I,J,L) * 100._r8 ! hPa to Pa 
      physta%pmid(n,l) = State_Met%PMID(I,J,L) * 100._r8
      physta%cld(n,l)= State_Met%CLDF(I,J,L)
      physta%relhum(n,l)= State_Met%RH(I,J,L)
      physta%qv(n,l) = State_Met%SPHU(I,J,L) * 1.E-3_r8 ! in kg/kgair  Caution here make sure
      physta%zm(n,l) =  sum(State_Met%BXHEIGHT(I,J,1:l))- 0.5 * State_Met%BXHEIGHT(I,J,l)

      if (L==1) physta%pblh(n) = State_Met%PBLH(I,J)

! load aerosol state fo
!IMPORTANT
! for now do not update with GC tracers  
!       state%q(n,L,:) =
! load ptend
!      ptend%lq    = dotend
!      ptend%q     = dqdt
! A FAIRE ZM
    END DO
    END DO
    END DO

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
! apply tendencies ! note ptend.lq can be modified by calcsize 
      do l = 1, pcnst
         if ( .not. ptend%lq(l) ) cycle
         do k = 1, pver
         do i = 1, pcols 
            physta%q(i,k,l) = physta%q(i,k,l) + ptend%q(i,k,l)*deltat  !fab: perhaps an exponential form ?
            physta%q(i,k,l) = max( physta%q(i,k,l), 0.0_r8 )
         end do
         end do
      end do
      
      physta%lchnk = lchnk      
      physta%ncol = pcols      

! WATER UPTAKE    
     call load_pbuf( pbuf, lchnk, pcols, &
        physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens )
     call modal_aero_wateruptake_dr( physta, pbuf )

     call unload_pbuf( pbuf, lchnk, pcols, &
         physta%cld, physta%qqcw, physta%dgncur_a, physta%dgncur_awet,  physta%qaerwat, physta%wetdens )

!-------------------------------------------------------------------------------- 
! switch from q & qqcw mass mixing ratios to volume mixing ratios  vmr and vmrcw
! only adress the gas/aerosol variables in q 
      vmr = 0.0_r8
      vmrcw = 0.0_r8
      do l = imozart, pcnst
         l2 = l - loffset
         vmr(  1:pcols,1:pver,l2) =physta%q(  1:pcols,1:pver,l2)*mwdry/adv_mass(l2)
         vmrcw(1:pcols,1:pver,l2) =physta%qqcw(1:pcols,1:pver,l2)*mwdry/adv_mass(l2)
      end do
!-------------------------------
! GASCHEM 
! not sure if h2so4 formed when aerosol is off in GC fullchem  ? GIVE MORE THOUGHTS 
! for now : consider simple gas chem identic to box model. The gaschem routine is directly
! inserted in the driver module rather than as part of the mam dir.
      vmr_svaa   = vmr  !save before gas chem

! global avg ~= 13 d = 1.12e6 s, daytime avg ~= 5.6e5, noontime peak ~= 3.7e5
      tau_gaschem_simple = 3.0e5  ! so2 gas-rxn timescale (s)

      if (mdo_gaschem > 0) then
         call gaschem_simple_sub(                      &
            lchnk,                                     &
            vmr,                tau_gaschem_simple      )
      end if 

      vmr_svbb = vmr    ! save before cloud chem
      vmrcw_svbb = vmrcw! save before cloud chem

!CLOUDCHEM modifies vmr and vmrcw.I skip for now,  question :  are cloud born aerosol also transported or
!can they be considered as diag variables  
      
 
! call the mam microphysics driver  
!------------------------------
  nstep = 1 ! improve later
  if(.true.) then ! .and. masterproc ) then
     call modal_aero_amicphys_intr(               &
         mdo_gasaerexch,     mdo_rename,          &
         mdo_newnuc,         mdo_coag,            &
         lchnk,    1,     nstep   ,           &
         loffset,  deltat,                        &
         latndx,   lonndx,                        &
         physta%t,   physta%pmid, physta%pdel,    &
         physta%zm,  physta%pblh,                 &
         physta%qv,  physta%cld ,                 &
         vmr,                vmrcw,               &   ! after  cloud chem
         vmr_svaa,                                &   ! before gas chem
         vmr_svbb,           vmrcw_svbb,          &   ! before cloud chem!
!         nqtendbb,           nqqcwtendbb,         &  ! ifdef cambox not enabled for now  
!         dvmrdt_bb,          dvmrcwdt_bb,         &  ! in the interface
         physta%dgncur_a,     physta%dgncur_awet,  &
         physta%wetdens,      physta%qaerwat              )
    end if
 





    IF(masterproc) THEN
      print*,'FAB DRIVER'   

      print*, 'vmr', vmr(10,2,:)  
      
      !      DO N=1,State_Chm%nSpecies 
      !      END DO
    END IF









    !FAB just keep qan example of communication with species
!    IF(IFDOISRP==1.or.IFDOISRP==4.or.IFDOISRP==5)THEN
!       !$OMP PARALLEL DO       &
!       !$OMP DEFAULT( SHARED ) &
!       !$OMP PRIVATE( I, J, L) &
!       !$OMP SCHEDULE( DYNAMIC )
!       DO L = 1, State_Grid%NZ
!       DO J = 1, State_Grid%NY
!       DO I = 1, State_Grid%NX
!          Spc(APMIDS%id_NH3)%Conc(I,J,L) = MASSISRP(I,J,L,1)
!          Spc(APMIDS%id_NH4)%Conc(I,J,L) = MASSISRP(I,J,L,2)
!          Spc(APMIDS%id_HNO3)%Conc(I,J,L)= MASSISRP(I,J,L,3)
!          Spc(APMIDS%id_NIT)%Conc(I,J,L) = MASSISRP(I,J,L,4)
!       ENDDO
!       ENDDO
!       ENDDO
!       !$OMP END PARALLEL DO
!    ENDIF

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

!    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
!    type(physics_state)  :: physta 

    integer :: species_class(pcnst) = -1
!-------initialise namelist parameters

   mdo_gaschem=1

   mdo_gasaerexch=1
   mdo_rename=1
   mdo_newnuc=1
   mdo_coag=1

   masterproc = Input_Opt%amIRoot

! Chemistry timestep [s]
   deltat  = GET_TS_CHEM()

! from regcm
!do j
!do i
!n = (j-jci1)+(i-ici1)*nj+1

!nj = (jci2-jci1+1)
!      npoints = nj*(ici2-ici1+1)


pcols = State_Grid%NX *  State_Grid%Ny
pver  = State_Grid%NZ
plev = pver

!pcols pver pourraient just passer par module plutot que par argument
call MAM_init_basics(  )


call MAM_ALLOCATE ( )

call MAM_init_run ()




END SUBROUTINE MAM_INIT        

!---------------------------------------------------------------
SUBROUTINE MAM_init_basics( )
! equivqlent to the cambox_init_basics

use precision_mod, only :r8 => f8
use constituents, only:   cnst_name, species_class  
use chem_mods, only: adv_mass, gas_pcnst, imozart
use mam_utils, only: solsym, endrun,  iulog, begchunk
use physics_buffer, only: physics_buffer_desc, pbuf_initialize,\
                          pbuf_init_time, pbuf_add_field, pbuf_get_chunk
use buffer, only: dtype_r8
use modal_aero_data, only: nbc, npoa, nsoa, nsoag
use modal_aero_initialize_data
use modal_aero_amicphys, only: mosaic
use modal_aero_calcsize, only: modal_aero_calcsize_reg
use modal_aero_wateruptake, only: modal_aero_wateruptake_reg, modal_aero_wateruptake_init


!type(physics_buffer_desc), pointer :: pbuf2d(:,:)

integer :: l, l2, n,idx

!-----------------------------------------------------------------------------

if (masterproc) then 
  OPEN( unit=iulog , file='mam.log',  status='replace', &
        action='write') 
end if

! configure siulation type
! now only MODAL_AERO_4MODE is enabled 

#if ( ( defined MODAL_AERO_7MODE ) && ( defined MOSAIC_SPECIES ) )
      n = 60
#elif ( defined MODAL_AERO_7MODE ) 
      n = 42
#elif ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined RAIN_EVAP_TO_COARSE_AERO ) ) 
      n = 35
#elif ( defined MODAL_AERO_4MODE_MOM ) 
      n = 31
#elif ( defined MODAL_AERO_4MODE ) 
      n = 28
#elif ( defined MODAL_AERO_3MODE ) 
      n = 25
#else
      call endrun( 'MODAL_AERO_3/4/4MOM/7MODE are all undefined' )
#endif


! perhaps simplify this since we are not in the mozart framework, or adap it to GC 
      n = n + 2*(nbc-1) + 2*(npoa-1) + 2*(nsoa-1)
      l = n - (imozart-1)


     if ( masterproc) write(iulog,'( a,3i5/)') 'pcnst, gas_pcnst, imozart =', pcnst, gas_pcnst, imozart
     if (pcnst /= gas_pcnst+imozart-1) call endrun( '*** bad pcnst aa' )
     if (pcnst /= n                  ) call endrun( '*** bad pcnst bb' )

#if ( defined MODAL_AERO_7MODE )
      if (nbc==1 .and. npoa==1 .and. nsoa==1) then

#if ( defined MOSAIC_SPECIES )
      solsym(:l) = &
      (/ 'H2O2    ','H2SO4   ','SO2     ','DMS     ','NH3     ', &
         'SOAG    ','HNO3    ','HCL     ',                       &
         'so4_a1  ','nh4_a1  ','pom_a1  ','soa_a1  ','bc_a1   ', &
         'ncl_a1  ','no3_a1  ','cl_a1   ','num_a1  ',            &
         'so4_a2  ','nh4_a2  ','soa_a2  ','ncl_a2  ','no3_a2  ', &
         'cl_a2   ','num_a2  ',                                  &
         'pom_a3  ','bc_a3   ','num_a3  ',                       &
         'ncl_a4  ','so4_a4  ','nh4_a4  ','no3_a4  ','cl_a4   ', &
         'num_a4  ',                                             &
         'dst_a5  ','so4_a5  ','nh4_a5  ','no3_a5  ','cl_a5   ', &
         'ca_a5   ','co3_a5  ','num_a5  ',                       &
         'ncl_a6  ','so4_a6  ','nh4_a6  ','no3_a6  ','cl_a6   ', &
         'num_a6  ',                                             &
         'dst_a7  ','so4_a7  ','nh4_a7  ','no3_a7  ','cl_a7   ', &
         'ca_a7   ','co3_a7  ','num_a7  '                        /)
      adv_mass(:l) = &
      (/ 34.0135994_r8, 98.0783997_r8, 64.0647964_r8, 62.1324005_r8, 17.0289402_r8, &
         12.0109997_r8, 63.0123400_r8, 36.4601000_r8,                               &
         96.0635986_r8, 18.0363407_r8, 12.0109997_r8, 12.0109997_r8, 12.0109997_r8, &
         22.9897667_r8, 62.0049400_r8, 35.4527000_r8, 1.00740004_r8,                &
         96.0635986_r8, 18.0363407_r8, 12.0109997_r8, 22.9897667_r8, 62.0049400_r8, &
         35.4527000_r8, 1.00740004_r8,                                              &
         12.0109997_r8, 12.0109997_r8, 1.00740004_r8,                               &
         22.9897667_r8, 96.0635986_r8, 18.0363407_r8, 62.0049400_r8, 35.4527000_r8, &
         1.00740004_r8,                                                             &
         135.064041_r8, 96.0635986_r8, 18.0363407_r8, 62.0049400_r8, 35.4527000_r8, &
         40.0780000_r8, 60.0092000_r8, 1.00740004_r8,                               &
         22.9897667_r8, 96.0635986_r8, 18.0363407_r8, 62.0049400_r8, 35.4527000_r8, &
         1.00740004_r8,                                                             &
         135.064041_r8, 96.0635986_r8, 18.0363407_r8, 62.0049400_r8, 35.4527000_r8, &
         40.0780000_r8, 60.0092000_r8, 1.00740004_r8                                /)
! nacl  58.4424667
! cl    35.4527000
! na    22.9897667
! hcl   36.4601000
! hno3  63.0123400
! no3   62.0049400
! ca    40.0780000
! co3   60.0092000


#else
      solsym(:l) = &
      (/ 'H2O2    ','H2SO4   ','SO2     ','DMS     ','NH3     ', &
         'SOAG    ','so4_a1  ','nh4_a1  ','pom_a1  ','soa_a1  ', &
         'bc_a1   ','ncl_a1  ','num_a1  ','so4_a2  ','nh4_a2  ', &
         'soa_a2  ','ncl_a2  ','num_a2  ','pom_a3  ','bc_a3   ', &
         'num_a3  ','ncl_a4  ','so4_a4  ','nh4_a4  ','num_a4  ', &
         'dst_a5  ','so4_a5  ','nh4_a5  ','num_a5  ','ncl_a6  ', &
         'so4_a6  ','nh4_a6  ','num_a6  ','dst_a7  ','so4_a7  ', &
         'nh4_a7  ','num_a7  ' /)
      adv_mass(:l) = &
      (/ 34.0135994_r8, 98.0783997_r8, 64.0647964_r8, 62.1324005_r8, 17.0289402_r8, &
         12.0109997_r8, 96.0635986_r8, 18.0363407_r8, 12.0109997_r8, 12.0109997_r8, &
         12.0109997_r8, 58.4424667_r8, 1.00740004_r8, 96.0635986_r8, 18.0363407_r8, &
         12.0109997_r8, 58.4424667_r8, 1.00740004_r8, 12.0109997_r8, 12.0109997_r8, &
         1.00740004_r8, 58.4424667_r8, 96.0635986_r8, 18.0363407_r8, 1.00740004_r8, &
         135.064041_r8, 96.0635986_r8, 18.0363407_r8, 1.00740004_r8, 58.4424667_r8, &
         96.0635986_r8, 18.0363407_r8, 1.00740004_r8, 135.064041_r8, 96.0635986_r8, &
         18.0363407_r8, 1.00740004_r8 /)

#endif

      else if (nbc==2 .and. npoa==2 .and. nsoa==1) then
      ! nbc=npoa=2 not fully implemented yet
      call endrun( '*** bad nbc and/or npoa and/or nsoa' )

      solsym(:l) = &
      (/ 'H2O2    ','H2SO4   ','SO2     ','DMS     ','NH3     ', &
         'SOAG    ','so4_a1  ','nh4_a1  ','poma_a1 ', &
                                          'pomb_a1 ','soa_a1  ', &
         'bca_a1  ', &
         'bcb_a1  ','ncl_a1  ','num_a1  ','so4_a2  ','nh4_a2  ', &
         'soa_a2  ','ncl_a2  ','num_a2  ','poma_a3 ','pomb_a3 ', &
                                          'bca_a3  ','bcb_a3  ', &
         'num_a3  ','ncl_a4  ','so4_a4  ','nh4_a4  ','num_a4  ', &
         'dst_a5  ','so4_a5  ','nh4_a5  ','num_a5  ','ncl_a6  ', &
         'so4_a6  ','nh4_a6  ','num_a6  ','dst_a7  ','so4_a7  ', &
         'nh4_a7  ','num_a7  ' /)
      adv_mass(:l) = &
      (/ 34.0135994_r8, 98.0783997_r8, 64.0647964_r8, 62.1324005_r8, 17.0289402_r8, &
         12.0109997_r8, 96.0635986_r8, 18.0363407_r8, 12.0109997_r8, 12.0109997_r8, 12.0109997_r8, &
         12.0109997_r8, 12.0109997_r8, 58.4424667_r8, 1.00740004_r8, 96.0635986_r8, 18.0363407_r8, &
         12.0109997_r8, 58.4424667_r8, 1.00740004_r8, 12.0109997_r8,12.0109997_r8,  12.0109997_r8, 12.0109997_r8, &
         1.00740004_r8, 58.4424667_r8, 96.0635986_r8, 18.0363407_r8, 1.00740004_r8, &
         135.064041_r8, 96.0635986_r8, 18.0363407_r8, 1.00740004_r8, 58.4424667_r8, &
         96.0635986_r8, 18.0363407_r8, 1.00740004_r8, 135.064041_r8, 96.0635986_r8, &
         18.0363407_r8, 1.00740004_r8 /)

      else
         call endrun( '*** bad nbc and/or npoa and/or nsoa' )
      end if

#elif ( defined MODAL_AERO_4MODE_MOM )
      if (nbc==1 .and. npoa==1 .and. nsoa==1 .and. nsoag==1) then
#if ( defined RAIN_EVAP_TO_COARSE_AERO )
      solsym(:l) = &
      (/ 'H2O2          ', 'H2SO4         ', 'SO2           ', 'DMS           ', 'SOAG          ', &
         'so4_a1        ', 'pom_a1        ', 'soa_a1        ', 'bc_a1         ', 'dst_a1        ', &
         'ncl_a1        ', 'mom_a1        ', 'num_a1        ', 'so4_a2        ', 'soa_a2        ', &
         'ncl_a2        ', 'mom_a2        ', 'num_a2        ', 'dst_a3        ', 'ncl_a3        ', &
         'so4_a3        ', 'bc_a3         ', 'pom_a3        ', 'soa_a3        ', 'mom_a3        ', &
         'num_a3        ', 'pom_a4        ', 'bc_a4         ', 'mom_a4        ', 'num_a4        ' /)
      adv_mass(:l) = &
      (/     34.013600_r8,     98.078400_r8,     64.064800_r8,     62.132400_r8,     12.011000_r8, &
            115.107340_r8,     12.011000_r8,     12.011000_r8,     12.011000_r8,    135.064039_r8, &
             58.442468_r8, 250092.672000_r8,      1.007400_r8,    115.107340_r8,     12.011000_r8, &
             58.442468_r8, 250092.672000_r8,      1.007400_r8,    135.064039_r8,     58.442468_r8, &
            115.107340_r8,     12.011000_r8,     12.011000_r8,     12.011000_r8, 250092.672000_r8, &
              1.007400_r8,     12.011000_r8,     12.011000_r8, 250092.672000_r8,      1.007400_r8 /)
#else
      solsym(:l) = &
      (/ 'H2O2    ', 'H2SO4   ', 'SO2     ', 'DMS     ',             &
         'SOAG    ', 'so4_a1  ',             'pom_a1  ', 'soa_a1  ', &
         'bc_a1   ', 'ncl_a1  ', 'dst_a1  ', 'mom_a1  ', 'num_a1  ', &
         'so4_a2  ', 'soa_a2  ', 'ncl_a2  ', 'mom_a2  ', 'num_a2  ', &
         'dst_a3  ', 'ncl_a3  ', 'so4_a3  ', 'num_a3  ',             &
         'pom_a4  ', 'bc_a4   ', 'mom_a4  ', 'num_a4  ' /)
      adv_mass(:l) = &
      (/ 34.0135994_r8, 98.0783997_r8, 64.0647964_r8, 62.1324005_r8,                &
         12.0109997_r8, 115.107340_r8,                12.0109997_r8, 12.0109997_r8, &
         12.0109997_r8, 58.4424667_r8, 135.064041_r8, 250092.672_r8, 1.00740004_r8, &
         115.107340_r8, 12.0109997_r8, 58.4424667_r8, 250092.672_r8, 1.00740004_r8, &
         135.064041_r8, 58.4424667_r8, 115.107340_r8, 1.00740004_r8,                &
         12.0109997_r8, 12.0109997_r8, 250092.672_r8, 1.00740004_r8 /)
#endif
      else
         call endrun( '*** bad nbc and/or npoa and/or nsoa' )
      end if


#elif ( defined MODAL_AERO_4MODE )
      if (nbc==1 .and. npoa==1 .and. nsoa==1 .and. nsoag==1) then

      solsym(:l) = &
      (/ 'H2O2    ', 'H2SO4   ', 'SO2     ', 'DMS     ',             &
         'SOAG    ', 'so4_a1  ',             'pom_a1  ', 'soa_a1  ', &
         'bc_a1   ', 'ncl_a1  ', 'dst_a1  ', 'num_a1  ', 'so4_a2  ', &
         'soa_a2  ', 'ncl_a2  ', 'num_a2  ',                         &
         'dst_a3  ', 'ncl_a3  ', 'so4_a3  ', 'num_a3  ',             &
         'pom_a4  ', 'bc_a4   ', 'num_a4  ' /)
      adv_mass(:l) = &
      (/ 34.0135994_r8, 98.0783997_r8, 64.0647964_r8, 62.1324005_r8,                &
         12.0109997_r8, 115.107340_r8,                12.0109997_r8, 12.0109997_r8, &
         12.0109997_r8, 58.4424667_r8, 135.064041_r8, 1.00740004_r8, 115.107340_r8, &
         12.0109997_r8, 58.4424667_r8, 1.00740004_r8,                               &
         135.064041_r8, 58.4424667_r8, 115.107340_r8, 1.00740004_r8,                &
         12.0109997_r8, 12.0109997_r8, 1.00740004_r8 /)
      else
         call endrun( '*** bad nbc and/or npoa and/or nsoa' )
      end if

#else
!if ( defined MODAL_AERO_3MODE )
      if (nbc==1 .and. npoa==1 .and. nsoa==1 .and. nsoag==1) then

      solsym(:l) = &
      (/ 'H2O2    ', 'H2SO4   ', 'SO2     ', 'DMS     ',             &
         'SOAG    ', 'so4_a1  ',             'pom_a1  ', 'soa_a1  ', &
         'bc_a1   ', 'ncl_a1  ', 'dst_a1  ', 'num_a1  ', 'so4_a2  ', &
         'soa_a2  ', 'ncl_a2  ', 'num_a2  ',                         &
         'dst_a3  ', 'ncl_a3  ', 'so4_a3  ', 'num_a3  ' /)
      adv_mass(:l) = &
      (/ 34.0135994_r8, 98.0783997_r8, 64.0647964_r8, 62.1324005_r8,                &
         12.0109997_r8, 115.107340_r8,                12.0109997_r8, 12.0109997_r8, &
         12.0109997_r8, 58.4424667_r8, 135.064041_r8, 1.00740004_r8, 115.107340_r8, &
         12.0109997_r8, 58.4424667_r8, 1.00740004_r8,                               &
         135.064041_r8, 58.4424667_r8, 115.107340_r8, 1.00740004_r8 /)

      else
         call endrun( '*** bad nbc and/or npoa and/or nsoa' )
      end if

#endif


      cnst_name(1) = 'QVAPOR'
      cnst_name(2) = 'CLDLIQ'
      cnst_name(3) = 'CLDICE'
      cnst_name(4) = 'NUMLIQ'
      cnst_name(5) = 'NUMICE'
      cnst_name(imozart:pcnst) = solsym(1:gas_pcnst)

      mosaic = .false.

     
     IF (masterproc) THEN
     
      write(iulog,'(/a)') &
         'l, l2, cnst_name(l), solsym(l2), adv_mass(l2)'
      do l = 1, pcnst
         if (l < imozart) then
            write(iulog,'(i4,6x,a)') l, cnst_name(l)
         else
            l2 = l - imozart + 1
            if (adv_mass(l2) < 1.0e5_r8) then
               write(iulog,'(2i4,2x,2a,f9.3)') l, l2, cnst_name(l), solsym(l2), adv_mass(l2)
            else
               write(iulog,'(2i4,2x,2a,1pe16.8)') l, l2, cnst_name(l), solsym(l2), adv_mass(l2)
            end if
         end if
      end do
     END IF 

     species_class = -1       
     call modal_aero_register(species_class)
     call modal_aero_calcsize_reg()
     call modal_aero_wateruptake_reg()

      call pbuf_init_time()
      call pbuf_add_field( 'CLD',  'global', dtype_r8, (/pcols, pver/), idx )
      call pbuf_initialize( pbuf2d)

      call modal_aero_initialize(pbuf2d, imozart, species_class )
      call modal_aero_wateruptake_init( pbuf2d )

      lchnk = begchunk
      loffset = imozart -1 
      pbuf => pbuf_get_chunk( pbuf2d, lchnk)
 


END SUBROUTINE MAM_init_basics         
!---------------------------------------------------------------------------------
SUBROUTINE MAM_ALLOCATE () 

use physics_types, only : physics_state 
use modal_aero_data, only : ntot_amode
! FAB peut etre remplacer physics_type par un MAM type .. 

!
! FAB: for now use parameter defined in mod_mam_utils        

!TYPE(physics_state), intent(out)  :: physta

integer ::  as 

allocate (physta%pblh(pcols) , stat=as)

allocate (physta%t(pcols,pver),stat=as)
allocate (physta%pmid(pcols,pver),stat=as)
allocate (physta%pdel(pcols,pver),stat=as)
allocate (physta%zm(pcols,pver), stat=as) 
allocate (physta%cld(pcols,pver), stat=as) 
allocate (physta%relhum(pcols,pver) , stat=as)
allocate (physta%qv(pcols,pver) , stat=as)


allocate (physta%q(pcols,pver,pcnst),stat=as)
allocate (physta%qqcw(pcols,pver,pcnst),stat=as) 

allocate (physta%dgncur_a(pcols,pver,ntot_amode),stat=as)
allocate(physta%dgncur_awet(pcols,pver,ntot_amode),stat=as)
allocate(physta%qaerwat(pcols,pver,ntot_amode),stat=as)
allocate(physta%wetdens(pcols,pver,ntot_amode),stat=as)

allocate(ptend%q(pcols,pver,pcnst))
allocate(ptend%lq(pcnst))



physta%pblh=0._r8 
physta%t=0._r8
physta%pmid=0._r8
physta%pdel=0._r8
physta%zm=0._r8
physta%cld=0._r8
physta%relhum =0._r8
physta%qv =0._r8
physta%q =0._r8
physta%qqcw =0._r8
physta%dgncur_a =0._r8
physta%dgncur_awet =0._r8
physta%qaerwat =0._r8
physta%wetdens =0._r8

ptend%lq=.false.
ptend%q=0._r8

END SUBROUTINE MAM_ALLOCATE         
 
!-----------------------------------------------------------------------------

SUBROUTINE MAM_init_run ()!

use physconst, only: pi, mwdry 
use mam_utils, only: pcols,pver, endrun

use modal_aero_amicphys, only :gaexch_h2so4_uptake_optaa, newnuc_h2so4_conc_optaa,&
           dens_aer, iaer_bc, iaer_pom, iaer_so4, iaer_soa, iaer_ncl, &
           iaer_mom, iaer_dst        

use modal_aero_data
use physics_types, only : physics_state
use constituents, only:  cnst_name, cnst_get_ind

!initial composition for q  
! should go on a namelist or initialized somehow from geos 
real(r8) :: numc1, numc2, numc3, numc4,                     &
                  mfso41, mfpom1, mfsoa1, mfbc1, mfdst1, mfncl1,  &
                  mfso42, mfsoa2, mfncl2,                         &
                  mfdst3, mfncl3, mfso43, mfbc3, mfpom3,  mfsoa3, &
                  mfpom4, mfbc4,                                  &
                  qso2, qh2so4, qsoag
real(r8) :: tmpfso4, tmpfnh4, tmpfsoa, tmpfpom, &
                  tmpfbcx, tmpfncl, tmpfdst, tmpfmom
real(r8) :: tmpfno3, tmpfclx, tmpfcax, tmpfco3

real(r8) :: tmpdens, tmpvol, tmpmass, sx


real(r8), allocatable :: aircon(:,:) 
real(r8), pointer :: q(:,:,:), dgncur_a(:,:,:)

integer :: i,k,n


!------------------------------------------------------------------------------

numc1          = 1.e8    ! unit: #/m3
numc2          = 1.e9
numc3          = 1.e5
numc4          = 2.e8

mfso41         = 0.3_r8
mfpom1         = 0._r8
mfsoa1         = 0.3_r8
mfbc1          = 0._r8
mfdst1         = 0._r8
mfncl1         = 0.4_r8

mfso42         = 0.3_r8
mfsoa2         = 0.3_r8
mfncl2         = 0.4_r8

mfdst3         = 0._r8
mfncl3         = 0.4_r8
mfso43         = 0.3_r8
mfbc3          = 0._r8
mfpom3         = 0._r8
mfsoa3         = 0.3_r8

mfpom4         = 0._r8
mfbc4          = 1._r8

      ! check if mass fraction is larger than one
      if (mfso41+mfpom1+mfsoa1+mfbc1+mfdst1+mfncl1 .gt. 1._r8) then
          print *, "The summed mass fraction is > 1 in mode 1"
          stop
      end if
      if (mfso42+mfsoa2+mfncl2 .gt. 1._r8) then
          print *, "The summed mass fraction is > 1 in mode 2"
          stop
      end if
      if (mfdst3+mfncl3+mfso43+mfbc3+mfpom3+mfsoa3 .gt. 1._r8) then
          print *, "The summed mass fraction is > 1 in mode 3"
          stop
      end if
      if (mfpom4+mfbc4 .gt. 1._r8) then
          print *, "The summed mass fraction is > 1 in mode 4"
          stop
      end if

gaexch_h2so4_uptake_optaa =  2 
newnuc_h2so4_conc_optaa   =  2 

!FAB quick and dirty , their should not be any met state related operation here
! fix this but perhaps this will go in a first call and not restart statement at the first integrartions tep  
allocate(aircon(pcols,pver))
aircon(:,:) = 0.0345_r8 * 0.8 ! kmol/m3 for a density of 0.8 kg/m3

q => physta%q
dgncur_a => physta%dgncur_a


! initialize the aerosol/number mixing ratio
! adapted to mam4 for now 
      do k = 1, pver
         do i = 1, pcols 
            do  n = 1, ntot_amode

                sx = log( sigmag_amode(n) )

                if      (n == 1) then
                   dgncur_a(i,k,n) = dgnum_amode(n)  ! 0.20e-6_r8 ! m
                   tmpfsoa      = mfsoa1
                   tmpfso4      = mfso41
                   tmpfncl      = mfncl1
                   tmpfdst      = mfdst1
                   tmpfpom      = mfpom1
                   tmpfbcx      = mfbc1
                   tmpfmom      = 1._r8 - tmpfsoa - tmpfso4 - &
                                  tmpfncl - tmpfdst - tmpfpom - tmpfbcx
                else if (n == 2) then
                   dgncur_a(i,k,n) = dgnum_amode(n)  ! 0.04e-6_r8
                   tmpfsoa      = mfsoa2
                   tmpfso4      = mfso42
                   tmpfncl      = mfncl2
                   tmpfdst      = 0._r8
                   tmpfpom      = 0._r8
                   tmpfbcx      = 0._r8
                   tmpfmom      = 1._r8 - tmpfsoa - tmpfso4 - &
                                  tmpfncl - tmpfdst - tmpfpom - tmpfbcx
                else if (n == 3) then
                   dgncur_a(i,k,n) = dgnum_amode(n)  ! 2.00e-6_r8
                   tmpfsoa      = mfsoa3
                   tmpfso4      = mfso43
                   tmpfncl      = mfncl3
                   tmpfdst      = mfdst3
                   tmpfpom      = mfpom3
                   tmpfbcx      = mfbc3
                   tmpfmom      = 1._r8 - tmpfsoa - tmpfso4 - &
                                  tmpfncl - tmpfdst - tmpfpom - tmpfbcx
                else if (n == 4) then
                   dgncur_a(i,k,n) = dgnum_amode(n)  ! 0.08e-6_r8
                   tmpfsoa      = 0._r8
                   tmpfso4      = 0._r8
                   tmpfncl      = 0._r8
                   tmpfdst      = 0._r8
                   tmpfpom      = mfpom4
                   tmpfbcx      = mfbc4
                   tmpfmom      = 1._r8 - tmpfsoa - tmpfso4 - &
                                  tmpfncl - tmpfdst - tmpfpom - tmpfbcx
                end if
                ! q(i,k,numptr_amode(n)) = #/kg-air
                if (n == modeptr_aitken) then
                   q(i,k,numptr_amode(n)) = numc2 / aircon(i,k) / mwdry
                   l_num_a2 = numptr_amode(n)
                   l_so4_a2 = lptr_so4_a_amode(n)
                else if (n == modeptr_accum) then
                   q(i,k,numptr_amode(n)) = numc1 / aircon(i,k) / mwdry
                   l_num_a1 = numptr_amode(n)
                   l_so4_a1 = lptr_so4_a_amode(n)
                else if (n == modeptr_pcarbon) then
                   q(i,k,numptr_amode(n)) = numc4 / aircon(i,k) / mwdry
                else
                   q(i,k,numptr_amode(n)) = numc3 / aircon(i,k) / mwdry
                end if

                ! tmpvol: m3-dry-aerosol/kg-air
                tmpvol  = q(i,k,numptr_amode(n)) * &
                          (dgncur_a(i,k,n)**3) * &
                          (pi/6.0_r8) * exp(4.5_r8*sx*sx)
                tmpdens = 1.0_r8 /                           &
                          ( (tmpfsoa / dens_aer(iaer_soa)) + &
                            (tmpfso4 / dens_aer(iaer_so4)) + &
                            (tmpfbcx / dens_aer(iaer_bc )) + &
                            (tmpfpom / dens_aer(iaer_pom)) + &
                            (tmpfncl / dens_aer(iaer_ncl)) + &
                            (tmpfdst / dens_aer(iaer_dst)) )
!                            (tmpfmom / dens_aer(iaer_mom))   )
                tmpmass = tmpvol*tmpdens   ! kg-dry-aerosol/kg-air
                l_so4a = lptr_so4_a_amode(n)
                l_nh4a = -1
                l_soaa = lptr_soa_a_amode(n)
                l_poma = lptr_pom_a_amode(n)
!                if (npoa == 2) l_poma = lptr_poma_a_amode(n)
                l_bcxa = lptr_bc_a_amode(n)
!                if (nbc  == 2) l_bcxa = lptr_bca_a_amode(n)
                l_ncla = lptr_nacl_a_amode(n)
                l_dsta = lptr_dust_a_amode(n)
                l_moma = lptr_mom_a_amode(n)
#if ( defined MOSAIC_SPECIES )
                l_no3a = lptr_no3_a_amode(n)
                l_clxa = lptr_cl_a_amode(n)
                l_caxa = lptr_ca_a_amode(n)
                l_co3a = lptr_co3_a_amode(n)
#else
                l_no3a = -1
                l_clxa = -1
                l_caxa = -1
                l_co3a = -1
#endif
                ! q array return kg-aer/kg-air
                if (l_so4a > 0) q(i,k,l_so4a) = tmpmass*tmpfso4
                if (l_nh4a > 0) q(i,k,l_nh4a) = tmpmass*tmpfnh4
                if (l_soaa > 0) q(i,k,l_soaa) = tmpmass*tmpfsoa
                if (l_poma > 0) q(i,k,l_poma) = tmpmass*tmpfpom
                if (l_bcxa > 0) q(i,k,l_bcxa) = tmpmass*tmpfbcx
                if (l_dsta > 0) q(i,k,l_dsta) = tmpmass*tmpfdst
                if (l_ncla > 0) q(i,k,l_ncla) = tmpmass*tmpfncl
                if (l_moma > 0) q(i,k,l_moma) = tmpmass*tmpfmom
                if (l_no3a > 0) q(i,k,l_no3a) = tmpmass*tmpfno3
                if (l_clxa > 0) q(i,k,l_clxa) = tmpmass*tmpfclx
                if (l_caxa > 0) q(i,k,l_caxa) = tmpmass*tmpfcax
                if (l_co3a > 0) q(i,k,l_co3a) = tmpmass*tmpfco3

            end do ! n
         end do ! i
      end do ! k   

!TEMPORAIRE EN ATTENDANT COUPLAGE  initialize the gas mixing ratio value from box namelist
! set trace gases
      call cnst_get_ind( 'SOAG',  l_soag,   .false. )
      call cnst_get_ind( 'SO2',   l_so2g,   .false. )
      call cnst_get_ind( 'NH3',   l_nh3g,   .false. )
      call cnst_get_ind( 'HNO3',  l_hno3g,  .false. )
      call cnst_get_ind( 'HCL',   l_hclg,   .false. )
      call cnst_get_ind( 'H2SO4', l_h2so4g, .false. )

      q(:,:,l_so2g)   = 1.e-4
      q(:,:,l_soag)   = 5.e-10
      q(:,:,l_h2so4g) = 1.e-13


END SUBROUTINE MAM_init_run


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

!----------------------------------------------------------------------
subroutine gaschem_simple_sub(                     &
   lchnk,                        &
   x,  tau_gaschem_simple      )

! !USES:
!use modal_aero_data

use constituents,      only:  cnst_name, cnst_get_ind
use chem_mods, only: gas_pcnst
implicit none

! !PARAMETERS:
   integer,  intent(in)    :: lchnk                ! chunk identifier
   real(r8), intent(in)    :: tau_gaschem_simple(pcols,pver)
   real(r8), intent(inout) :: x(pcols,pver,gas_pcnst) ! tracer mixing ratio (TMR) array
                                                   ! *** MUST BE  #/kmol-air for number
! local variables
   integer, parameter :: method_soa = 2
!     method_soa=0 is no uptake
!     method_soa=1 is irreversible uptake done like h2so4 uptake
!     method_soa=2 is reversible uptake using subr modal_aero_soaexch
! FAB not treated here think about it !! 
   integer :: i
   integer :: k
   integer :: i_h2so4g, i_so2g

   real (r8) :: tmpa, tmpb

   ! set gas species indices

   i_h2so4g = l_h2so4g - loffset
   i_so2g = l_so2g - loffset

   do k = 1, pver
   do i = 1, pcols
      tmpa = x(i,k,i_so2g)*exp( -deltat/tau_gaschem_simple(i,k) )
      tmpb = x(i,k,i_so2g) - tmpa
      x(i,k,i_so2g) = tmpa
      x(i,k,i_h2so4g) = tmpb
   end do
   end do

   return
   end subroutine gaschem_simple_sub


END MODULE MAM_DRIV_MOD

!FAB#endif

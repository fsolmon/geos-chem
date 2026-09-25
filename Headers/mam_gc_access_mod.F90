!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mam_gc_access_mod.F90
!
! !DESCRIPTION: Small read-only accessors for MAM aerosol mass in GEOS-Chem
!  units, for code OUTSIDE the KPP integrator that would otherwise read the
!  standard-aerosol tracers.
!
!  FAB (MAM-decouple-std, Step 5). Rationale and the rule this module exists to
!  serve (devnotes 11.22):
!
!    * Values consumed INSIDE the KPP integrator -- the C() array as seen by
!      Integrate and by the rate-law functions it calls -- must come from the
!      SHADOW (MAM_KPP_Shadow_In/Out in mam_driv_mod). Those functions are
!      re-evaluated many times within one step while the species are evolving,
!      so they have to see the integrated value, not a frozen pre-step copy.
!      Do NOT replace those with the accessors here.
!
!    * Values read OUTSIDE the integration, i.e. anything doing
!      State_Chm%Species(id_<STD>)%Conc(I,J,L) before or after KPP, should use
!      MAM explicitly. There is no consistency constraint there, and the STD
!      tracer is at best a copy of MAM (NIT/NITs since Step 3b) and at worst not
!      a copy at all (SALA/SALC are never shadowed -- they are a parallel legacy
!      sea-salt population with their own emissions and removal).
!
!  Deliberately a LEAF module: it uses only State_Chm_Mod (for Ind_) and
!  precision_mod, resolving species by NAME rather than through
!  modal_aero_data. That keeps it free of the MAM library build guards, lets it
!  compile unconditionally, and avoids the dependency cycles that would follow
!  from pulling mam_driv_mod (physics_buffer, physics_types, constituents) into
!  light modules such as photolysis_mod.
!
!  Lives in Headers/ (moved from GeosCore/ for the SET_SO2 cloud-pH fix,
!  devnotes 11.24.5) so that the KPP library, which cannot depend on GeosCore,
!  can use it too.
!
!  UNITS: whatever State_Chm%Species(:)%Conc currently holds -- these routines
!  only sum existing species, they never convert. Ratios of two accessors are
!  therefore always unit-safe, which is how they are meant to be used.
!
!  MODE CONVENTION: the trailing digit of the species name is the MAM mode,
!  1 = accumulation, 2 = Aitken, 3 = coarse. MAM_INIT asserts this ordering
!  against modeptr_accum/aitken/coarse and stops the run if it ever changes,
!  so the name-based lookup here is safe.
!
!  "Fine" = accumulation + Aitken, matching exactly what feeds the KPP NIT
!  shadow; "coarse" = mode 3, matching NITs.
!
! !INTERFACE:
!
MODULE Mam_Gc_Access_Mod
!
! !USES:
!
  USE Precision_Mod
  USE State_Chm_Mod, ONLY : ChmState, Ind_

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: MAM_Init_Access
  PUBLIC :: MAM_NO3_Fine,     MAM_NO3_Coarse
  PUBLIC :: MAM_SeaSalt_Fine, MAM_SeaSalt_Coarse
  PUBLIC :: MAM_SO4_Fine,     MAM_SO4_Coarse
  PUBLIC :: MAM_NH4_Fine,     MAM_NH4_Coarse
  PUBLIC :: MAM_Na_Fine,      MAM_Na_Coarse
  PUBLIC :: MAM_Cl_Fine,      MAM_Cl_Coarse
  PUBLIC :: MAM_SeaSaltMW_Coarse
  PUBLIC :: MAM_Access_Ready
!
! !PRIVATE DATA:
!
  ! Cached GC species indices per MAM mode (-1 = species not carried).
  INTEGER, SAVE :: id_NO3 (3) = -1   ! MAMNO31/2/3   nitrate
  INTEGER, SAVE :: id_SSLT(3) = -1   ! MAMSSLT1/2/3  sea salt (Na+ only in the
                                     !               MOSAIC build)
  INTEGER, SAVE :: id_CL  (3) = -1   ! MAMCL1/2/3    chloride
  INTEGER, SAVE :: id_SO4 (3) = -1   ! MAMSO41/2/3   sulfate (incl. primary
                                     !               sea-salt SO4 when emitted)
  INTEGER, SAVE :: id_NH4 (3) = -1   ! MAMNH41/2/3   ammonium
  LOGICAL, SAVE :: isReady    = .FALSE.

CONTAINS

!------------------------------------------------------------------------------
! Resolve and cache the species indices. Idempotent. Called explicitly from
! MAM_INIT, and lazily by the accessors, so no call-ordering dependency exists.
!------------------------------------------------------------------------------
  SUBROUTINE MAM_Init_Access()

    IF ( isReady ) RETURN

    id_NO3 (1) = Ind_('MAMNO31' ) ; id_NO3 (2) = Ind_('MAMNO32' )
    id_NO3 (3) = Ind_('MAMNO33' )
    id_SSLT(1) = Ind_('MAMSSLT1') ; id_SSLT(2) = Ind_('MAMSSLT2')
    id_SSLT(3) = Ind_('MAMSSLT3')
    id_CL  (1) = Ind_('MAMCL1'  ) ; id_CL  (2) = Ind_('MAMCL2'  )
    id_CL  (3) = Ind_('MAMCL3'  )
    id_SO4 (1) = Ind_('MAMSO41' ) ; id_SO4 (2) = Ind_('MAMSO42' )
    id_SO4 (3) = Ind_('MAMSO43' )
    id_NH4 (1) = Ind_('MAMNH41' ) ; id_NH4 (2) = Ind_('MAMNH42' )
    id_NH4 (3) = Ind_('MAMNH43' )

    isReady = .TRUE.

  END SUBROUTINE MAM_Init_Access

!------------------------------------------------------------------------------
  FUNCTION MAM_Access_Ready() RESULT( yes )
    LOGICAL :: yes
    yes = isReady
  END FUNCTION MAM_Access_Ready

!------------------------------------------------------------------------------
! Sum of the non-negative concentrations of ids(:) in box (I,J,L).
! Absent species (id <= 0) contribute nothing.
!------------------------------------------------------------------------------
  FUNCTION SumPos( State_Chm, ids, I, J, L ) RESULT( c )

    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: ids(:), I, J, L
    REAL(fp)                   :: c
    INTEGER                    :: k

    c = 0.0_fp
    DO k = 1, SIZE( ids )
       IF ( ids(k) <= 0 ) CYCLE
       c = c + MAX( State_Chm%Species(ids(k))%Conc(I,J,L), 0.0_fp )
    END DO

  END FUNCTION SumPos

!------------------------------------------------------------------------------
! MAM fine-mode (accumulation + Aitken) nitrate. Analogue of STD NIT.
!------------------------------------------------------------------------------
  FUNCTION MAM_NO3_Fine( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_NO3(1:2), I, J, L )
  END FUNCTION MAM_NO3_Fine

!------------------------------------------------------------------------------
! MAM coarse-mode nitrate. Analogue of STD NITs.
!------------------------------------------------------------------------------
  FUNCTION MAM_NO3_Coarse( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_NO3(3:3), I, J, L )
  END FUNCTION MAM_NO3_Coarse

!------------------------------------------------------------------------------
! MAM fine-mode sea salt as Na+ + Cl-, the analogue of STD SALA.
!
! Why Na+ + Cl- and not Na+ alone or a reconstructed total mass: in the MOSAIC
! build MAM's 'seasalt' species carries Na+ ONLY (MAM4 naming heritage), with
! the chloride in a separate MAMCL species. STD SALA is a single lumped
! sea-salt species with MW_g = 31.4, and since NaCl is 58.44 g/mol, SALA
! expressed as a number density is ~58.44/31.4 = 1.86 ~ 2 times the NaCl
! formula-unit count -- i.e. STD SALA already counts sea salt roughly per ION.
! Na+ + Cl- is therefore the like-for-like MAM quantity, and unlike
! MAMSSLT/SS_NA_MF it needs no hard-coded mass fraction.
!------------------------------------------------------------------------------
  FUNCTION MAM_SeaSalt_Fine( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_SSLT(1:2), I, J, L ) +                          &
        SumPos( State_Chm, id_CL  (1:2), I, J, L )
  END FUNCTION MAM_SeaSalt_Fine

!------------------------------------------------------------------------------
! MAM coarse-mode sea salt as Na+ + Cl-, the analogue of STD SALC.
!------------------------------------------------------------------------------
  FUNCTION MAM_SeaSalt_Coarse( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_SSLT(3:3), I, J, L ) +                          &
        SumPos( State_Chm, id_CL  (3:3), I, J, L )
  END FUNCTION MAM_SeaSalt_Coarse

!------------------------------------------------------------------------------
! Single-ion accessors (fine = accumulation + Aitken, coarse = mode 3).
!
! Added for the cloud-pH inputs of SET_SO2 (devnotes 11.24.5), which need each
! ion on its own rather than lumped. The MAM species carry the molecular weight
! of the ion itself (Na+ 22.99, NH4+ 18.04, SO4= 96, NO3- 62, Cl- 35.45), so in
! number-density units these return ion counts directly.
!------------------------------------------------------------------------------
  FUNCTION MAM_SO4_Fine( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_SO4(1:2), I, J, L )
  END FUNCTION MAM_SO4_Fine

  FUNCTION MAM_SO4_Coarse( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_SO4(3:3), I, J, L )
  END FUNCTION MAM_SO4_Coarse

  FUNCTION MAM_NH4_Fine( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_NH4(1:2), I, J, L )
  END FUNCTION MAM_NH4_Fine

  FUNCTION MAM_NH4_Coarse( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_NH4(3:3), I, J, L )
  END FUNCTION MAM_NH4_Coarse

  ! Na+ alone (MAMSSLT is Na+ only in the MOSAIC build). Contrast with
  ! MAM_SeaSalt_*, which is Na+ + Cl- as the analogue of STD SALA/SALC.
  FUNCTION MAM_Na_Fine( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_SSLT(1:2), I, J, L )
  END FUNCTION MAM_Na_Fine

  FUNCTION MAM_Na_Coarse( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_SSLT(3:3), I, J, L )
  END FUNCTION MAM_Na_Coarse

  FUNCTION MAM_Cl_Fine( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_CL(1:2), I, J, L )
  END FUNCTION MAM_Cl_Fine

  FUNCTION MAM_Cl_Coarse( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = SumPos( State_Chm, id_CL(3:3), I, J, L )
  END FUNCTION MAM_Cl_Coarse

!------------------------------------------------------------------------------
! MAM coarse sea salt as a MASS-like quantity: sum over Na+ and Cl- of
! Conc * MW_g, i.e. [Conc units] x [g/mol]. Divide by AIRMW (and apply the
! usual CVF/AD/AIRVOL factors) to get kg/m3. Replaces STD SALC * MW_SALC where
! a sea-salt MASS is needed (fullchem_HetDropChem coarse number, 11.24.6).
!------------------------------------------------------------------------------
  FUNCTION MAM_SeaSaltMW_Coarse( State_Chm, I, J, L ) RESULT( c )
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp)                   :: c
    IF ( .NOT. isReady ) CALL MAM_Init_Access()
    c = 0.0_fp
    IF ( id_SSLT(3) > 0 ) c = c + MAX( State_Chm%Species(id_SSLT(3))%Conc(I,J,L), 0.0_fp ) &
                                * State_Chm%SpcData(id_SSLT(3))%Info%MW_g
    IF ( id_CL(3)   > 0 ) c = c + MAX( State_Chm%Species(id_CL(3))%Conc(I,J,L),   0.0_fp ) &
                                * State_Chm%SpcData(id_CL(3))%Info%MW_g
  END FUNCTION MAM_SeaSaltMW_Coarse

END MODULE Mam_Gc_Access_Mod

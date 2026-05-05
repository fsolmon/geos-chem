!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: species_mod.F90
!
! !DESCRIPTION: Module SPECIES\_MOD contains types and routines to define
!  the GEOS-Chem species object.
!\\
!\\
! !INTERFACE:
!
MODULE Mam_container_Mod
!
! USES:
!
  USE Precision_Mod
  USE ErrCode_Mod
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_Mam_Container
  PUBLIC :: Cleanup_Mam_Container
!
! !PUBLIC TYPES:
!
  ! Type for single MAM modes
  !=========================================================================
  TYPE, PUBLIC :: MAMContainer 
     real(fp)          :: vol2num        ! converts aerosol volume in m3 to number 
     REAL(fp), POINTER :: dryrad (:,:,:) ! vol geo mean dry radius
     REAL(fp), POINTER :: wetrad (:,:,:) ! ------------ wet radius 
     REAL(fp), POINTER :: nudryrad (:,:,:) !num geo mean dry radius 
     REAL(fp), POINTER :: nuwetrad (:,:,:) !------------ wet radius 
     REAL(fp), POINTER :: aerdens(:,:,:) ! aerosol effective density  
     REAL(fp), POINTER :: hygro(:,:,:) ! aerosol hygroscopicity (volume average)
     REAL(fp), POINTER :: pH(:,:,:)   ! aerosol pH [-3, 14]

     REAL(fp), POINTER :: aerwat(:,:,:)!mam water mass concentration
     REAL(fp), POINTER :: so4(:,:,:) ! mam so4 mass concentration 
     REAL(fp), POINTER :: bc(:,:,:) ! mam  mass concentration 
     REAL(fp), POINTER :: pom(:,:,:) ! mam  mass concentration 
     REAL(fp), POINTER :: soa(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: sslt(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: dust(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: nh4(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: no3(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: ca(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: co3(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: cl(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: mom(:,:,:) ! mam  mass concentration
   
     REAL(fp), POINTER :: nu(:,:,:) ! mam number concentration 

     ! Per-mode SW optical properties (all bands).
     ! 4th dimension is nswbands; the mode index is the GCMAM array subscript.
     REAL(fp), POINTER :: tauxar(:,:,:,:) ! (NX,NY,NZ,NSWBANDS) aerosol optical depth
     REAL(fp), POINTER :: ssa   (:,:,:,:) ! (NX,NY,NZ,NSWBANDS) single-scattering albedo
     REAL(fp), POINTER :: g     (:,:,:,:) ! (NX,NY,NZ,NSWBANDS) asymmetry parameter

     LOGICAL           :: lso4, lbc, lpom, lsoa, lsslt, ldust
     LOGICAL           :: lnh4, lno3, lca, lco3, lcl, lmom


  END TYPE MAMContainer 

!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SpcData_Init
!
! !DESCRIPTION: Routine SpcData\_Init initializes species database object.
!  This is an array where each element is of type Species.  This object holds
!  the metadata for each species (name, molecular weight, Henry's law
!  constants, drydep info, wetdep info, etc.
!\\
!\\
! !INTERFACE:
!
!EOC
 SUBROUTINE Init_MAM_Container( Input_Opt, State_Grid,SpcLocData, GCMAM , RC )
!
! !USES:
!
    USE CMN_Size_Mod,   ONLY : NAER
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE Species_Mod, ONLY :SpcPtr
    !
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(GrdState),      INTENT(IN)  :: State_Grid ! Grid object

    TYPE(SpcPtr),       INTENT(IN)  :: SpcLocData    (:) ! GC Species database
    !
! !INPUT/OUTPUT PARAMETERS:
!
    Type(MAMcontainer) ,   POINTER   :: GCMAM(:)       ! mam object 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,             INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: errMsg, thisLoc
    INTEGER            :: NX, NY, NZ, n, s
    integer, parameter :: nswbands= 14 ! It has to match exactly the RRTMG , think about making this cleaner
    !======================================================================
    ! Init_AerMass_Container starts here
    !======================================================================
! Initialize local variables
    RC = GC_SUCCESS
    thisLoc = ' -> at Init_MAM_Container (in module Headers/aermass_container_mod.F90)'
    NX = State_Grid%NX
    NY = State_Grid%NY
    NZ = State_Grid%NZ

    ! Exit immediately if this is a dry-run
    IF ( Input_Opt%DryRun ) RETURN

    !======================================================================
    ! Initialize arrays
    !======================================================================
    
    
    ! please define a MAM mode number , jeez..
    do n= 1, size(GCMAM) 

     GCMAM(n)%lso4 = .false.
     GCMAM(n)%lbc = .false.
     GCMAM(n)%lpom = .false.
     GCMAM(n)%lsoa = .false.
     GCMAM(n)%lsslt = .false.
     GCMAM(n)%ldust = .false.
     GCMAM(n)%lnh4 = .false.
     GCMAM(n)%lno3 = .false.
     GCMAM(n)%lca = .false.
     GCMAM(n)%lco3 = .false.
     GCMAM(n)%lcl = .false.
     GCMAM(n)%lmom = .false.

    ! modal geo dry radius number
    ALLOCATE(GCMAM(n)%nudryrad( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'NUDRYRAD', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array NUDRYRAD!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%nudryrad = 0.1E-6_fp

    ALLOCATE(GCMAM(n)%nuwetrad( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'NUWETRAD', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array NUWETRAD!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%nuwetrad = 0.1E-6_fp

    ALLOCATE(GCMAM(n)%dryrad( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'DRYRAD', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array DRYRAD!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%DRYRAD = 0.1E-6_fp

    ! wet radius
    ALLOCATE( GCMAM(n)%wetrad( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'WETRAD', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array WETRAD!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%WETRAD = 0.1E-6_fp

    ! 
    ALLOCATE( GCMAM(n)%aerdens( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'AERDENS', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array AERDENS!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%AERDENS = 1500._fp ! default needs to be fixes for first time step

    ALLOCATE( GCMAM(n)%hygro( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'HYGRO', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array AERDENS!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%hygro = 0.2_fp ! default needs to be fixes for first time step

    ALLOCATE( GCMAM(n)%pH( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'PH', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array PH!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%pH = 7._fp

   ALLOCATE( GCMAM(n)%aerwat( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'AERWAT', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array AERWAT!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%aerwat = 0._fp ! default needs to be fixes for first time step
    
    ! Per-mode SW optical properties (all bands)
    ALLOCATE( GCMAM(n)%tauxar( NX, NY, NZ, NSWBANDS ), STAT=RC )
    CALL GC_CheckVar( 'TAUXAR', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array TAUXAR!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%tauxar = 0._fp

    ALLOCATE( GCMAM(n)%ssa( NX, NY, NZ, NSWBANDS ), STAT=RC )
    CALL GC_CheckVar( 'SSA', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array SSA!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%ssa = 1._fp

    ALLOCATE( GCMAM(n)%g( NX, NY, NZ, NSWBANDS ), STAT=RC )
    CALL GC_CheckVar( 'G', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array G!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%g = 0.5_fp
    
    
! now define the mass concentration per mode. They will be used for diagnostics and for 
! input to e.g. optical depth, Fast-J etc . Note that for each MAM mode, 
! not all species are relevant and corresponding mass array are not allocated to save mem !
! Always think making an appropriate est when using GCMAM(n)%spec(:,:,:) elsewhere in the code.
! The info on relevant species per mode is accessible through the species_data.yml  
! 
    do s = 1,size(SpcLocData)
    if(SpcLocData(s)%info%MamModId == n ) then ! test if mam species

          if (SpcLocData(s)%info%name(4:6) == 'SO4') then  
              GCMAM(n)%lso4 = .true.
              ALLOCATE( GCMAM(n)%so4( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'SO4', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:5) == 'BC') then 
              GCMAM(n)%lbc = .true.
              ALLOCATE( GCMAM(n)%bc( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'BC', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:6) == 'POM') then 
              GCMAM(n)%lpom = .true.
              ALLOCATE( GCMAM(n)%pom( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'POM', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:6) == 'SOA') then 
              GCMAM(n)%lsoa = .true.
              ALLOCATE( GCMAM(n)%soa( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'SOA', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:7) == 'SSLT') then 
              GCMAM(n)%lsslt = .true.
              ALLOCATE( GCMAM(n)%sslt( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'SSLT', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:7) == 'DUST') then 
              GCMAM(n)%ldust = .true.
              ALLOCATE( GCMAM(n)%dust( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'DUST', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:6) == 'NH4') then 
              GCMAM(n)%lnh4 = .true.
              ALLOCATE( GCMAM(n)%nh4( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'NH4', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:6) == 'NO3') then 
              GCMAM(n)%lno3 = .true.
              ALLOCATE( GCMAM(n)%no3( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'NO3', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:5) == 'CA') then 
              GCMAM(n)%lca = .true.
              ALLOCATE( GCMAM(n)%ca( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'CA', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:6) == 'CO3') then 
              GCMAM(n)%lco3 = .true.
              ALLOCATE( GCMAM(n)%co3( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'CO3', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:5) == 'CL') then 
              GCMAM(n)%lcl = .true.
              ALLOCATE( GCMAM(n)%cl( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'CL', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:6) == 'MOM') then 
              GCMAM(n)%lmom = .true.
              ALLOCATE( GCMAM(n)%mom( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'MOM', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
          if (SpcLocData(s)%info%name(4:5) == 'Nu') then 
              ALLOCATE( GCMAM(n)%Nu( NX, NY, NZ ), STAT=RC )
              CALL GC_CheckVar( 'Nu', 0, RC )
              IF ( RC /= GC_SUCCESS ) THEN
                errMsg = 'Error allocating array MAM !'
                CALL GC_Error( errMsg, RC, thisLoc )
                RETURN
              ENDIF
          end if
! UPDATE WHEN NEW SPECIES WILL BE INTRODUCED
      end if ! its a mam species mode n   
    end do ! loop species 

! initialize ratios by considering the default MAM modal diameters and standard dev
!
    if (n==1)   GCMAM(n)%vol2num =  3.0312191595848506E+020_fp
    if (n==2)   GCMAM(n)%vol2num =  4.0212792866476384E+022_fp
    if (n==3)   GCMAM(n)%vol2num =  50431908767592952_fp
    if (n==4)   GCMAM(n)%vol2num =  5.6542403793695120E+021
    ! update and test for future cases where mode > 4 ( MAM7 etc)    
    end do
 END SUBROUTINE Init_MAM_Container

 !DESCRIPTION: Subroutine CLEANUP\_AER\_CONTAINER deallocates all fields
!  of the aer container object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Mam_Container( GCMAM, RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MAMContainer), POINTER :: GCMAM(:)   ! MAM data container
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: RC    ! Return code
    INTEGER                         :: n
!
! !REVISION HISTORY:
!  28 Mar 2023 - E. Lundgren- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! Cleanup_AerMass_Container starts here
    !======================================================================

    ! Assume success
    RC = GC_SUCCESS

    
    DO n= 1,4
    ! Deallocate arrays and nullify pointer
    IF ( ASSOCIATED( GCMAM(n)%nuwetrad ) ) THEN
       DEALLOCATE( GCMAM(n)%nuwetRAD, STAT=RC )
       CALL GC_CheckVar( 'MAM%nuwetrad', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%nuwetrad => NULL()
    ENDIF

    IF ( ASSOCIATED( GCMAM(n)%nudryrad ) ) THEN
       DEALLOCATE( GCMAM(n)%nudryrad, STAT=RC )
       CALL GC_CheckVar( 'MAM%nudryrad', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%nudryrad => NULL()
    ENDIF

    IF ( ASSOCIATED( GCMAM(n)%wetrad ) ) THEN
       DEALLOCATE( GCMAM(n)%wetRAD, STAT=RC )
       CALL GC_CheckVar( 'MAM%wetrad', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%wetrad => NULL()
    ENDIF

    IF ( ASSOCIATED( GCMAM(n)%dryrad ) ) THEN
       DEALLOCATE( GCMAM(n)%dryrad, STAT=RC )
       CALL GC_CheckVar( 'MAM%dryrad', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%dryrad => NULL()
    ENDIF

     IF ( ASSOCIATED(GCMAM(n)%aerdens ) ) THEN
       DEALLOCATE( GCMAM(n)%aerdens, STAT=RC )
       CALL GC_CheckVar( 'MAM%aerdens', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%aerdens => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%hygro ) ) THEN
       DEALLOCATE( GCMAM(n)%hygro, STAT=RC )
       CALL GC_CheckVar( 'MAM%hygro', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%hygro => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%pH ) ) THEN
       DEALLOCATE( GCMAM(n)%pH, STAT=RC )
       CALL GC_CheckVar( 'MAM%pH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%pH => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%aerwat ) ) THEN
       DEALLOCATE( GCMAM(n)%aerwat, STAT=RC )
       CALL GC_CheckVar( 'MAM%aerwat', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%aerwat => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%so4 ) ) THEN
       DEALLOCATE( GCMAM(n)%so4, STAT=RC )
       CALL GC_CheckVar( 'MAM%so4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%so4 => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%bc ) ) THEN
       DEALLOCATE( GCMAM(n)%bc, STAT=RC )
       CALL GC_CheckVar( 'MAM%bc', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%bc => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%pom ) ) THEN
       DEALLOCATE( GCMAM(n)%pom, STAT=RC )
       CALL GC_CheckVar( 'MAM%pom', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%pom => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%soa ) ) THEN
       DEALLOCATE( GCMAM(n)%soa, STAT=RC )
       CALL GC_CheckVar( 'MAM%soa', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%soa => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%sslt ) ) THEN
       DEALLOCATE( GCMAM(n)%sslt, STAT=RC )
       CALL GC_CheckVar( 'MAM%sslt', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%sslt => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%dust ) ) THEN
       DEALLOCATE( GCMAM(n)%dust, STAT=RC )
       CALL GC_CheckVar( 'MAM%dust', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%dust => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%nh4 ) ) THEN
       DEALLOCATE( GCMAM(n)%nh4, STAT=RC )
       CALL GC_CheckVar( 'MAM%nh4', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%nh4 => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%no3 ) ) THEN
       DEALLOCATE( GCMAM(n)%no3, STAT=RC )
       CALL GC_CheckVar( 'MAM%no3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%no3 => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%ca ) ) THEN
       DEALLOCATE( GCMAM(n)%ca, STAT=RC )
       CALL GC_CheckVar( 'MAM%ca', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%ca => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%co3 ) ) THEN
       DEALLOCATE( GCMAM(n)%co3, STAT=RC )
       CALL GC_CheckVar( 'MAM%co3', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%co3 => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%cl ) ) THEN
       DEALLOCATE( GCMAM(n)%cl, STAT=RC )
       CALL GC_CheckVar( 'MAM%cl', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%cl => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%mom ) ) THEN
       DEALLOCATE( GCMAM(n)%mom, STAT=RC )
       CALL GC_CheckVar( 'MAM%mom', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%mom => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%nu ) ) THEN
       DEALLOCATE( GCMAM(n)%nu, STAT=RC )
       CALL GC_CheckVar( 'MAM%nu', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%nu => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%tauxar ) ) THEN
       DEALLOCATE( GCMAM(n)%tauxar, STAT=RC )
       CALL GC_CheckVar( 'MAM%tauxar', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%tauxar => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%ssa ) ) THEN
       DEALLOCATE( GCMAM(n)%ssa, STAT=RC )
       CALL GC_CheckVar( 'MAM%ssa', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%ssa => NULL()
     ENDIF

     IF ( ASSOCIATED(GCMAM(n)%g ) ) THEN
       DEALLOCATE( GCMAM(n)%g, STAT=RC )
       CALL GC_CheckVar( 'MAM%g', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%g => NULL()
     ENDIF

    END DO 
   
  END SUBROUTINE Cleanup_MAM_Container
!EOC

END MODULE Mam_container_Mod

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
 
     REAL(fp), POINTER :: so4(:,:,:) ! mam so4 mass concentration 
     REAL(fp), POINTER :: bc(:,:,:) ! mam  mass concentration 
     REAL(fp), POINTER :: pom(:,:,:) ! mam  mass concentration 
     REAL(fp), POINTER :: soa(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: sslt(:,:,:) ! mam  mass concentration
     REAL(fp), POINTER :: dust(:,:,:) ! mam  mass concentration
   
     REAL(fp), POINTER :: nu(:,:,:) ! mam number concentration 

     LOGICAL           :: lso4, lbc, lpom, lsoa, lsslt,ldust


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


! now define the mass concentration per mode. They will be used for diagnostics and for 
! input to e.g. optical depth, Fast-J etc . Note that for each MAM mode, 
! not all species are relevant and corresponding mass array are not allocated to save mem !
! Always think making an if allocated test when using GCMAM(n)%spec(:,:,:) elsewhere in the code.
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

    print*, 'FAB mam container lso4' , n, GCMAM(n)%lso4

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
    END DO 
   
  END SUBROUTINE Cleanup_MAM_Container
!EOC

END MODULE Mam_container_Mod



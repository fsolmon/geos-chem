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
     REAL(fp), POINTER :: dryrad (:,:,:) ! dry radius
     REAL(fp), POINTER :: wetrad (:,:,:) ! wet radius 
     REAL(fp), POINTER :: aerdens(:,:,:) ! aerosol effective density  
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
 SUBROUTINE Init_MAM_Container( Input_Opt, State_Grid,GCMAM , RC )
!
! !USES:
!
    USE CMN_Size_Mod,   ONLY : NAER
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    !
! !INPUT PARAMETERS:
!
    TYPE(OptInput),      INTENT(IN)  :: Input_Opt  ! Input Options object
    TYPE(GrdState),      INTENT(IN)  :: State_Grid ! Grid object

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
    INTEGER            :: NX, NY, NZ, n

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
    
    ! dry radius
    do n= 1, 4 
    ALLOCATE(GCMAM(n)%dryrad( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'WETRAD', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array DRYRAD!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%DRYRAD = 0.0_fp
    ! wet radius
    ALLOCATE( GCMAM(n)%wetrad( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'WETRAD', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array WETRAD!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%WETRAD = 0.0_fp

    ! wet radius
    ALLOCATE( GCMAM(n)%aerdens( NX, NY, NZ ), STAT=RC )
    CALL GC_CheckVar( 'AERDENS', 0, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error allocating array AERDENS!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF
    GCMAM(n)%AERDENS = 0.0_fp
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
       CALL GC_CheckVar( 'MAM%wetrad', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GCMAM(n)%aerdens => NULL()
    ENDIF
    END DO 
   
  END SUBROUTINE Cleanup_MAM_Container
!EOC

END MODULE Mam_container_Mod



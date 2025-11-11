! constituents.F90
!    This F90 module file is a special version of the equivalent ACME (and CAM5) module.
!    It provides the functionality needed by the cambox offline code
!    that is used for development and testing of the modal aerosol module (MAM),
!    but (in most cases) not all the functionality of the equivalent ACME module.
!    Also, it may have been taken from a version of CAM5 that was older
!    than ACME-V0 (i.e., pre 2014).

      module constituents

!      use mam_utils, only:  endrun, iulog

      implicit none

      public
!FAB integer, parameter :: pcnst = PCNST 
!Initially PCNST is provided as a precompiling key. 
!Here add precompiling hard-coded options for limiting errorsin dimensioning 
!mam arrays (no dynamical allocations) 
!Now only retains the options considered for Geos-Chem.
!A consistency test has also been added in modal_aero_initialize (inspired from box model).  

#if ( ( defined MODAL_AERO_4MODE_MOM ) && ( defined RAIN_EVAP_TO_COARSE_AERO ) && ( defined MOSAIC_SPECIES ) )
      integer, parameter :: pcnst = 54 
#elif ( defined MODAL_AERO_4MODE )
      integer, parameter :: pcnst = 28
#endif
      character(len=16) :: cnst_name(pcnst)     ! constituent names
      integer           :: species_class(pcnst)
      contains

!==============================================================================

  subroutine cnst_get_ind (name, ind, abort)
!----------------------------------------------------------------------- 
! 
! Purpose: Get the index of a constituent 
! 
! Author:  B.A. Boville
! 
!-----------------------------Arguments---------------------------------
!
    character(len=*),  intent(in)  :: name  ! constituent name
    integer,           intent(out) :: ind   ! global constituent index (in q array)
    logical, optional, intent(in)  :: abort ! optional flag controlling abort

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
    logical :: abort_on_error
!-----------------------------------------------------------------------

! Find tracer name in list
    do m = 1, pcnst
       if (name == cnst_name(m)) then
          ind  = m
          return
       end if
    end do

! Unrecognized name
    abort_on_error = .true.
    if ( present(abort) ) abort_on_error = abort

    if ( abort_on_error ) then
       write(*,*) 'CNST_GET_IND, name:', name, &
                      ' not found in list:', cnst_name(:)
    end if

! error return
    ind = -1

  end subroutine cnst_get_ind

!==============================================================================================

      end module constituents

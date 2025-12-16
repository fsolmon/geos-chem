MODULE MAM_UTILS

  use precision_mod, only: r8 => f8
  use chem_mods, only:  gas_pcnst

  implicit none
  !private
  public
  save
  integer,parameter :: i4 = selected_int_kind ( 6) ! 4 byte integer MOVE that in precision later
  integer :: iulog = 51 ! comme le pastis 
  integer, parameter :: fieldname_len = 64 
  integer, parameter :: phys_decomp = 1   ! *** this is a kludge to avoid the pio file

  integer, public :: ncol_for_outfld

  character(len=10), public :: horiz_only = "horiz_only"

  integer :: ptimelevels = 2 

  character(len=16) :: solsym(gas_pcnst) = '????????????????'



  integer, parameter :: plat = 1
  integer, parameter :: plon = 1
 

  integer :: mdo_mambox, mdo_gaschem, mdo_cloudchem, mdo_coldstart
  integer :: mdo_gasaerexch, mdo_rename, mdo_newnuc, mdo_coag

  integer :: pcols 
  integer :: pver 
  integer, parameter :: psubcols = 1
  integer  :: plev 

  integer, parameter :: begchunk = 1
  integer, parameter :: endchunk = 2


  integer, parameter :: trop_cloud_top_lev = 1
  integer, parameter :: clim_modal_aero_top_lev = 1
  
  logical, public    :: masterproc 
  integer, parameter :: iam = 0
  
  logical, public :: is_first_step_save = .true.

  integer, public :: l_h2so4g, l_soag, l_hno3g, l_so2g, l_hclg, l_nh3g

  integer, parameter, public :: n_ocean_data = 4
#if (defined MODAL_AERO_9MODE || defined MODAL_AERO_4MODE_MOM)
  logical, parameter, public :: has_mam_mom = .true.
#else
  logical, parameter, public :: has_mam_mom = .false.
#endif


#if ( ( defined __PGI ) || ( defined CPRPGI ) )
  ! quiet nan for portland group compilers
  real(r8), parameter :: inf = O'0777600000000000000000'
  real(r8), parameter :: nan = O'0777700000000000000000'

#elif ( defined CPRNAG )
  ! *** THIS NEEDS MORE WORK                                                ***
  ! *** the 2 lines in the "else" block produce errors with nagfor compiler ***
  ! *** the 1.0e300 lines here are a temporary bandaid                      ***
  ! *** NOTE ALSO that infnan is only used in the amicphys-mosaic interface ***
  ! ***      and mosaic is not incorporated in E3SM or cmdv-se verification ***
  real(r8), parameter :: inf =  1.0e300_r8
  real(r8), parameter :: nan = -1.0e300_r8

#else
  ! signaling nan otherwise
  real(r8), parameter :: inf = O'0777600000000000000000'
  real(r8), parameter :: nan = O'0777610000000000000000'
#endif

  integer,  parameter :: bigint = O'17777777777'           ! largest possible 32-bit integer
  real(r8), parameter :: uninit_r8 = inf                   ! uninitialized floating point number

  !PUBLIC MEMBER FUNCTIONS:

  public :: get_spc_ndx

  
  
CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine endrun( msg )
    character(len=*), optional, intent(in) :: msg
    integer :: lunout
    lunout = 6
    if ( present( msg ) ) then
       if (msg /= ' ') then
          write(lunout,'(/a/a/)') &
               '*** stopping in ENDRUN with msg =', trim(msg)
       else
          write(lunout,'(/a/)') &
               '*** stopping in ENDRUN with msg = blank'
       end if
    else
       write(lunout,'(/a/)') &
            '*** stopping in ENDRUN with no msg'
    end if
    stop
  end subroutine endrun



  function aerodep_flx_prescribed()
    logical :: aerodep_flx_prescribed
    aerodep_flx_prescribed = .false.
    return
  end function aerodep_flx_prescribed



  subroutine addfldv1 (fname, units, numlev, avgflag, long_name, &
       decomp_type, flag_xyfill,sampling_seq, &
       begdim1, enddim1, &
       begdim2, enddim2,&
       begdim3, enddim3, mdimnames, fill_value )

    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Add a field to the master field list
    ! 
    ! Method: Put input arguments of field name, units, number of levels, averaging flag, and 
    !         long name into a type entry in the global master field list (masterlist).
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------

    !
    ! Arguments
    !
    character(len=*), intent(in) :: fname      ! field name--should be "max_fieldname_len" characters long
    ! or less
    character(len=*), intent(in) :: units      ! units of fname--should be 8 chars
    character(len=1), intent(in) :: avgflag    ! averaging flag
    character(len=*), intent(in) :: long_name  ! long name of field

    integer, intent(in) :: numlev              ! number of vertical levels (dimension and loop)
    integer, intent(in) :: decomp_type         ! decomposition type

    logical, intent(in), optional :: flag_xyfill ! non-applicable xy points flagged with fillvalue

    character(len=*), intent(in), optional :: sampling_seq ! sampling sequence - if not every timestep, 
    ! how often field is sampled:  
    ! every other; only during LW/SW radiation calcs, etc.
    integer, intent(in), optional :: begdim1, enddim1
    integer, intent(in), optional :: begdim2, enddim2
    integer, intent(in), optional :: begdim3, enddim3

    character(len=*), intent(in), optional :: mdimnames(:)
    real(r8), intent(in), optional :: fill_value

    character(len=16) :: txtaa
    txtaa = fname
    write(95,'(2a,i6)') 'addfld - ', txtaa, numlev

    return
  end subroutine addfldv1

  !#######################################################################

  subroutine addfld (fname, numlev_txt, avgflag, units, long_name, &
       flag_xyfill, sampling_seq, &
       begdim1, enddim1, &
       begdim2, enddim2,&
       begdim3, enddim3, mdimnames, fill_value )

    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Add a field to the master field list
    ! 
    ! Method: Put input arguments of field name, units, number of levels, averaging flag, and 
    !         long name into a type entry in the global master field list (masterlist).
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------

    !      use ppgrid, only:  pver
    !
    ! Arguments
    !
    character(len=*), intent(in) :: fname      ! field name--should be "max_fieldname_len" characters long
    ! or less
    character(len=*), intent(in) :: units      ! units of fname--should be 8 chars
    character(len=1), intent(in) :: avgflag    ! averaging flag
    character(len=*), intent(in) :: long_name  ! long name of field

    character(len=*), intent(in) :: numlev_txt(*)              ! number of vertical levels (dimension and loop)

    logical, intent(in), optional :: flag_xyfill ! non-applicable xy points flagged with fillvalue

    character(len=*), intent(in), optional :: sampling_seq ! sampling sequence - if not every timestep, 
    ! how often field is sampled:  
    ! every other; only during LW/SW radiation calcs, etc.
    integer, intent(in), optional :: begdim1, enddim1
    integer, intent(in), optional :: begdim2, enddim2
    integer, intent(in), optional :: begdim3, enddim3

    character(len=*), intent(in), optional :: mdimnames(:)
    real(r8), intent(in), optional :: fill_value

    character(len=16) :: txtaa
    txtaa = fname
    if (numlev_txt(1) == 'horiz_only') then
       write(95,'(2a,i6)') 'addfld - ', txtaa, 1
    else
       write(95,'(2a,i6)') 'addfld - ', txtaa, pver
    end if

    return
  end subroutine addfld

  !#######################################################################

  subroutine add_default (name, tindex, flag)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Add a field to the default "on" list for a given history file
    ! 
    ! Method: 
    ! 
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    character(len=*), intent(in) :: name  ! field name
    character(len=1), intent(in) :: flag  ! averaging flag

    integer, intent(in) :: tindex         ! history tape index

    character(len=16) :: txtaa
    txtaa = name
    write(95,'(2a,i6)') 'adddef - ', txtaa, tindex

    return
  end subroutine add_default

  !#######################################################################

  subroutine outfld (fname, field, idim, c)
    !
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Accumulate (or take min, max, etc. as appropriate) input field
    !          into its history buffer for appropriate tapes
    ! 
    ! Method: Check 'masterlist' whether the requested field 'fname' is active
    !         on one or more history tapes, and if so do the accumulation.
    !         If not found, return silently.
    ! 
    ! Author: CCM Core Group
    ! 
    !-----------------------------------------------------------------------
    !
    !      use ppgrid, only: pcols, pver

    ! Arguments
    !
    character(len=*), intent(in) :: fname ! Field name--should be 8 chars long

    integer, intent(in) :: idim           ! Longitude dimension of field array
    integer, intent(in) :: c              ! chunk (physics) or latitude (dynamics) index

    real(r8), intent(in) :: field(idim,*) ! Array containing field values

    integer :: k
    character(len=16) :: txtaa

    txtaa = fname
    write(90,'(/a,1p,20e11.3)') txtaa, field(1:ncol_for_outfld,1)
!FAB    print*, 'ahha' , ncol_for_outfld, pver,idim  
!    if ( txtaa == 'SOAG_sfgaex3d   ' .or. &
!         txtaa == 'num_a2_nuc1     ' .or. &
!         txtaa == 'num_a2_nuc2     ' ) then
!       do k = 2, pver
!          write(90,'(a,i2,10x,1p,20e11.3)') '  k=', k, field(1:ncol_for_outfld,k)
 !      end do
 !   end if


  end subroutine outfld




  integer function get_spc_ndx( spc_name )
    !-----------------------------------------------------------------------
    !     ... return overall species index associated with spc_name
    !-----------------------------------------------------------------------

    use chem_mods,     only : gas_pcnst
    !    use mo_tracname,   only : tracnam => solsym


    !-----------------------------------------------------------------------
    !     ... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: spc_name

    !-----------------------------------------------------------------------
    !     ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    get_spc_ndx = -1
    do m = 1,gas_pcnst
       if( trim( spc_name ) == trim( solsym(m) ) ) then
          get_spc_ndx = m
          exit
       end if
    end do

  end function get_spc_ndx


  subroutine seasalt_init()
    return
  end subroutine seasalt_init



  !    than ACME-V0 (i.e., pre 2014).

  !!moduleunits

  function getunit( )
    integer getunit
    getunit = 6
    return
  end function getunit

  !!end module units

   function is_first_step()
     logical :: is_first_step
     is_first_step = is_first_step_save
     return
    end function is_first_step



END MODULE MAM_UTILS


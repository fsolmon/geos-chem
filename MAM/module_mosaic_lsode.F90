module module_mosaic_lsode
  !
  ! double precision lsodes solver
  !
  implicit none
  private

  public :: MOSAIC_LSODE

contains
  subroutine MOSAIC_LSODE(dtchem)
    use precision_mod, only: r8 => f8

    implicit none


    ! subr arguments
    real(r8), intent(in) :: dtchem

    write(*,'(//a//)') &
         '*** error - mosaic_lsode has been deactivated ***'
    stop

    return
  end subroutine MOSAIC_LSODE

end module module_mosaic_lsode

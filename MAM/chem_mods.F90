 module chem_mods
!--------------------------------------------------------------
! ... Basic chemistry parameters and arrays
!--------------------------------------------------------------

      use precision_mod, only : r8 => f8
      use constituents, only:  pcnst
      implicit none
      
      public
      integer, parameter :: imozart = 6
      integer, parameter :: gas_pcnst = pcnst - (imozart - 1)
      real(r8) :: adv_mass(gas_pcnst) = 0._r8

      end module chem_mods

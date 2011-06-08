module types_numeriques

!*************************************************************
!** Variables permettant de definir la precision des reels ***
!** Version 1.0 - juin 2008
!*************************************************************

integer, parameter :: simple_precision = selected_real_kind(6,37)
integer, parameter :: double_precision = selected_real_kind(15,307)

end module types_numeriques

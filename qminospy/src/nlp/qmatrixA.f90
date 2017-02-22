module qmatrixA

  !------------------------------------------------------------------
  ! qmatrixA stores a global copy of the tinyME matrix A,
  ! which is needed by funcon to define the function and Jacobian
  ! for the nonlinear constraints.
  !
  ! 30 Dec 2014: First version of qmatrixA.f90 to go with qsolveME.f90.
  !              Ding Ma and Michael Saunders, SOL, Stanford University.
  !------------------------------------------------------------------

  implicit none

  public

  integer, parameter        :: ip = 4, dp = 8, qp = 16
!  integer, parameter        :: maxnnCon = 40000
!  integer, parameter        :: maxnnJac = 80000
!  integer, parameter        :: maxneJac = 80000

!  integer(ip)               :: haa(maxneJac), kaa(maxnnJac+1)
!  real(qp)                  ::  aa(maxneJac)
  integer(ip),  allocatable :: haa(:), kaa(:)
  real(qp), allocatable     :: aa(:)

end module qmatrixA

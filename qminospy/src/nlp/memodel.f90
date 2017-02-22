!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File memodel.f90
!
! module memodel
!
! Contains data matrices required by qMEfuns.f90. 
! Enables passing all ME model data to qsolveME in RAM, instead of flatfiles
! Thus, more or less moved all allocatable variables needed by minoss from 
! qMEfuns.f90 to here.
! NOTE: f2py has basic support for module allocatable arrays
!
! Usage example from python (after compiling using f2py):
!   import qsolveME
!   qsolveME
!   qsolveME.memodel.hs = hs
!   qsolveME.memodel.ha = ha
!   qsolveME.memodel.ka = ka
!   qsolveME.memodel.a = a
!   qsolveME.memodel.bl = bl
!   qsolveME.memodel.bu = bu
!   qsolveME.memodel.cObj = cObj
!   qsolveME.memodel.xn = xn
!   qsolveME.memodel.pi = pi
!   qsolveME.memodel.rc = rc
!   qsolveME.memodel.name1= = name1
!   qsolveME.memodel.name2 = name2
!
!   qsolveME.qsolveme(nb, mu0, m, n, ne, nnCon, nnJac, neJac,
!       hs, ha, ka, a, bl, bu, cObj, xn, pi, rc, name1, name2)
!
! 
! Laurence Yang, SBRG, UCSD
! 
! 22 Jul 2015: [LY] first version to be used with qMEfuns.f90 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module memodel

  implicit none

  public

  integer, parameter        :: ip = 4, dp = 8, qp = 16

  ! integer(ip),  allocatable, dimension(:) :: ha, ka !ha(:), ka(:)
  ! integer,  allocatable, dimension(:) :: ha, ka !ha(:), ka(:)
  ! real(qp),     allocatable, dimension(:) :: a !a(:)
  ! real(qp),     allocatable, dimension(:) :: bl, bu
  ! User interface provides double
  ! real(dp),     allocatable, dimension(:) :: ad
  ! real(dp),     allocatable, dimension(:) :: bld, bud
  ! real,     allocatable, dimension(:) :: ad
  real,     allocatable, dimension(:) :: bld, bud

end module memodel

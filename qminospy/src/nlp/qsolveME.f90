!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File qMEfuns.f90: modified from qsolveME.f90 by Ding Ma and Michael Saunders
!
! subroutine qsolveME (changed from program, for f2py)
!     contains
!         subroutine qdq
!         subroutine qd
! subroutine funcon
! subroutine funobj
! subroutine matmod (dummy)
!
! Modified code from Ding Ma and Michael Saunders for compatibility with f2py
! and calling from python
! 
! Laurence Yang, SBRG, UCSD
!
! How to generate python extension module via f2py, and use from python:
!   0) [prerequisite] make qminos using -fPIC compiler flag.
!      Otherwise, will get error at step 2) about requiring the -fPIC flag

!   1) Generate signature file, qsolveME.pyf:
!      f2py -m qsolveME -h qsolveME.pyf qMEfuns.f90 qmatrixA.f90
!
!   2) Compile with qminos library linked, and generate .so file to import into python:
!      f2py -c qsolveME.pyf qMEfuns.f90 qmatrixA.f90 -L/home/laurence/Software/qminos1114/qminos56/lib -lquadminos
!
!   3) From python: 
!      import qsolveME
!      import numpy as np
!      inform = np.array(0)     # gets modified by qsolveme: in/output exit flag
!      mu0 = 0.1
!      x = qsolveME.qsolveme(inform, mu0, probname, M, nnCon, nnJac, neJac, ha,
!      ka, ad, bld, bud, nb, N, ne)
!
! 19 Jul 2015: [LY] first version. 
! 20 Jul 2015: [LY] reproduced fortran solution from python--still fort.8 input
! file, though
! 30 Jul 2015: [LY] return inform (exit status) as intent(inout)
! 14 Aug 2015: [LY] qsolveME is now warm-startable with intent(inout) basis
!                   if warm=True
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! 20 Jul 2015: [LY] modified such that solution is returned to python as double,
!                   rounded from quad.
! 22 Jul 2015: [LY] mu0 is an input argument
!                   Also specified intent(in) and intent(out) for nb and xnd,
!                   so the .pyf file from py2f need not be manually edited
subroutine qsolveME(xnd, inform, nb, mu0d, Probname, m, n, ne, nnCon, nnJac, neJac, &
    ha, ka, ad, bld, bud, hs, warm, &
    nStrOpts, nIntOpts, nRealOpts, stropts, intopts, realopts, intvals, realvalsd)
  ! From qmatrixA, use only: ...
  ! use qmatrixA,  only: maxnnCon, maxnnJac, maxneJac, haa, kaa, aa
  ! 24 Mar 2016: allocatable haa, kaa, aa
  use qmatrixA,  only: haa, kaa, aa

  implicit none

  integer, parameter        :: ip = 4, dp = 8, qp = 16
  ! Pass options directly instead of using specs file
  integer(ip), intent(in)                               :: nStrOpts
  integer(ip), intent(in)                               :: nIntOpts
  integer(ip), intent(in)                               :: nRealOpts
  character(len=72), dimension(nStrOpts), intent(in)    :: stropts
  character(len=55), dimension(nIntOpts), intent(in)    :: intopts
  character(len=55), dimension(nRealOpts),intent(in)    :: realopts
  integer(ip), dimension(nIntOpts), intent(in)          :: intvals
  ! Real options input as real(8) but converted to quad when passing to minos
  real(dp), dimension(nRealOpts), intent(in)            :: realvalsd
  real(qp), dimension(nRealOpts)                        :: realvals

  ! qMINOS common block
  real(16)        :: dparm
  integer         :: iparm
  common    /m2parm/ dparm(30),iparm(30)
  ! 22 Jul 2015: [LY] ensure these params in the right order...
  !              provide these params as arguments instead of reading from file
  integer(ip),  intent(in)  :: nb
  real(dp),     intent(in)  :: mu0d
  character(8), intent(in)  :: Probname
  integer(ip),  intent(in)  :: m, n, ne, nnCon, nnJac, neJac
  integer(ip), dimension(ne),  intent(in)  :: ha
  integer(ip), dimension(n+1),  intent(in)  :: ka
  real(dp), dimension(ne), intent(in)    :: ad
  real(dp), dimension(nb), intent(in)    :: bld, bud
  logical, intent(in)       :: warm

  ! Local variables for MINOS
  integer(ip)               :: INFO, iExit
  !  integer(ip)               :: iSpecs, nnCon, nnJac, neJac, nnObj, nwcore
  integer(ip)               :: iSpecs, nnObj, nwcore
  ! integer(ip)               :: iObj, iPrint, iSumm, inform, nout
  ! 30 Jul 2015: [LY] want to return inform (exit status) to user
  integer(ip)               :: iObj, iPrint, iSumm, nout
  integer(ip),intent(inout) :: inform
  ! integer(ip)               :: m, n, nb, ne, mincor, nInf, nname, nS
  integer(ip)               :: mincor, nInf, nname, nS
  real(qp)                  :: ObjAdd, Obj, sInf
  ! LY: allocatable: at run time, allocate() used to assign memory for storage during execution
  !integer(ip),  allocatable :: hs(:)
  integer(ip),intent(inout) :: hs(nb)
  ! integer(ip),  allocatable :: ha(:), ka(:)

  ! real(qp),     allocatable :: a(:)
  real(qp)                  :: a(ne)
  ! 22 Jul 2015: [LY] no need to allocate dynamically
  ! real(dp),     intent(in)  :: bld(nb), bud(nb)
  real(qp)                  :: bl(nb), bu(nb)
  !real(qp),     allocatable :: xn(:), pi(:), rc(:)
  real(qp),     allocatable :: pi(:), rc(:)

  ! LY: make double version of xn, which will be intent(output) to python
  real(dp),   intent(inout) :: xnd(nb)      ! double-rounded version of xn, which is quad 
  real(qp)                  :: xn(nb)       ! quad  

  integer(ip),  allocatable :: name1(:), name2(:)
  !character(8)              :: names(5), Probname
  character(8)              :: names(5)
  ! LY: as argument


  integer(ip),  parameter   :: lenz = 100000000 !5000000   ! As big as you like.
  real(qp)                  :: z(lenz)          ! This is the MINOS workspace.

  ! Local variables for the problem
  integer(ip)               :: in, i, j, k, l, nnCon0, nnObj0, nnJac0
  real(qp)                  :: mu0
  real(qp),     parameter   :: zero = 0.0_qp
  !------------------------------------------------------------------
  ! LY: user inputs mu0d as double, convert to quad for qMINOS
  mu0 = real(mu0d, qp)

  ! LY: instead of reading a file, pass model data directly from python
  ! Assign nnCon, nnJac, etc. directly from python
  ! Or make a internal procedure to interface the required variable 
  ! assignments clearly to python?

   in     = 8         ! File from ~/matlab/tinyME/dumpME_NLP.m (e.g. tinyME.txt)
                      ! cp tinyME.txt fort.10
   ! open(in, status='old')
   ! read(in,*) Probname! Problem name, up to 8 characters
   ! read(in,*) nnCon   ! No of rows in Jacobian, the matrix [Abar*xbar mu*Abar]
   ! read(in,*) nnJac   ! No of cols in Jacobian
   ! read(in,*) neJac   ! No of nonzeros in Jacobian
   ! read(in,*) m       ! No of rows in A
   ! read(in,*) n       ! No of cols in A
   ! read(in,*) ne      ! No of nonzeros in A

!  if      (nnCon > maxnnCon) then
!     write(*,*) 'maxnnCon in qmatrixA  is smaller than', nnCon
!     stop
!  else if (nnJac > maxnnJac) then
!     write(*,*) 'maxnnJac in qmatrixA  is smaller than', nnJac
!     stop
!  else if (neJac > maxneJac) then
!     write(*,*) 'maxneJac in qmatrixA  is smaller than', neJac
!     stop
!  end if

  ! 24 Mar 2016 [LY]
  ! Allocate A
  write(*,*) 'Allocating haA (neJac):', neJac
  write(*,*) 'Allocating kaA (nnJac+1):', nnJac+1
  write(*,*) 'Allocating aA (neJac): ', neJac
  allocate( haA(neJac), kaA(nnJac+1))
  allocate( aA(neJac))

  ! Allocate A as pointers, row indices, and Aij
  ! allocate( ka(n+1), a(ne), ha(ne) )
  ! read(in,*) ka      ! Pointers for A
  ! read(in,*) ha      ! Row indices
  ! read(in,*) a       ! Aij
  ! 22 Jul 2015: [LY] allocate quad Aij
  ! allocate(a(ne))
  a = real(ad, qp)

  ! nb       = n + m
  nname    = 1

  print*, 'Allocating name1, name2'
  allocate( name1(nname), name2(nname) )

!  print*, 'Allocating pi, rc, xn'
!  allocate( pi(m), rc(nb), xn(nb) )
  print*, 'Allocating pi, rc '
  allocate( pi(m), rc(nb) )

  print *, 'Allocated pi, rc, xn'

  ! 22 Jul 2015: [LY] allocate quad bl, bu. Populate from user-provided dp
  ! allocate( bl(nb), bu(nb) )
  bl = real(bld, qp)
  bu = real(bud, qp)

  ! read(in,*) bl        ! Lower bounds for x and slacks
  ! read(in,*) bu        ! Upper bounds for x and slacks
  ! close(in)

  ! The data (a,bl,bu) just read into Quad arrays was really just Double.
  ! For example, a(1) may now be [a11 a12] in Quad.
  ! Round from Quad to Double [a11] and then back to Quad as [a11 0]

  call qdq(  a, ne )
  call qdq( bl, nb )
  call qdq( bu, nb )

  ! Save the Jacobian part of A in the arrays of module qmatrixA.
  ! (This is the f90 alternative to a common block.)
  ! LY: comments from qminos56/src/mi02lib.f
  ! a(ne)   constraint matrix (Jacobian), stored column-wise
  ! ha(ne)  list of row indices for each nonzero in a(*)
  ! ka(n+1) set of pointers to beginning of each column of the constraint matrix
  !         within a(*) and ha(*). Must have ka(1) = 1 and ka(n+1) = ne+1.

  l = 0
  do j = 1, nnJac
     kaA(j) = l+1
     do k = ka(j), ka(j+1)-1
        i = ha(k)
        if (i > nnCon) exit
        l = l+1
        haA(l) = i
         aA(l) = a(k)
     end do
  end do

  kaA(nnJac+1) = l+1   ! Last pointer

  ! Assign various names.
  ! These are relics from the days of MPS files.
  ! They appear in the MINOS Print file and/or Solution file.

  names(1) = Probname
  names(2) = 'c       '  ! Objective name
  names(3) = 'b       '  ! RHS name
  names(4) = '        '  ! Ranges name (bounds on slacks)
  names(5) = 'bounds  '  ! Bounds name (bounds on variables)

  ! Specify file numbers for MINOS.  (Others may be in the SPECS file.)
  ! 0 means that there should be no file.

!  ispecs = 4   ! The MINOS SPECS   file.
  ispecs = 0   ! The MINOS SPECS   file.
! ispecs2= 5   ! The second SPECS  file (for warm start)
  iprint = 9   ! The MINOS PRINT   file.
  isumm  = 6   ! The MINOS SUMMARY file.
  nout   = 6   ! Local output file (6 = screen).

  ! Now we may open any number of files ourselves
  ! (perhaps to give them sensible names).
  ! For example:

  ! open( iprint, file='TMA_ME.out', status='UNKNOWN')
  ! open( ispecs, file='TMA_ME.spc', status='OLD')

  ! Alternatively, we may let mistart and minoss open them,
  ! using the method selected in subroutine m1open
  ! in the mi10*.f file that was used to build MINOS.
  !
  ! RULE OF THUMB:
  ! MINOS won't open file ispec (for example)
  ! if that unit is already open, or if ispec = 0.
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! mistart MUST BE CALLED BEFORE ANY OTHER MINOS ROUTINE.
  !------------------------------------------------------------------
  call mistart( iprint, isumm, ispecs )  ! Initialize MINOS and open
                                         ! the specified files.

  ! call mispec( ispecs, inform )          ! Read the SPECS file
                                         ! (if ispecs > 0).
  !------------------------------------------------------------------
  ! Specify options directly
  ! Options written as one string
  do i=1, nStrOpts
    write(*,*) 'Calling miopt to set option: ', trim(stropts(i))
    call miopt(trim(stropts(i)), iprint, isumm, inform)
  end do

  ! Integer valued options
  do i=1, nIntOpts
    write(*,*) 'Calling miopti to set option ', trim(intopts(i)), ' to ', intvals(i)
    call miopti(trim(intopts(i)), intvals(i), iprint, isumm, inform)
  end do

  ! Real-valued options
  ! Convert double to quad options
  call dq( realvals, realvalsd, nRealOpts )

  do i=1, nRealOpts
    write(*,*) 'Calling mioptr to set option ', trim(realopts(i)), ' to ', realvalsd(i)
    call mioptr(trim(realopts(i)), realvals(i), iprint, isumm, inform)
  end do

  !------------------------------------------------------------------
  if (inform >= 2) then
     write(nout, *) 'Error: ispecs > 0 but no SPECS file found'
     stop
  end if

  !----------------------------------------------------------------------
  ! Solve the problem.
  ! iobj   = 0    means there is no linear objective row in a(*).
  ! objadd = zero means there is no constant to be added to the objective.
  ! nname  = 1    means there are no meaningful names for the
  !               variables and constraints inside name1(*) and name2(*).
  !               MINOS will print default names.
  !----------------------------------------------------------------------
  ! [LY] 22 Jul 2015: provide mu0 as argument
  ! mu0 = dparm(20)
  write(*, "(/ ' mu0d =', f4.1)") mu0d
  write(*, "(/ ' mu0  =', f4.1)") mu0

  write(iprint, "(/ ' mu0  =', f4.1)") mu0
  write(isumm , "(/ ' mu0  =', f4.1)") mu0

  !hs(1:nb) = 0
  !xn(1:nb) = zero
  !xn(1)    = mu0    ! Initialize mu
  pi(1:m)  = zero
  iobj     = 0
  objadd   = zero
  nwcore   = lenz

  ! dumpLP.m generates slack bounds for SQOPT: Ax - s = 0
  ! Minos has Ax + s = 0, so we have to flip the bounds on s.

  z(1:m)     = -bl(n+1:nb)
  bl(n+1:nb) = -bu(n+1:nb)
  bu(n+1:nb) =   z(1:m)

  ! make quad-precision initial solution based on double
  call dq( xn, xnd, nb )

!!!!!!!!!!!!!!!!!!!! This was trying to solve an LP first to get feasible

!   nnCon0   = 0
!   nnObj0   = 0
!   nnJac0   = 0

!   call minoss( 'Cold', m, n, nb, ne, nname,         &
!                nnCon0, nnObj0, nnJac0,              &
!                iobj, objadd, names,                 &
!                a, ha, ka, bl, bu, name1, name2,     &
!                hs, xn, pi, rc,                      &
!                inform, mincor, ns, ninf, sinf, obj, &
!                z, nwcore )

!   write(nout, *) ' '
!   write(nout, *) 'Quad MINOS finished.'
!   write(nout, *) 'inform =', inform
!   write(nout, *) 'ninf   =', ninf
!   write(nout, *) 'sinf   =', sinf
!   write(nout, *) 'obj    =', obj

!   !======================================================
!   ! Warm start with 2nd SPECS file
!   !======================================================
!   call mispec( ispecs, inform )

!   if (inform .ge. 2) then
!      write(nout, *) 'No second SPECS file'
!      go to 900
!   end if

!   ! 20 Apr 2014:
!   ! Step 2 normally uses scaling (and then unscaling).
!   ! The unscaling might fill in the 2nd half of each quad word.
!   ! Round from Quad to Double and back to Quad again.

!   call qdq(  a, ne )
!   call qdq( bl, nb )
!   call qdq( bu, nb )

  iObj  = 0
  nnObj = 1
  bl(1) = mu0  ! For now, mu0 must be BELOW the OPTIMAL (maximized) value

  if (warm) then
      print *, 'Calling minoss. Warm start with provided basis (hs)'

      call minoss( 'Warm', m, n, nb, ne, nname,        &
          nncon, nnobj, nnjac,                &
          iobj, objadd, names,                &
          a, ha, ka, bl, bu, name1, name2,    &
          hs, xn, pi, rc,                     &
          inform, mincor, ns, ninf, sinf, obj, &
          z, nwcore )

      write(nout, *) ' '
      write(nout, *) 'Quad MINOS finished.'
      write(nout, *) 'inform =', inform
      write(nout, *) 'ninf   =', ninf
      write(nout, *) 'sinf   =', sinf
      write(nout, *) 'obj    =', obj

      ! Return xn rounded to double precision
      call qd( xn, xnd, nb)
  else
      print *, 'Calling minoss. Cold start' 
      print *, 'Setting xn to zero and xn(1) to mu0'
      xn(1:nb) = zero
      xn(1)    = mu0    ! Initialize mu

      call minoss( 'Cold', m, n, nb, ne, nname,        &
          nncon, nnobj, nnjac,                &
          iobj, objadd, names,                &
          a, ha, ka, bl, bu, name1, name2,    &
          hs, xn, pi, rc,                     &
          inform, mincor, ns, ninf, sinf, obj, &
          z, nwcore )

      write(nout, *) ' '
      write(nout, *) 'Quad MINOS finished.'
      write(nout, *) 'inform =', inform
      write(nout, *) 'ninf   =', ninf
      write(nout, *) 'sinf   =', sinf
      write(nout, *) 'obj    =', obj

      ! Return xn rounded to double precision
      call qd( xn, xnd, nb)
  end if

900 close( iprint )
  close( ispecs )

  ! Return xn rounded to double precision
  call qd( xn, xnd, nb)

  !deallocate( pi, rc, xn )
  deallocate( pi, rc )
  deallocate( name1, name2 )
  deallocate( haA, kaA, aA)
  ! 22 Jul 2015: [LY] arrays below are now inputs
  !deallocate( bl, bu, ka, a, ha )

contains

  subroutine qdq( a, n )

    integer(ip), intent(in)    :: n
    real(qp),    intent(inout) :: a(n)

    ! 15 Apr 2014: First version of qdq.
    !              Round from Quad to Double and back to Quad.

    real(dp)                  :: da(n)   ! local array

    write(*,'(z32)') a(n)
    da = real( a,dp)
    a  = real(da,qp)
    write(*,'(z17)') da(n)
    write(*,'(z32)') a(n)
  end subroutine qdq

  ! LY, 20 Jul 2015: First version of qd.
  !                  Round quad vector to double for python
  subroutine qd( aq, ad, n )

      integer(ip),  intent(in)  :: n
      real(qp),     intent(in)  :: aq(n)
      real(dp),     intent(inout) :: ad(n)

      write(*,'(z32)') aq(n)
      ad = real(aq, dp)
      write(*,'(z17)') ad(n)
  end subroutine qd

  ! 15 Aug 2015: double to quad
  subroutine dq( aq, ad, n )
      integer(ip),  intent(in)  :: n
      real(dp),     intent(in)  :: ad(n)
      real(qp),    intent(inout):: aq(n)

      write(*,'(z17)') ad(n)
      aq = real(ad, qp)
      write(*,'(z32)') aq(n)
  end subroutine dq

end subroutine qsolveME


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine funobj( mode, nnObj, x, fobj, gobj, nstate, nprob, z, lenz )

  implicit none
  integer,  intent(inout) :: mode
  integer,  intent(in)    :: nnObj, nstate, nprob, lenz
  real(16), intent(in)    :: x(nnObj)
  real(16), intent(out)   :: fobj, gobj(nnObj)
  real(16), intent(inout) :: z(lenz)

  !------------------------------------------------------------------
  ! This is funobj for tinyME-NLP.
  ! nnObj = 1.  The only nonlinear objective variable is max mu = x(1).
  !------------------------------------------------------------------

  fobj    = x(1)
  gobj(1) = 1.0
    
end subroutine funobj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine funcon( mode, nnCon, nnJac, neJac, x, &
                   f, gcon, nstate, nprob, z, lenz )

  use qmatrixA,  only: haA, kaA, aA

  implicit none
  integer,  intent(inout) :: mode
  integer,  intent(in)    :: nnCon, nnJac, neJac, nstate, nprob, lenz
  real(16), intent(in)    :: x(nnJac)
  real(16), intent(out)   :: f(nnCon), gcon(neJac)
  real(16), intent(inout) :: z(lenz)

  !------------------------------------------------------------------
  ! This is funcon for tinyME-NLP.
  ! Module qmatrixA holds a copy of A in arrays haA, kaA, aA.
  !------------------------------------------------------------------

  integer             :: i, j, k
  real(16)            :: Ai, mu, xj
  real(16), parameter :: zero = 0.0

  mu = x(1)

  if (nstate == 1) then   ! First call of funcon
     write(*,*) 'Calling funcon.  mu =', mu
     write(*,*) 'nnCon, nnJac, neJac', nncon, nnJac, neJac
     write(9,*) 'Calling funcon.  mu =', mu
  end if

  if (nstate >= 2) then   ! Last call of funcon
     write(*,*) 'Final value of   mu =', mu
     write(9,*) 'Final value of   mu =', mu
  end if

  gcon(1:nnCon) = zero    ! Initialize first column of Jacobian

  do j = 2, nnJac         ! f(1:nnCon) and all other columns of Jacobian
     xj = x(j)
   ! write(*,*) "j=", j
     do k = kaA(j), kaA(j+1)-1
        i       = haA(k)
        Ai      =  aA(k)
        gcon(i) = gcon(i) + Ai*xj  ! Accumulate Abar*xbar in col 1
        gcon(k) = mu*Ai
      ! write(*,*) "k=", k, "   i=", i
     end do
  end do

  f = mu*gcon(1:nnCon)
    
end subroutine funcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine matmod
end subroutine matmod

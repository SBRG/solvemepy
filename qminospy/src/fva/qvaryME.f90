!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File qvaryME.f90: adapted from qsolveME to solve FVA in quad-precision with warm-start
!
! subroutine qvaryME : adapted from qsolveME.f90 and qsolveLP.f90 by Ding Ma
!                      and Michael A. Saunders, to interface with Python
!     contains
!         subroutine qdq
!         subroutine qd
!         subroutine dq
!
! Laurence Yang, SBRG, UCSD
!
! How to generate python extension module via f2py, and use from python:
!   0) [prerequisite] make qminos using -fPIC compiler flag.
!      Otherwise, will get error at step 2) about requiring the -fPIC flag

!   1) Generate signature file, qvaryME.pyf:
!      f2py -m qvaryME -h qvaryME.pyf qvaryME.f90
!
!   2) Compile with qminos library linked, and generate .so file to import into python:
!      f2py -c qvaryME.pyf qvaryME.f90 -L/home/laurence/Software/qminos1114/qminos56/lib -lquadminos
!
!   3) From python: 
!      import qvaryME
!      import numpy as np
!      inform = np.array(0)     # gets modified by qsolveme: in/output exit flag
!      mu0 = 0.1
!      x = qvaryME.qvaryme(inform, mu0, probname, M, nncon, nnJac, neJac, ha,
!      ka, ad, bld, bud, nb, N, ne)
!
! 10 Aug 2015: first version. 
! 11 Aug 2015: put warm-start in loop
! 16 May 2016: allow first LP to be warm-started. Provide minos options
!              via function interface instead of fort.14.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


subroutine qvaryME(xnd, informs, nb, Probname, m, n, ne, &
    ha, ka, ad, bld, bud, hs, warm_first, &
    obj_inds, obj_coeffs, obj_vals, nVary, &
    nStrOpts, nIntOpts, nRealOpts, stropts, intopts, realopts, intvals, realvalsd)
  ! all allocatable arrays passed to minoss are stored here 
  ! Need to provide ha, ka, ad for all obj_inds, so that we can
  ! efficiently update Jacobian without having to shift columns each time

  implicit none

  integer, parameter        :: ip = 4, dp = 8, qp = 16

  !----------------------------------------------------------
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
  !----------------------------------------------------------

  ! qMINOS common block
  real(16)        :: dparm
  integer         :: iparm
  common    /m2parm/ dparm(30),iparm(30)
  ! 22 Jul 2015: [LY] ensure these params in the right order...
  !              provide these params as arguments instead of reading from file
  integer(ip),  intent(in)  :: nb
  character(8), intent(in)  :: Probname
  integer(ip),  intent(in)  :: m, n, ne
  integer(ip)               :: nnCon, nnJac, neJac
  integer(ip), dimension(ne),  intent(in)  :: ha
  integer(ip), dimension(n+1),  intent(in)  :: ka
  real(dp), dimension(ne), intent(in)    :: ad
  real(dp), dimension(nb), intent(in)    :: bld, bud
  logical, intent(in)       :: warm_first

  ! 11 Aug 2015: [LY] list of obj inds, obj coeffs, obj vals
  integer(ip),  intent(in)  :: nVary
  integer(ip),  intent(in)  :: obj_inds(nVary)
  real(dp),     intent(in)  :: obj_coeffs(nVary)
  real(qp)                  :: obj_coeffsq(nVary)
  real(dp),  intent(inout)  :: obj_vals(nVary)
  integer(ip),intent(inout) :: informs(nVary)

  ! Local variables for MINOS
  integer(ip)               :: INFO, iExit
  integer(ip)               :: iSpecs, nnObj, nwcore
  ! 30 Jul 2015: [LY] want to return inform (exit status) to user
  integer(ip)               :: iObj, iPrint, iSumm, nout
  integer(ip)               :: inform
  integer(ip)               :: mincor, nInf, nname, nS
  real(qp)                  :: ObjAdd, obj, sInf
  ! LY: allocatable: at run time, allocate() used to assign memory for storage during execution
  ! integer(ip),  allocatable :: hs(:)
  integer(ip),intent(inout) :: hs(nb)
  real(qp)                  :: a(ne)
  ! 22 Jul 2015: [LY] no need to allocate dynamically
  real(qp)                  :: bl(nb), bu(nb)
  real(qp),     allocatable :: cObj(:)
  real(qp),     allocatable :: xn(:), pi(:), rc(:)

  ! LY: make double version of xn, which will be intent(output) to python
  real(dp),     intent(out) :: xnd(nb)         ! double-rounded version of xn, which is quad 

  integer(ip),  allocatable :: name1(:), name2(:)
  character(8)              :: names(5)

  integer(ip),  parameter   :: lenz = 100000000 !5000000   ! As big as you like.
  real(qp)                  :: z(lenz)          ! This is the MINOS workspace.

  ! Local variables for the problem
  real(qp),     parameter   :: zero = 0.0_qp
  integer(ip)               :: i, l, col, ptr, row

  !------------------------------------------------------------------

  a = real(ad, qp)


  print*, 'Allocating name1, name2'
  allocate( name1(nname), name2(nname) )

!  print*, 'Allocating hs, pi, rc, xn'
!  allocate( hs(nb), pi(m), rc(nb), xn(nb) )
  print*, 'Allocating pi, rc, xn'
  allocate( pi(m), rc(nb), xn(nb) )

  print *, 'Allocated hs, pi, rc, xn'

  print*, 'No Jacobian matrix since no nonlinear obj or constraints'
  nnCon     = 0
  nnJac     = 0
  nnObj     = 0

  ! 22 Jul 2015: [LY] allocate quad bl, bu. Populate from user-provided dp
  bl = real(bld, qp)
  bu = real(bud, qp)

  ! The data (a,bl,bu) just read into Quad arrays was really just Double.
  ! For example, a(1) may now be [a11 a12] in Quad.
  ! Round from Quad to Double [a11] and then back to Quad as [a11 0]

  call qdq(  a, ne )
  call qdq( bl, nb )
  call qdq( bu, nb )

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

!  ispecs = 14   ! The MINOS SPECS   file.
  ispecs = 0
! ispecs2= 5   ! The second SPECS  file (for warm start)
  iprint = 9   ! The MINOS PRINT   file.
  isumm  = 6   ! The MINOS SUMMARY file.
  nout   = 6   ! Local output file (6 = screen).

  ! Initialize inform
  inform = 0

  ! RULE OF THUMB:
  ! MINOS won't open file ispec (for example)
  ! if that unit is already open, or if ispec = 0.
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! mistart MUST BE CALLED BEFORE ANY OTHER MINOS ROUTINE.
  !------------------------------------------------------------------
  call mistart( iprint, isumm, ispecs )  ! Initialize MINOS and open
                                         ! the specified files.

  !call mispec( ispecs, inform )          ! Read the SPECS file
                                         ! (if ispecs > 0).

  !----------------------------------------------------------------------
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

  !----------------------------------------------------------------------
  if (inform >= 2) then
     write(nout, *) 'inform =', inform
     write(nout, *) 'Error: ispecs > 0 but no SPECS file found'
     stop
  end if

  !======================================================================
  ! Solve the first problem: sets hs, which is used for subsequent warm-
  ! starts
  !======================================================================
  ! Initialize obj vals
  obj_vals(1:nVary) = 0
  ! And informs
  informs(1:nVary) = 0
  !----------------------------------------------------------------------
  ! Solve the problem.
  ! iobj   = 0    means there is no linear objective row in a(*).
  ! objadd = zero means there is no constant to be added to the objective.
  ! nname  = 1    means there are no meaningful names for the
  !               variables and constraints inside name1(*) and name2(*).
  !               MINOS will print default names.
  !----------------------------------------------------------------------
  ! hs(1:nb) = 0  ! Now provided
  xn(1:nb) = zero
  pi(1:m)  = zero
! iobj    says which row of A is a free row containing a linear
!         objective vector  c  (iobj = 0 if none).
!         iobj = 0  or  nncon < iobj le m.
  iObj     = m      ! Always the last row, appended on constraint matrix
  objadd   = zero
  nwcore   = lenz
  nname    = 1
  ObjAdd   = zero

  ! dumpLP.m generates slack bounds for SQOPT: Ax - s = 0
  ! Minos has Ax + s = 0, so we have to flip the bounds on s.

  z(1:m)     = -bl(n+1:nb)
  bl(n+1:nb) = -bu(n+1:nb)
  bu(n+1:nb) =   z(1:m)


  !======================================================================
  ! Solve all subsequent problems using hs from previous run
  ! Thus, no need to use Old basis file
  !======================================================================
  ! Warm start with 2nd SPECS file
  !======================================================================
  ! call mispec( ispecs, inform )
  ! if (inform .ge. 2) then
  !    write(nout, *) 'No second SPECS file'
  !    go to 900
  ! end if

  ! Get quad-precision for obj_coeffs
  obj_coeffsq(1:nVary) = zero
  call dq(obj_coeffsq, obj_coeffs, nVary)

  do i = 1, nVary
      !--------------------------------------------------------
      ! Update the objective function:
      ! A: last row; bl & bu objective becomes 
      ! A is pointers, row indices, and Aij
      ! ka: pointers (n+1). Column j is in ka[j]:(ka[j+1]-1)
      ! ha: row indices (ne)
      ! a : Aij (ne)
      !--------------------------------------------------------
      ! LY: comments from qminos56/src/mi02lib.f
      ! a(ne)   constraint matrix (Jacobian), stored column-wise
      ! ha(ne)  list of row indices for each nonzero in a(*)
      ! ka(n+1) set of pointers to beginning of each column of the constraint matrix
      !         within a(*) and ha(*). Must have ka(1) = 1 and ka(n+1) = ne+1.
      !--------------------------------------------------------
      ! No need to change bl & bu:
      ! bl(nb) = -bigbnd; bu(nb) = bigbnd, since that is the slack on obj row

      ! Reset coeffs to zero for all other columns, and current obj column to coeff
      l = 0
      do col = 1, n
        do ptr = ka(col), ka(col+1)-1
            l = l+1
            row = ha(ptr)
            if (row == m) then
                if (col == obj_inds(i)) then 
                    a(l) = obj_coeffsq(i) 
                else
                    a(l) = zero
                end if
            end if
        end do
      end do

      if (i == 1) then
          call qdq(  a, ne )
          call qdq( bl, nb )
          call qdq( bu, nb )

          ! Was basis provided to warm-start even the very first LP?
          if (warm_first) then
              print *, 'Calling minoss. Warm start with provided basis (hs)'

              call minoss( 'Warm', m, n, nb, ne, nname,        &
                  nncon, nnobj, nnjac,                &
                  iobj, objadd, names,                &
                  a, ha, ka, bl, bu, name1, name2,    &
                  hs, xn, pi, rc,                     &
                  inform, mincor, ns, ninf, sinf, obj, &
                  z, nwcore )
          else
              print *, 'Calling minoss. Cold start for first run'
              call minoss( 'Cold', m, n, nb, ne, nname,        &
                  nncon, nnobj, nnjac,                &
                  iobj, objadd, names,                &
                  a, ha, ka, bl, bu, name1, name2,    &
                  hs, xn, pi, rc,                     &
                  inform, mincor, ns, ninf, sinf, obj, &
                  z, nwcore )
          end if

          write(nout, *) ' '
          write(nout, *) 'Quad MINOS finished.'
          write(nout, *) 'inform =', inform
          write(nout, *) 'ninf   =', ninf
          write(nout, *) 'sinf   =', sinf
          write(nout, *) 'obj    =', obj

          ! Return xn rounded to double precision
          call qd( xn, xnd, nb)
          ! Save result
          obj_vals(i) = xnd(obj_inds(i))
          ! Save inform
          informs(i) = inform
      else
          ! Do we need to call qdq again?
          print *, 'Calling minoss. Warm start for subsequent run'

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
          ! Save result
          obj_vals(i) = xnd(obj_inds(i))
          ! Save inform
          informs(i) = inform
      end if
  end do

  900 close( iprint )
  close( ispecs )

  ! Return xn rounded to double precision
  ! call qd( xn, xnd, nb)

  ! deallocate( hs, pi, rc, xn )
  deallocate( pi, rc, xn )
  deallocate( name1, name2 )
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

  ! 11 Aug 2015: double to quad
  subroutine dq( aq, ad, n )
      integer(ip),  intent(in)  :: n
      real(dp),     intent(in)  :: ad(n)
      real(qp),    intent(inout):: aq(n)

      write(*,'(z17)') ad(n)
      aq = real(ad, qp)
      write(*,'(z32)') aq(n)
  end subroutine dq

!  ! 11 Aug 2015: first version of get_nz_ind.
!  !              Return index of A (2D array flattened column-wise)
!  !              corresponding to nonzero row/col index
!  subroutine get_nz_ind( a, n, row, col )
!      integer(ip), intent(in)   :: n, row, col
!      real(qp),    intent(in)   :: a(n)
!
!
!  end subroutine get_nz_ind

end subroutine qvaryME

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine funcon
end subroutine funcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine funobj
end subroutine funobj

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine matmod
end subroutine matmod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

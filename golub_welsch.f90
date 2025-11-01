module golub_welsch

!-------------------------------------------------------------------------------
!
! inputs: n_points - number of nodes (including -1 and 1)
! outputs: x - nodes (zeros of orthogonal polynomial, plus endpoints)
!          w - weights
! abstract:     computes the nodes (zeros) and weights for Lobatto quadrature
!               using the Golub-Welsch algorithm
!               using the lapack routine dstev
!---------------written by: Kawaii_Bokuchin_Puyuyu(2025/10/28)-----------------------------

  implicit none
  private
  public :: lobatto_nodes_weights, get_lobatto, dbg_enabled, set_dbg, dbg_dump_recurrence

  logical :: dbg_enabled = .true.

  ! Simple one-entry cache: remembers last requested n and arrays
  integer :: cached_n = 0
  real(8), allocatable :: cached_xi(:)
  real(8), allocatable :: cached_w(:)

contains

  subroutine set_dbg(flag)
    logical, intent(in) :: flag
    dbg_enabled = flag
  end subroutine set_dbg

  !-------------------------------------------------------------------
  ! Debugging utility: dump recurrence and tridiagonal diagnostics
  ! Controlled by module-level `dbg_enabled` flag. Prints a compact view
  ! (first 10 entries) of recurrence arrays and optional tridiag arrays.
  !-------------------------------------------------------------------
  subroutine dbg_dump_recurrence(a, b, c, d, diag, offdiag, eigvals, w_nodes, N, label)
    real(8), intent(in) :: a(:), b(:), c(:), d(:)
    real(8), intent(in), optional :: diag(:), offdiag(:), eigvals(:), w_nodes(:)
    integer, intent(in) :: N
    character(len=*), intent(in), optional :: label
    integer :: i, nprint
    if (.not. dbg_enabled) return
    nprint = min(10, N)
    if (present(label)) then
      print '(A)', 'DBG(' // trim(label) // '): recurrence (first 10)'
    else
      print '(A)', 'DBG: recurrence (first 10)'
    end if
    do i = 1, nprint
      print '(I3,2X,4(1P,E12.4))', i, a(i), b(i), c(i), d(i)
    end do
    if (present(diag)) then
      print '(A)', ' DBG: diag, offdiag (first 10)'
      do i = 1, nprint
        if (i <= size(offdiag)) then
          print '(I3,2X,2(1P,E12.4))', i, diag(i), offdiag(i)
        else
          print '(I3,2X,1P,E12.4)', i, diag(i)
        end if
      end do
    end if
    if (present(eigvals)) then
      print '(A)', ' DBG: eigvals (first 10)'
      do i = 1, nprint
        ! Detect NaN/Inf/oversize and print a clear marker instead of field overflow stars
        if (eigvals(i) /= eigvals(i) .or. abs(eigvals(i)) > 1.0d300) then
          write(*,'(I3,2X,A)') i, ' <invalid>'
        else
          write(*,'(I3,2X,1P,E16.8)') i, eigvals(i)
        end if
      end do
    end if
    if (present(w_nodes)) then
      print '(A)', ' DBG: raw w_nodes before endpoint adjustment (first 10)'
      do i = 1, nprint
        if (w_nodes(i) /= w_nodes(i) .or. abs(w_nodes(i)) > 1.0d300) then
          write(*,'(I3,2X,A)') i, ' <invalid>'
        else
          write(*,'(I3,2X,1P,E16.8)') i, w_nodes(i)
        end if
      end do
    end if
  end subroutine dbg_dump_recurrence

  !-------------------------------------------------------------------
  ! build_recurrence
  ! Fill recurrence coefficients a,b,c,d for the interior
  ! this is especially for Lobatto nodes (Jacobi(1,1))
  !-------------------------------------------------------------------
  subroutine build_recurrence(N_internal, a, b, c, d)
    integer, intent(in) :: N_internal
    real(8), intent(out) :: a(:), b(:), c(:), d(:)
    integer :: i, n
    ! Standardized form: a_n p_{n+1} = (b_n + c_n x) p_n - d_n p_{n-1}
    ! We adopt a_n = 1 and return arrays with n = i-1 mapping.
    ! For Lobatto (Jacobi alpha=1,beta=1) the DLMF-specialized formulas are:
    !   A_n = (2n+3)(2n+4) / (2 (n+1)(n+3))  (this will be c_n)
    !   B_n = 0.0                                (this will be b_n)
    !   C_n = (n+2)/(n+3)                       (this will be d_n)
    do i = 1, N_internal
      n = i - 1
      a(i) = 1.0d0
      b(i) = 0.0d0
      ! c_n = A_n = (2n+3)(2n+4) / (2 (n+1)(n+3))
      c(i) = dble( (2*n + 3) * (2*n + 4) ) / ( 2.0d0 * dble(n + 1) * dble(n + 3) )
      ! d_n = C_n = (n+2)/(n+3)
      d(i) = dble(n + 2) / dble(n + 3)
    end do
  end subroutine build_recurrence

  !-------------------------------------------------------------------
  ! build_jacobi_from_recurrence
  ! Convert three-term recurrence coefficients to symmetric Jacobi
  ! tridiagonal (diag, offdiag). 
  !-------------------------------------------------------------------
  subroutine build_jacobi_from_recurrence(N, a, b, c, d, diag, offdiag, info)
    integer, intent(in) :: N
    real(8), intent(in) :: a(:), b(:), c(:), d(:)
    real(8), intent(out) :: diag(:), offdiag(:)
    integer, intent(out) :: info
    integer :: i

    info = 0
    if (N < 1) then
      info = -1; return
    end if

    do i = 1, N
      diag(i) = b(i)
    end do

    if (N > 1) then
      do i = 1, N-1
        if (c(i) == 0d0 .or. c(i+1) == 0d0) then
          offdiag(i) = 0d0
        else
          offdiag(i) = sqrt( max(0d0, a(i) * d(i+1) / ( c(i) * c(i+1) ) ) )
        end if
      end do
    end if

  end subroutine build_jacobi_from_recurrence

  !-------------------------------------------------------------------
  ! diagonalize_tridiag
  ! Diagonalize symmetric tridiagonal using LAPACK dstev (compute vecs)
  !-------------------------------------------------------------------
  subroutine diagonalize_tridiag(N, diag_in, offdiag_in, eigvals, eigvecs, info)
    integer, intent(in) :: N
    real(8), intent(in) :: diag_in(:), offdiag_in(:)
    real(8), intent(out) :: eigvals(:), eigvecs(:,:)
    integer, intent(out) :: info
    real(8), allocatable :: D(:), E(:), Z(:,:), WORK(:)
  integer :: lwork, i
    character(len=1) :: jobz
    integer :: info_lap

    info = 0
    if (N < 1) then
      info = -1
      return
    end if

    allocate(D(N))
    if (N > 1) then
      allocate(E(N-1))
    else
      allocate(E(1)); E(1) = 0d0
    end if

    D = diag_in
    if (N > 1) then
      E(1:N-1) = offdiag_in(1:N-1)
    end if

    allocate(Z(N,N))
    jobz = 'V'
    lwork = max(1, 2*N-1)
    allocate(WORK(lwork))

    call dstev(jobz, N, D, E, Z, N, WORK, info_lap)
    if (info_lap /= 0) then
      ! Debug print to help diagnose LAPACK failure
      print '(A,I0)', 'DBG golub_welsch: dstev failed, info_lap=', info_lap
      print '(A,I0)', ' DBG: N=', N
      info = info_lap
      return
    end if

    ! Basic sanity check on returned eigenvalues (D) to catch corrupted results
    do i = 1, N
      if (D(i) /= D(i) .or. abs(D(i)) > 1.0d300) then
        print '(A)', 'DBG golub_welsch: dstev returned invalid eigenvalues (NaN/Inf/overflow)'
        info = -999
        deallocate(D, E, Z, WORK)
        return
      end if
    end do

    eigvals(1:N) = D(1:N)
    eigvecs(1:N,1:N) = Z(1:N,1:N)

    deallocate(D, E, Z, WORK)

  end subroutine diagonalize_tridiag

  !-------------------------------------------------------------------
  ! lobatto_nodes_weights
  !-------------------------------------------------------------------
  subroutine lobatto_nodes_weights(n_points, x, w, info)
    integer, intent(in) :: n_points
    integer, intent(out) :: info
    real(8), intent(out) :: x(:), w(:)

    integer :: N
    real(8), allocatable :: diag(:), offdiag(:)
    real(8), allocatable :: eigvals(:), eigvecs(:,:)
    real(8), allocatable :: w_nodes(:)
    real(8), allocatable :: a(:), b(:), c(:), d(:)
    integer :: i, ierr, j
    real(8) :: mu0

    info = 0
    if (n_points < 2) then
      info = -1; return
    end if

    if (n_points == 2) then
      if (size(x) < 2 .or. size(w) < 2) then
        info = -2; return
      end if
      x(1) = -1d0; x(2) = 1d0
      w(1) = 1d0; w(2) = 1d0
      return
    end if

    N = n_points - 2

    allocate(a(N), b(N), c(N), d(N))
    call build_recurrence(N, a, b, c, d)

    allocate(diag(N), offdiag(max(1,N-1)))
    call build_jacobi_from_recurrence(N, a, b, c, d, diag, offdiag, ierr)
    if (ierr /= 0) then
      info = ierr; deallocate(a,b,c,d); return
    end if

    allocate(eigvals(N), eigvecs(N,N))
    call diagonalize_tridiag(N, diag, offdiag, eigvals, eigvecs, ierr)
    if (ierr /= 0) then
      info = ierr; deallocate(a,b,c,d,diag,offdiag); return
    end if

    allocate(w_nodes(N))
    mu0 = 4.0d0 / 3.0d0  ! Zeroth moment for Jacobi(1,1): integral_{-1}^1 (1 - x^2) dx
    do j = 1, N
      w_nodes(j) = mu0 * eigvecs(1, j)**2
    end do

    ! Debug: optional dump of recurrence and tridiagonal pieces
    if (dbg_enabled) then
      call dbg_dump_recurrence(a, b, c, d, diag, offdiag, eigvals, w_nodes, N, 'lobatto')
    end if

    ! Adjust to Lobatto weights
    do i = 1, N
      w_nodes(i) = w_nodes(i) / (1.0d0 - eigvals(i)**2)
    end do

    if (size(x) < n_points .or. size(w) < n_points) then
      info = -3; deallocate(diag,offdiag,eigvals,eigvecs,w_nodes,a,b,c,d); return
    end if

    x(1) = -1d0
    w(1) = 2.0d0 / ( dble(n_points) * dble(n_points - 1) )
    do i = 1, N
      x(i+1) = eigvals(i)
      w(i+1) = w_nodes(i)
    end do
    x(n_points) = 1d0
    w(n_points) = w(1)

    deallocate(diag, offdiag, eigvals, eigvecs, w_nodes, a, b, c, d)

  end subroutine lobatto_nodes_weights

  subroutine get_lobatto(n_points, xi_ref, w_ref, info)
    ! Public wrapper that returns Lobatto nodes/weights on reference [-1,1]
    ! Uses a simple cache for the most recent n to avoid repeated diagonalizations.
    integer, intent(in) :: n_points
    real(8), intent(out) :: xi_ref(:), w_ref(:)
    integer, intent(out), optional :: info

    integer :: ierr

    ! Check sizes
    if (size(xi_ref) < n_points .or. size(w_ref) < n_points) then
      ierr = -3
      if (present(info)) info = ierr
      return
    end if

    ! If cache hit, copy and return
    if (cached_n == n_points) then
      xi_ref(1:n_points) = cached_xi(1:n_points)
      w_ref(1:n_points) = cached_w(1:n_points)
      if (present(info)) info = 0
      return
    end if

    ! Cache miss: compute via existing routine
    call lobatto_nodes_weights(n_points, xi_ref, w_ref, ierr)
    if (ierr /= 0) then
      if (present(info)) info = ierr
      return
    end if

    ! Store into cache (replace previous)
    if (allocated(cached_xi)) deallocate(cached_xi)
    if (allocated(cached_w)) deallocate(cached_w)
    allocate(cached_xi(n_points))
    allocate(cached_w(n_points))
    cached_xi = xi_ref(1:n_points)
    cached_w = w_ref(1:n_points)
    cached_n = n_points

    if (present(info)) info = 0

  end subroutine get_lobatto

end module golub_welsch

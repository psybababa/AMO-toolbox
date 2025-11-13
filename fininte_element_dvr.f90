module finite_element_dvr

    !-------------------------------------------------------------------
    !
    ! Finite Element Discrete Variable Representation (FEDVR) module
    ! This module provides below functionalities and self-contained API:
    !
    ! - compute_barycentric_weights : help function to compute lagrange polynomial via barycentric formula
    ! - compute_derivative_matrix : compute derivative matrix D_ref on reference element
    ! - build_element_matrices : compute M_e (diag) and K_e (dense) for element
    ! - build_potential_diag : build a potential vector for each element in hamiltonianian
    ! - build_hamiltonian_matrix : assemble H = T + V from K and V_diag
    ! - assemble_global_matrices : assemble global M and K matrices from element matrices
    ! - assemble_global_hamiltonian : assemble global H from element H matrices
    ! - build_global_grid : create global node coordinates, weights, and maps
    !
    !--------------------written by: Kawaii_Bokuchin_Puyuyu(2025/10/29)----------------------

    use, intrinsic :: iso_fortran_env, only: wp => real64
    use golub_welsch
    implicit none
    private
    public ::  build_element_matrices, compute_barycentric_weights, compute_differentiation_matrix
    public ::  assemble_global_matrices, build_hamiltonian_matrix, assemble_global_hamiltonian
    public ::  build_global_grid, build_potential_diag

    contains

    !-------------------------------------------------------------------
    ! compute_barycentric_weights
    ! Compute barycentric weights for Lagrange polynomials at given nodes
    ! this barycentric weights is different from the weights for quadrature
    ! thus it is named as b_weights
    ! reference: Berrut, J. P., & Trefethen, L. N. (2004). Barycentric lagrange interpolation. SIAM review, 46(3), 501-517.
    !-------------------------------------------------------------------

subroutine compute_barycentric_weights(n, x_node, b_weights)
    integer, intent(in) :: n
    real(wp), intent(in) :: x_node(n)
    real(wp), intent(out) :: b_weights(n)
    integer :: i, j
    real(wp) :: d, sumlog
    integer :: sign_prod
    real(wp), parameter :: collision_tol = 1.0e-15_wp

    ! Numerically stable barycentric weights using log-sum of absolute differences.
    ! w_i = 1 / prod_{j != i} (x_i - x_j) = sign * exp(- sum_j log(|x_i-x_j|) )
    do i = 1, n
        sumlog = 0.0_wp
        sign_prod = 1
        do j = 1, n
            if (j == i) cycle
            d = x_node(i) - x_node(j)
            if (abs(d) < collision_tol) then
                b_weights = 0.0_wp
                return
            end if
            if (d < 0.0_wp) sign_prod = -sign_prod
            sumlog = sumlog + log(abs(d))
        end do
        b_weights(i) = real(sign_prod, wp) * exp(-sumlog)
    end do
end subroutine compute_barycentric_weights

!--------------------------------------------------------------
!
! compute Derivative Matrix D_ref on [-1,1]
! this subroutine should only work once to make global D matrix
! D matrix is used to construct kinetic energy matrix
! use compute_barycentric_weights first to get b_weights
! D(i,j) = (w_j / w_i) / (x_i - x_j)  for i != j
! D(i,i) = - sum_{j != i} D(i,j)
! barycentric weights and formula reference: Berrut, J. P., & Trefethen, L. N. (2004). Barycentric lagrange interpolation. SIAM review, 46(3), 501-517.
!--------------------------------------------------------------

subroutine compute_differentiation_matrix(n, x_node, b_weights, D_ref)
    integer, intent(in) :: n
    real(wp), intent(in) :: x_node(n), b_weights(n)
    real(wp), intent(out) :: D_ref(n, n)

    integer :: i, j
    real(wp) :: denom, tmp1, tmp2
    real(wp), allocatable :: diag(:)

    allocate(diag(n))
    diag = 0.0_wp
    D_ref = 0.0_wp

    ! Pairwise loop: compute off-diagonals for (i,j) and (j,i) in one pass
    ! and accumulate diagonal contributions using identity
    !   D(ii) = sum_{k != i} 1/(x_i - x_k)
    do i = 1, n-1
        do j = i+1, n
            denom = x_node(i) - x_node(j)
            if (abs(denom) < 1.0e-12_wp) then  ! Near-colliding nodes -> invalid; zero output and return
                D_ref = 0.0_wp
                diag = 0.0_wp
                deallocate(diag)
                return
            end if

            tmp1 = b_weights(j) / ( b_weights(i) * denom )    ! D(i,j)
            tmp2 = - b_weights(i) / ( b_weights(j) * denom )  ! D(j,i) = negative reciprocal adjusted

            D_ref(i,j) = tmp1
            D_ref(j,i) = tmp2

            ! accumulate diagonal sums: diag(i) += 1/(x_i-x_j), diag(j) += 1/(x_j-x_i) = -1/(x_i-x_j)
            diag(i) = diag(i) + 1.0_wp / denom
            diag(j) = diag(j) - 1.0_wp / denom
        end do
    end do

    ! Set diagonals
    do i = 1, n
        D_ref(i,i) = diag(i)
    end do

    deallocate(diag)
end subroutine compute_differentiation_matrix

    !------------------------------------------------------------------
    ! build_element_matrices
    !
    ! Build local element mass (M_e) and stiffness/kinetic (K_e) matrices for
    ! an element with physical endpoints a_e,b_e.  The reference derivative
    ! matrix D_ref and reference weights w_ref are inputs.
    !
    ! Formulas used:
    !   J = (b_e - a_e)/2 (jacobian  corresponds to [-1,1] -> [a_e,b_e])
    !   x(ξ) = mid + J * ξ
    !   M_e(ii) = w_ref(i) * J        (mass diagonal entries)
    !   K_e = (1/J) * ( D_ref^T * diag(w_ref) * D_ref )
    !------------------------------------------------------------------

subroutine build_element_matrices(N, D_ref, w_ref, a_e, b_e, M_e, K_e, info)
        integer, intent(in) :: N
        real(wp), intent(in) :: D_ref(:,:), w_ref(:)
        real(wp), intent(in) :: a_e, b_e
        real(wp), intent(out) :: M_e(:), K_e(:,:)
        integer, intent(out), optional :: info
        real(wp) :: J
        real(wp), allocatable :: W(:,:), temp(:,:)
        integer :: ierr, i

        ierr = 0
        if (size(D_ref,1) < N .or. size(D_ref,2) < N .or. size(w_ref) < N) then
            ierr = -1
            if (present(info)) info = ierr
            return
        end if
        if (size(M_e) < N .or. size(K_e,1) < N .or. size(K_e,2) < N) then
            ierr = -2
            if (present(info)) info = ierr
            return
        end if

        J = 0.5_wp * (b_e - a_e)

        ! M_e diagonal
        do i = 1, N
            M_e(i) = w_ref(i) * J
        end do

        ! K_e = (1/J) * D_ref^T * diag(w_ref) * D_ref
        allocate(W(N,N), temp(N,N))
        W = 0.0_wp
        do i = 1, N
            W(i,i) = w_ref(i)
        end do
        temp = matmul(W, D_ref)
        K_e = matmul(transpose(D_ref), temp)
        K_e = K_e / J

        deallocate(W, temp)

        if (present(info)) info = ierr
end subroutine build_element_matrices

    !------------------------------------------------------------------
    !
    ! build_potential_diag
    ! currently the potential is hardcoded as V(x) = -Z_eff / x
    ! Z_eff = 1 + (Z_nuclei - 1) * WS  (WS Woods-Saxon screening, residual +1 at infinity)
    ! WS = 1 / (1 + (\eta/kia) * (exp(kia*r) - 1) )
    ! Z_nuclei : nuclear charge 54 for Xe
    ! eta : screening parameter 5.197
    ! kia : 1.048
    !
    ! key_soi = 1
    ! dV/dr = - Z_eff / r^2 - DZ_eff / r
    ! Dz_eff = -(Z_nuclei - 1) * (eta * exp(kia*r)) * WS^2
    ! V_soi = alpha^2 / (4 * r) * dV/dr * (J^2 - L^2 - 3/4 ) 
    ! 
    ! this will be modified to take various potential forms later 
    !
    !------------------------------------------------------------------

subroutine build_potential_diag(n, x, v, key_soi, j_tot, L, info)
    integer, intent(in) :: n                 ! number of global nodes
    integer, intent(in) :: key_soi           ! include spin-orbit if 1
    real(wp), intent(in) :: j_tot            ! total angular momentum j (can be half-integer)
    integer, intent(in) :: L                 ! orbital angular momentum quantum number
    real(wp), intent(in) :: x(:)             ! partition node coordinates (r)
    real(wp), intent(out) :: v(:)            ! potential diagonal entries
    integer, intent(out), optional :: info

    real(wp) :: Z_nuclei, eta, kia, alpha, z_core
    real(wp) :: r, WS, Z_eff, dZ_eff, dV, V_soi
    integer :: i
    integer :: ierr
    real(wp), parameter :: rmin = 1.0e-12_wp  ! small-radius clip

    ierr = 0
    if (size(x) < n .or. size(v) < n) then
        ierr = -1
        if (present(info)) info = ierr
        return
    end if

    ! ---- model constants (Data extracted from TDSE simulator) ----
    Z_nuclei = 54.0_wp
    z_core = Z_nuclei - 1.0_wp
    eta = 5.197_wp
    kia = 1.048_wp
    ! fine-structure constant alpha ~ 1/137
    alpha = 1.0_wp / 137.035999084_wp

    do i = 1, n
        r = x(i)
        if (r < rmin) r = rmin

        WS = 1.0_wp / (1.0_wp + (eta/kia) * (exp(kia*r) - 1.0_wp) )
        ! Ensure Z_eff -> Z_nuclei near the origin and -> 1 for large r (residual ion)
        Z_eff = 1.0_wp + z_core * WS

        if (key_soi == 1) then
            ! derivative of Z_eff wrt r
            dZ_eff = -z_core * (eta * exp(kia*r)) * WS**2
            ! derivative of V = -Z_eff(r)/r  => dV/dr = - dZ_eff/r + Z_eff/r^2
            dV = - dZ_eff / r + Z_eff / (r**2)

            ! angular factor: j(j+1) - l(l+1) - s(s+1) with s=1/2
            V_soi = (alpha**2) / (4.0_wp * r) * dV * ( j_tot*(j_tot+1.0_wp) - real(L*(L+1),wp) - 3.0_wp/4.0_wp ) 
            v(i) = - Z_eff / r + V_soi
        else
            v(i) = - Z_eff / r
        end if
    end do

    if (present(info)) info = ierr
end subroutine build_potential_diag

!------------------------------------------------------------------
!
! build_hamiltonian_matrix
! Assemble H = T + V in atomic units from global stiffness K_global
! and potential values V_diag and global mass diagonal M_global.
! In atomic units: T = 1/2 * K_global  (K defined as integral phi' phi' dx)
!
! Inputs:
!   K_global(:,:) - global stiffness matrix
!   V_diag(:)     - potential values at global nodes
!   M_global(:)   - global mass diagonal (integration weights)
!
! Output:
!   H_global(:,:) - assembled Hamiltonian (returned)
!
!------------------------------------------------------------------

subroutine build_hamiltonian_matrix(K_global, V_diag, M_global, H_global, info)
    real(wp), intent(in) :: K_global(:,:), V_diag(:), M_global(:)
    real(wp), intent(out) :: H_global(size(K_global,1), size(K_global,2))
    integer, intent(out), optional :: info

    integer :: N, ierr, i, j
    real(wp) :: sqrt_mi, sqrt_mj, sqrt_mimj

    ierr = 0
    N = size(K_global,1)


    H_global = 0.0_wp

    ! Orthogonalized kinetic: T_ik = 0.5 * K_ik / sqrt(M_i * M_k)
    do i = 1, N
        if (M_global(i) <= 0.0_wp) then  ! Avoid div by zero
            ierr = -3
            cycle
        end if
        sqrt_mi = sqrt(M_global(i))
        do j = 1, N
            if (M_global(j) <= 0.0_wp) cycle
            sqrt_mj = sqrt(M_global(j))
            sqrt_mimj = sqrt_mi * sqrt_mj
            H_global(i,j) = 0.5_wp * K_global(i,j) / sqrt_mimj
        end do
    end do

    ! Potential: V_ii = V_diag(i) (no M multiplication)
    do i = 1, N
        H_global(i,i) = H_global(i,i) + V_diag(i)
    end do

    if (present(info)) info = ierr
end subroutine build_hamiltonian_matrix

!------------------------------------------------------------------
! assemble_global_matrices
! Assemble global mass (diagonal) and stiffness matrices from element
! matrices using elem_map(e,i) -> global index.
!------------------------------------------------------------------

subroutine assemble_global_matrices(nelems, N, elem_map, M_elems, K_elems, N_global, M_global, K_global, info)
    integer, intent(in) :: nelems, N
    integer, intent(in) :: elem_map(nelems, N)
    real(wp), intent(in) :: M_elems(N, nelems)
    real(wp), intent(in) :: K_elems(N, N, nelems)
    integer, intent(inout) :: N_global
    real(wp), intent(inout) :: M_global(:)
    real(wp), intent(inout) :: K_global(:,:)
    integer, intent(out), optional :: info

    integer :: e, i, j, gi, gj, ierr

    ierr = 0
    if (nelems < 1 .or. N < 1) then
        ierr = -1
        if (present(info)) info = ierr
        return
    end if

    if (N_global <= 0) then
        N_global = 0
        do e = 1, nelems
            do i = 1, N
                if (elem_map(e,i) > N_global) N_global = elem_map(e,i)
            end do
        end do
    end if

    if (size(M_global) < N_global .or. size(K_global,1) < N_global .or. size(K_global,2) < N_global) then
        ierr = -2
        if (present(info)) info = ierr
        return
    end if

    M_global = 0.0_wp
    K_global = 0.0_wp

    do e = 1, nelems
        do i = 1, N
            gi = elem_map(e,i)
            if (gi < 1 .or. gi > N_global) cycle
            M_global(gi) = M_global(gi) + M_elems(i,e)
            do j = 1, N
                gj = elem_map(e,j)
                if (gj < 1 .or. gj > N_global) cycle
                K_global(gi,gj) = K_global(gi,gj) + K_elems(i,j,e)
            end do
        end do
    end do

    if (present(info)) info = ierr
end subroutine assemble_global_matrices

!------------------------------------------------------------------
!
! assemble_global_hamiltonian
! to assemble global hamiltonian matrix H = T + V not in partition but globally
! so this matrix is what we finally need to diagonalize and get eigenvalues/vectors
!
! Hamiltonian is expressed in DVR basis and atomic units are used throughout
! 
!------------------------------------------------------------------

    subroutine assemble_global_hamiltonian(key_soi, j_tot, L, N_global, x_global, M_global, H_global, K_global, info)
        integer, intent(in) :: key_soi          ! include spin-orbit if 1
        real(wp), intent(in) :: j_tot           ! total angular momentum j (can be half-integer)
        integer, intent(in) :: L                 ! orbital angular momentum quantum number
        integer, intent(inout) :: N_global      ! total number of global DVR points (may be inferred)
        real(wp), intent(in) :: x_global(:)       ! global node coordinates
        real(wp), intent(in) :: M_global(:)     ! global weights (mass diagonal entries)
        real(wp), intent(out) :: H_global(:,:)  ! global hamiltonian matrix
        real(wp), intent(in), optional :: K_global(:,:) ! optional global stiffness matrix
        integer, intent(out), optional :: info

        integer :: ierr, ierr2
        real(wp), allocatable :: V_diag(:)
        integer :: i, N

        ierr = 0
        N = N_global
        if (N <= 0) N = size(x_global)

        if (size(x_global) < N .or. size(M_global) < N) then
            ierr = -1
            if (present(info)) info = ierr
            return
        end if

        allocate(V_diag(N))
        call build_potential_diag(N, x_global, V_diag, key_soi, j_tot, L, ierr)
        if (ierr /= 0) then
            if (present(info)) info = ierr
            deallocate(V_diag)
            return
        end if

        ! Initialize H
        H_global = 0.0_wp

        ! If stiffness provided, use full build_hamiltonian_matrix
        if (present(K_global)) then
            call build_hamiltonian_matrix(K_global, V_diag, M_global, H_global, ierr2)
            if (ierr2 /= 0) then
                ierr = ierr2
                if (present(info)) info = ierr
                deallocate(V_diag)
                return
            end if
        else
            ! Only potential part (mass-weighted diagonal)
            do i = 1, N
                H_global(i,i) = H_global(i,i) + V_diag(i) * M_global(i)
            end do
        end if

        if (present(info)) info = ierr
        deallocate(V_diag)
    end subroutine assemble_global_hamiltonian

    !------------------------------------------------------------------
    !
    ! build_global_grid
    ! Create global node coordinates, weights, and element-to-global mapping
    ! using refence nodes and weights on [-1,1] and element boundaries.
    ! Inputs:
    !   nelems        - number of all elements
    !   N             - number of nodes per partition
    !   a, b          - global domain [a,b]
    ! Outputs:
    !   x_global(:)   - global node coordinates
    !   N_global      - total number of global nodes
    !   w_global(:)   - global weights (mass diagonal entries)
    !   elem_map(:,:) - element to global node index mapping
    ! ------------------------------------------------------------------
    
    subroutine build_global_grid(nelems, N, x_ref, w_ref, &
                                    a, b, N_global, &
                                    x_global, M_global, elem_map, info)
        integer, intent(in) :: nelems, N
        real(wp), intent(in) :: a, b
        real(wp), intent(out) :: x_ref(:), w_ref(:)
        integer, intent(out) :: N_global
    real(wp), intent(out) :: x_global(:), M_global(:)
        integer, intent(out) :: elem_map(:,:)
        integer, intent(out), optional :: info

    integer :: e, i, ierr, required_N, gi
        real(wp) :: x_left, x_right
        ! temporaries for first-element reference
        real(wp), allocatable :: x_ref_first(:), w_ref_first(:)

        ierr = 0
        ! Diagnostic prints are guarded by the shared dbg flag from golub_welsch.
        if (dbg_enabled) then
            ! Debug: dump entry-state
            print '(A)', 'DBG(build_global_grid) ENTRY'
            print '(A,I0)', ' DBG: nelems=', nelems
            print '(A,I0)', ' DBG: N=', N
            print '(A,I0)', ' DBG: size(x_ref)=', size(x_ref)
            print '(A,I0)', ' DBG: size(w_ref)=', size(w_ref)
            print '(A,I0)', ' DBG: size(x_ref)=', size(x_ref)
            print '(A,I0)', ' DBG: size(w_ref)=', size(w_ref)
            if (size(x_global) > 0) then
                print '(A,I0)', ' DBG: size(x_global)=', size(x_global)
            else
                print '(A)', ' DBG: x_global size is 0'
            end if
            if (size(M_global) > 0) then
                print '(A,I0)', ' DBG: size(M_global)=', size(M_global)
            else
                print '(A)', ' DBG: M_global size is 0'
            end if
            if (size(elem_map,1) > 0 .and. size(elem_map,2) > 0) then
                print '(A,I0,A,I0)', ' DBG: size(elem_map)=', size(elem_map,1), ',', size(elem_map,2)
            else
                print '(A)', ' DBG: elem_map size is zero in at least one dimension'
            end if
        end if

    required_N = nelems*(N-1) + 1
    N_global = required_N

        ! Basic size checks for output arrays
        if (size(x_ref) < N .or. size(w_ref) < N) then
            ierr = -2
            if (present(info)) info = ierr
            if (dbg_enabled) then
                print '(A)', 'DBG build_global_grid: FAIL size(x_ref) or size(w_ref) < N'
                print '(A,I0)', ' DBG: size(x_ref)=', size(x_ref)
                print '(A,I0)', ' DBG: size(w_ref)=', size(w_ref)
                print '(A,I0)', ' DBG: N=', N
            end if
            return
        end if
        if (size(x_global) < required_N .or. size(M_global) < required_N) then
            ierr = -3
            if (present(info)) info = ierr
            if (dbg_enabled) then
                print '(A)', 'DBG build_global_grid: FAIL size(x_global) or size(M_global) < required_N'
                print '(A,I0)', ' DBG: size(x_global)=', size(x_global)
                print '(A,I0)', ' DBG: size(M_global)=', size(M_global)
                print '(A,I0)', ' DBG: required_N=', required_N
            end if
            return
        end if
        if (size(elem_map,1) < nelems .or. size(elem_map,2) < N) then
            ierr = -4
            if (present(info)) info = ierr
            if (dbg_enabled) then
                print '(A)', 'DBG build_global_grid: FAIL size(elem_map) too small'
                print '(A,I0)', ' DBG: size(elem_map,1)=', size(elem_map,1)
                print '(A,I0)', ' DBG: size(elem_map,2)=', size(elem_map,2)
                print '(A,I0)', ' DBG: nelems=', nelems
                print '(A,I0)', ' DBG: N=', N
            end if
            return
        end if

        ! construct reference nodes and weights on [-1,1]
        ! use public wrapper in golub_welsch
        call get_lobatto(N, x_ref, w_ref, ierr)
        if (ierr /= 0) then
            if (present(info)) info = ierr
            print '(A,I0)', 'DBG build_global_grid: get_lobatto failed, ierr=', ierr
            return
        end if

    ! Initialize global weights to zero 
    M_global = 0.0_wp

        ! For simplicity use Lobatto reference nodes for all elements.
        ! Keep x_ref_first/w_ref_first as copies of the reference Lobatto rule.
    allocate(x_ref_first(N), w_ref_first(N))
    x_ref_first = x_ref(1:N)
    w_ref_first = w_ref(1:N)

        ! Map reference nodes to global coordinates and accumulate weights
        do e = 1, nelems
            x_left = a + (e-1)*(b-a)/nelems
            x_right = a + e*(b-a)/nelems
            if (e == 1) then
                do i = 1, N
                    gi = (e-1)*(N-1) + i
                    x_global(gi) = 0.5_wp*(x_left + x_right) + 0.5_wp*(x_right - x_left)*x_ref_first(i)
                    M_global(gi) = M_global(gi) + 0.5_wp*(x_right - x_left)*w_ref_first(i)
                    elem_map(e,i) = gi
                end do
            else
                do i = 1, N
                    gi = (e-1)*(N-1) + i
                    x_global(gi) = 0.5_wp*(x_left + x_right) + 0.5_wp*(x_right - x_left)*x_ref(i)
                    M_global(gi) = M_global(gi) + 0.5_wp*(x_right - x_left)*w_ref(i)
                    elem_map(e,i) = gi
                end do
            end if
        end do

        deallocate(x_ref_first, w_ref_first)

        if (present(info)) info = ierr
    end subroutine build_global_grid


end module finite_element_dvr

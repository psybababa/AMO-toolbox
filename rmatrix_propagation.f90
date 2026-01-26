module rmatrix_propagation
    !------------------------------------------------------------------
    ! This module implements the R-matrix propagation method for solving
    ! coupled-channel equations in quantum scattering problems.
    ! It handles sector basis initialization, Hamiltonian construction,
    ! and R-matrix propagation steps.
    !------------------------------------------------------------------
    ! Assisted by GitHub Copilot
    !------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use finite_element_dvr
    use aeigen_solver
    use golub_welsch, only: get_lobatto
    implicit none
    private
    public :: sector_basis_t, init_sector_basis, build_and_solve_sector
    public :: propagate_step_outward, propagate_step_inward
    public :: export_coupling_matrices

    ! Structure to hold sector data
    type :: sector_basis_t
        real(wp) :: r_min, r_max, r_center
        integer :: n_channels ! Number of coupled channels (adiabatic states)
        integer :: n_basis    ! Number of DVR basis functions per channel in this sector
        
        ! Sector Hamiltonian Eigenpairs
        ! We solve (H_sector + L) u = E u
        ! eigenvalues: E_k
        ! surface amplitudes: u_k(r_min) and u_k(r_max)
        ! Dimensions: [n_channels * n_basis]
        complex(wp), allocatable :: eigvals(:) 
        
        ! Surface amplitudes for each eigenstate k
        ! Dimensions: [n_channels, n_total_states] where n_total_states = n_channels * n_basis
        ! For each eigenstate k, we store the channel vector at the boundaries.
        complex(wp), allocatable :: surface_amp_L(:,:) ! [channel, k] value at r_min
        complex(wp), allocatable :: surface_amp_R(:,:) ! [channel, k] value at r_max
        
        ! Matrices P and W at sector center (stored for reference/debugging)
        complex(wp), allocatable :: P_mat(:,:)      ! Weighted Coupling (SOI_DPHI)
        complex(wp), allocatable :: W_mat(:,:)      ! Potential + SOI (SOI_PHI)
        complex(wp), allocatable :: P_nac_mat(:,:)  ! Standard NAC (NAC)
    end type sector_basis_t

contains

    !----------------------------------------------------------------
    ! init_sector_basis
    ! Initializes a sector structure
    !----------------------------------------------------------------
    subroutine init_sector_basis(sector, r_min, r_max, n_channels, n_basis)
        type(sector_basis_t), intent(out) :: sector
        real(wp), intent(in) :: r_min, r_max
        integer, intent(in) :: n_channels, n_basis
        
        sector%r_min = r_min
        sector%r_max = r_max
        sector%r_center = 0.5_wp * (r_min + r_max)
        sector%n_channels = n_channels
        sector%n_basis = n_basis
        
        ! Allocations will happen in build_and_solve_sector
    end subroutine init_sector_basis

    !----------------------------------------------------------------
    ! build_and_solve_sector
    !
    ! Constructs the sector Hamiltonian and diagonalizes it.
    ! The equation solved is: d^2/deta^2 + A(eta) d/deta + B(eta)
    ! A and B matrices are evaluated at the sector center.
    ! The sector Hamiltonian is constructed using a local DVR basis.
    !
    ! Inputs:
    !   sector: The sector structure to be populated
    !   nelems_ad, N_per_elem_ad, x_max_ad: Parameters for the adiabatic solver
    !   Mj, F_val: Physical parameters (Angular momentum, Field strength)
    !----------------------------------------------------------------
    subroutine build_and_solve_sector(sector, nelems_ad, N_per_elem_ad, x_max_ad, Mj, F_val, info)
        type(sector_basis_t), intent(inout) :: sector
        integer, intent(in) :: nelems_ad, N_per_elem_ad
        real(wp), intent(in) :: x_max_ad, Mj, F_val
        ! E_target is removed from arguments as it is handled during propagation
        integer, intent(out), optional :: info
        
        integer :: ierr, n_tot, i, k, mu, nu, bi, bj, idx_i, idx_j
        complex(wp), allocatable :: P_center(:,:), W_center(:,:), P_nac_center(:,:)
        complex(wp), allocatable :: W_up_down(:,:), W_down_up(:,:)
        complex(wp), allocatable :: evals_ad(:)
        complex(wp), allocatable :: A_mat(:,:), B_mat(:,:)
        complex(wp), allocatable :: H_sector(:,:)
        real(wp), allocatable :: D_local(:,:)
        real(wp), allocatable :: b_weights(:), x_ref(:)
        real(wp) :: eta, h_width
        real(wp) :: val_L, val_R, w_phys_L, w_phys_R
        complex(wp) :: term_val
        complex(wp), allocatable :: evals_sector(:), evecs_sector(:,:)
        complex(wp), allocatable :: WORK(:), dummy_VL(:,:)
        real(wp), allocatable :: RWORK(:)
        integer :: LWORK
        
        ierr = 0
        n_tot = sector%n_channels * sector%n_basis
        eta = sector%r_center
        h_width = sector%r_max - sector%r_min
        
        ! 1. Get Raw Matrices at sector center
        allocate(P_center(sector%n_channels, sector%n_channels))
        allocate(W_center(sector%n_channels, sector%n_channels))
        allocate(P_nac_center(sector%n_channels, sector%n_channels))
        allocate(W_up_down(sector%n_channels, sector%n_channels))
        allocate(W_down_up(sector%n_channels, sector%n_channels))
        allocate(evals_ad(sector%n_channels))
        
        call get_matrices_at_eta(nelems_ad, N_per_elem_ad, x_max_ad, eta, Mj, F_val, &
                                 sector%n_channels, P_center, W_center, P_nac_center, evals_ad, ierr, &
                                 W_up_down, W_down_up)
        
        if (ierr /= 0) then
            if (present(info)) info = ierr
            return
        end if
        
        ! Store for reference
        if (allocated(sector%P_mat)) deallocate(sector%P_mat)
        allocate(sector%P_mat(sector%n_channels, sector%n_channels))
        sector%P_mat = P_center ! Note: This is the raw "Weighted Coupling" from aeigen_solver
        
        if (allocated(sector%W_mat)) deallocate(sector%W_mat)
        allocate(sector%W_mat(sector%n_channels, sector%n_channels))
        sector%W_mat = W_center

        if (allocated(sector%P_nac_mat)) deallocate(sector%P_nac_mat)
        allocate(sector%P_nac_mat(sector%n_channels, sector%n_channels))
        sector%P_nac_mat = P_nac_center

        ! 2. Construct A and B Matrices
        allocate(A_mat(sector%n_channels, sector%n_channels))
        allocate(B_mat(sector%n_channels, sector%n_channels))
        
        ! A = 2*P_nac + 1/(2eta) * (W_down_up - W_up_down)
        A_mat = 2.0_wp * P_nac_center + (1.0_wp / (2.0_wp * eta)) * (W_down_up - W_up_down)
        
        ! B = B_diag + B_coup
        ! B_diag = [F*eta/4 + epsilon_nu/eta + 1/(4*eta^2)] * delta_nu_mu
        ! Note: The energy term (E/2) is handled as a shift during propagation.
        B_mat = (0.0_wp, 0.0_wp)
        do i = 1, sector%n_channels
            term_val = (F_val * eta) / 4.0_wp + &
                       evals_ad(i) / eta + 1.0_wp / (4.0_wp * eta**2)
            B_mat(i, i) = term_val
        end do
        
        ! B_coup = -1/(4*eta^2) * (W_down_up - W_up_down) 
        !          + 1/(2*eta) * (W_down_up * P_nac - W_up_down * P_nac)
        B_mat = B_mat - (1.0_wp / (4.0_wp * eta**2)) * (W_down_up - W_up_down)
        
        ! Calculate (W_down_up - W_up_down) * P_nac
        P_center = matmul((W_down_up - W_up_down), P_nac_center)
        B_mat = B_mat + (1.0_wp / (2.0_wp * eta)) * P_center
        
        ! 3. Build Sector Hamiltonian in DVR basis
        ! We use Lobatto DVR for the sector coordinate r (mapped to [-1, 1])
        
        allocate(x_ref(sector%n_basis))
        allocate(b_weights(sector%n_basis))
        allocate(D_local(sector%n_basis, sector%n_basis))
        
        call get_lobatto_points_local(sector%n_basis, x_ref, b_weights)
        call compute_differentiation_matrix(sector%n_basis, x_ref, b_weights, D_local)
        
        allocate(H_sector(n_tot, n_tot))
        H_sector = (0.0_wp, 0.0_wp)
        
        ! Loop over channels
        do mu = 1, sector%n_channels
            do nu = 1, sector%n_channels
                
                do bi = 1, sector%n_basis
                    do bj = 1, sector%n_basis
                        idx_i = (mu - 1) * sector%n_basis + bi
                        idx_j = (nu - 1) * sector%n_basis + bj
                        
                        term_val = (0.0_wp, 0.0_wp)
                        
                        ! 1. Kinetic Energy: d^2/dr^2 (Only diagonal in channels)
                        if (mu == nu) then
                            do k = 1, sector%n_basis
                                term_val = term_val + (2.0_wp/h_width)**2 * D_local(bi, k) * D_local(k, bj)
                            end do
                        end if
                        
                        ! 2. A term: A * d/deta
                        if (abs(A_mat(mu, nu)) > 1.0e-14_wp) then
                            term_val = term_val + A_mat(mu, nu) * (2.0_wp/h_width) * D_local(bi, bj)
                        end if
                        
                        ! 3. B term: B
                        if (bi == bj) then
                            term_val = term_val + B_mat(mu, nu)
                        end if
                        
                        ! Apply -1/2 factor to convert from Operator O to Hamiltonian H
                        H_sector(idx_i, idx_j) = -0.5_wp * term_val
                    end do
                end do
            end do
        end do
        
        ! 4. Add Bloch Operator to restore Hermiticity
        ! L1s = 1/2 [ delta(r-r_max) d/dr - delta(r-r_min) d/dr ]
        ! L2s = 1/4 A [ delta(r-r_max) - delta(r-r_min) ]
        
        w_phys_L = b_weights(1) * h_width / 2.0_wp
        w_phys_R = b_weights(sector%n_basis) * h_width / 2.0_wp
        val_L = 1.0_wp / sqrt(w_phys_L)
        val_R = 1.0_wp / sqrt(w_phys_R)
        
        do mu = 1, sector%n_channels
            do nu = 1, sector%n_channels
                
                ! L2 Term (Proportional to A matrix)
                if (abs(A_mat(mu, nu)) > 1.0e-14_wp) then
                    ! At r_max (Index n_basis)
                    idx_i = (mu - 1) * sector%n_basis + sector%n_basis
                    idx_j = (nu - 1) * sector%n_basis + sector%n_basis
                    H_sector(idx_i, idx_j) = H_sector(idx_i, idx_j) + &
                                             0.25_wp * A_mat(mu, nu) * (val_R * val_R)
                                             
                    ! At r_min (Index 1) - Note the minus sign in L2 definition
                    idx_i = (mu - 1) * sector%n_basis + 1
                    idx_j = (nu - 1) * sector%n_basis + 1
                    H_sector(idx_i, idx_j) = H_sector(idx_i, idx_j) - &
                                             0.25_wp * A_mat(mu, nu) * (val_L * val_L)
                end if
                
                ! L1 Term (Derivative)
                if (mu == nu) then
                    ! At r_max
                    idx_i = (mu - 1) * sector%n_basis + sector%n_basis
                    do bj = 1, sector%n_basis
                        idx_j = (nu - 1) * sector%n_basis + bj
                        H_sector(idx_i, idx_j) = H_sector(idx_i, idx_j) + &
                                                 0.5_wp * val_R * (D_local(sector%n_basis, bj) * val_R) * (2.0_wp/h_width)
                    end do
                    
                    ! At r_min
                    idx_i = (mu - 1) * sector%n_basis + 1
                    do bj = 1, sector%n_basis
                        idx_j = (nu - 1) * sector%n_basis + bj
                        H_sector(idx_i, idx_j) = H_sector(idx_i, idx_j) - &
                                                 0.5_wp * val_L * (D_local(1, bj) * val_L) * (2.0_wp/h_width)
                    end do
                end if
                
            end do
        end do

        ! Diagonalize H_sector
        allocate(evals_sector(n_tot))
        allocate(evecs_sector(n_tot, n_tot))
        allocate(dummy_VL(1, n_tot))
        
        ! Use ZGEEV for general complex matrix
        allocate(WORK(1))
        allocate(RWORK(2*n_tot))
        
        call ZGEEV('N', 'V', n_tot, H_sector, n_tot, evals_sector, &
                   dummy_VL, 1, evecs_sector, n_tot, WORK, -1, RWORK, ierr)
                   
        LWORK = int(real(WORK(1)))
        deallocate(WORK)
        allocate(WORK(LWORK))
        
        call ZGEEV('N', 'V', n_tot, H_sector, n_tot, evals_sector, &
                   dummy_VL, 1, evecs_sector, n_tot, WORK, LWORK, RWORK, ierr)
                   
        if (ierr /= 0) then
            if (present(info)) info = ierr
            return
        end if
        
        ! Store results in sector structure
        if (allocated(sector%eigvals)) deallocate(sector%eigvals)
        allocate(sector%eigvals(n_tot))
        sector%eigvals = evals_sector
        
        ! 5. Extract Surface Amplitudes
        ! u_k(r_min) and u_k(r_max)
        ! In DVR, the first and last basis functions correspond to the boundaries.
        ! Basis 1 is at r_min, Basis N is at r_max.
        ! The amplitude is just the coefficient of the first/last basis function
        ! normalized by the weight?
        ! In Lobatto DVR, basis functions are cardinal: f_i(x_j) = delta_ij / sqrt(w_j)?
        ! Or f_i(x_j) = delta_ij.
        ! Standard FEM-DVR basis: u(x) = sum c_i f_i(x).
        ! At boundary x_1, u(x_1) = c_1 * f_1(x_1).
        ! If f_1(x_1) = 1/sqrt(w_1), then u(x_1) = c_1 / sqrt(w_1).
        ! We need to check the normalization used in `finite_element_dvr`.
        ! Usually orthonormal basis: phi_i(x_j) = delta_ij / sqrt(w_j).
        ! So amplitude = eigenvector_component / sqrt(w_boundary).
        
        allocate(sector%surface_amp_L(sector%n_channels, n_tot))
        allocate(sector%surface_amp_R(sector%n_channels, n_tot))
        
        do k = 1, n_tot
            do mu = 1, sector%n_channels
                ! Index of first basis function for channel mu
                idx_i = (mu - 1) * sector%n_basis + 1
                ! Index of last basis function for channel mu
                idx_j = (mu - 1) * sector%n_basis + sector%n_basis
                
                ! Amplitude at Left Boundary (Basis 1)
                ! We need weights. b_weights are for [-1, 1].
                ! Physical weights w_phys = b_weights * (h_width/2).
                ! Amp = Evec(idx, k) / sqrt(w_phys)
                
                sector%surface_amp_L(mu, k) = evecs_sector(idx_i, k) / sqrt(b_weights(1) * h_width / 2.0_wp)
                sector%surface_amp_R(mu, k) = evecs_sector(idx_j, k) / sqrt(b_weights(sector%n_basis) * h_width / 2.0_wp)
                
            end do
        end do
        
        deallocate(P_center, W_center, P_nac_center, W_up_down, W_down_up, evals_ad)
        deallocate(A_mat, B_mat, H_sector, x_ref, b_weights, D_local)
        deallocate(evals_sector, evecs_sector, dummy_VL)

        if (present(info)) info = ierr
    end subroutine build_and_solve_sector

    ! Helper to get Lobatto points (simple implementation for [-1, 1])
    subroutine get_lobatto_points_local(n, x, w)
        integer, intent(in) :: n
        real(wp), intent(out) :: x(:), w(:)
        integer :: ierr
        
        call get_lobatto(n, x, w, ierr)
    end subroutine get_lobatto_points_local

    !----------------------------------------------------------------
    ! propagate_step_outward
    !
    ! Propagates the R-matrix from sector boundary a to b (outward).
    ! Uses the stable R-matrix propagation formula:
    ! R(b) = R_bb - R_ba * (R(a) + R_aa)^(-1) * R_ab
    !
    ! Inputs:
    !   R_curr: R-matrix at the inner boundary (r_min)
    !   sector: The current sector data
    !   E_energy: The energy at which to propagate
    !
    ! Outputs:
    !   R_next: R-matrix at the outer boundary (r_max)
    !----------------------------------------------------------------
    subroutine propagate_step_outward(R_curr, sector, E_energy, R_next, info)
        complex(wp), intent(in) :: R_curr(:,:) ! R matrix at r_min (R(a))
        type(sector_basis_t), intent(in) :: sector
        complex(wp), intent(in) :: E_energy
        complex(wp), intent(out) :: R_next(:,:) ! R matrix at r_max (R(b))
        integer, intent(out), optional :: info
        
        integer :: ierr, n_ch, n_tot, k, i
        complex(wp), allocatable :: R_aa(:,:), R_bb(:,:), R_ab(:,:), R_ba(:,:)
        complex(wp), allocatable :: Mat_to_inv(:,:), RHS_mat(:,:)
        complex(wp), allocatable :: D_vec(:)
        complex(wp), allocatable :: U_L_scaled(:,:), U_R_scaled(:,:)
        integer, allocatable :: IPIV(:)
        complex(wp) :: denom
        
        ierr = 0
        n_ch = sector%n_channels
        n_tot = size(sector%eigvals)
        
        allocate(R_aa(n_ch, n_ch), R_bb(n_ch, n_ch))
        allocate(R_ab(n_ch, n_ch), R_ba(n_ch, n_ch))
        allocate(D_vec(n_tot))
        allocate(U_L_scaled(n_ch, n_tot), U_R_scaled(n_ch, n_tot))
        
        ! 1. Construct diagonal term 1/(e_k - E_eff) and scaled amplitudes
        ! The equation is (H_sector - E/4) Psi = 0
        ! So the effective energy shift is -E/4.
        ! denom = e_k - E/4
        
        do k = 1, n_tot
            denom = sector%eigvals(k) - E_energy / 4.0_wp
            
            if (abs(denom) > 1.0e-14_wp) then
                D_vec(k) = 0.5_wp / denom
            else
                D_vec(k) = 0.5_wp / 1.0e-14_wp ! Avoid division by zero
            end if
            
            ! Scale columns of U_L and U_R by D_vec
            do i = 1, n_ch
                U_L_scaled(i, k) = sector%surface_amp_L(i, k) * D_vec(k)
                U_R_scaled(i, k) = sector%surface_amp_R(i, k) * D_vec(k)
            end do
        end do
        
        ! 2. Compute Sector R-matrices
        ! R_xy = U_x * D * U_y^T = U_x_scaled * U_y^T
        
        ! R_aa = U_L_scaled * U_L^T
        R_aa = matmul(U_L_scaled, transpose(sector%surface_amp_L))
        
        ! R_bb = U_R_scaled * U_R^T
        R_bb = matmul(U_R_scaled, transpose(sector%surface_amp_R))
        
        ! R_ab = U_L_scaled * U_R^T
        R_ab = matmul(U_L_scaled, transpose(sector%surface_amp_R))
        
        ! R_ba = U_R_scaled * U_L^T
        R_ba = matmul(U_R_scaled, transpose(sector%surface_amp_L))
        
        ! 3. Solve linear system (R_curr + R_aa) X = R_ab
        allocate(Mat_to_inv(n_ch, n_ch))
        allocate(RHS_mat(n_ch, n_ch))
        allocate(IPIV(n_ch))
        
        Mat_to_inv = R_curr + R_aa
        RHS_mat = R_ab
        
        ! ZGESV computes the solution to system A * X = B
        ! On exit, RHS_mat contains the solution X
        call ZGESV(n_ch, n_ch, Mat_to_inv, n_ch, IPIV, RHS_mat, n_ch, ierr)
        
        if (ierr /= 0) then
            if (present(info)) info = ierr
            return
        end if
        
        ! 4. Compute R_next = R_bb - R_ba * X
        R_next = R_bb - matmul(R_ba, RHS_mat)
        
        if (present(info)) info = ierr
    end subroutine propagate_step_outward

    !----------------------------------------------------------------
    ! propagate_step_inward
    !
    ! Propagates the R-matrix from sector boundary b to a (inward).
    ! Uses the stable R-matrix propagation formula:
    ! R(a) = R_aa - R_ab * (R(b) + R_bb)^(-1) * R_ba
    !
    ! Inputs:
    !   R_curr: R-matrix at the outer boundary (r_max)
    !   sector: The current sector data
    !   E_energy: The energy at which to propagate
    !
    ! Outputs:
    !   R_prev: R-matrix at the inner boundary (r_min)
    !----------------------------------------------------------------
    subroutine propagate_step_inward(R_curr, sector, E_energy, R_prev, info)
        complex(wp), intent(in) :: R_curr(:,:) ! R matrix at r_max (R(b))
        type(sector_basis_t), intent(in) :: sector
        complex(wp), intent(in) :: E_energy
        complex(wp), intent(out) :: R_prev(:,:) ! R matrix at r_min (R(a))
        integer, intent(out), optional :: info
        
        integer :: ierr, n_ch, n_tot, k, i
        complex(wp), allocatable :: R_aa(:,:), R_bb(:,:), R_ab(:,:), R_ba(:,:)
        complex(wp), allocatable :: Mat_to_inv(:,:), RHS_mat(:,:)
        complex(wp), allocatable :: D_vec(:)
        complex(wp), allocatable :: U_L_scaled(:,:), U_R_scaled(:,:)
        integer, allocatable :: IPIV(:)
        complex(wp) :: denom
        
        ierr = 0
        n_ch = sector%n_channels
        n_tot = size(sector%eigvals)
        
        allocate(R_aa(n_ch, n_ch), R_bb(n_ch, n_ch))
        allocate(R_ab(n_ch, n_ch), R_ba(n_ch, n_ch))
        allocate(D_vec(n_tot))
        allocate(U_L_scaled(n_ch, n_tot), U_R_scaled(n_ch, n_tot))
        
        ! 1. Construct diagonal term 1/(e_k - E_eff) and scaled amplitudes
        ! denom = e_k - E/4
        do k = 1, n_tot
            denom = sector%eigvals(k) - E_energy / 4.0_wp
            
            if (abs(denom) > 1.0e-14_wp) then
                D_vec(k) = 0.5_wp / denom
            else
                D_vec(k) = 0.5_wp / 1.0e-14_wp
            end if
            
            do i = 1, n_ch
                U_L_scaled(i, k) = sector%surface_amp_L(i, k) * D_vec(k)
                U_R_scaled(i, k) = sector%surface_amp_R(i, k) * D_vec(k)
            end do
        end do
        
        ! 2. Compute Sector R-matrices
        R_aa = matmul(U_L_scaled, transpose(sector%surface_amp_L))
        R_bb = matmul(U_R_scaled, transpose(sector%surface_amp_R))
        R_ab = matmul(U_L_scaled, transpose(sector%surface_amp_R))
        R_ba = matmul(U_R_scaled, transpose(sector%surface_amp_L))
        
        ! 3. Solve linear system (R_curr + R_bb) X = R_ba
        allocate(Mat_to_inv(n_ch, n_ch))
        allocate(RHS_mat(n_ch, n_ch))
        allocate(IPIV(n_ch))
        
        Mat_to_inv = R_curr + R_bb
        RHS_mat = R_ba
        
        call ZGESV(n_ch, n_ch, Mat_to_inv, n_ch, IPIV, RHS_mat, n_ch, ierr)
        
        if (ierr /= 0) then
            if (present(info)) info = ierr
            return
        end if
        
        ! 4. Compute R_prev = R_aa - R_ab * X
        R_prev = R_aa - matmul(R_ab, RHS_mat)
        
        if (present(info)) info = ierr
    end subroutine propagate_step_inward

    !----------------------------------------------------------------
    ! export_coupling_matrices
    !
    ! Exports the coupling matrices (P_nac, W) of a sector to files.
    ! This allows visualizing the "dynamic selection rules" (checkerboard pattern).
    !
    ! Inputs:
    !   sector: The sector to export
    !   prefix: Filename prefix (e.g., "sector_050")
    !----------------------------------------------------------------
    subroutine export_coupling_matrices(sector, prefix)
        type(sector_basis_t), intent(in) :: sector
        character(len=*), intent(in) :: prefix
        
        integer :: i, j, unit_p, unit_w
        character(len=100) :: filename_p, filename_w
        
        if (.not. allocated(sector%P_nac_mat)) return
        
        filename_p = trim(prefix) // "_P_nac.dat"
        filename_w = trim(prefix) // "_W_mat.dat"
        
        unit_p = 30
        unit_w = 31
        
        open(unit=unit_p, file=trim(filename_p), status='replace')
        open(unit=unit_w, file=trim(filename_w), status='replace')
        
        write(unit_p, '(A)') '# Row Col Re(P_nac) Im(P_nac) Abs(P_nac)'
        write(unit_w, '(A)') '# Row Col Re(W) Im(W) Abs(W)'
        
        do i = 1, sector%n_channels
            do j = 1, sector%n_channels
                write(unit_p, '(I4, I4, ES16.8, ES16.8, ES16.8)') &
                    i, j, real(sector%P_nac_mat(i,j)), aimag(sector%P_nac_mat(i,j)), abs(sector%P_nac_mat(i,j))
                    
                write(unit_w, '(I4, I4, ES16.8, ES16.8, ES16.8)') &
                    i, j, real(sector%W_mat(i,j)), aimag(sector%W_mat(i,j)), abs(sector%W_mat(i,j))
            end do
            write(unit_p, *) ! Blank line for gnuplot pm3d
            write(unit_w, *)
        end do
        
        close(unit_p)
        close(unit_w)
        
    end subroutine export_coupling_matrices

end module rmatrix_propagation
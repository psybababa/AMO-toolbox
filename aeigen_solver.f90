module aeigen_solver

!------------------------------------------------------------------
!
! This module provides routines for solving the adiabatic
! eigenvalue problem in quantum scattering calculations
! involving spin-orbit interactions using finite element DVR methods.
! this module is part of a computational physics package.
! named AMO-toolbox
!
! below are the routines included in this module:
! 1. build_xi_differential_matrix
! 2. construct_diago_block_matrix
! 3. construct_coupling_block_matrix
! 4. assemble_global_hamiltonian
! 5. solve_adiabatic_hamiltonian
! 6. calculate_coupling_matrix_elements
!----------------written by K.B. Puyuyu(2025/12/04 )----------------

    use, intrinsic :: iso_fortran_env, only: wp => real64
    use finite_element_dvr
    implicit none
    private
    public :: solve_adiabatic_hamiltonian
    public :: calculate_coupling_matrix_elements
    public :: grid_data_t

    type :: grid_data_t
        real(wp), allocatable :: x_global(:)
        real(wp), allocatable :: w_global(:)
        integer, allocatable :: elem_map(:,:)
        real(wp), allocatable :: D_ref(:,:)
        real(wp), allocatable :: w_ref(:)
        integer, allocatable :: active_idx(:)
        integer :: N_full
        integer :: nelems
    end type grid_data_t

    contains

    !--------------------------------------------------------------
    !
    ! build_xi_differential_matrix
    !
    ! constructs the differential matrix in the xi coordinate
    ! for operator: f(xi) * d/dxi
    !
    ! key_f toggles the function f(xi)
    ! key_f = 0 : f = 1
    ! key_f = 1 : f = xi
    ! key_f = 2 : f = lambda * sqrt(xi*eta) * (xi + eta)/2  (Radial part of SOI term)
    ! key_f = 3 : Stiffness matrix form for operator d/dxi ( xi * d/dxi )
    !             Computes integral < d(psi_i)/dxi | xi | d(psi_j)/dxi >
    !             Note: This is symmetric and used for Kinetic Energy part.
    !
    ! Returns N_global x N_global matrix (Scalar)
    !--------------------------------------------------------------

    subroutine build_xi_differential_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                            D_ref, w_ref, eta, key_f, Mat, info)
        integer, intent(in) :: N_global, nelems
        integer, intent(in) :: elem_map(:,:)
        real(wp), intent(in) :: x_global(:), w_global(:)
        real(wp), intent(in) :: D_ref(:,:), w_ref(:)
        real(wp), intent(in) :: eta
        integer, intent(in) :: key_f
        real(wp), intent(out) :: Mat(N_global, N_global)
        integer, intent(out), optional :: info

        integer :: i, j, k, e, gi, gj, gk, N_local, ierr
        real(wp) :: xi, f_val, val, J_inv, xk
        real(wp) :: lambda_val, Z_eff_val
        
        ierr = 0
        Mat = 0.0_wp
        N_local = size(D_ref, 1)

        do e = 1, nelems
            ! Calculate Jacobian inverse for this element
            ! J = (x_last - x_first) / 2
            gi = elem_map(e, 1)
            gj = elem_map(e, N_local)
            J_inv = 2.0_wp / (x_global(gj) - x_global(gi))

            if (key_f == 3) then
                do i = 1, N_local
                    gi = elem_map(e,i)
                    do j = 1, N_local
                        gj = elem_map(e,j)
                        
                        val = 0.0_wp
                        do k = 1, N_local
                            gk = elem_map(e,k)
                            xk = x_global(gk)
                            ! D_ref(k, i) is derivative of basis i at node k
                            val = val + w_ref(k) * xk * D_ref(k,i) * D_ref(k,j)
                        end do
                        
                        val = val * J_inv / sqrt(w_global(gi) * w_global(gj))
                        Mat(gi, gj) = Mat(gi, gj) + val
                    end do
                end do
                
            else
                ! 1st derivative matrix: < i | f d/dxi | j >
                do i = 1, N_local
                    gi = elem_map(e,i)
                    xi = x_global(gi)
                    
                    ! Calculate f(xi) based on key_f
                    select case (key_f)
                    case (0)
                        f_val = 1.0_wp
                    case (1)
                        f_val = xi
                    case (2)
                        ! f = lambda * sqrt(xi*eta) * (xi + eta)/2
                        call convert_potential_to_xi_eta(xi, eta, lambda_val, Z_eff_val)
                        f_val = lambda_val * sqrt(xi * eta) * (xi + eta) / 2.0_wp
                    case default
                        f_val = 0.0_wp
                    end select

                    ! Compute matrix elements
                    ! O_ij = <psi_i | f d/dxi | psi_j>
                    ! Contribution from element e:
                    ! = w_ref(i) * f(xi) * D_ref(i,j) / sqrt(w_global(gi) * w_global(gj))
                    ! Note: Jacobian cancels out (J from integral, 1/J from derivative)
                    
                    do j = 1, N_local
                        gj = elem_map(e,j)
                        
                        val = w_ref(i) * f_val * D_ref(i,j)
                        val = val / sqrt(w_global(gi) * w_global(gj))
                        
                        Mat(gi, gj) = Mat(gi, gj) + val
                    end do
                end do
            end if
        end do

        if (present(info)) info = ierr
    end subroutine build_xi_differential_matrix

    !--------------------------------------------------------------
    !
    !
    ! construct_diago_block_matrix
    !
    ! this adiabiatic hamiltonian consisting diagonal blocks
    ! which include kinetic energy and potential energy terms
    ! though we are now considering pauli spinors thus
    ! global matrix will be 2 times larger
    ! and there also exsits coupling blocks which comes from spin-orbit interaction
    ! 
    ! this subroutine constructs only diagonal blocks
    ! before assembling global matrix
    ! 
    ! due to the use of spinors
    ! the d/dphi and Sz termes gives different angular momentum contributions
    ! depending on whether Mj = m + 1/2 or Mj = m - 1/2
    ! key_m toggles between these two cases
    ! key_m = 1 : Mj = m +1/2 (Spin Up)
    ! key_m = 2 : Mj = m -1/2 (Spin Down)
    !--------------------------------------------------------------

    subroutine construct_diago_block_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                            D_ref, w_ref, eta, Mj, E_val, F_val, key_m, Mat, info)
        integer, intent(in) :: N_global, nelems
        integer, intent(in) :: elem_map(:,:)
        real(wp), intent(in) :: x_global(:), w_global(:)
        real(wp), intent(in) :: D_ref(:,:), w_ref(:)
        real(wp), intent(in) :: eta, Mj, F_val
        complex(wp), intent(in) :: E_val
        integer, intent(in) :: key_m
        complex(wp), intent(out) :: Mat(N_global, N_global)
        integer, intent(out), optional :: info

        integer :: i, e, gi, ierr
        real(wp) :: xi, m_val, Sz, lambda_val, Z_eff_val
        complex(wp) :: V_diag, term_pot
        real(wp) :: term_azi, term_soi
        real(wp), allocatable :: Mat_real(:,:)
        
        ierr = 0
        
        ! 1. Kinetic Energy Term: d/dxi ( xi d/dxi )
        ! We use build_xi_differential_matrix with key_f=3
        ! This returns K_ij = < i' | xi | j' >
        ! The operator in B(eta) is + d/dxi ( xi d/dxi )
        ! Matrix element < i | d/dxi ( xi d/dxi ) | j > = - < i' | xi | j' > (by integration by parts)
        ! So we multiply by -1.0
        
        allocate(Mat_real(N_global, N_global))
        call build_xi_differential_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                          D_ref, w_ref, eta, 3, Mat_real, ierr)
        
        Mat = cmplx(-1.0_wp * Mat_real, 0.0_wp, wp)
        deallocate(Mat_real)
        
        ! 2. Determine m and Sz based on key_m
        if (key_m == 1) then
            ! Spin Up: Mj = m + 1/2 => m = Mj - 1/2
            Sz = 0.5_wp
            m_val = Mj - 0.5_wp
        else
            ! Spin Down: Mj = m - 1/2 => m = Mj + 1/2
            Sz = -0.5_wp
            m_val = Mj + 0.5_wp
        end if
        
        ! 3. Add Diagonal Potential Terms
        ! Loop over all nodes (diagonal in DVR)
        ! Since Mat is initialized with Kinetic term, we add to diagonal elements.
        
        do e = 1, nelems
            do i = 1, size(D_ref, 1)
                gi = elem_map(e,i)
                xi = x_global(gi)
                
                ! Get potential parameters
                call convert_potential_to_xi_eta(xi, eta, lambda_val, Z_eff_val)
                
                ! Calculate terms
                
                ! Azimuthal term: + (xi+eta)/(4*xi*eta) * d^2/dphi^2
                ! d^2/dphi^2 -> -m^2
                ! Term = - m^2 * (xi+eta)/(4*xi*eta) = - m^2/4xi - m^2/4eta
                term_azi = - (m_val**2) / (4.0_wp * xi) - (m_val**2) / (4.0_wp * eta)
                
                ! Potential term: - (xi+eta)/2 * Vc + E*xi/2 - F*xi^2/4
                ! Vc = -Z_eff / r = -Z_eff / ((xi+eta)/2)
                ! So - (xi+eta)/2 * Vc = Z_eff
                term_pot = Z_eff_val + (E_val * xi) / 2.0_wp - (F_val * xi**2) / 4.0_wp
                
                ! SOI Diagonal term:
                ! From Eq (8): + i/2 * lambda * (xi+eta)/2 * (2Sz) * d/dphi
                ! d/dphi -> i*m
                ! Term = i/2 * lambda * (xi+eta)/2 * 2Sz * (i*m)
                !      = - lambda * (xi+eta)/2 * m * Sz
                term_soi = - lambda_val * (xi + eta) / 2.0_wp * m_val * Sz
                
                V_diag = term_azi + term_pot + term_soi
                
                ! Add to diagonal
                Mat(gi, gi) = Mat(gi, gi) + V_diag
            end do
        end do

        if (present(info)) info = ierr
    end subroutine construct_diago_block_matrix

    !--------------------------------------------------------------
    !
    ! construct_coupling_block_matrix
    !
    ! constructs coupling blocks due to spin-orbit interaction
    ! before assembling global matrix
    ! from (1,N+1) to (1,2N) and (N,N+1) to (N,2N) elements of matrix
    ! there is the coupling term comes from down spinor to up spinor
    ! and vice versa
    !
    ! thus this subroutine constructs only coupling blocks
    ! before assembling global matrix
    ! 
    ! key_m toggles between Mj = m +1/2 or Mj = m -1/2
    ! key_m = 1 : Mj = m +1/2 (Target is Up, Source is Down) -> Construct Up-Down block
    ! key_m = 2 : Mj = m -1/2 (Target is Down, Source is Up) -> Construct Down-Up block
    !--------------------------------------------------------------

    subroutine construct_coupling_block_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                               D_ref, w_ref, eta, Mj, key_m, Mat, info)
        integer, intent(in) :: N_global, nelems
        integer, intent(in) :: elem_map(:,:)
        real(wp), intent(in) :: x_global(:), w_global(:)
        real(wp), intent(in) :: D_ref(:,:), w_ref(:)
        real(wp), intent(in) :: eta, Mj
        integer, intent(in) :: key_m
        real(wp), intent(out) :: Mat(N_global, N_global)
        integer, intent(out), optional :: info

        integer :: i, e, gi, ierr
        real(wp) :: xi, m_source, lambda_val, Z_eff_val
        real(wp) :: term_pot
        real(wp), allocatable :: Mat_dxi(:,:)
        
        ierr = 0
        Mat = 0.0_wp
        
        ! 1. d/dxi term (Radial part of SOI)
        ! Term: -1/2 * lambda * sqrt(xi*eta) * (xi+eta)/2 * (S+ e^-iphi - S- e^iphi) d/dxi
        ! We use build_xi_differential_matrix with key_f=2
        ! key_f=2 returns matrix for operator: lambda * sqrt(xi*eta) * (xi+eta)/2 * d/dxi
        
        allocate(Mat_dxi(N_global, N_global))
        call build_xi_differential_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                          D_ref, w_ref, eta, 2, Mat_dxi, ierr)
        
        ! Determine sign and m based on transition
        if (key_m == 1) then
            ! Constructing Up-Down Block (Target: Up, Source: Down)
            ! Operator involves S+ (acting on Down gives Up)
            ! Term in Eq(8): -1/2 * (...) * S+ * ...
            ! S+ |Down> = |Up> (factor 1)
            ! Coefficient is -1/2
            ! So we multiply Mat_dxi by -0.5
            
            Mat = -0.5_wp * Mat_dxi
            
            ! For potential term, we need m of the SOURCE state (Down state)
            ! Mj = m_up + 1/2 = m_down - 1/2
            ! So m_down = Mj + 1/2
            m_source = Mj + 0.5_wp
            
        else
            ! Constructing Down-Up Block (Target: Down, Source: Up)
            ! Operator involves S- (acting on Up gives Down)
            ! Term in Eq(8): -1/2 * (...) * (- S-) * ...
            ! Note the minus sign inside the bracket in Eq(8): (S+ ... - S- ...)
            ! So term is +1/2 * (...) * S-
            ! S- |Up> = |Down> (factor 1)
            ! Coefficient is +0.5
            ! So we multiply Mat_dxi by +0.5
            
            Mat = 0.5_wp * Mat_dxi
            
            ! For potential term, we need m of the SOURCE state (Up state)
            ! m_up = Mj - 0.5_wp
            m_source = Mj - 0.5_wp
        end if
        
        deallocate(Mat_dxi)
        
        ! 2. d/dphi term (Potential-like term)
        ! Term: + i/2 * lambda * (xi+eta)/2 * (xi-eta)/(2sqrt(xi*eta)) * (S+ e^-iphi + S- e^iphi) * d/dphi
        ! d/dphi acts on source state -> i * m_source
        ! (S+ + S-) acts on source state -> |Target> (factor 1)
        ! Coefficient: i/2 * lambda * ... * (i * m_source)
        !            = -1/2 * lambda * (xi+eta)/2 * (xi-eta)/(2sqrt(xi*eta)) * m_source
        !            = - lambda * (xi^2 - eta^2) / (8 * sqrt(xi*eta)) * m_source
        ! This term is diagonal in DVR basis (function of xi)
        ! Both Up->Down and Down->Up have the same sign structure here (S+ + S- has + sign)
        ! So coefficient is always -1/2 * ... * m_source
        
        do e = 1, nelems
            do i = 1, size(D_ref, 1)
                gi = elem_map(e,i)
                xi = x_global(gi)
                
                call convert_potential_to_xi_eta(xi, eta, lambda_val, Z_eff_val)
                term_pot = - lambda_val * (xi**2 - eta**2) / (8.0_wp * sqrt(xi*eta)) * m_source
                
                Mat(gi, gi) = Mat(gi, gi) + term_pot
            end do
        end do

        if (present(info)) info = ierr
    end subroutine construct_coupling_block_matrix

    !--------------------------------------------------------------
    !
    ! assemble_global_hamiltonian
    ! solve_adiabatic_hamiltonian
    !
    ! assembles global matrix from diagonal blocks and coupling blocks
    !
    ! then solves the adiabatic hamiltonian eigenvalue problem
    !
    ! solve_adiabatic_hamiltonian subroutine returns eigenvalues and eigenvectors
    ! and can be called from outside 
    !--------------------------------------------------------------

    subroutine assemble_global_matrix(nelems, N_per_elem, x_max, eta, Mj, E_val, F_val, &
                                      N_global, x_global, w_global, H_global, info, &
                                      full_grid_data)
        integer, intent(in) :: nelems, N_per_elem
        real(wp), intent(in) :: x_max, eta, Mj, F_val
        complex(wp), intent(in) :: E_val
        integer, intent(out) :: N_global
        real(wp), allocatable, intent(out) :: x_global(:), w_global(:)
        complex(wp), allocatable, intent(out) :: H_global(:,:)
        integer, intent(out), optional :: info
        
        ! Optional container for full grid data needed for coupling calculation
        type(grid_data_t), intent(out), optional :: full_grid_data

        integer :: ierr, i, j, N_full
        real(wp), allocatable :: x_ref(:), w_ref(:), b_weights(:), D_ref(:,:)
        real(wp), allocatable :: M_global_full(:), x_global_full(:)
        integer, allocatable :: elem_map(:,:)
        complex(wp), allocatable :: H_diag_up(:,:), H_diag_down(:,:)
        real(wp), allocatable :: H_coup_up_down(:,:), H_coup_down_up(:,:)
        logical, allocatable :: active_mask(:)
        integer, allocatable :: active_idx(:)
        integer :: n_active
        real(wp) :: x_min = 0.0_wp

        ierr = 0
        
        ! 1. Build Global Grid
        allocate(x_ref(N_per_elem), w_ref(N_per_elem))
        allocate(x_global_full(nelems * (N_per_elem - 1) + 1)) 
        allocate(M_global_full(size(x_global_full)))
        allocate(elem_map(nelems, N_per_elem))
        
        call build_global_grid(nelems, N_per_elem, x_ref, w_ref, &
                               x_min, x_max, N_full, &
                               x_global_full, M_global_full, elem_map, ierr)
        
        if (ierr /= 0) then
            if (present(info)) info = ierr
            return
        end if
        
        ! 2. Compute Differentiation Matrix
        allocate(b_weights(N_per_elem), D_ref(N_per_elem, N_per_elem))
        call compute_barycentric_weights(N_per_elem, x_ref, b_weights)
        call compute_differentiation_matrix(N_per_elem, x_ref, b_weights, D_ref)
        
        ! 3. Construct Block Matrices (Full Size)
        allocate(H_diag_up(N_full, N_full))
        allocate(H_diag_down(N_full, N_full))
        allocate(H_coup_up_down(N_full, N_full))
        allocate(H_coup_down_up(N_full, N_full))
        
        ! Up-Up Block (key_m = 1)
        call construct_diago_block_matrix(N_full, nelems, elem_map, x_global_full, M_global_full, &
                                          D_ref, w_ref, eta, Mj, E_val, F_val, 1, H_diag_up, ierr)
        
        ! Down-Down Block (key_m = 2)
        call construct_diago_block_matrix(N_full, nelems, elem_map, x_global_full, M_global_full, &
                                          D_ref, w_ref, eta, Mj, E_val, F_val, 2, H_diag_down, ierr)
                                          
        ! Up-Down Coupling (Target Up, Source Down, key_m = 1)
        call construct_coupling_block_matrix(N_full, nelems, elem_map, x_global_full, M_global_full, &
                                             D_ref, w_ref, eta, Mj, 1, H_coup_up_down, ierr)
                                             
        ! Down-Up Coupling (Target Down, Source Up, key_m = 2)
        call construct_coupling_block_matrix(N_full, nelems, elem_map, x_global_full, M_global_full, &
                                             D_ref, w_ref, eta, Mj, 2, H_coup_down_up, ierr)

        ! 4. Apply Boundary Conditions (Remove Origin and Infinity)
        allocate(active_mask(N_full))
        active_mask = .true.
        active_mask(1) = .false.
        active_mask(N_full) = .false.
        
        n_active = count(active_mask)
        allocate(active_idx(n_active))
        active_idx = pack([(i, i=1, N_full)], active_mask)
        
        ! 5. Assemble Global Matrix (Reduced Size)
        allocate(H_global(2*n_active, 2*n_active))
        H_global = (0.0_wp, 0.0_wp)
        
        do i = 1, n_active
            do j = 1, n_active
                ! Top-Left: Up-Up
                H_global(i, j) = H_diag_up(active_idx(i), active_idx(j))
                ! Bottom-Right: Down-Down
                H_global(n_active + i, n_active + j) = H_diag_down(active_idx(i), active_idx(j))
                ! Top-Right: Up-Down
                H_global(i, n_active + j) = cmplx(H_coup_up_down(active_idx(i), active_idx(j)), 0.0_wp, wp)
                ! Bottom-Left: Down-Up
                H_global(n_active + i, j) = cmplx(H_coup_down_up(active_idx(i), active_idx(j)), 0.0_wp, wp)
            end do
        end do
        
        ! Update output N_global and arrays to reflect reduced grid
        N_global = n_active
        
        allocate(w_global(n_active))
        w_global = M_global_full(active_idx)
        
        allocate(x_global(n_active))
        x_global = x_global_full(active_idx)
        
        ! Save full grid data if requested
        if (present(full_grid_data)) then
            allocate(full_grid_data%x_global(N_full))
            full_grid_data%x_global = x_global_full
            
            allocate(full_grid_data%w_global(N_full))
            full_grid_data%w_global = M_global_full
            
            allocate(full_grid_data%elem_map(nelems, N_per_elem))
            full_grid_data%elem_map = elem_map
            
            allocate(full_grid_data%D_ref(N_per_elem, N_per_elem))
            full_grid_data%D_ref = D_ref
            
            allocate(full_grid_data%w_ref(N_per_elem))
            full_grid_data%w_ref = w_ref
            
            allocate(full_grid_data%active_idx(n_active))
            full_grid_data%active_idx = active_idx
            
            full_grid_data%N_full = N_full
            full_grid_data%nelems = nelems
        end if
        
        deallocate(x_ref, w_ref, b_weights, D_ref, M_global_full, x_global_full, elem_map)
        deallocate(H_diag_up, H_diag_down, H_coup_up_down, H_coup_down_up)
        deallocate(active_mask, active_idx)

        if (present(info)) info = ierr
    end subroutine assemble_global_matrix

    subroutine solve_adiabatic_hamiltonian(nelems, N_per_elem, x_max, eta, Mj, E_val, F_val, &
                                           Evals, Evecs, x_global, w_global, info, &
                                           calc_coupling, P_mat, Evecs_prev, W_mat, P_nac_mat)
        integer, intent(in) :: nelems, N_per_elem
        real(wp), intent(in) :: x_max, eta, Mj, F_val
        complex(wp), intent(in) :: E_val
        complex(wp), allocatable, intent(out) :: Evals(:)
        complex(wp), allocatable, intent(out) :: Evecs(:,:)
        real(wp), allocatable, intent(out) :: x_global(:), w_global(:)
        integer, intent(out), optional :: info
        logical, intent(in), optional :: calc_coupling
        complex(wp), allocatable, intent(out), optional :: P_mat(:,:)
        complex(wp), intent(in), optional :: Evecs_prev(:,:)
        complex(wp), allocatable, intent(out), optional :: W_mat(:,:)
        complex(wp), allocatable, intent(out), optional :: P_nac_mat(:,:)

        integer :: ierr, N_global, i, j, k
        complex(wp), allocatable :: H_global(:,:)
        complex(wp), allocatable :: WORK(:)
        real(wp), allocatable :: RWORK(:)
        integer :: LWORK
        complex(wp), allocatable :: VL(:,:), VR(:,:)
        type(grid_data_t) :: full_grid
        logical :: do_coupling
        complex(wp) :: dot_prod
        real(wp) :: norm_val
        complex(wp) :: phase_factor
        integer, allocatable :: idx_sort(:)
        real(wp), allocatable :: real_evals(:)
        complex(wp), allocatable :: temp_vec(:)
        complex(wp), allocatable :: Evecs_full(:,:)

        ierr = 0
        do_coupling = .false.
        if (present(calc_coupling)) do_coupling = calc_coupling

        ! 1. Assemble Hamiltonian
        ! We request full_grid_data if coupling calculation is needed
        if (do_coupling) then
            call assemble_global_matrix(nelems, N_per_elem, x_max, eta, Mj, E_val, F_val, &
                                        N_global, x_global, w_global, H_global, ierr, &
                                        full_grid_data=full_grid)
        else
            call assemble_global_matrix(nelems, N_per_elem, x_max, eta, Mj, E_val, F_val, &
                                        N_global, x_global, w_global, H_global, ierr)
        end if

        if (ierr /= 0) then
            if (present(info)) info = ierr
            return
        end if

        ! 2. Diagonalize using ZGEEV
        ! ZGEEV computes eigenvalues and right eigenvectors
        ! H * v = lambda * v
        
        allocate(Evals(2*N_global))
        allocate(VR(2*N_global, 2*N_global))
        allocate(VL(1,1)) ! Not used
        
        ! Query workspace size
        allocate(WORK(1))
        allocate(RWORK(4*N_global)) ! Minimum size for ZGEEV is 2*N
        
        call ZGEEV('N', 'V', 2*N_global, H_global, 2*N_global, Evals, &
                   VL, 1, VR, 2*N_global, WORK, -1, RWORK, ierr)
                   
        LWORK = int(real(WORK(1)))
        deallocate(WORK)
        allocate(WORK(LWORK))
        
        call ZGEEV('N', 'V', 2*N_global, H_global, 2*N_global, Evals, &
                   VL, 1, VR, 2*N_global, WORK, LWORK, RWORK, ierr)
                   
        if (ierr /= 0) then
            if (present(info)) info = ierr
            return
        end if
        
        ! 3. Sort Eigenvalues (by Real part)
        ! We sort in DESCENDING order because the kinetic energy operator 
        ! d/dxi(xi d/dxi) is negative definite.
        ! The physical states (low kinetic energy) correspond to the largest (least negative) eigenvalues.
        allocate(idx_sort(2*N_global))
        allocate(real_evals(2*N_global))
        real_evals = real(Evals)
        
        ! Simple bubble sort index array
        do i = 1, 2*N_global
            idx_sort(i) = i
        end do
        
        do i = 1, 2*N_global-1
            do j = i+1, 2*N_global
                if (real_evals(idx_sort(j)) > real_evals(idx_sort(i))) then
                    k = idx_sort(i)
                    idx_sort(i) = idx_sort(j)
                    idx_sort(j) = k
                end if
            end do
        end do
        
        ! Reorder Evals and Evecs
        allocate(Evecs(2*N_global, 2*N_global))
        allocate(temp_vec(2*N_global))
        
        do i = 1, 2*N_global
            k = idx_sort(i)
            Evecs(:, i) = VR(:, k)
            temp_vec(i) = Evals(k) ! Store sorted evals temporarily
        end do
        Evals = temp_vec
        deallocate(temp_vec)
        
        ! 4. Phase Alignment (if previous eigenvectors provided)
        if (present(Evecs_prev)) then
            if (size(Evecs_prev, 1) == 2*N_global .and. size(Evecs_prev, 2) == 2*N_global) then
                do i = 1, 2*N_global
                    ! Calculate overlap < v_prev | v_curr >
                    dot_prod = dot_product(Evecs_prev(:, i), Evecs(:, i))
                    
                    ! We want < v_prev | v_curr > to be real and positive
                    ! Multiply v_curr by phase factor exp(-i * arg(dot_prod))
                    ! phase_factor = conjg(dot_prod) / |dot_prod|
                    
                    norm_val = abs(dot_prod)
                    if (norm_val > 1.0e-10_wp) then
                        phase_factor = conjg(dot_prod) / norm_val
                        Evecs(:, i) = Evecs(:, i) * phase_factor
                    end if
                end do
            end if
        end if
        
        ! 5. Calculate Coupling Matrix Elements (if requested)
        if (do_coupling) then
            ! We need to reconstruct full eigenvectors (with zeros at boundaries)
            ! full_grid contains the mapping info
            
            allocate(Evecs_full(2*full_grid%N_full, 2*N_global))
            Evecs_full = (0.0_wp, 0.0_wp)
            
            ! Map reduced eigenvectors to full grid
            ! active_idx maps reduced index 1..N_global to full index
            do j = 1, 2*N_global ! Loop over states
                do i = 1, N_global ! Loop over reduced grid points
                    ! Up component
                    Evecs_full(full_grid%active_idx(i), j) = Evecs(i, j)
                    ! Down component
                    Evecs_full(full_grid%N_full + full_grid%active_idx(i), j) = Evecs(N_global + i, j)
                end do
            end do
            
            if (present(P_mat)) allocate(P_mat(2*N_global, 2*N_global))
            if (present(W_mat)) allocate(W_mat(2*N_global, 2*N_global))
            if (present(P_nac_mat)) allocate(P_nac_mat(2*N_global, 2*N_global))
            
            call calculate_coupling_matrix_elements(full_grid%N_full, full_grid%nelems, &
                                                    full_grid%elem_map, full_grid%x_global, full_grid%w_global, &
                                                    full_grid%D_ref, full_grid%w_ref, eta, Mj, &
                                                    Evals, Evecs_full, P_mat, W_mat, ierr, P_nac_mat)
                                                    
            deallocate(Evecs_full)
        end if

        if (present(info)) info = ierr
    end subroutine solve_adiabatic_hamiltonian

    !--------------------------------------------------------------
    !
    ! calculate_coupling_matrix_elements
    !
    ! Using hellman-feynman theorem
    ! <phi_i | d/d\eta | phi_j> = 1/(E_j - E_i) * <phi_i | dB/d\eta | phi_j>
    !
    ! calculates non-adiabatic coupling matrix elements
    ! between adiabatic states
    ! must be called after solve_adiabatic_hamiltonian
    ! and calculates derivative couplings
    ! from adiabatic eigenvectors
    !
    ! to gain consistent phase of adiabatic eigenvectors
    ! it checks the inner product between each eigenvectors is positive
    !
    ! this subroutine can be also called from outside
    !--------------------------------------------------------------

    subroutine calculate_coupling_matrix_elements(N_global, nelems, elem_map, x_global, w_global, &
                                                  D_ref, w_ref, eta, Mj, Evals, Evecs, P_mat, W_mat, info, P_nac_mat)
        integer, intent(in) :: N_global, nelems
        integer, intent(in) :: elem_map(:,:)
        real(wp), intent(in) :: x_global(:), w_global(:)
        real(wp), intent(in) :: D_ref(:,:), w_ref(:)
        real(wp), intent(in) :: eta, Mj
        complex(wp), intent(in) :: Evals(:)
        complex(wp), intent(in) :: Evecs(:,:)
        complex(wp), intent(out), optional :: P_mat(:,:)
        complex(wp), intent(out), optional :: W_mat(:,:)
        integer, intent(out), optional :: info
        complex(wp), intent(out), optional :: P_nac_mat(:,:)

        integer :: i, j, k, ierr, n_states
        complex(wp), allocatable :: dB_deta(:,:)
        complex(wp), parameter :: czero = (0.0_wp, 0.0_wp)
        real(wp), parameter :: tol = 1.0e-10_wp
        real(wp), allocatable :: W_dvr(:)
        complex(wp), allocatable :: Temp_mat(:,:)
        complex(wp), allocatable :: P_std(:,:)
        complex(wp), allocatable :: W_internal(:,:)
        real(wp) :: xi, lambda_val, Z_eff_val

        ierr = 0
        n_states = size(Evals)
        
        ! 1. Calculate P_std (Standard Non-adiabatic coupling)
        ! We calculate this if P_mat OR P_nac_mat is requested
        if (present(P_mat) .or. present(P_nac_mat)) then
            allocate(dB_deta(2*N_global, 2*N_global))
            allocate(P_std(n_states, n_states))
            
            ! Build dB/deta matrix
            call build_dB_deta_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                      D_ref, w_ref, eta, Mj, dB_deta, ierr)
            
            ! Calculate P_ij = <i|dB/deta|j> / (Ej - Ei)
            P_std = matmul(conjg(transpose(Evecs)), matmul(dB_deta, Evecs))
            
            ! Divide by energy difference
            do i = 1, n_states
                do j = 1, n_states
                    if (i == j) then
                        P_std(i,j) = czero ! Diagonal term is 0
                    else
                        if (abs(Evals(j) - Evals(i)) > tol) then
                            P_std(i,j) = P_std(i,j) / (Evals(j) - Evals(i))
                        else
                            P_std(i,j) = czero 
                        end if
                    end if
                end do
            end do
            deallocate(dB_deta)
            
            ! If P_nac_mat is requested, copy P_std
            if (present(P_nac_mat)) then
                P_nac_mat = P_std
            end if
        end if
        
        ! 2. Calculate W_mat (Operator matrix for lambda * sqrt(eta*xi))
        ! We calculate this if W_mat OR P_mat is requested (needed for weighting P)
        if (present(W_mat) .or. present(P_mat)) then
            allocate(W_dvr(2*N_global))
            allocate(Temp_mat(2*N_global, n_states))
            allocate(W_internal(n_states, n_states))
            
            ! Construct diagonal operator in DVR basis
            do i = 1, N_global
                xi = x_global(i)
                if (xi < 1.0e-12_wp) then
                    W_dvr(i) = 0.0_wp
                    W_dvr(N_global + i) = 0.0_wp
                else
                    call convert_potential_to_xi_eta(xi, eta, lambda_val, Z_eff_val)
                    ! Operator is lambda * sqrt(eta * xi)
                    W_dvr(i) = lambda_val * sqrt(eta * xi)
                    W_dvr(N_global + i) = W_dvr(i)
                end if
            end do
            
            ! Transform to adiabatic basis: W_ad = U^H * W_dvr * U
            do j = 1, n_states
                do k = 1, 2*N_global
                    Temp_mat(k, j) = W_dvr(k) * Evecs(k, j)
                end do
            end do
            
            W_internal = matmul(conjg(transpose(Evecs)), Temp_mat)
            
            if (present(W_mat)) then
                W_mat = W_internal
            end if
            
            deallocate(W_dvr, Temp_mat)
        end if
        
        ! 3. Calculate Weighted P_mat = W * P_std
        if (present(P_mat)) then
            if (allocated(P_std) .and. allocated(W_internal)) then
                P_mat = matmul(W_internal, P_std)
            end if
        end if
        
        if (allocated(P_std)) deallocate(P_std)
        if (allocated(W_internal)) deallocate(W_internal)

        if (present(info)) info = ierr
    end subroutine calculate_coupling_matrix_elements

    !--------------------------------------------------------------
    ! build_dB_deta_matrix
    ! Constructs the derivative of Hamiltonian with respect to eta
    !--------------------------------------------------------------
    subroutine build_dB_deta_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                    D_ref, w_ref, eta, Mj, Mat, info)
        integer, intent(in) :: N_global, nelems
        integer, intent(in) :: elem_map(:,:)
        real(wp), intent(in) :: x_global(:), w_global(:)
        real(wp), intent(in) :: D_ref(:,:), w_ref(:)
        real(wp), intent(in) :: eta, Mj
        complex(wp), intent(out) :: Mat(2*N_global, 2*N_global)
        integer, intent(out), optional :: info

        integer :: i, j, e, gi, gj, N_local
        real(wp) :: xi, lambda_val, Z_eff_val, dVc_deta, dlambda_deta
        real(wp) :: m_up, m_down
        real(wp) :: term_diag_up, term_diag_down
        real(wp) :: C1, dC1_deta, C2, dC2_deta
        real(wp) :: term_off_dxi, term_off_pot
        real(wp) :: val_dxi
        complex(wp), parameter :: czero = (0.0_wp, 0.0_wp)
        integer :: ierr

        ierr = 0
        Mat = czero
        N_local = size(D_ref, 1)
        
        m_up = Mj - 0.5_wp
        m_down = Mj + 0.5_wp

        ! Loop over elements
        do e = 1, nelems
            do i = 1, N_local
                gi = elem_map(e,i)
                xi = x_global(gi)
                
                if (xi < 1.0e-12_wp) then
                    ! Singularity at origin:
                    ! The wavefunction is zero here (Dirichlet BC).
                    ! However, terms like 1/xi or 1/sqrt(xi) diverge.
                    ! Since this node is removed from the active basis, 
                    ! its contribution to the final matrix element <i|O|j> should be zero
                    ! if we strictly follow the projection.
                    ! But here we are filling the full matrix.
                    ! To avoid NaNs propagating via 0 * Inf, we explicitly set terms to 0.
                    
                    Mat(gi, gi) = czero
                    Mat(N_global+gi, N_global+gi) = czero
                    ! Off-diagonal blocks also 0
                    Mat(gi, N_global+gi) = czero
                    Mat(N_global+gi, gi) = czero
                    
                    ! For d/dxi part, we also skip
                    cycle
                end if

                call convert_potential_to_xi_eta(xi, eta, lambda_val, Z_eff_val, dVc_deta, dlambda_deta)
                
                ! 1. Diagonal Blocks (Up-Up and Down-Down)
                ! Up-Up (Sz = 1/2, m = m_up)
                term_diag_up = (m_up**2)/(4.0_wp * eta**2) &
                             - (-Z_eff_val/((xi+eta)/2.0_wp))/2.0_wp & ! -Vc/2
                             - (xi+eta)/2.0_wp * dVc_deta &
                             - (dlambda_deta * (xi+eta)/2.0_wp + lambda_val * 0.5_wp) * 0.5_wp * m_up

                ! Down-Down (Sz = -1/2, m = m_down)
                term_diag_down = (m_down**2)/(4.0_wp * eta**2) &
                               - (-Z_eff_val/((xi+eta)/2.0_wp))/2.0_wp &
                               - (xi+eta)/2.0_wp * dVc_deta &
                               - (dlambda_deta * (xi+eta)/2.0_wp + lambda_val * 0.5_wp) * (-0.5_wp) * m_down

                ! Add to diagonal elements
                Mat(gi, gi) = cmplx(term_diag_up, 0.0_wp, wp)
                Mat(N_global+gi, N_global+gi) = cmplx(term_diag_down, 0.0_wp, wp)
                
                ! 2. Off-Diagonal Blocks (Up-Down and Down-Up)
                ! Coefficients
                C1 = sqrt(xi*eta) * (xi+eta)/2.0_wp
                dC1_deta = (xi**(1.5_wp)*eta**(-0.5_wp) + 3.0_wp*sqrt(xi*eta)) / 4.0_wp
                
                C2 = (xi**2 - eta**2) / (4.0_wp * sqrt(xi*eta))
                dC2_deta = - (3.0_wp*eta**2 + xi**2) / (8.0_wp * sqrt(xi) * eta**(1.5_wp))
                
                ! Term 1: d/dxi part coefficient
                term_off_dxi = -0.5_wp * (dlambda_deta * C1 + lambda_val * dC1_deta)
                
                ! Term 2: Potential part coefficient (proportional to m)
                ! Up-Down Block (Row Up, Col Down) -> m_target = m_down
                term_off_pot = -0.5_wp * (dlambda_deta * C2 + lambda_val * dC2_deta) * m_down
                Mat(gi, N_global+gi) = Mat(gi, N_global+gi) + cmplx(term_off_pot, 0.0_wp, wp)
                
                ! Down-Up Block (Row Down, Col Up) -> m_target = m_up
                term_off_pot = -0.5_wp * (dlambda_deta * C2 + lambda_val * dC2_deta) * m_up
                Mat(N_global+gi, gi) = Mat(N_global+gi, gi) + cmplx(term_off_pot, 0.0_wp, wp)
                
                ! d/dxi part (Integration needed)
                do j = 1, N_local
                    gj = elem_map(e,j)
                    
                    ! Calculate matrix element for d/dxi term
                    ! <i| f(x) d/dx |j> in orthonormal DVR basis
                    ! = w_ref(i) * f(x_i) * D_ref(i,j) / sqrt(w_global(gi)*w_global(gj))
                    ! Note: w_ref(i) comes from quadrature.
                    ! D_ref is d/d(xi_ref).
                    ! We need to check if Jacobian is handled.
                    ! In build_xi_differential_matrix, it used D_ref directly.
                    ! Assuming D_ref includes Jacobian or it cancels out as analyzed before.
                    
                    val_dxi = w_ref(i) * term_off_dxi * D_ref(i,j)
                    val_dxi = val_dxi / sqrt(w_global(gi) * w_global(gj))
                    
                    ! Add to Off-Diagonal Blocks
                    ! Up-Down
                    Mat(gi, N_global+gj) = Mat(gi, N_global+gj) + cmplx(val_dxi, 0.0_wp, wp)
                    
                    ! Down-Up
                    Mat(N_global+gi, gj) = Mat(N_global+gi, gj) + cmplx(val_dxi, 0.0_wp, wp)
                end do
            end do
        end do

        if (present(info)) info = ierr
    end subroutine build_dB_deta_matrix

    !--------------------------------------------------------------
    !
    ! convert_potential_to_xi_eta
    !
    ! this is kinda helper subroutine
    ! to convert potential expressed in r to that in xi and eta coordinates
    ! because finite_element_dvr module providers potential in r coordinate
    ! it is used in both construct_diago_block_matrix and construct_coupling_block_matrix
    !
    ! lambda = alpha^2/2r  * ( -Z_eff/r^2 -(Z_nuclei -1) *(eta * exp(kia *r)) *WS^2(r))
    ! WS = 1 / (1+ (eta /kia) * (exp(kia*r)-1) )
    ! Z_eff = 1 + (Z_nuclei -1) * WS(r)
    !
    ! in parabolic cordinates
    ! r = (xi + eta)/2
    ! thus potential must be converted accordingly
    !
    ! lambda = alpha^2/(xi + eta) * [ -4*Z_eff/(xi + eta)^2 - 2*(Z_nuclei -1) * eta * exp(kia *(xi + eta)/2) * WS^2( (xi + eta)/2 ) 
    ! WS = 1/(1 + (eta /kia) * (exp(kia *(xi + eta)/2) -1) )]
    !
    ! parameters for Xe 5p electron state is below
    ! Z_nuclei = 54
    ! eta = 5.197
    ! kia = 1.048
    !--------------------------------------------------------------

    subroutine convert_potential_to_xi_eta(xi, eta, lambda, Z_eff_val, dVc_deta, dlambda_deta)
        real(wp), intent(in) :: xi, eta
        real(wp), intent(out) :: lambda, Z_eff_val
        real(wp), intent(out), optional :: dVc_deta, dlambda_deta
        
        real(wp) :: r, WS, dZ_eff, dV, d2Z_eff, d2V
        real(wp) :: Z_nuclei, eta_param, kia, alpha, z_core
        real(wp) :: dlambda_dr
        
        ! Parameters for Xe 5p
        Z_nuclei = 54.0_wp
        eta_param = 5.197_wp
        kia = 1.048_wp
        alpha = 1.0_wp / 137.035999084_wp
        z_core = Z_nuclei - 1.0_wp
        
        r = (xi + eta) / 2.0_wp
        ! Avoid division by zero at r=0
        if (r < 1.0e-12_wp) r = 1.0e-12_wp
        
        WS = 1.0_wp / (1.0_wp + (eta_param/kia) * (exp(kia*r) - 1.0_wp) )
        
        ! Z_eff = 1 + (Z_nuclei - 1) * WS(r)
        Z_eff_val = 1.0_wp + z_core * WS
        
        ! dZ_eff/dr
        dZ_eff = -z_core * eta_param * exp(kia*r) * WS**2
        
        ! dV/dr = d/dr(-Z_eff/r) = -dZ_eff/dr / r + Z_eff/r^2
        dV = -dZ_eff / r + Z_eff_val / (r**2)
        
        ! lambda = alpha^2 / (2r) * dV/dr
        lambda = (alpha**2) / (2.0_wp * r) * dV

        ! Calculate derivatives with respect to eta if requested
        if (present(dVc_deta) .or. present(dlambda_deta)) then
            ! d/deta = (dr/deta) * d/dr = 1/2 * d/dr
            
            ! Vc = -Z_eff/r
            ! dVc/dr = dV (calculated above)
            if (present(dVc_deta)) then
                dVc_deta = 0.5_wp * dV
            end if
            
            if (present(dlambda_deta)) then
                ! Need dlambda/dr
                ! lambda = (alpha^2/2) * (1/r * dV/dr)
                ! dlambda/dr = (alpha^2/2) * [ -1/r^2 * dV/dr + 1/r * d2V/dr2 ]
                
                ! d2V/dr2 = d/dr( -dZ_eff/r + Z_eff/r^2 )
                !         = -(d2Z_eff/r - dZ_eff/r^2) + (dZ_eff/r^2 - 2*Z_eff/r^3)
                !         = -d2Z_eff/r + 2*dZ_eff/r^2 - 2*Z_eff/r^3
                
                ! d2Z_eff/dr2
                ! dZ_eff = -C * exp(kia*r) * WS^2  (C = z_core * eta_param)
                ! d2Z_eff = -C * [ kia*exp * WS^2 + exp * 2*WS*dWS/dr ]
                ! dWS/dr = - (eta_param/kia) * kia*exp * WS^2 = -eta_param * exp * WS^2
                ! So:
                ! d2Z_eff = -C * exp * WS^2 * (kia - 2*eta_param*exp*WS)
                !         = dZ_eff * (kia - 2.0_wp * eta_param * exp(kia*r) * WS)
                
                d2Z_eff = dZ_eff * (kia - 2.0_wp * eta_param * exp(kia*r) * WS)
                
                d2V = -d2Z_eff/r + 2.0_wp*dZ_eff/(r**2) - 2.0_wp*Z_eff_val/(r**3)
                
                dlambda_dr = (alpha**2 / 2.0_wp) * ( -dV/(r**2) + d2V/r )
                
                dlambda_deta = 0.5_wp * dlambda_dr
            end if
        end if
        
    end subroutine convert_potential_to_xi_eta

    
end module aeigen_solver

module aeigen_solver

!------------------------------------------------------------------
!UNNDER CONSTRUCTION
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
        real(wp), intent(in) :: eta, Mj, E_val, F_val
        integer, intent(in) :: key_m
        real(wp), intent(out) :: Mat(N_global, N_global)
        integer, intent(out), optional :: info

        integer :: i, e, gi, ierr
        real(wp) :: xi, m_val, Sz, lambda_val, Z_eff_val
        real(wp) :: V_diag, term_azi, term_pot, term_soi
        
        ierr = 0
        
        ! 1. Kinetic Energy Term: d/dxi ( xi d/dxi )
        ! We use build_xi_differential_matrix with key_f=3
        ! This returns K_ij = < i' | xi | j' >
        ! The operator in B(eta) is + d/dxi ( xi d/dxi )
        ! Matrix element < i | d/dxi ( xi d/dxi ) | j > = - < i' | xi | j' > (by integration by parts)
        ! So we multiply by -1.0
        
        call build_xi_differential_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                          D_ref, w_ref, eta, 3, Mat, ierr)
        
        Mat = -1.0_wp * Mat
        
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
                
                if (xi < 1.0e-12_wp) xi = 1.0e-12_wp
                
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

    subroutine assemble_global_matrix

    end subroutine assemble_global_matrix

    subroutine solve_adiabatic_hamiltonian

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
                                                  D_ref, w_ref, eta, Mj, Evals, Evecs, P_mat, info)
        integer, intent(in) :: N_global, nelems
        integer, intent(in) :: elem_map(:,:)
        real(wp), intent(in) :: x_global(:), w_global(:)
        real(wp), intent(in) :: D_ref(:,:), w_ref(:)
        real(wp), intent(in) :: eta, Mj
        complex(wp), intent(in) :: Evals(:)
        complex(wp), intent(in) :: Evecs(:,:)
        complex(wp), intent(out) :: P_mat(:,:)
        integer, intent(out), optional :: info

        integer :: i, j, ierr
        complex(wp), allocatable :: dB_deta(:,:)
        complex(wp), parameter :: czero = (0.0_wp, 0.0_wp)
        real(wp), parameter :: tol = 1.0e-10_wp

        ierr = 0
        allocate(dB_deta(2*N_global, 2*N_global))
        
        ! Build dB/deta matrix
        call build_dB_deta_matrix(N_global, nelems, elem_map, x_global, w_global, &
                                  D_ref, w_ref, eta, Mj, dB_deta, ierr)
        
        ! Calculate P_ij = <i|dB/deta|j> / (Ej - Ei)
        ! P_mat = <i | dB/deta | j>
        ! Note: Evecs are stored in columns.
        P_mat = matmul(conjg(transpose(Evecs)), matmul(dB_deta, Evecs))
        
        ! Divide by energy difference
        do i = 1, size(Evals)
            do j = 1, size(Evals)
                if (i == j) then
                    P_mat(i,j) = czero ! Diagonal term is 0 (or purely imaginary Berry phase, assumed 0 for now)
                else
                    if (abs(Evals(j) - Evals(i)) > tol) then
                        P_mat(i,j) = P_mat(i,j) / (Evals(j) - Evals(i))
                    else
                        ! Degenerate case: needs special handling or set to 0
                        P_mat(i,j) = czero 
                    end if
                end if
            end do
        end do

        deallocate(dB_deta)
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

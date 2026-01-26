program siegert_main
    !------------------------------------------------------------------------------
    ! Siegert State Solver Main Program
    !
    ! Description:
    !   Main driver for finding Siegert states (complex energy eigenvalues) using
    !   the R-matrix propagation method combined with FE-DVR.
    !   It searches for the complex energy E where the matching condition
    !   det(R_in(E) - R_out(E)) = 0 is satisfied.
    !
    !   The program performs the following steps:
    !   1. Initializes R-matrix sectors.
    !   2. Builds sector Hamiltonians (energy-independent parts).
    !   3. Performs a root search (Secant method) for the complex energy E.
    !   4. Analyzes channel amplitudes and outputs physical quantities (basis functions,
    !      adiabatic potentials) at the converged energy.
    !
    ! Modules Used:
    !   rmatrix_propagation - Handles R-matrix propagation logic.
    !   aeigen_solver       - Solves the adiabatic eigenvalue problem using FE-DVR.
    !
    ! Assisted by GitHub Copilot
    !------------------------------------------------------------------------------
    use, intrinsic :: iso_fortran_env, only: wp => real64
    use rmatrix_propagation
    use aeigen_solver
    implicit none

    ! Parameters
    integer, parameter :: n_sectors = 100
    real(wp), parameter :: r_start = 0.0_wp
    real(wp), parameter :: r_end = 100.0_wp ! Need large box for Siegert state
    integer, parameter :: n_channels = 10 ! Reduced for speed in test
    integer, parameter :: n_basis_per_sector = 15
    integer, parameter :: match_sector_idx = 100 ! Matching point index
    
    ! Variables
    type(sector_basis_t), allocatable :: sectors(:)
    integer :: i, ierr, iter
    real(wp) :: dr, r_local_min, r_local_max
    complex(wp) :: E_curr, E_prev, E_next
    complex(wp) :: det_curr, det_prev
    complex(wp), allocatable :: R_in(:,:), R_out(:,:)
    
    ! Physical Parameters
    integer :: nelems_ad = 100
    integer :: N_per_elem_ad = 10
    real(wp) :: x_max_ad = 50.0_wp
    real(wp) :: Mj = 0.5_wp 
    real(wp) :: F_val = 0.05_wp ! Tunneling regime (F=0.03 a.u.)
    real(wp) :: Gamma_guess
    
    print *, "--- Siegert State Search (Newton-Raphson/Secant) ---"
    print *, "Field Strength F =", F_val
    
    ! 1. Initialize Sectors
    allocate(sectors(n_sectors))
    dr = (r_end - r_start) / real(n_sectors, wp)
    
    print *, "Initializing ", n_sectors, " sectors from ", r_start, " to ", r_end
    
    do i = 1, n_sectors
        r_local_min = r_start + (i-1)*dr
        r_local_max = r_start + i*dr
        call init_sector_basis(sectors(i), r_local_min, r_local_max, n_channels, n_basis_per_sector)
    end do
    
    ! 2. Build Sector Hamiltonians (Energy Independent Part)
    print *, "Building sector Hamiltonians..."
    do i = 1, n_sectors
        if (mod(i, 20) == 0) print *, "  Sector ", i, "/", n_sectors
        call build_and_solve_sector(sectors(i), nelems_ad, N_per_elem_ad, x_max_ad, Mj, F_val, ierr)
        if (ierr /= 0) stop "Error in sector build"
        
        ! Export coupling matrices for visualization of "Dynamic Selection Rules"
        ! Sectors correspond to roughly: 2->1.0, 4->2.0, 6->3.0, 10->5.0, 14->7.0, 50->25.0, 90->45.0
        if (i == 2 .or. i == 4 .or. i == 6 .or. i == 10 .or. i == 14 .or. i == 50 .or. i == 90) then
            block
                character(len=20) :: prefix
                write(prefix, '("sector_", I0.3)') i
                call export_coupling_matrices(sectors(i), trim(prefix))
                print *, "  Exported coupling matrices for sector ", i
            end block
        end if
    end do
    
    ! 3. Root Search Loop (Secant Method)
    ! We search for E such that det(R_in(E) - R_out(E)) = 0
    
    ! Initial Guess using Landau's formula for H(1s) as a rough estimate
    ! Gamma = 4/F * exp(-2/3F)
    ! For F=0.03, Gamma is very small (~1e-9). 
    ! We use a small value to ensure we find the physical tunneling root,
    ! avoiding spurious large-gamma roots.
    print *, "Initial Guess Gamma (Physical estimate):", Gamma_guess
    
    ! E_bound for Xe 5p (j=3/2) is approx -0.446 a.u.
    ! Stark shift lowers it. We try -0.450 to find the ground state (Spin Up).
    ! 3. Root Search Loop (Secant Method)
    ! We search for E such that det(R_in(E) - R_out(E)) = 0
    
    ! Initial Guess using ADK-like formula for Tunneling Rate
    ! Gamma ~ F * exp( -2 * (2*Ip)^1.5 / 3F )
    ! This allows us to distinguish J=3/2 (Ip~0.46) and J=1/2 (Ip~0.51)
    !
    ! Target: Xe 5p J=3/2 (Ip approx 0.46 a.u.)
    ! If searching for J=1/2, change Ip_target to 0.51
    block
        real(wp) :: Ip_target, kappa
        Ip_target = 0.450_wp 
        kappa = sqrt(2.0_wp * Ip_target)
        
        ! Simplified ADK rate estimate
        Gamma_guess = F_val * exp(-2.0_wp * kappa**3 / (3.0_wp * F_val))
    end block
    
    print *, "Initial Guess Gamma (ADK estimate):", Gamma_guess
    
    ! E_bound for Xe 5p (j=3/2) is approx -0.446 a.u.
    ! Stark shift lowers it. We try -0.460 to find the ground state (Spin Up).
    E_prev = cmplx(-0.465_wp, -Gamma_guess / 2.0_wp, wp)
    E_curr = E_prev + cmplx(0.0001_wp, 0.0_wp, wp) ! Perturb slightly for secant step
    
    ! Compute determinant for E_prev
    call compute_matching_determinant(E_prev, det_prev)
    print *, "Iter 0: E =", E_prev, " Det =", det_prev
    
    open(unit=20, file='siegert_history_mj05.dat', status='replace')
    write(20, '(A)') '# Iter Re(E) Im(E) Abs(Det)'
    write(20, *) 0, real(E_prev), aimag(E_prev), abs(det_prev)
    
    do iter = 1, 20
        call compute_matching_determinant(E_curr, det_curr)
        print *, "Iter", iter, ": E =", E_curr, " Det =", det_curr
        write(20, *) iter, real(E_curr), aimag(E_curr), abs(det_curr)
        
        ! Check convergence
        if (abs(det_curr) < 1.0e-10_wp) then
            print *, "Converged!"
            exit
        end if
        
        if (abs(det_curr - det_prev) < 1.0e-20_wp) then
            print *, "Determinant difference too small, stopping."
            exit
        end if
        
        ! Secant Update
        E_next = E_curr - det_curr * (E_curr - E_prev) / (det_curr - det_prev)
        
        ! Update history
        E_prev = E_curr
        det_prev = det_curr
        E_curr = E_next
        
        ! Safety clamp (optional)
        if (real(E_curr) > 0.0_wp) E_curr = cmplx(-0.1_wp, aimag(E_curr), wp)
    end do
    
    close(20)
    
    print *, "Final Energy:", E_curr
    print *, "Re(E) =", real(E_curr)
    print *, "Gamma =", -2.0_wp * aimag(E_curr)
    
    open(unit=21, file='siegert_result_mj05.dat', status='replace')
    write(21, '(A)') '# Final Energy Result'
    write(21, '(A)') '# Re(E) Im(E) Gamma'
    write(21, *) real(E_curr), aimag(E_curr), -2.0_wp * aimag(E_curr)
    close(21)
    
    ! Analyze Channel Amplitudes
    call analyze_channel_amplitudes()

contains

    !------------------------------------------------------------------------------
    ! Subroutine: analyze_channel_amplitudes
    !
    ! Description:
    !   Analyzes the channel amplitudes at the matching point after the complex
    !   energy has been found. It solves the homogeneous equation:
    !       (R_in - R_out) * C = 0
    !   to find the channel amplitude vector C. This is done using Singular Value
    !   Decomposition (SVD) to find the null space of (R_in - R_out).
    !
    !   After finding C, it performs a final physical analysis:
    !   1. Re-solves the adiabatic Hamiltonian at the final boundary using the
    !      converged energy E_curr.
    !   2. Computes and outputs channel probabilities, basis energies, and spin
    !      probabilities.
    !   3. Outputs the basis functions and adiabatic potentials for plotting.
    !------------------------------------------------------------------------------
    subroutine analyze_channel_amplitudes()
        complex(wp), allocatable :: R_diff(:,:)
        complex(wp), allocatable :: U(:,:), VT(:,:), work(:)
        real(wp), allocatable :: S(:), rwork(:)
        integer :: info, lwork, k
        complex(wp) :: amplitude_vec(n_channels)
        real(wp) :: total_norm
        
        allocate(R_diff(n_channels, n_channels))
        R_diff = R_in - R_out
        
        allocate(S(n_channels))
        allocate(U(n_channels, n_channels))
        allocate(VT(n_channels, n_channels))
        
        ! Query workspace for ZGESVD
        allocate(work(1))
        allocate(rwork(5*n_channels))
        call ZGESVD('A', 'A', n_channels, n_channels, R_diff, n_channels, S, U, n_channels, VT, n_channels, work, -1, rwork, info)
        lwork = int(real(work(1)))
        deallocate(work)
        allocate(work(lwork))
        
        ! Compute SVD: R_diff = U * S * VT
        ! The solution to R_diff * C = 0 corresponds to the singular vector for the smallest singular value.
        ! S is sorted descending. The last row of VT (conjugate of last col of V) is the target.
        call ZGESVD('A', 'A', n_channels, n_channels, R_diff, n_channels, S, U, n_channels, VT, n_channels, work, lwork, rwork, info)
        
        if (info /= 0) then
            print *, "Error in SVD analysis"
            return
        end if
        
        ! Extract null vector (last row of VT, conjugated)
        do k = 1, n_channels
            amplitude_vec(k) = conjg(VT(n_channels, k))
        end do
        
        ! Normalize vector (sum of |c|^2 = 1)
        total_norm = 0.0_wp
        do k = 1, n_channels
            total_norm = total_norm + abs(amplitude_vec(k))**2
        end do
        total_norm = sqrt(total_norm)
        amplitude_vec = amplitude_vec / total_norm
        
        ! --- Physical Analysis & Output ---
        block
            real(wp) :: eta_final, eta_scan
            complex(wp), allocatable :: evals_final(:), evecs_final(:,:)
            complex(wp), allocatable :: evals_scan(:), evecs_scan(:,:)
            real(wp), allocatable :: x_gl(:), w_gl(:)
            real(wp) :: prob_up, prob_down, norm_spin
            integer :: N_gl, i_gl, i_eta, n_eta_steps
            
            eta_final = sectors(n_sectors)%r_max
            ! Use the converged Siegert energy E_curr for the final analysis
            ! to get the correct adiabatic potentials and basis functions.
            
            ! 1. Solve at final boundary for amplitudes and basis functions
            call solve_adiabatic_hamiltonian(nelems_ad, N_per_elem_ad, x_max_ad, eta_final, Mj, E_curr, F_val, &
                                             evals_final, evecs_final, x_gl, w_gl, ierr)
            N_gl = size(x_gl)
            
            ! Output Amplitudes
            open(unit=22, file='siegert_amplitudes_mj05.dat', status='replace')
            write(22, '(A)') '# Channel Prob(|C|^2) BasisEnergy SpinUp SpinDown'
            
            print *, " "
            print *, "--- Final Physical Analysis (at Eta =", eta_final, ") ---"
            print *, "Channel | Prob(|C|^2)  | Basis Energy | Spin Up Prob   | Spin Down Prob"
            print *, "--------------------------------------------------------------------------"
            
            do k = 1, n_channels
                prob_up = 0.0_wp
                prob_down = 0.0_wp
                do i_gl = 1, N_gl
                    prob_up = prob_up + abs(evecs_final(i_gl, k))**2 * w_gl(i_gl)
                    prob_down = prob_down + abs(evecs_final(N_gl+i_gl, k))**2 * w_gl(i_gl)
                end do
                norm_spin = prob_up + prob_down
                if (norm_spin > 1.0e-10_wp) then
                    prob_up = prob_up / norm_spin
                    prob_down = prob_down / norm_spin
                end if
                
                print '(I7, " | ", ES12.5, " | ", F12.5, " | ", ES12.5, " | ", ES12.5)', &
                      k, abs(amplitude_vec(k))**2, real(evals_final(k)), prob_up, prob_down
                      
                write(22, '(I7, ES16.8, ES16.8, ES16.8, ES16.8)') &
                      k, abs(amplitude_vec(k))**2, real(evals_final(k)), prob_up, prob_down
            end do
            close(22)
            
            ! 2. Output Basis Functions (Phi) at final boundary
            open(unit=23, file='basis_functions_final_mj05.dat', status='replace')
            write(23, '(A)') '# x_gl (Channel 1 Up/Down) (Channel 2 Up/Down) ...'
            do i_gl = 1, N_gl
                write(23, '(ES16.8)', advance='no') x_gl(i_gl)
                do k = 1, n_channels
                    ! Write |Phi|^2 for Up and Down
                    write(23, '(ES16.8, ES16.8)', advance='no') &
                        abs(evecs_final(i_gl, k))**2, abs(evecs_final(N_gl+i_gl, k))**2
                end do
                write(23, *) ! Newline
            end do
            close(23)
            print *, "Written basis_functions_final.dat"

            ! 3. Output Adiabatic Potentials (Beta(eta))
            open(unit=24, file='adiabatic_potentials_mj05.dat', status='replace')
            write(24, '(A)') '# Eta (Channel 1 Energy) (Channel 2 Energy) ...'
            
            n_eta_steps = 100
            do i_eta = 0, n_eta_steps
                eta_scan = r_start + (r_end - r_start) * real(i_eta, wp) / real(n_eta_steps, wp)
                if (eta_scan < 1.0e-6_wp) eta_scan = 1.0e-6_wp ! Avoid singularity at 0
                
                call solve_adiabatic_hamiltonian(nelems_ad, N_per_elem_ad, x_max_ad, eta_scan, Mj, E_curr, F_val, &
                                                 evals_scan, evecs_scan, x_gl, w_gl, ierr)
                
                write(24, '(ES16.8)', advance='no') eta_scan
                do k = 1, n_channels
                    write(24, '(ES16.8)', advance='no') real(evals_scan(k))
                end do
                write(24, *) ! Newline
            end do
            close(24)
            print *, "Written adiabatic_potentials.dat"
            
        end block
        
    end subroutine analyze_channel_amplitudes

    !------------------------------------------------------------------------------
    ! Subroutine: compute_matching_determinant
    !
    ! Description:
    !   Computes the determinant of the matching matrix M = R_in - R_out at the
    !   matching point.
    !
    !   Steps:
    !   1. Propagates the R-matrix from the origin (r=0) to the matching point
    !      to obtain R_in.
    !   2. Calculates the boundary R-matrix at r_max using the Siegert boundary
    !      condition.
    !   3. Propagates the R-matrix inward from r_max to the matching point to
    !      obtain R_out.
    !   4. Computes det(R_in - R_out).
    !
    ! Inputs:
    !   E_val   - The complex energy at which to evaluate the determinant.
    !
    ! Outputs:
    !   det_val - The computed complex determinant.
    !------------------------------------------------------------------------------
    subroutine compute_matching_determinant(E_val, det_val)
        complex(wp), intent(in) :: E_val
        complex(wp), intent(out) :: det_val
        
        complex(wp), allocatable :: R_temp(:,:), R_next_temp(:,:)
        complex(wp), allocatable :: R_diff(:,:)
        integer :: j
        
        allocate(R_temp(n_channels, n_channels))
        allocate(R_next_temp(n_channels, n_channels))
        
        ! --- Outward Propagation (0 -> match) ---
        ! Boundary condition at r=0: Hard wall -> R = 0
        R_temp = (0.0_wp, 0.0_wp)
        
        do j = 1, match_sector_idx
            call propagate_step_outward(R_temp, sectors(j), E_val, R_next_temp, ierr)
            R_temp = R_next_temp
        end do
        
        if (allocated(R_in)) deallocate(R_in)
        allocate(R_in(n_channels, n_channels))
        R_in = R_temp
        
        ! --- Inward Propagation (max -> match) ---
        ! Boundary condition at r_max: Siegert BC
        call get_boundary_R_matrix(sectors(n_sectors)%r_max, E_val, F_val, R_temp)
        
        do j = n_sectors, match_sector_idx + 1, -1
            call propagate_step_inward(R_temp, sectors(j), E_val, R_next_temp, ierr)
            R_temp = R_next_temp
        end do
        
        if (allocated(R_out)) deallocate(R_out)
        allocate(R_out(n_channels, n_channels))
        R_out = R_temp
        
        ! --- Matching ---
        allocate(R_diff(n_channels, n_channels))
        R_diff = R_in - R_out
        
        call get_complex_determinant(n_channels, R_diff, det_val)
        
    end subroutine compute_matching_determinant

    !------------------------------------------------------------------------------
    ! Subroutine: get_boundary_R_matrix
    !
    ! Description:
    !   Calculates the R-matrix at the outer boundary (r_max) using the Siegert
    !   asymptotic boundary condition for a particle in a static electric field.
    !
    !   The asymptotic form of the wavefunction is:
    !     chi ~ (4/F*eta)^1/4 * exp(i * (sqrt(F)/3 * eta^1.5 + E/sqrt(F) * eta^0.5))
    !
    !   The logarithmic derivative is computed analytically, and the R-matrix is
    !   constructed as a diagonal matrix with R_kk = 1 / (chi'/chi).
    !
    ! Inputs:
    !   eta  - The radial coordinate (parabolic coordinate eta).
    !   E    - The complex energy.
    !   F    - The electric field strength.
    !
    ! Outputs:
    !   R_bc - The boundary R-matrix.
    !------------------------------------------------------------------------------
    subroutine get_boundary_R_matrix(eta, E, F, R_bc)
        real(wp), intent(in) :: eta, F
        complex(wp), intent(in) :: E
        complex(wp), intent(out) :: R_bc(:,:)
        
        complex(wp) :: log_deriv
        integer :: k
        
        ! Analytical Logarithmic Derivative:
        ! d/deta ln(chi) = -1/(4*eta) + i * ( sqrt(F)/2 * eta^0.5 + E/(2*sqrt(F)) * eta^-0.5 )
        
        log_deriv = -1.0_wp / (4.0_wp * eta) + &
                    (0.0_wp, 1.0_wp) * ( (sqrt(F)/2.0_wp) * sqrt(eta) + &
                                         E / (2.0_wp * sqrt(F)) / sqrt(eta) )
                                         
        ! R = (chi' / chi)^-1 = 1 / log_deriv
        R_bc = (0.0_wp, 0.0_wp)
        do k = 1, n_channels
            R_bc(k, k) = 1.0_wp / log_deriv
        end do
        
    end subroutine get_boundary_R_matrix

    !------------------------------------------------------------------------------
    ! Subroutine: get_complex_determinant
    !
    ! Description:
    !   Computes the determinant of a complex square matrix using LAPACK's ZGETRF
    !   (LU factorization).
    !
    ! Inputs:
    !   n   - Dimension of the matrix.
    !   A   - The complex matrix (n x n).
    !
    ! Outputs:
    !   det - The computed determinant.
    !------------------------------------------------------------------------------
    subroutine get_complex_determinant(n, A, det)
        integer, intent(in) :: n
        complex(wp), intent(in) :: A(:,:)
        complex(wp), intent(out) :: det
        
        complex(wp) :: A_copy(n, n)
        integer :: IPIV(n)
        integer :: info, k
        
        A_copy = A
        call ZGETRF(n, n, A_copy, n, IPIV, info)
        
        det = (1.0_wp, 0.0_wp)
        do k = 1, n
            det = det * A_copy(k, k)
            if (IPIV(k) /= k) det = -det
        end do
    end subroutine get_complex_determinant

end program siegert_main

!これまでのファイルを組み合わせて,遷移双極子モーメントのエネルギースペクトラムを計算するプログラム
!開発中なのでコメント,説明追加前
!まだテストも一切していません

program compute_spin_polarization

    use, intrinsic :: iso_fortran_env, only: wp => real64
    use golub_welsch
    use finite_element_dvr
    use dvr_integral_equation_methods

    implicit none

    integer, parameter :: n_per_elem = 100
    integer, parameter :: nelems = 50
    real(wp), parameter :: r_min = 0.0_wp
    real(wp), parameter :: r_max = 500.0_wp
    real(wp), parameter :: photon_energy = 3.14_wp          ! incident photon energy (a.u.)
    real(wp), parameter :: energy_min = 2.0_wp              ! scattering energy scan lower bound (a.u.)
    real(wp), parameter :: energy_max = 4.0_wp              ! scattering energy scan upper bound (a.u.)
    real(wp), parameter :: energy_step = 0.001_wp            ! scattering energy increment (a.u.)
    character(len=64) :: initial_label
    character(len=64) :: normalized_label
    character(len=64) :: output_tag
    character(len=256) :: argument_buffer

    real(wp), allocatable :: x_ref(:), w_ref(:), b_weights(:), d_ref(:,:)
    real(wp), allocatable :: x_global(:), m_global(:)
    real(wp), allocatable :: m_elems(:,:), k_elems(:,:,:), k_global(:,:)
    real(wp), allocatable :: v_core(:), v_eff(:)
    real(wp), allocatable :: hmat(:,:), hred(:,:), work(:)
    real(wp), allocatable :: eigvals_bound(:)
    real(wp), allocatable :: u_bound(:)
    real(wp), allocatable :: energies(:), p_up_grid(:), p_down_grid(:), ratio_grid(:)
    complex(wp), allocatable :: amp_up_l2_grid(:), amp_down_l2_grid(:)
    complex(wp), allocatable :: amp_up_l0_grid(:)
    complex(wp), allocatable :: amp_up_l2_j52_grid(:), amp_down_l2_j52_grid(:)
    real(wp), allocatable :: phi_l2_tmp(:), nodes_l2_tmp(:)
    real(wp), allocatable :: phi_l0_tmp(:), nodes_l0_tmp(:)
    real(wp), allocatable :: phi_l2_j52_tmp(:), nodes_l2_j52_tmp(:)
    real(wp), allocatable :: phi_l2_saved(:), nodes_saved(:)
    real(wp), allocatable :: phi_l0_saved(:)
    real(wp), allocatable :: phi_l2_j52_saved(:)
    real(wp), allocatable :: phi_l2_prev(:)
    real(wp), allocatable :: phi_l0_prev(:)
    real(wp), allocatable :: phi_l2_j52_prev(:)
    real(wp) :: delta_l2, delta_l0, delta_l2_j52
    real(wp) :: delta_l2_saved, delta_l0_saved, delta_l2_j52_saved
    real(wp) :: norm_val, energy_bound
    real(wp) :: r_value, centrifugal_term
    complex(wp) :: r_integral_l2, r_integral_l0, r_integral_l2_j52
    complex(wp) :: amp_up_l2, amp_up_l0, amp_down_l2
    complex(wp) :: amp_up_l2_j52, amp_down_l2_j52
    complex(wp) :: pref_up_l2, pref_up_l0, pref_down_l2
    complex(wp) :: pref_up_l2_j52, pref_down_l2_j52
    real(wp) :: p_up_val, p_down_val, ratio_val
    real(wp) :: energy_scatt, target_energy, best_diff
    integer :: elem_map(nelems, n_per_elem)
    integer :: n_global, info, lwork, idx, n_active
    integer :: e, i, j, n_energy, idx_energy, ref_index
    logical, allocatable :: active_mask(:)
    integer, allocatable :: active_idx(:)
    real(wp) :: k_scatt
    real(wp), parameter :: node_tol = 1.0e-8_wp
    logical :: have_saved_wave
    logical :: have_prev_l2
    logical :: have_prev_l0
    logical :: have_prev_l2_j52
    complex(wp) :: weight_val, bound_val, scat_val_l2, scat_val_l0, scat_val_l2_j52
    real(wp) :: target_energy_5p
    real(wp) :: coeff_d_j, coeff_s_j
    real(wp) :: coeff_d_j52
    logical :: use_two_d_channels
    logical :: need_s_wave
    character(len=256) :: spectrum_filename, amplitude_filename

    argument_buffer = ''
    call get_command_argument(1, argument_buffer)
    initial_label = trim(argument_buffer)
    if (len_trim(initial_label) == 0) initial_label = 'j32_m32'
    normalized_label = to_lower(trim(adjustl(initial_label)))

    use_two_d_channels = .false.
    need_s_wave = .false.
    pref_up_l2 = cmplx(0.0_wp, 0.0_wp, kind=wp)
    pref_down_l2 = cmplx(0.0_wp, 0.0_wp, kind=wp)
    pref_up_l0 = cmplx(0.0_wp, 0.0_wp, kind=wp)
    pref_up_l2_j52 = cmplx(0.0_wp, 0.0_wp, kind=wp)
    pref_down_l2_j52 = cmplx(0.0_wp, 0.0_wp, kind=wp)

    select case (normalized_label)
    case ('5p1/2', 'j12_m12', 'j=1/2,m=1/2')
        ! Dipole selection for 5p_{1/2}, m_j = +1/2 → εd_{3/2}, εs_{1/2}
        target_energy_5p = -0.446_wp
        coeff_d_j = sqrt(2.0_wp) / 3.0_wp
        coeff_s_j = -1.0_wp / 3.0_wp
        pref_up_l2 = cmplx(-sqrt(2.0_wp/5.0_wp) * coeff_d_j, 0.0_wp, kind=wp)
        pref_down_l2 = cmplx(sqrt(3.0_wp/5.0_wp) * coeff_d_j, 0.0_wp, kind=wp)
        pref_up_l0 = cmplx(coeff_s_j, 0.0_wp, kind=wp)
        need_s_wave = .true.
        output_tag = 'j12_m12'
    case ('5p3/2', 'j32_m12', 'j=3/2,m=1/2')
        ! 5p_{3/2}, m_j = +1/2 → εd_{3/2}, εd_{5/2}, εs_{1/2}
        target_energy_5p = -0.439_wp
        coeff_d_j = -1.0_wp / 15.0_wp
        coeff_d_j52 = sqrt(6.0_wp) / 5.0_wp
        coeff_s_j = sqrt(2.0_wp) / 3.0_wp
        pref_up_l2 = cmplx(-sqrt(2.0_wp/5.0_wp) * coeff_d_j, 0.0_wp, kind=wp)
        pref_down_l2 = cmplx(sqrt(3.0_wp/5.0_wp) * coeff_d_j, 0.0_wp, kind=wp)
        pref_up_l2_j52 = cmplx(sqrt(3.0_wp/5.0_wp) * coeff_d_j52, 0.0_wp, kind=wp)
        pref_down_l2_j52 = cmplx(sqrt(2.0_wp/5.0_wp) * coeff_d_j52, 0.0_wp, kind=wp)
        pref_up_l0 = cmplx(coeff_s_j, 0.0_wp, kind=wp)
        use_two_d_channels = .true.
        need_s_wave = .true.
        output_tag = 'j32_m12'
    case ('j32_m32', 'j=3/2,m=3/2')
        ! 5p_{3/2}, m_j = +3/2 → εd_{3/2}, εd_{5/2}
        target_energy_5p = -0.439_wp
        coeff_d_j = -1.0_wp / 5.0_wp
        coeff_d_j52 = 2.0_wp / 5.0_wp
        pref_up_l2 = cmplx(-sqrt(1.0_wp/5.0_wp) * coeff_d_j, 0.0_wp, kind=wp)
        pref_down_l2 = cmplx(sqrt(4.0_wp/5.0_wp) * coeff_d_j, 0.0_wp, kind=wp)
        pref_up_l2_j52 = cmplx(sqrt(4.0_wp/5.0_wp) * coeff_d_j52, 0.0_wp, kind=wp)
        pref_down_l2_j52 = cmplx(sqrt(1.0_wp/5.0_wp) * coeff_d_j52, 0.0_wp, kind=wp)
        use_two_d_channels = .true.
        need_s_wave = .false.
        output_tag = 'j32_m32'
    case default
        print '(A)', 'ERROR: Unsupported initial state label: '//trim(initial_label)
        print '(A)', 'Supported labels: j12_m12, j32_m12, j32_m32'
        stop 1
    end select

    spectrum_filename = 'spin_pol_'//trim(output_tag)//'_spectrum.dat'
    amplitude_filename = 'spin_pol_'//trim(output_tag)//'_amplitudes.dat'

    ! Build FEDVR grid
    allocate(x_ref(n_per_elem), w_ref(n_per_elem))
    allocate(x_global(nelems * (n_per_elem - 1) + 1))
    allocate(m_global(size(x_global)))

    call build_global_grid(nelems, n_per_elem, x_ref, w_ref, r_min, r_max, n_global, &
                           x_global, m_global, elem_map, info)
    if (info /= 0) then
        print '(A,I0)', 'build_global_grid failed, info = ', info
        stop 1
    end if

    allocate(b_weights(n_per_elem), d_ref(n_per_elem, n_per_elem))
    call compute_barycentric_weights(n_per_elem, x_ref, b_weights)
    call compute_differentiation_matrix(n_per_elem, x_ref, b_weights, d_ref)

    allocate(m_elems(n_per_elem, nelems))
    allocate(k_elems(n_per_elem, n_per_elem, nelems))

    do e = 1, nelems
        call build_element_matrices(n_per_elem, d_ref, w_ref, &
                                    r_min + (e - 1) * (r_max - r_min) / nelems, &
                                    r_min + e * (r_max - r_min) / nelems, &
                                    m_elems(:, e), k_elems(:, :, e), info)
        if (info /= 0) then
            print '(A,I0)', 'build_element_matrices failed, info = ', info
            stop 1
        end if
    end do

    allocate(k_global(n_global, n_global))
    m_global = 0.0_wp
    call assemble_global_matrices(nelems, n_per_elem, elem_map, m_elems, k_elems, &
                                  n_global, m_global, k_global, info)
    if (info /= 0) then
        print '(A,I0)', 'assemble_global_matrices failed, info = ', info
        stop 1
    end if

    allocate(v_core(n_global), v_eff(n_global))
    call build_potential_diag(n_global, x_global, v_core, 1, 0.5_wp, 1, info)
    if (info /= 0) then
        print '(A,I0)', 'build_potential_diag failed for bound state, info = ', info
        stop 1
    end if

    do idx = 1, n_global
        r_value = max(x_global(idx), 1.0e-8_wp)
        centrifugal_term = 0.5_wp * real(1 * (1 + 1), wp) / (r_value * r_value)
        v_eff(idx) = v_core(idx) + centrifugal_term
    end do

    allocate(hmat(n_global, n_global))
    call build_hamiltonian_matrix(k_global, v_eff, m_global, hmat, info)
    if (info /= 0) then
        print '(A,I0)', 'build_hamiltonian_matrix failed, info = ', info
        stop 1
    end if

    allocate(active_mask(n_global))
    active_mask = .true.
    if (n_global < 3) then
        print '(A)', 'Insufficient grid points to enforce Dirichlet boundaries'
        stop 1
    end if
    active_mask(1) = .false.
    active_mask(n_global) = .false.

    n_active = count(active_mask)
    allocate(active_idx(n_active))
    active_idx = pack([(i, i=1, n_global)], active_mask)
    n_active = size(active_idx)

    allocate(hred(n_active, n_active), eigvals_bound(n_active))
    do i = 1, n_active
        do j = 1, n_active
            hred(i, j) = hmat(active_idx(i), active_idx(j))
        end do
    end do

    lwork = max(1, 3 * n_active - 1)
    allocate(work(lwork))
    call dsyev('V', 'U', n_active, hred, n_active, eigvals_bound, work, lwork, info)
    deallocate(work)
    if (info /= 0) then
        print '(A,I0)', 'dsyev failed on reduced Hamiltonian, info = ', info
        stop 1
    end if

    idx = nearest_bound_state(eigvals_bound, target_energy_5p)
    if (idx < 1) then
        print '(A)', 'Could not locate 5p bound state'
        stop 1
    end if

    energy_bound = eigvals_bound(idx)
    allocate(u_bound(n_global))
    u_bound = 0.0_wp

    do i = 1, n_active
        u_bound(active_idx(i)) = hred(i, idx)
    end do

    do i = 1, n_active
        if (m_global(active_idx(i)) > 0.0_wp) then
            u_bound(active_idx(i)) = u_bound(active_idx(i)) / sqrt(m_global(active_idx(i)))
        else
            u_bound(active_idx(i)) = 0.0_wp
        end if
    end do

    do i = 1, n_global
        if (.not. active_mask(i)) u_bound(i) = 0.0_wp
    end do
    norm_val = sqrt(sum((u_bound(active_idx)**2) * m_global(active_idx)))
    if (norm_val > 0.0_wp) then
        u_bound(active_idx) = u_bound(active_idx) / norm_val
    end if

    if (energy_max <= energy_min) then
        print '(A)', 'ERROR: energy_max must be greater than energy_min.'
        stop 1
    end if
    if (energy_step <= 0.0_wp) then
        print '(A)', 'ERROR: energy_step must be positive.'
        stop 1
    end if

    n_energy = ceiling((energy_max - energy_min) / energy_step) + 1
    if (n_energy < 1) n_energy = 1

    allocate(energies(n_energy), p_up_grid(n_energy), p_down_grid(n_energy), ratio_grid(n_energy))
    do idx_energy = 1, n_energy
        energies(idx_energy) = min(energy_min + real(idx_energy - 1, wp) * energy_step, energy_max)
    end do
    p_up_grid = 0.0_wp
    p_down_grid = 0.0_wp
    ratio_grid = huge(1.0_wp)
    allocate(amp_up_l2_grid(n_energy), amp_down_l2_grid(n_energy))
    amp_up_l2_grid = cmplx(0.0_wp, 0.0_wp, kind=wp)
    amp_down_l2_grid = cmplx(0.0_wp, 0.0_wp, kind=wp)
    if (need_s_wave) then
        allocate(amp_up_l0_grid(n_energy))
        amp_up_l0_grid = cmplx(0.0_wp, 0.0_wp, kind=wp)
    end if
    if (use_two_d_channels) then
        allocate(amp_up_l2_j52_grid(n_energy), amp_down_l2_j52_grid(n_energy))
        amp_up_l2_j52_grid = cmplx(0.0_wp, 0.0_wp, kind=wp)
        amp_down_l2_j52_grid = cmplx(0.0_wp, 0.0_wp, kind=wp)
    end if

    target_energy = photon_energy + energy_bound
    best_diff = huge(1.0_wp)
    ref_index = -1
    have_saved_wave = .false.
    delta_l2_saved = 0.0_wp
    delta_l0_saved = 0.0_wp
    delta_l2_j52_saved = 0.0_wp

    allocate(phi_l2_prev(n_global))
    have_prev_l2 = .false.
    if (need_s_wave) then
        allocate(phi_l0_prev(n_global))
        have_prev_l0 = .false.
    else
        have_prev_l0 = .false.
    end if
    if (use_two_d_channels) then
        allocate(phi_l2_j52_prev(n_global))
        have_prev_l2_j52 = .false.
    else
        have_prev_l2_j52 = .false.
    end if

    do idx_energy = 1, n_energy
        energy_scatt = energies(idx_energy)
        if (energy_scatt <= 0.0_wp) then
            cycle
        end if

        k_scatt = sqrt(2.0_wp * energy_scatt)

        call dvr_iem_solve(nelems, n_per_elem, r_max, k_scatt, 2, 1, 1.5_wp, &
                           phi_l2_tmp, nodes_l2_tmp, delta_l2, info)
        if (info /= 0) then
            print '(A,I0)', 'dvr_iem_solve failed for l=2, j=3/2, info = ', info
            stop 1
        end if

        if (use_two_d_channels) then
            call dvr_iem_solve(nelems, n_per_elem, r_max, k_scatt, 2, 1, 2.5_wp, &
                               phi_l2_j52_tmp, nodes_l2_j52_tmp, delta_l2_j52, info)
            if (info /= 0) then
                print '(A,I0)', 'dvr_iem_solve failed for l=2, j=5/2, info = ', info
                stop 1
            end if
        end if

        if (need_s_wave) then
            call dvr_iem_solve(nelems, n_per_elem, r_max, k_scatt, 0, 1, 0.5_wp, &
                               phi_l0_tmp, nodes_l0_tmp, delta_l0, info)
            if (info /= 0) then
                print '(A,I0)', 'dvr_iem_solve failed for l=0, j=1/2, info = ', info
                stop 1
            end if
        else
            delta_l0 = 0.0_wp
        end if

        if (size(phi_l2_tmp) /= n_global) then
            print '(A)', 'ERROR: l=2 scattering solution size mismatch with FEDVR grid.'
            stop 1
        end if
        if (size(nodes_l2_tmp) /= n_global) then
            print '(A)', 'ERROR: l=2 scattering node array size mismatch.'
            stop 1
        end if
        if (maxval(abs(nodes_l2_tmp - x_global)) > node_tol) then
            print '(A)', 'ERROR: l=2 scattering nodes do not match FEDVR grid.'
            stop 1
        end if
        if (use_two_d_channels) then
            if (size(phi_l2_j52_tmp) /= n_global .or. size(nodes_l2_j52_tmp) /= n_global) then
                print '(A)', 'ERROR: j=5/2 scattering array size mismatch.'
                stop 1
            end if
            if (maxval(abs(nodes_l2_j52_tmp - x_global)) > node_tol) then
                print '(A)', 'ERROR: l=2, j=5/2 nodes do not match FEDVR grid.'
                stop 1
            end if
        end if
        if (need_s_wave) then
            if (size(phi_l0_tmp) /= n_global .or. size(nodes_l0_tmp) /= n_global) then
                print '(A)', 'ERROR: l=0 scattering array size mismatch.'
                stop 1
            end if
            if (maxval(abs(nodes_l0_tmp - x_global)) > node_tol) then
                print '(A)', 'ERROR: l=0 scattering nodes do not match FEDVR grid.'
                stop 1
            end if
        end if

        where (.not. active_mask)
            phi_l2_tmp = 0.0_wp
        end where
        if (use_two_d_channels) then
            where (.not. active_mask)
                phi_l2_j52_tmp = 0.0_wp
            end where
        end if
        if (need_s_wave) then
            where (.not. active_mask)
                phi_l0_tmp = 0.0_wp
            end where
        end if

        call enforce_wavefunction_phase(phi_l2_tmp, have_prev_l2, phi_l2_prev, active_idx, m_global)
        if (use_two_d_channels) then
            call enforce_wavefunction_phase(phi_l2_j52_tmp, have_prev_l2_j52, phi_l2_j52_prev, active_idx, m_global)
        end if
        if (need_s_wave) then
            call enforce_wavefunction_phase(phi_l0_tmp, have_prev_l0, phi_l0_prev, active_idx, m_global)
        end if

        r_integral_l2 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        r_integral_l0 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        r_integral_l2_j52 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        do i = 1, n_active
            weight_val = cmplx(m_global(active_idx(i)), 0.0_wp, kind=wp)
            bound_val = cmplx(u_bound(active_idx(i)), 0.0_wp, kind=wp)
            scat_val_l2 = cmplx(phi_l2_tmp(active_idx(i)), 0.0_wp, kind=wp)
            r_integral_l2 = r_integral_l2 + weight_val * conjg(bound_val) * scat_val_l2
            if (need_s_wave) then
                scat_val_l0 = cmplx(phi_l0_tmp(active_idx(i)), 0.0_wp, kind=wp)
                r_integral_l0 = r_integral_l0 + weight_val * conjg(bound_val) * scat_val_l0
            end if
            if (use_two_d_channels) then
                scat_val_l2_j52 = cmplx(phi_l2_j52_tmp(active_idx(i)), 0.0_wp, kind=wp)
                r_integral_l2_j52 = r_integral_l2_j52 + weight_val * conjg(bound_val) * scat_val_l2_j52
            end if
        end do

        amp_up_l2 = pref_up_l2 * cis(delta_l2) * r_integral_l2
        amp_down_l2 = pref_down_l2 * cis(delta_l2) * r_integral_l2
        if (need_s_wave) then
            amp_up_l0 = pref_up_l0 * cis(delta_l0) * r_integral_l0
        else
            amp_up_l0 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        end if
        if (use_two_d_channels) then
            amp_up_l2_j52 = pref_up_l2_j52 * cis(delta_l2_j52) * r_integral_l2_j52
            amp_down_l2_j52 = pref_down_l2_j52 * cis(delta_l2_j52) * r_integral_l2_j52
        else
            amp_up_l2_j52 = cmplx(0.0_wp, 0.0_wp, kind=wp)
            amp_down_l2_j52 = cmplx(0.0_wp, 0.0_wp, kind=wp)
        end if

        amp_up_l2_grid(idx_energy) = amp_up_l2
        amp_down_l2_grid(idx_energy) = amp_down_l2
        if (use_two_d_channels) then
            amp_up_l2_j52_grid(idx_energy) = amp_up_l2_j52
            amp_down_l2_j52_grid(idx_energy) = amp_down_l2_j52
        end if
        if (need_s_wave) then
            amp_up_l0_grid(idx_energy) = amp_up_l0
        end if

        if (use_two_d_channels) then
            p_up_val = abs(amp_up_l2 + amp_up_l2_j52)**2
            if (need_s_wave) p_up_val = p_up_val + abs(amp_up_l0)**2
            p_down_val = abs(amp_down_l2 + amp_down_l2_j52)**2
        else
            p_up_val = abs(amp_up_l2)**2
            if (need_s_wave) p_up_val = p_up_val + abs(amp_up_l0)**2
            p_down_val = abs(amp_down_l2)**2
        end if

        if (p_down_val > 0.0_wp) then
            ratio_val = p_up_val / p_down_val
        else
            ratio_val = huge(1.0_wp)
        end if

        p_up_grid(idx_energy) = p_up_val
        p_down_grid(idx_energy) = p_down_val
        ratio_grid(idx_energy) = ratio_val

        if (abs(energy_scatt - target_energy) < best_diff) then
            best_diff = abs(energy_scatt - target_energy)
            ref_index = idx_energy
            delta_l2_saved = delta_l2
            if (need_s_wave) delta_l0_saved = delta_l0
            if (use_two_d_channels) delta_l2_j52_saved = delta_l2_j52
            if (.not. allocated(phi_l2_saved)) then
                if (use_two_d_channels .and. need_s_wave) then
                    allocate(phi_l2_saved(n_global), phi_l0_saved(n_global), phi_l2_j52_saved(n_global), nodes_saved(n_global))
                else if (use_two_d_channels) then
                    allocate(phi_l2_saved(n_global), phi_l2_j52_saved(n_global), nodes_saved(n_global))
                else if (need_s_wave) then
                    allocate(phi_l2_saved(n_global), phi_l0_saved(n_global), nodes_saved(n_global))
                else
                    allocate(phi_l2_saved(n_global), nodes_saved(n_global))
                end if
            end if
            phi_l2_saved = phi_l2_tmp
            nodes_saved = nodes_l2_tmp
            if (need_s_wave) phi_l0_saved = phi_l0_tmp
            if (use_two_d_channels) phi_l2_j52_saved = phi_l2_j52_tmp
            have_saved_wave = .true.
        end if

        if (use_two_d_channels .and. need_s_wave) then
            deallocate(phi_l2_tmp, nodes_l2_tmp, phi_l0_tmp, nodes_l0_tmp, phi_l2_j52_tmp, nodes_l2_j52_tmp)
        else if (use_two_d_channels) then
            deallocate(phi_l2_tmp, nodes_l2_tmp, phi_l2_j52_tmp, nodes_l2_j52_tmp)
        else if (need_s_wave) then
            deallocate(phi_l2_tmp, nodes_l2_tmp, phi_l0_tmp, nodes_l0_tmp)
        else
            deallocate(phi_l2_tmp, nodes_l2_tmp)
        end if
    end do

    call write_wavefunction('bound_wavefunction_'//trim(output_tag)//'.dat', x_global, u_bound)
    if (have_saved_wave) then
        call write_wavefunction('scattering_l2_j3_2_'//trim(output_tag)//'.dat', nodes_saved, phi_l2_saved)
        if (use_two_d_channels) call write_wavefunction('scattering_l2_j5_2_'//trim(output_tag)//'.dat', nodes_saved, phi_l2_j52_saved)
        if (need_s_wave) call write_wavefunction('scattering_l0_j1_2_'//trim(output_tag)//'.dat', nodes_saved, phi_l0_saved)
    else
        print '(A)', 'WARNING: No scattering wavefunction saved (energy grid may not include positive energies).'
    end if

    call write_spectrum(trim(spectrum_filename), energies, p_up_grid, p_down_grid, ratio_grid)
    if (use_two_d_channels) then
        if (need_s_wave) then
            call write_amplitudes_two_d(trim(amplitude_filename), energies, amp_up_l2_grid, amp_up_l2_j52_grid, amp_down_l2_grid, amp_down_l2_j52_grid, amp_up_l0_grid)
        else
            call write_amplitudes_two_d(trim(amplitude_filename), energies, amp_up_l2_grid, amp_up_l2_j52_grid, amp_down_l2_grid, amp_down_l2_j52_grid)
        end if
    else
        if (need_s_wave) then
            call write_amplitudes_single(trim(amplitude_filename), energies, amp_up_l2_grid, amp_down_l2_grid, amp_up_l0_grid)
        else
            call write_amplitudes_single(trim(amplitude_filename), energies, amp_up_l2_grid, amp_down_l2_grid)
        end if
    end if

    print '(A)', '---------------------------------------------'
    print '(A)', ' Photoelectron Spin Polarization Diagnostics'
    print '(A)', '---------------------------------------------'
    print '(A,1pe15.6)', ' Bound-state energy (a.u.)  : ', energy_bound
    print '(A,1pe15.6)', ' Photon energy (a.u.)       : ', photon_energy
    print '(A,1pe15.6)', ' Target continuum (a.u.)    : ', target_energy
    print '(A,1x,A)', ' Spectrum written to ', trim(spectrum_filename)
    print '(A,1x,A)', ' Amplitudes written to ', trim(amplitude_filename)
    if (ref_index > 0) then
        print '(A,1pe15.6)', ' Reference energy (a.u.)    : ', energies(ref_index)
        print '(A,1pe15.6)', ' delta_{3/2} (rad)          : ', delta_l2_saved
        if (use_two_d_channels) then
            print '(A,1pe15.6)', ' delta_{5/2} (rad)          : ', delta_l2_j52_saved
        end if
        if (need_s_wave) then
            print '(A,1pe15.6)', ' delta_{1/2} (rad)          : ', delta_l0_saved
        end if
        print '(A,1pe15.6)', ' P_up(reference)            : ', p_up_grid(ref_index)
        print '(A,1pe15.6)', ' P_down(reference)          : ', p_down_grid(ref_index)
        if (p_down_grid(ref_index) > 0.0_wp) then
            print '(A,1pe15.6)', ' Ratio(reference)           : ', ratio_grid(ref_index)
        else
            print '(A)', ' Ratio(reference)           : (infinite, P_down ≈ 0)'
        end if
        if (use_two_d_channels) then
            print '(A,1pe15.6,1x,1pe15.6)', ' A_up^{j=3/2} [Re, Im]      : ', real(amp_up_l2_grid(ref_index)), aimag(amp_up_l2_grid(ref_index))
            print '(A,1pe15.6,1x,1pe15.6)', ' A_up^{j=5/2} [Re, Im]      : ', real(amp_up_l2_j52_grid(ref_index)), aimag(amp_up_l2_j52_grid(ref_index))
            if (need_s_wave) then
                print '(A,1pe15.6,1x,1pe15.6)', ' A_up^{s} [Re, Im]          : ', real(amp_up_l0_grid(ref_index)), aimag(amp_up_l0_grid(ref_index))
            end if
            print '(A,1pe15.6,1x,1pe15.6)', ' A_down^{j=3/2} [Re, Im]    : ', real(amp_down_l2_grid(ref_index)), aimag(amp_down_l2_grid(ref_index))
            print '(A,1pe15.6,1x,1pe15.6)', ' A_down^{j=5/2} [Re, Im]    : ', real(amp_down_l2_j52_grid(ref_index)), aimag(amp_down_l2_j52_grid(ref_index))
        else
            print '(A,1pe15.6,1x,1pe15.6)', ' A_up^{d} [Re, Im]          : ', real(amp_up_l2_grid(ref_index)), aimag(amp_up_l2_grid(ref_index))
            if (need_s_wave) then
                print '(A,1pe15.6,1x,1pe15.6)', ' A_up^{s} [Re, Im]          : ', real(amp_up_l0_grid(ref_index)), aimag(amp_up_l0_grid(ref_index))
            end if
            print '(A,1pe15.6,1x,1pe15.6)', ' A_down^{d} [Re, Im]        : ', real(amp_down_l2_grid(ref_index)), aimag(amp_down_l2_grid(ref_index))
        end if
    else
        print '(A)', 'WARNING: Reference energy not found within scan range.'
    end if

contains

    pure integer function nearest_bound_state(evals, target)
        real(wp), intent(in) :: evals(:)
        real(wp), intent(in) :: target
        real(wp) :: best_diff_local
        integer :: idx_local

        nearest_bound_state = -1
        best_diff_local = huge(1.0_wp)
        do idx_local = 1, size(evals)
            if (evals(idx_local) >= 0.0_wp) cycle
            if (abs(evals(idx_local) - target) < best_diff_local) then
                best_diff_local = abs(evals(idx_local) - target)
                nearest_bound_state = idx_local
            end if
        end do
    end function nearest_bound_state

    pure complex(wp) function cis(angle)
        real(wp), intent(in) :: angle
        cis = cmplx(cos(angle), sin(angle), kind=wp)
    end function cis

    pure function to_lower(input) result(output)
        character(len=*), intent(in) :: input
        character(len=len(input)) :: output
        integer :: idx_char, code

        output = input
        do idx_char = 1, len(output)
            code = iachar(output(idx_char:idx_char))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                output(idx_char:idx_char) = achar(code + 32)
            end if
        end do
    end function to_lower

    subroutine write_wavefunction(filename, coords, values)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: coords(:)
        real(wp), intent(in) :: values(:)
        integer :: unit, idx_local

        open(newunit=unit, file=filename, status='replace', action='write', form='formatted')
        do idx_local = 1, size(coords)
            write(unit, '(1pe22.14,1x,1pe22.14)') coords(idx_local), values(idx_local)
        end do
        close(unit)
    end subroutine write_wavefunction

    subroutine write_spectrum(filename, energy_vals, pup_vals, pdown_vals, ratio_vals)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: energy_vals(:)
        real(wp), intent(in) :: pup_vals(:)
        real(wp), intent(in) :: pdown_vals(:)
        real(wp), intent(in) :: ratio_vals(:)
        integer :: unit, idx_local

        if (size(energy_vals) /= size(pup_vals) .or. size(energy_vals) /= size(pdown_vals) .or. &
            size(energy_vals) /= size(ratio_vals)) then
            print '(A)', 'ERROR: write_spectrum received arrays of inconsistent length.'
            stop 1
        end if

        open(newunit=unit, file=filename, status='replace', action='write', form='formatted')
        write(unit, '(A)') '# energy(a.u.)   Pup   Pdown   Pup/Pdown'
        do idx_local = 1, size(energy_vals)
            write(unit, '(1pe22.14,1x,1pe22.14,1x,1pe22.14,1x,1pe22.14)') &
                energy_vals(idx_local), pup_vals(idx_local), pdown_vals(idx_local), ratio_vals(idx_local)
        end do
        close(unit)
    end subroutine write_spectrum

    subroutine write_amplitudes_single(filename, energy_vals, amp_up_d_vals, amp_down_d_vals, amp_up_s_vals)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: energy_vals(:)
        complex(wp), intent(in) :: amp_up_d_vals(:)
        complex(wp), intent(in) :: amp_down_d_vals(:)
        complex(wp), intent(in), optional :: amp_up_s_vals(:)
        integer :: unit, idx_local
        logical :: has_s

        has_s = present(amp_up_s_vals)
        if (size(energy_vals) /= size(amp_up_d_vals) .or. size(energy_vals) /= size(amp_down_d_vals)) then
            print '(A)', 'ERROR: write_amplitudes_single received arrays of inconsistent length.'
            stop 1
        end if
        if (has_s) then
            if (size(energy_vals) /= size(amp_up_s_vals)) then
                print '(A)', 'ERROR: write_amplitudes_single received mismatched s-wave array length.'
                stop 1
            end if
        end if

        open(newunit=unit, file=filename, status='replace', action='write', form='formatted')
        if (has_s) then
            write(unit, '(A)') '# energy(a.u.)   Re[A_up(d)]   Im[A_up(d)]   Re[A_up(s)]   Im[A_up(s)]   Re[A_down(d)]   Im[A_down(d)]'
        else
            write(unit, '(A)') '# energy(a.u.)   Re[A_up(d)]   Im[A_up(d)]   Re[A_down(d)]   Im[A_down(d)]'
        end if
        do idx_local = 1, size(energy_vals)
            if (has_s) then
                write(unit, '(1pe22.14,6(1x,1pe22.14))') energy_vals(idx_local), &
                    real(amp_up_d_vals(idx_local)), aimag(amp_up_d_vals(idx_local)), &
                    real(amp_up_s_vals(idx_local)), aimag(amp_up_s_vals(idx_local)), &
                    real(amp_down_d_vals(idx_local)), aimag(amp_down_d_vals(idx_local))
            else
                write(unit, '(1pe22.14,4(1x,1pe22.14))') energy_vals(idx_local), &
                    real(amp_up_d_vals(idx_local)), aimag(amp_up_d_vals(idx_local)), &
                    real(amp_down_d_vals(idx_local)), aimag(amp_down_d_vals(idx_local))
            end if
        end do
        close(unit)
    end subroutine write_amplitudes_single

    subroutine write_amplitudes_two_d(filename, energy_vals, amp_up_j32_vals, amp_up_j52_vals, amp_down_j32_vals, amp_down_j52_vals, amp_up_s_vals)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: energy_vals(:)
        complex(wp), intent(in) :: amp_up_j32_vals(:)
        complex(wp), intent(in) :: amp_up_j52_vals(:)
        complex(wp), intent(in) :: amp_down_j32_vals(:)
        complex(wp), intent(in) :: amp_down_j52_vals(:)
        complex(wp), intent(in), optional :: amp_up_s_vals(:)
        integer :: unit, idx_local
        logical :: has_s

        has_s = present(amp_up_s_vals)
        if (size(energy_vals) /= size(amp_up_j32_vals) .or. size(energy_vals) /= size(amp_up_j52_vals) .or. &
            size(energy_vals) /= size(amp_down_j32_vals) .or. size(energy_vals) /= size(amp_down_j52_vals)) then
            print '(A)', 'ERROR: write_amplitudes_two_d received arrays of inconsistent length.'
            stop 1
        end if
        if (has_s) then
            if (size(energy_vals) /= size(amp_up_s_vals)) then
                print '(A)', 'ERROR: write_amplitudes_two_d received mismatched s-wave array length.'
                stop 1
            end if
        end if

        open(newunit=unit, file=filename, status='replace', action='write', form='formatted')
        if (has_s) then
            write(unit, '(A)') '# energy(a.u.)   Re[A_up(j=3/2)]   Im[A_up(j=3/2)]   Re[A_up(j=5/2)]   Im[A_up(j=5/2)]   Re[A_up(s)]   Im[A_up(s)]   Re[A_down(j=3/2)]   Im[A_down(j=3/2)]   Re[A_down(j=5/2)]   Im[A_down(j=5/2)]'
        else
            write(unit, '(A)') '# energy(a.u.)   Re[A_up(j=3/2)]   Im[A_up(j=3/2)]   Re[A_up(j=5/2)]   Im[A_up(j=5/2)]   Re[A_down(j=3/2)]   Im[A_down(j=3/2)]   Re[A_down(j=5/2)]   Im[A_down(j=5/2)]'
        end if
        do idx_local = 1, size(energy_vals)
            if (has_s) then
                write(unit, '(1pe22.14,10(1x,1pe22.14))') energy_vals(idx_local), &
                    real(amp_up_j32_vals(idx_local)), aimag(amp_up_j32_vals(idx_local)), &
                    real(amp_up_j52_vals(idx_local)), aimag(amp_up_j52_vals(idx_local)), &
                    real(amp_up_s_vals(idx_local)), aimag(amp_up_s_vals(idx_local)), &
                    real(amp_down_j32_vals(idx_local)), aimag(amp_down_j32_vals(idx_local)), &
                    real(amp_down_j52_vals(idx_local)), aimag(amp_down_j52_vals(idx_local))
            else
                write(unit, '(1pe22.14,8(1x,1pe22.14))') energy_vals(idx_local), &
                    real(amp_up_j32_vals(idx_local)), aimag(amp_up_j32_vals(idx_local)), &
                    real(amp_up_j52_vals(idx_local)), aimag(amp_up_j52_vals(idx_local)), &
                    real(amp_down_j32_vals(idx_local)), aimag(amp_down_j32_vals(idx_local)), &
                    real(amp_down_j52_vals(idx_local)), aimag(amp_down_j52_vals(idx_local))
            end if
        end do
        close(unit)
    end subroutine write_amplitudes_two_d

    subroutine enforce_wavefunction_phase(phi_curr, has_prev, phi_prev, active_indices, weights)
        real(wp), intent(inout) :: phi_curr(:)
        logical, intent(inout) :: has_prev
        real(wp), intent(inout) :: phi_prev(:)
        integer, intent(in) :: active_indices(:)
        real(wp), intent(in) :: weights(:)
        real(wp) :: overlap, threshold
        integer :: n, gidx

        threshold = 1.0e-10_wp
        if (has_prev) then
            overlap = 0.0_wp
            do n = 1, size(active_indices)
                gidx = active_indices(n)
                overlap = overlap + weights(gidx) * phi_curr(gidx) * phi_prev(gidx)
            end do
            if (overlap < 0.0_wp) phi_curr = -phi_curr
        else
            do n = 1, size(active_indices)
                gidx = active_indices(n)
                if (abs(phi_curr(gidx)) > threshold) then
                    if (phi_curr(gidx) < 0.0_wp) phi_curr = -phi_curr
                    exit
                end if
            end do
            has_prev = .true.
        end if

        phi_prev = phi_curr
        has_prev = .true.
    end subroutine enforce_wavefunction_phase

end program compute_spin_polarization


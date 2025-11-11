module dvr_integral_equation_methods
    !--------------------------------------------------------------------------------------
    ! DVR Integral Equation Method module (Gonzales et al. formulation)
    !
    ! Method references:
    !   - R.A. Gonzales et al., J. Comput. Phys. 134 (1997) 134–149 (IEM)
    !   - Barnett COUL90 continued–fraction Coulomb functions (wrapped via module COULFG)
    !   - FEDVR infrastructure in this repository for Gauss–Lobatto nodes & weights
    !
    ! 概要 (Japanese summary):
    !   本モジュールは有限要素DVR(FEDVR)上で連続スペクトル向け放射シュレーディンガー方程式を
    !   Gonzales型積分方程式(IEM)として解き、散乱位相シフトをクーロン関数 (COUL90) とマッチングして
    !   抽出する。各要素で局所積分方程式 (y_i, z_i) を離散化し DGESV で解き、
    !   要素間結合係数 cy, cz, sy, sz を用いて (A_i, B_i) の 2×2 ブロック三重対角系を構築する。
    !   その後グローバル波動関数 φ(r) を再構成し、外端 r=T で Coulomb 関数 F_l, G_l による境界マッチング
    !   から位相シフト δ_l を計算する。
    !
    !    自由粒子,箱型ポテンシャルで 10^(-13)程度の位相のずれの精度をチェック済
    !    修正・改造はご自由にどうぞ
    !
    ! Design notes:
    !   - Local IE discretisation splits integrals at r_j relative ordering (j>i or j<i).
    !   - Diagonal correction enforces (I - K_i) structure for stable DGESV solve.
    !   - Block system currently uses a simplified pattern; future refinement can follow
    !     exact Eq.(3.10/3.11) row operations for a minimal tridiagonal form.
    !   - Phase shift uses tan δ = (φ_T F'_l - φ'_T F_l)/(φ'_T G_l - φ_T G'_l).
    !   - Potential includes centrifugal term l(l+1)/r^2 plus chosen atomic potential.
    !
    ! TODO / Future improvements:
    !   * Implement adaptive element sizing based on |V(r)| falloff to reduce nelems.
    !   * Provide optional iterative refinement of local (I-K_i) systems.
    !   * Support coupled-channel (matrix V) generalisation.
    !   * Cache trigonometric factors sin(k r_j), cos(k r_j) for performance.
    !   * Add unit tests verifying convergence vs known analytic Coulomb phase shifts.
    !---------------------------written by Kawaii_Bokuchin_Puyuyu(2025/11/11)----------------

    use, intrinsic :: iso_fortran_env, only : wp => real64
    use finite_element_dvr
    use COULFG

    implicit none
    private

    public :: dvr_iem_solve

    real(wp), parameter :: two_pi = 6.28318530717958647692_wp

    ! LAPACK (no explicit interface to avoid compiler-specific module deps)
    external :: dgesv

contains

    subroutine dvr_iem_solve(nelems, nloc, r_max, k_wave, ell, key_soi, j_tot, solution, nodes, phase_shift, info, potential_model)
        ! High-level driver for the Gonzales-style DVR-IEM solver.
        ! nelems   : number of FEDVR elements
        ! nloc     : local Lobatto points per element
        ! r_max    : outer truncation radius T
        ! k_wave   : scattering wavenumber k (>0)
        ! ell      : orbital quantum number l
        ! key_soi  : spin-orbit switch passed to finite_element_dvr potentials
        ! j_tot    : total angular momentum (only used when key_soi=1)
        ! solution : allocatable output array φ(r) on global DVR grid
        ! nodes    : allocatable output array of radial nodes (same length as solution)
        ! phase_shift : computed Coulomb phase shift δ_l (radians)
    ! info     : status flag (0 on success)
    ! potential_model (optional) : 0=atomic Coulomb core (default), 1=zero potential

        integer, intent(in) :: nelems, nloc, ell, key_soi
        integer, intent(in), optional :: potential_model
        real(wp), intent(in) :: r_max, k_wave, j_tot
        real(wp), allocatable, intent(out) :: solution(:), nodes(:)
        real(wp), intent(out) :: phase_shift
        integer, intent(out) :: info

        integer :: ierr, N_global, dim_block, lda, ldb
        integer :: e, i, j, p, rowA, rowB
        integer, allocatable :: elem_map(:,:), ipiv(:), ipiv_local(:)
        real(wp) :: r_left, r_right, h_elem, scale_deriv
        real(wp) :: weight_j, kr_i, kr_j, val
        real(wp) :: phi_T, phi_prime_T, eta_tail, rho_tail
        real(wp) :: centrifugal
        real(wp), allocatable :: x_ref(:), w_ref(:), bary_w(:)
        real(wp), allocatable :: D_ref(:,:), local_matrix(:,:), rhs(:,:), rhs_copy(:,:)
        real(wp), allocatable :: y_vals(:,:), z_vals(:,:)
        real(wp), allocatable :: cy(:), cz(:), sy(:), sz(:)
        real(wp), allocatable :: local_nodes(:), local_weights(:)
        real(wp), allocatable :: V_global(:), V_eff(:)
        real(wp), allocatable :: x_global(:), M_global(:)
        real(wp), allocatable :: block_matrix(:,:), block_rhs(:)
        real(wp), allocatable :: phi_local(:), counts(:)
        real(wp) :: proj_num, proj_den
        real(wp) :: numer, denom
        integer :: ifail
        integer :: pot_model
        real(wp), allocatable, target :: fc(:), gc(:), fcp(:), gcp(:)
        real(wp), pointer :: fc_lb0(:), gc_lb0(:), fcp_lb0(:), gcp_lb0(:)
        integer :: dirichlet_idx

        info = 0
        if (nelems < 1 .or. nloc < 2) then
            info = -1
            return
        end if
        if (k_wave <= 0.0_wp) then
            info = -2
            return
        end if

    !--------------------------------------------------------------
    ! 1. 参照節点 / 重みの確保 (Gauss-Lobatto from FEDVR util)
    !--------------------------------------------------------------
        allocate(x_ref(nloc), w_ref(nloc), stat=ierr)
        if (ierr /= 0) then
            info = -10
            return
        end if

    ! Global grid のサイズ: nelems*(nloc-1)+1 (端点共有)
        allocate(x_global(nelems * (nloc - 1) + 1), M_global(nelems * (nloc - 1) + 1), stat=ierr)
        if (ierr /= 0) then
            info = -11
            return
        end if

    ! elem_map(e,i) → global index
        allocate(elem_map(nelems, nloc), stat=ierr)
        if (ierr /= 0) then
            info = -12
            return
        end if

    ! 2. グローバル節点生成
        call build_global_grid(nelems, nloc, x_ref, w_ref, 0.0_wp, r_max, N_global, x_global, M_global, elem_map, ierr)
        if (ierr /= 0) then
            info = -13
            return
        end if

        allocate(bary_w(nloc), D_ref(nloc, nloc), stat=ierr)
        if (ierr /= 0) then
            info = -14
            return
        end if

    ! 3. 1要素参照上の補間用バリセントリック重みと微分行列
        call compute_barycentric_weights(nloc, x_ref, bary_w)
        call compute_differentiation_matrix(nloc, x_ref, bary_w, D_ref)

        allocate(V_global(N_global), V_eff(N_global), stat=ierr)
        if (ierr /= 0) then
            info = -15
            return
        end if

    ! 4. 原子＋SOI ポテンシャル (V_global) 構築 / モデル切替
        pot_model = 0
        if (present(potential_model)) pot_model = potential_model

        select case (pot_model)
        case (0)
            call build_potential_diag(N_global, x_global, V_global, key_soi, j_tot, ell, ierr)
            if (ierr /= 0) then
                info = -16
                return
            end if
        case (1)
            V_global = 0.0_wp
        case default
            info = -30
            return
        end select

    ! 5. 有効ポテンシャル V_eff = V + l(l+1)/r^2
        do i = 1, N_global
            val = x_global(i)
            if (val > 0.0_wp) then
                centrifugal = real(ell * (ell + 1), wp) / (val * val)
            else
                centrifugal = 0.0_wp
            end if
            V_eff(i) = V_global(i) + centrifugal
        end do

    ! 6. 出力配列確保
        allocate(solution(N_global), nodes(N_global), counts(N_global), stat=ierr)
        if (ierr /= 0) then
            info = -17
            return
        end if
        solution = 0.0_wp
        counts = 0.0_wp
        nodes = x_global

    ! y_i, z_i 局所解格納領域
        allocate(y_vals(nloc, nelems), z_vals(nloc, nelems), stat=ierr)
        if (ierr /= 0) then
            info = -18
            return
        end if

        allocate(local_nodes(nloc), local_weights(nloc), stat=ierr)
        if (ierr /= 0) then
            info = -19
            return
        end if

    ! (I - K_i) 行列 + RHS(sin, cos)
        allocate(local_matrix(nloc, nloc), rhs(nloc, 2), rhs_copy(nloc, 2), stat=ierr)
        if (ierr /= 0) then
            info = -20
            return
        end if

    ! 要素間結合係数配列
        allocate(cy(nelems), cz(nelems), sy(nelems), sz(nelems), stat=ierr)
        if (ierr /= 0) then
            info = -21
            return
        end if

        cy = 0.0_wp; cz = 0.0_wp; sy = 0.0_wp; sz = 0.0_wp

    ! DGESV pivot 配列
        allocate(ipiv_local(nloc), stat=ierr)
        if (ierr /= 0) then
            info = -22
            return
        end if

        dirichlet_idx = elem_map(1, 1)

        do e = 1, nelems
            do j = 1, nloc
                local_nodes(j) = x_global(elem_map(e, j))
            end do
            r_left = local_nodes(1)
            r_right = local_nodes(nloc)
            h_elem = r_right - r_left
            if (h_elem <= 0.0_wp) then
                info = -23
                return
            end if

            do j = 1, nloc
                local_weights(j) = 0.5_wp * h_elem * w_ref(j)
            end do

            ! 7. 局所 IE 離散化: K_i 構築 (sin/cos カーネル分割)
            local_matrix = 0.0_wp

            ! 行列の I - K_i 形式化 (対角・非対角の符号反転)
            do i = 1, nloc
                kr_i = modulo(k_wave * local_nodes(i), two_pi)
                do j = 1, nloc
                    if (e == 1 .and. (i == 1 .or. j == 1)) cycle
                    kr_j = modulo(k_wave * local_nodes(j), two_pi)
                    val = V_eff(elem_map(e, j))
                    weight_j = local_weights(j)
                    if (j > i) then
                        local_matrix(i, j) = local_matrix(i, j) + (sin(kr_i) * cos(kr_j) * val * weight_j) / k_wave
                    else if (j < i) then
                        local_matrix(i, j) = local_matrix(i, j) + (cos(kr_i) * sin(kr_j) * val * weight_j) / k_wave
                    else
                        local_matrix(i, j) = local_matrix(i, j) + &
                            ((sin(kr_i) * cos(kr_j) + cos(kr_i) * sin(kr_j)) * val * weight_j) / k_wave
                    end if
                end do
            end do

            do i = 1, nloc
                local_matrix(i, i) = 1.0_wp - local_matrix(i, i)
                do j = 1, nloc
                    if (j /= i) local_matrix(i, j) = -local_matrix(i, j)
                end do
            end do

            do i = 1, nloc
                kr_i = k_wave * local_nodes(i)
                rhs(i, 1) = sin(kr_i)
                rhs(i, 2) = cos(kr_i)
            end do

            if (e == 1) then
                local_matrix(1, :) = 0.0_wp
                local_matrix(:, 1) = 0.0_wp
                local_matrix(1, 1) = 1.0_wp
                rhs(1, :) = 0.0_wp
            end if

            ! 8. DGESV で (I-K_i) y_i = sin, (I-K_i) z_i = cos を同時解
            rhs_copy = rhs
            call dgesv(nloc, 2, local_matrix, nloc, ipiv_local, rhs_copy, nloc, ierr)
            if (ierr /= 0) then
                info = -24
                return
            end if
            y_vals(:, e) = rhs_copy(:, 1)
            z_vals(:, e) = rhs_copy(:, 2)

            if (e == 1) then
                y_vals(1, e) = 0.0_wp
                z_vals(1, e) = 0.0_wp
            end if

            ! 9. cy,cz,sy,sz の加重求積 (DVR 重み使用)
            do j = 1, nloc
                val = V_eff(elem_map(e, j)) * local_weights(j) / k_wave
                kr_j = k_wave * local_nodes(j)
                if (e == 1 .and. j == 1) cycle
                cy(e) = cy(e) + cos(kr_j) * val * y_vals(j, e)
                cz(e) = cz(e) + cos(kr_j) * val * z_vals(j, e)
                sy(e) = sy(e) + sin(kr_j) * val * y_vals(j, e)
                sz(e) = sz(e) + sin(kr_j) * val * z_vals(j, e)
            end do
        end do

        dim_block = 2 * nelems
        allocate(block_matrix(dim_block, dim_block), block_rhs(dim_block), ipiv(dim_block), stat=ierr)
        if (ierr /= 0) then
            info = -25
            return
        end if

    ! 10. 2×2 ブロック系の構築
    !     指定仕様: 
    !       A_i = 1 - Σ_{p=i+1}^{m} [ c_y(p) A_p + c_z(p) B_p ]
    !       B_i = - Σ_{p=1}^{i-1} [ s_y(p) A_p + s_z(p) B_p ]
    !     → 行列表現: A_i + Σ_{p=i+1}^{m} (...) = 1,  B_i + Σ_{p=1}^{i-1} (...) = 0
        block_matrix = 0.0_wp
        block_rhs = 0.0_wp

        do e = 1, nelems
            rowA = 2 * e - 1
            rowB = 2 * e

            ! Row for A_i
            block_matrix(rowA, rowA) = 1.0_wp
            do p = e + 1, nelems
                block_matrix(rowA, 2 * p - 1) = block_matrix(rowA, 2 * p - 1) + cy(p)
                block_matrix(rowA, 2 * p)     = block_matrix(rowA, 2 * p)     + cz(p)
            end do
            block_rhs(rowA) = 1.0_wp

            ! Row for B_i
            block_matrix(rowB, rowB) = 1.0_wp
            do p = 1, e - 1
                block_matrix(rowB, 2 * p - 1) = block_matrix(rowB, 2 * p - 1) + sy(p)
                block_matrix(rowB, 2 * p)     = block_matrix(rowB, 2 * p)     + sz(p)
            end do
            block_rhs(rowB) = 0.0_wp
        end do

        lda = dim_block
        ldb = dim_block
    ! 11. グローバルブロック系解法
        call dgesv(dim_block, 1, block_matrix, lda, ipiv, block_rhs, ldb, ierr)
        if (ierr /= 0) then
            info = -26
            return
        end if

        allocate(phi_local(nloc), stat=ierr)
        if (ierr /= 0) then
            info = -27
            return
        end if

    ! 12. φ(r) 再構成 (重複節点は平均化)
        do e = 1, nelems
            do j = 1, nloc
                phi_local(j) = block_rhs(2 * e - 1) * y_vals(j, e) + block_rhs(2 * e) * z_vals(j, e)
            end do
            do j = 1, nloc
                if (e == 1 .and. j == 1) cycle
                i = elem_map(e, j)
                solution(i) = solution(i) + phi_local(j)
                counts(i) = counts(i) + 1.0_wp
            end do
        end do

        solution(dirichlet_idx) = 0.0_wp
        counts(dirichlet_idx) = 1.0_wp

        do i = 1, N_global
            if (counts(i) > 0.0_wp) solution(i) = solution(i) / counts(i)
        end do

    ! 13. 外端値 φ(T) と導関数 φ'(T) 評価
        phi_T = solution(N_global)
        do j = 1, nloc
            phi_local(j) = block_rhs(2 * nelems - 1) * y_vals(j, nelems) + block_rhs(2 * nelems) * z_vals(j, nelems)
        end do
        r_left = local_nodes(1)
        r_right = local_nodes(nloc)
        h_elem = r_right - r_left
        scale_deriv = 2.0_wp / h_elem
        phi_prime_T = 0.0_wp
        do j = 1, nloc
            phi_prime_T = phi_prime_T + D_ref(nloc, j) * phi_local(j)
        end do
        phi_prime_T = phi_prime_T * scale_deriv

        rho_tail = k_wave * nodes(N_global)
        eta_tail = -V_global(N_global) / k_wave

        allocate(fc(0:ell), gc(0:ell), fcp(0:ell), gcp(0:ell), stat=ierr)
        if (ierr /= 0) then
            info = -29
            return
        end if

        fc_lb0(0:ell) => fc
        gc_lb0(0:ell) => gc
        fcp_lb0(0:ell) => fcp
        gcp_lb0(0:ell) => gcp

        call coul90(rho_tail, eta_tail, real(ell, wp), 0, fc_lb0, gc_lb0, fcp_lb0, gcp_lb0, 0, ifail)
        if (ifail /= 0) then
            info = -28
            return
        end if

        do i = 0, ell
            fcp(i) = k_wave * fcp(i)
            gcp(i) = k_wave * gcp(i)
        end do

        proj_num = phi_T * fcp(ell) - phi_prime_T * fc(ell)
        proj_den = phi_prime_T * gc(ell) - phi_T * gcp(ell)
        numer = proj_num
        denom = proj_den
    ! 14. 位相シフト δ_l = atan2(num, den)
        phase_shift = atan2(numer, denom)

        deallocate(phi_local)
        deallocate(block_matrix, block_rhs, ipiv)
        deallocate(cy, cz, sy, sz)
        deallocate(local_matrix, rhs, rhs_copy)
        deallocate(local_nodes, local_weights)
        deallocate(y_vals, z_vals)
        deallocate(ipiv_local)
        deallocate(counts)
        deallocate(fc, gc, fcp, gcp)
        deallocate(V_global, V_eff)
        deallocate(bary_w, D_ref)
        deallocate(elem_map)
        deallocate(x_global, M_global)
        deallocate(x_ref, w_ref)
    end subroutine dvr_iem_solve

end module dvr_integral_equation_methods

! ***********************************************************************
!
!   Copyright (C) 2025 Matthias Fabry & The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

module binary_roche_deformation

   ! computes the stellar structure correction factors fp/ft and the specific
   ! moment of inertia i_rot assuming shellularity in the Roche potential of
   ! a binary star. Follows Fabry, Marchant and Sana, 2022, A&A 661, A123.

   use ieee_arithmetic, only: ieee_is_nan
   use interp_2d_lib_db
   use auto_diff
   use star_def
   use const_def
   use binary_def
   use binary_utils, only : eval_rlobe


   implicit none

   real(dp), parameter :: nudge = 1d-4
   real(dp), pointer :: xvals(:), yvals(:), fpfunc1d(:), ftfunc1d(:), &
         irotfunc1d(:), otherrfunc1d(:), afunc1d(:), psifunc1d(:)
   logical :: inter_ok = .false., dbg = .false.
   integer :: num_xpts, num_ypts

   ! where a point sits in the (lq, ar) tables; see set_roche_stencil
   type roche_stencil
      integer :: ix = 0, jy = 0, jlo = 0
      real(dp) :: dx = 0d0, y = 0d0
      real(dp) :: ys(4) = 0d0
   end type roche_stencil

contains

   subroutine build_roche_interpolators
      real(dp) :: xtest, ytest, testval
      integer :: ierr
      character(len=strlen) :: upstairs
      include 'formats'
      if (.not. inter_ok) then
         upstairs = trim(mesa_data_dir) // '/roche_data/'  ! where fp/ft data lives
         if (dbg) write(*, 1) 'starting interpolator setup'
         call setup_interpolator(trim(upstairs) // 'fp_data.txt', xvals, num_xpts, yvals, &
               num_ypts, fpfunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'ft_data.txt', xvals, num_xpts, yvals, &
               num_ypts, ftfunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'irot_data.txt', xvals, num_xpts, &
               yvals, num_ypts, irotfunc1d, ierr)
         call setup_interpolator(trim(upstairs) // 'psi_data.txt', xvals, num_xpts, &
               yvals, num_ypts, psifunc1d, ierr)

         if (dbg) then
            xtest = -0.5
            ytest = 1.35
            write(*, 11) 'grid size', num_xpts, num_ypts
            write(*, 1) 'setup interpolators successful,'

            call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
                  fpfunc1d, num_xpts, testval, ierr)
            write(*, 1) 'fp   test gave should be close to 0.6', testval
            call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
                  ftfunc1d, num_xpts, testval, ierr)
            write(*, 1) 'ft   test gave should be close to 0.8', testval
            call interp_evbipm_db(xtest, ytest, xvals, num_xpts, yvals, num_ypts,&
                  irotfunc1d, num_xpts, testval, ierr)
            write(*, 1) 'irot test gave should be close to 1.4', testval
         end if
         inter_ok = .true.
      end if

      contains
      ! interpolator data reading
      subroutine setup_interpolator(filename, xs, num_xs, ys, num_ys, func1d, ierr)
         integer, intent(out) :: ierr, num_xs, num_ys
         character(len = *) :: filename
         integer :: k, iounit
         real(dp), pointer, intent(out) :: xs(:), ys(:), func1d(:)
         real(dp), pointer :: func(:,:,:)

         include 'formats'

         if (dbg) then
            write(*, '(a100)') 'loading ' // filename
         end if
         ! open data to interpolate
         open(newunit = iounit, file = trim(filename), status = 'old', action = 'read',&
               iostat = ierr)

         read(iounit, *, iostat = ierr) num_xs
         allocate(xs(num_xs))
         do k = 1, num_xs
            read(iounit, *, iostat = ierr) xs(k)
         end do

         read(iounit, *, iostat = ierr) num_ys
         allocate(ys(num_ys))
         do k = 1, num_ys
            read(iounit, *, iostat = ierr) ys(k)
         end do

         ! create a 1d array with all the data, point func to it
         allocate(func1d(4 * num_xs * num_ys))
         func(1:4, 1:num_xs, 1:num_ys) => func1d(1:4 * num_xs * num_ys)
         do k = 1, num_xs
            read(iounit, *, iostat = ierr) func(1, k, :)
         end do

         if (ierr /= 0) then
            close(iounit)
         end if
         ! create interpolator
         call interp_mkbipm_db(xs, num_xs, ys, num_ys, func1d, num_xs, ierr)
      end subroutine setup_interpolator

   end subroutine build_roche_interpolators

   real(dp) function eval_fp(lq, ar, ierr) result(fp)
      ! evaluates fp of the equipotential shell with fractional equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr

      include 'formats'

      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            fpfunc1d, num_xpts, fp, ierr)
      if (ierr /= 0) write(*, 1) "error in eval fp", ar, fp
   end function eval_fp

   real(dp) function eval_ft(lq, ar, ierr) result(ft)
      ! evaluates ft of the equipotential shell with fractional equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr

      include 'formats'

      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            ftfunc1d, num_xpts, ft, ierr)
      if (ierr /= 0) write(*, 1) "error in eval ft", ar, ft

   end function eval_ft

   real(dp) function eval_irot(lq, ar, ierr) result(irot)
      ! evaluates moment of inertia irot divided by the spherical equivalent moi (2/3 r_\psi^2)
      ! of the equipotential shell with fractional equivalent radius
      ! ar = r/rl and lq = log10(m(other)/m(this)), no units (dimensionless)
      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr

      include 'formats'

      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            irotfunc1d, num_xpts, irot, ierr)
      if (ierr /= 0) write(*, 1) "error in eval irot", ar, irot

   end function eval_irot

   real(dp) function eval_psi(lq, ar, ierr) result(psi)
      ! evaluates the dimensionless Roche potential of a shell of a given fractional equivalent radius ar = r/rl and
      ! lq = log10(m(other)/m(this)); result has no units. Scaling to dimensionfull units should be done by:
      ! Psi_dim_full = (Psi_dim_less - q^2/(2(1+q)) * G * M_this / separation, with q = m_other / m_this
      ! note there was an error in Fabry et al. (2022) regarding this scaling.

      real(dp), intent(in) :: lq
      real(dp) :: ar
      integer, intent(out) :: ierr

      include 'formats'

      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      call interp_evbipm_db(lq, ar, xvals, num_xpts, yvals, num_ypts,&
            psifunc1d, num_xpts, psi, ierr)
      if (ierr /= 0) write(*, 1) "error in eval psi", ar, psi
   end function eval_psi

   subroutine roche_psi(id, r, psi, ierr)
      integer, intent(in) :: id
      real(dp), intent(in) :: r
      integer, intent(out) :: ierr
      real(dp), intent(out) :: psi

      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: j, this_star, other_star
      real(dp) :: m1, m2, a, r_roche, lq

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr /= 0) return

      m1 = b% m(this_star)
      m2 = b% m(other_star)
      a = b% separation
      r_roche = b% rl(this_star)
      ! if masses or sep are not yet set at the binary level, use these
      if (m1 <= 0) m1 = b% m1 * Msun
      if (m2 <= 0) m2 = b% m2 * Msun
      if (a <= 0) a = pow(standard_cgrav * (m1 + m2) * &
            pow((b% initial_period_in_days) * 86400, 2) / (4 * pi2), one_third)
      if (r_roche <= 0) r_roche = eval_rlobe(m1, m2, a)
      lq = log10(m2 / m1)

      psi = eval_psi(lq, r / r_roche, ierr)
      if (ieee_is_nan(psi)) then
         psi = -9d99  ! let's put negative a lot for regions outside what we can interpolate
                      ! this likely happens in the core where r ~= 0, and so psi == -infty
      end if

   end subroutine roche_psi

   subroutine roche_fp_ft(id, r, fp, ft, r_polar, r_equatorial, report_ierr, ierr)
      use utils_lib, only: mesa_error
      integer, intent(in) :: id
      real(dp), intent(in) :: r
      logical, intent(in) :: report_ierr
      real(dp), intent(inout) :: r_polar, r_equatorial
      type (auto_diff_real_star_order1), intent(out) :: fp, ft
      integer, intent(out) :: ierr

      type (roche_stencil) :: st, st_nudged
      real(dp) :: r_roche, lq, ar, ar_nudged, fp_val, ft_val, fp_nudged, ft_nudged

      include 'formats'

      if (.not. inter_ok) then
         write(*, 1) "interpolators not setup, should happen in startup"
         call mesa_error(__FILE__, __LINE__)
      end if

      call get_roche_params(id, lq, r_roche, ierr)
      if (ierr /= 0) return

      ar = r / r_roche
      call clip_roche_ar(ar, ar_nudged)

      ! one bracket at (lq, ar) and one at (lq, ar + nudge), shared by the fp and ft tables
      call set_roche_stencil(lq, ar, st, ierr)
      if (ierr == 0) call move_roche_stencil(st, ar_nudged, st_nudged, ierr)
      if (ierr == 0) call eval_roche_table(fpfunc1d, st, st_nudged, fp_val, fp_nudged, ierr)
      if (ierr == 0) call eval_roche_table(ftfunc1d, st, st_nudged, ft_val, ft_nudged, ierr)
      if (ierr /= 0) then
         write(*, 1) "error in roche_fp_ft", lq, ar
         ! outside the tables: no tidal correction
         fp = 1d0
         ft = 1d0
      else
         ! set values
         fp = fp_val
         ft = ft_val
         ! set log derivatives
         fp% d1Array(i_lnR_00) = (fp_nudged - fp% val) / nudge * ar
         ft% d1Array(i_lnR_00) = (ft_nudged - ft% val) / nudge * ar
      end if
      ! fix these to the current radius, they're only used for some wind mass loss enhancement
      r_polar = r
      r_equatorial = r
   end subroutine roche_fp_ft

   subroutine roche_irot(id, r00, i_rot)
      integer, intent(in) :: id
      real(dp), intent(in) :: r00
      type (auto_diff_real_star_order1), intent(out) :: i_rot

      type (roche_stencil) :: st, st_nudged
      integer :: ierr
      real(dp) :: r_roche, ar, ar_nudged, lq, eval, eval_nudged, d_eval_d_ar

      include 'formats'

      call get_roche_params(id, lq, r_roche, ierr)
      if (ierr /= 0) return

      i_rot = 0d0

      if (inter_ok) then
         ar = r00 / r_roche
         call clip_roche_ar(ar, ar_nudged)
         ! one bracket at (lq, ar) and one at (lq, ar + nudge)
         call set_roche_stencil(lq, ar, st, ierr)
         if (ierr == 0) call move_roche_stencil(st, ar_nudged, st_nudged, ierr)
         if (ierr == 0) call eval_roche_table(irotfunc1d, st, st_nudged, eval, eval_nudged, ierr)
         if (ierr /= 0) then
            write(*, 1) "error in roche_irot", lq, ar
            ! outside the tables: spherical moment of inertia
            eval = 1d0
            eval_nudged = 1d0
         end if
         ! set value
         d_eval_d_ar = (eval_nudged - eval) / nudge
         ! scale value with spherical moi (see eval_irot)
         i_rot = eval * (two_thirds * r00 * r00)
         ! set radius derivative
         i_rot% d1Array(i_lnR_00) = two_thirds * r00 * r00 * (ar * d_eval_d_ar + 2 * eval)
      end if
   end subroutine roche_irot

   subroutine get_roche_params(id, lq, r_roche, ierr)
      ! the quantities shared by every Roche-table lookup for star id:
      ! lq = log10(m(other)/m(this)) and the Roche lobe radius r_roche
      integer, intent(in) :: id
      real(dp), intent(out) :: lq, r_roche
      integer, intent(out) :: ierr

      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer :: this_star, other_star
      real(dp) :: m1, m2, a

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr /= 0) return
      call which_star(s, b, this_star, other_star, ierr)
      if (ierr /= 0) return

      m1 = b% m(this_star)
      m2 = b% m(other_star)
      a = b% separation
      r_roche = b% rl(this_star)
      ! if masses or sep are not yet set at the binary level, use these
      if (m1 <= 0d0) m1 = b% m1 * Msun
      if (m2 <= 0d0) m2 = b% m2 * Msun
      if (a <= 0d0) a = pow(standard_cgrav * (m1 + m2) * &
            pow((b% initial_period_in_days) * secday, 2) / (4d0 * pi2), one_third)
      if (r_roche <= 0d0) r_roche = eval_rlobe(m1, m2, a)
      lq = log10(m2 / m1)
   end subroutine get_roche_params

   subroutine clip_roche_ar(ar, ar_nudged)
      ! keeps ar, and the point ar + nudge used for the finite-difference derivative,
      ! below the top of the tables. Same clipping as eval_fp, eval_ft and eval_irot.
      real(dp), intent(inout) :: ar
      real(dp), intent(out) :: ar_nudged
      if (ar >= yvals(num_ypts)) ar = yvals(num_ypts) - 2 * nudge
      ar_nudged = ar + nudge
      if (ar_nudged >= yvals(num_ypts)) ar_nudged = yvals(num_ypts) - 2 * nudge
   end subroutine clip_roche_ar

   ! The routines below split interp_evbipm_db (see interp_2d/private/bipm_db.f90) into
   ! a bracket search, a fit in y, and an evaluation of that fit. Tables on the same
   ! (lq, ar) grid can then share one bracket, and the value at ar + nudge can reuse the
   ! fit made for ar when both use the same four y nodes. The arithmetic is the same as
   ! in do_evbipm_db, so the results are identical.

   subroutine set_roche_stencil(xget, yget, st, ierr)
      use num_lib, only: binary_search
      real(dp), intent(in) :: xget, yget
      type (roche_stencil), intent(out) :: st
      integer, intent(out) :: ierr

      type (roche_stencil) :: st_x

      ierr = 0
      st_x% ix = binary_search(num_xpts, xvals, 0, xget)  ! x(ix) <= xget < x(ix+1)
      if (st_x% ix < 1 .or. st_x% ix >= num_xpts) then
         ierr = -1
         return
      end if
      st_x% dx = xget - xvals(st_x% ix)
      st_x% jy = 0
      call move_roche_stencil(st_x, yget, st, ierr)
   end subroutine set_roche_stencil

   subroutine move_roche_stencil(st_in, yget, st, ierr)
      ! stencil at the same x as st_in, at a new y
      use num_lib, only: binary_search
      type (roche_stencil), intent(in) :: st_in
      real(dp), intent(in) :: yget
      type (roche_stencil), intent(out) :: st
      integer, intent(out) :: ierr
      integer :: i, jy

      ierr = 0
      jy = binary_search(num_ypts, yvals, st_in% jy, yget)  ! y(jy) <= yget < y(jy+1)
      if (jy < 1 .or. jy >= num_ypts) then
         ierr = -1
         return
      end if
      st% ix = st_in% ix
      st% dx = st_in% dx
      st% jy = jy
      st% y = yget
      if (jy == 1) then
         st% jlo = 1
      else if (jy >= num_ypts - 1) then
         st% jlo = jy - 2
      else
         st% jlo = jy - 1
      end if
      do i = 1, 4
         st% ys(i) = yvals(st% jlo + i - 1)
      end do
   end subroutine move_roche_stencil

   subroutine fit_roche_table(fin1, st, ff1, ierr)
      ! evaluates the table in x at the four y nodes of st, then fits a
      ! monotonic cubic in y through those four values
      use interp_1d_def, only: pm_work_size
      use interp_1d_lib, only: interp_pm
      real(dp), intent(in), pointer :: fin1(:)  ! (4, num_xpts, num_ypts)
      type (roche_stencil), intent(in) :: st
      real(dp), intent(inout), pointer :: ff1(:)  ! (4, 4)
      integer, intent(out) :: ierr

      integer, parameter :: nwork = pm_work_size
      real(dp), target :: work_ary(4*nwork)
      real(dp), pointer :: work(:), fin(:,:,:), ff(:,:)
      real(dp) :: dx
      integer :: i, j, ix

      work => work_ary
      ff(1:4, 1:4) => ff1(1:4*4)
      fin(1:4, 1:num_xpts, 1:num_ypts) => fin1(1:4*num_xpts*num_ypts)
      ix = st% ix
      dx = st% dx
      do i = 1, 4
         j = st% jlo + i - 1
         ff(1, i) = fin(1, ix, j) + dx*(fin(2, ix, j) + dx*(fin(3, ix, j) + dx*fin(4, ix, j)))
      end do
      call interp_pm(st% ys, 4, ff1, nwork, work, 'fit_roche_table', ierr)
   end subroutine fit_roche_table

   subroutine eval_roche_table(fin1, st, st_nudged, val, val_nudged, ierr)
      ! value of one table at st and at st_nudged, which must have the same x
      use interp_1d_lib, only: interp_values
      real(dp), intent(in), pointer :: fin1(:)
      type (roche_stencil), intent(in) :: st, st_nudged
      real(dp), intent(out) :: val, val_nudged
      integer, intent(out) :: ierr

      real(dp), target :: ff_ary(4*4)
      real(dp), pointer :: ff1(:)
      real(dp) :: ynew(1), v(1)

      ff1 => ff_ary
      call fit_roche_table(fin1, st, ff1, ierr)
      if (ierr /= 0) return
      ynew(1) = st% y
      call interp_values(st% ys, 4, ff1, 1, ynew, v, ierr)
      if (ierr /= 0) return
      val = v(1)

      ! the fit only depends on the four y nodes, so reuse it when they are the same
      if (st_nudged% jlo /= st% jlo) then
         call fit_roche_table(fin1, st_nudged, ff1, ierr)
         if (ierr /= 0) return
      end if
      ynew(1) = st_nudged% y
      call interp_values(st_nudged% ys, 4, ff1, 1, ynew, v, ierr)
      if (ierr /= 0) return
      val_nudged = v(1)
   end subroutine eval_roche_table

   subroutine assign_stars(id, this_star, other_star, ierr)
      ! determine which star is which in the binary
      integer, intent(in) :: id
      integer, intent(out) :: this_star, other_star, ierr
      type (binary_info), pointer :: b
      type (star_info), pointer :: s

      ierr = 0
      call star_ptr(id, s, ierr)
      if (ierr /= 0) return
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr /= 0) return
      call which_star(s, b, this_star, other_star, ierr)
   end subroutine assign_stars

   subroutine which_star(s, b, this_star, other_star, ierr)
      ! same as assign_stars, for callers that already have s and b
      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      integer, intent(out) :: this_star, other_star, ierr

      ierr = 0
      this_star = 0
      other_star = 0
      if (b% s1% id == s% id .and. b% have_star_1) then
         this_star = 1
         other_star = 2
      else if (b% s2% id == s% id .and. b% have_star_2) then
         this_star = 2
         other_star = 1
      else
         ierr = 1
      end if
   end subroutine which_star

   ! determines by what fraction the tidal deformation corrections should be used
   ! eg fp = f_switch * fp_tidal + (1-f_switch) * fp_single
   ! it uses the synchronicity parameter to estimate this. When the shell is quite synchronous,
   ! you probably want to use the tidal corrections, so f_switch -> 1, when not synchronous, the single
   ! star rotation deformation is likely more accurate, so f_switch -> 0 in that case.
   subroutine synchronicity(id, k, omega_in, f_switch, ierr)
      integer, intent(in) :: id, k
      real(dp), intent(in) :: omega_in
      integer, intent(out) :: ierr
      real(dp), intent(out) :: f_switch

      type (star_info), pointer :: s
      type (binary_info), pointer :: b
      real(dp) :: omega, omega_sync, f_sync, p

      include 'formats'

      call star_ptr(id, s, ierr)
      call binary_ptr(s% binary_id, b, ierr)
      if (ierr /= 0) return

      if (b% use_other_tidal_deformation_switch_function) then
         call b% other_tidal_deformation_switch_function(id, k, omega_in, f_switch, ierr)
         return
      end if

      p = b% period
      if (p <= 0d0) p = b% initial_period_in_days * secday
      omega_sync = 2*pi/p
      if (ieee_is_nan(omega_in)) then
         omega = 0d0
      else if (omega_in < 0d0) then
         omega = abs(omega_in)
      else
         omega = omega_in
      end if
      f_sync = omega / omega_sync
      f_sync = min(f_sync, 1d0 / f_sync)  ! we could be super or sub synchronous

      f_switch = 1d0 / (1d0 + exp(-b% f_sync_switch_width * (f_sync - b% f_sync_switch_from_rot_defor)))
      ! apply limiting values if f_switch is far up (or down) the sigmoid
      if (f_switch < b% f_sync_switch_lim) then
         f_switch = 0d0
      else if (1d0 - f_switch < b% f_sync_switch_lim) then
         f_switch = 1d0
      end if

      if (ieee_is_nan(f_switch) .and. k == 1) then
         write(*, 1) "error in synchronicity", f_switch, f_sync, omega, omega_sync, p
         ierr = 1
      end if

   end subroutine synchronicity
end module binary_roche_deformation

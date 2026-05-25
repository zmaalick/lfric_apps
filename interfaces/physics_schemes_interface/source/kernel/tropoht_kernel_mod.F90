!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Kernel to compute tropopause height, temperature, pressure and ICAO height.

module tropoht_kernel_mod

  use argument_mod,         only: arg_type,                  &
                                  GH_FIELD,                  &
                                  GH_READ, GH_WRITE,         &
                                  GH_REAL, GH_INTEGER,       &
                                  CELL_COLUMN,               &
                                  ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,        only: r_def, i_def, r_um
  use fs_continuity_mod,    only: Wtheta
  use kernel_mod,           only: kernel_type

  implicit none

  private

  !> PSyclone kernel metadata
  type, public, extends(kernel_type) :: tropoht_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                         &
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  Wtheta),                      & ! theta
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  Wtheta),                      & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL,    GH_READ,  Wtheta),                      & ! height_wth
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),   & ! trop_level
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! trop_ht
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! trop_temp
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   & ! trop_press
         arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)    & ! trop_icao_ht
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tropoht_code
  end type tropoht_kernel_type

  public :: tropoht_code

contains

  !> @details Computes tropopause height, temperature, pressure and ICAO height
  !> from the tropopause level found by locate_tropopause_kernel.
  !> Uses the WMO lapse-rate intersection method from the UM pws_tropoht_mod.
  !>
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     theta         Potential temperature (K)
  !> @param[in]     exner_in_wth  Exner pressure on Wtheta levels
  !> @param[in]     height_wth    Height above surface on Wtheta levels (m)
  !> @param[in]     trop_level    Tropopause level index from locate_tropopause_kernel
  !> @param[out]    trop_ht       Tropopause height (m)
  !> @param[out]    trop_temp     Tropopause temperature (K)
  !> @param[out]    trop_press    Tropopause pressure (Pa)
  !> @param[out]    trop_icao_ht  Tropopause ICAO height (kft)
  !> @param[in]     ndf_wth       Number of DOFs per cell for Wtheta space
  !> @param[in]     undf_wth      Number of unique DOFs for Wtheta space
  !> @param[in]     map_wth       Dofmap for Wtheta space column base cell
  !> @param[in]     ndf_2d        Number of DOFs per cell for 2D space
  !> @param[in]     undf_2d       Number of unique DOFs for 2D space
  !> @param[in]     map_2d        Dofmap for 2D space column base cell
  subroutine tropoht_code(nlayers,                         &
                          theta,                           &
                          exner_in_wth,                    &
                          height_wth,                      &
                          trop_level,                      &
                          trop_ht,                         &
                          trop_temp,                       &
                          trop_press,                      &
                          trop_icao_ht,                    &
                          ndf_wth, undf_wth, map_wth,      &
                          ndf_2d, undf_2d, map_2d)

    use planet_constants_mod, only: p_zero, kappa, g_over_r

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    real(kind=r_def),    intent(in),    dimension(undf_wth) :: theta
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: exner_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: height_wth
    integer(kind=i_def), intent(in),    dimension(undf_2d)  :: trop_level
    real(kind=r_def),    intent(inout), dimension(undf_2d)  :: trop_ht
    real(kind=r_def),    intent(inout), dimension(undf_2d)  :: trop_temp
    real(kind=r_def),    intent(inout), dimension(undf_2d)  :: trop_press
    real(kind=r_def),    intent(inout), dimension(undf_2d)  :: trop_icao_ht

    integer(kind=i_def) :: k
    real(kind=r_def) :: t_km2, t_km1, t_k, t_kp1
    real(kind=r_def) :: h_km2, h_km1, h_k, h_kp1
    real(kind=r_def) :: lapselwr, lapseupr, delta_lapse
    real(kind=r_def) :: trop_ht_val, trop_temp_val, trop_press_val, p_km1
    real(kind=r_def), parameter :: vsmall = 1.0e-6_r_def

    k = trop_level(map_2d(1))

    ! Guard: locate_tropopause_kernel guarantees [3, nlayers-1] for the
    ! lapse-rate path; cold-point fallback can produce k=1. Both the
    ! [k-2, k-1] and [k, k+1] intervals must be within bounds.
    if (k < 3 .or. k > nlayers - 1) then
      trop_ht(map_2d(1))      = 0.0_r_def
      trop_temp(map_2d(1))    = 0.0_r_def
      trop_press(map_2d(1))   = 0.0_r_def
      trop_icao_ht(map_2d(1)) = 0.0_r_def
      return
    end if

    ! Temperatures at levels k-2, k-1, k, k+1
    t_km2 = theta(map_wth(1)+k-2) * exner_in_wth(map_wth(1)+k-2)
    t_km1 = theta(map_wth(1)+k-1) * exner_in_wth(map_wth(1)+k-1)
    t_k   = theta(map_wth(1)+k)   * exner_in_wth(map_wth(1)+k)
    t_kp1 = theta(map_wth(1)+k+1) * exner_in_wth(map_wth(1)+k+1)

    h_km2 = height_wth(map_wth(1)+k-2)
    h_km1 = height_wth(map_wth(1)+k-1)
    h_k   = height_wth(map_wth(1)+k)
    h_kp1 = height_wth(map_wth(1)+k+1)

    ! locate_tropopause_kernel returns k where the lapse rate in [k-1, k]
    ! first drops below the WMO threshold. The lapse rate below that
    ! transition interval is in [k-2, k-1] (lapselwr, tropospheric), and
    ! the lapse rate above is in [k, k+1] (lapseupr, stratospheric).
    ! This matches the UM pws_tropoht_mod convention with UM level m
    ! corresponding to LFRic level k-1.
    lapselwr = (t_km2 - t_km1) / (h_km1 - h_km2)
    lapseupr = (t_k   - t_kp1) / (h_kp1 - h_k)

    delta_lapse = lapselwr - lapseupr
    if (abs(delta_lapse) < vsmall) then
      if (delta_lapse >= 0.0_r_def) delta_lapse =  vsmall
      if (delta_lapse <  0.0_r_def) delta_lapse = -vsmall
    end if

    ! Exact tropopause height from intersection of lapse-rate lines
    trop_ht_val = ((t_km1 + lapselwr*h_km1) - (t_k + lapseupr*h_k)) / delta_lapse

    ! Clamp to the transition interval
    trop_ht_val = max(trop_ht_val, h_km1)
    trop_ht_val = min(trop_ht_val, h_k)
    trop_ht(map_2d(1)) = trop_ht_val

    ! Temperature at tropopause
    trop_temp_val = t_km1 - lapselwr * (trop_ht_val - h_km1)
    trop_temp(map_2d(1)) = trop_temp_val

    ! Guard lapselwr before use as exponent denominator
    if (abs(lapselwr) < vsmall) then
      if (lapselwr >= 0.0_r_def) lapselwr =  vsmall
      if (lapselwr <  0.0_r_def) lapselwr = -vsmall
    end if

    ! Pressure at tropopause via hypsometric equation, using level k-1 as reference
    p_km1 = p_zero * exner_in_wth(map_wth(1)+k-1)**(1.0_r_def/kappa)
    trop_press_val = p_km1 * (trop_temp_val / t_km1)**(real(g_over_r, r_def) / lapselwr)
    trop_press(map_2d(1)) = trop_press_val

    ! ICAO height (kft) from tropopause pressure
    trop_icao_ht(map_2d(1)) = icao_height_kft(real(trop_press_val, r_um), &
                                               real(g_over_r, r_um))

  end subroutine tropoht_code

  !> @brief Convert pressure (Pa) to ICAO standard atmosphere height (kft).
  !> Ported from UM icao_ht_fc.F90.
  pure function icao_height_kft(pressure_pa, g_over_r_val) result(height_kft)

    implicit none

    real(kind=r_um), intent(in) :: pressure_pa
    real(kind=r_um), intent(in) :: g_over_r_val
    real(kind=r_um)             :: height_kft

    real(kind=r_um), parameter :: lapse_rate_l = 6.5e-3_r_um   ! K/m, troposphere
    real(kind=r_um), parameter :: lapse_rate_u = -1.0e-3_r_um  ! K/m, above 20 km
    real(kind=r_um), parameter :: press_bot    = 101325.0_r_um  ! Pa, surface
    real(kind=r_um), parameter :: press_mid    = 22632.0_r_um   ! Pa, 11000 gpm
    real(kind=r_um), parameter :: press_top    = 5474.87_r_um   ! Pa, 20000 gpm
    real(kind=r_um), parameter :: temp_bot     = 288.15_r_um    ! K, surface
    real(kind=r_um), parameter :: temp_top     = 216.65_r_um    ! K, isothermal layer
    real(kind=r_um), parameter :: gpm1         = 11000.0_r_um   ! m
    real(kind=r_um), parameter :: gpm2         = 20000.0_r_um   ! m
    real(kind=r_um), parameter :: mtokft       = 1.0_r_um / (0.3048_r_um * 1000.0_r_um)

    real(kind=r_um) :: p, zp1, zp2

    zp1 = lapse_rate_l / g_over_r_val
    zp2 = lapse_rate_u / g_over_r_val

    p = max(min(pressure_pa, press_bot), 1000.0_r_um)

    if (p > press_mid) then
      ! Troposphere: up to 11000 gpm
      height_kft = (1.0_r_um - (p/press_bot)**zp1) * temp_bot / lapse_rate_l * mtokft
    else if (p > press_top) then
      ! Isothermal layer: 11000-20000 gpm
      height_kft = (gpm1 + (-log(p/press_mid)) * temp_top / g_over_r_val) * mtokft
    else
      ! Above 20000 gpm
      height_kft = (gpm2 + (1.0_r_um - (p/press_top)**zp2) * temp_top / lapse_rate_u) * mtokft
    end if

  end function icao_height_kft

end module tropoht_kernel_mod

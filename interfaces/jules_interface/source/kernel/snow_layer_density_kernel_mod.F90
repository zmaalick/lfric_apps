!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute the bulk density of each snow layer on tiles.
!> @details The density of a snow layer is the total mass (frozen plus liquid)
!>          divided by its thickness. Layers with zero thickness (inactive
!>          layers) are set to a density of zero.
module snow_layer_density_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_REAL,         &
                           GH_READ, GH_READWRITE,     &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           CELL_COLUMN
  use constants_mod, only: r_def, i_def
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: snow_layer_density_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: snow_layer_density_code
  end type snow_layer_density_kernel_type

  public :: snow_layer_density_code

contains

  !> @param[in]     nlayers               The number of layers
  !> @param[in,out] snow_layer_density    Bulk density of snow layers (kg m-3)
  !> @param[in]     snow_layer_ice_mass   Mass of ice in snow layers (kg m-2)
  !> @param[in]     snow_layer_liq_mass   Mass of liquid in snow layers (kg m-2)
  !> @param[in]     snow_layer_thickness  Thickness of snow layers (m)
  !> @param[in]     ndf_snow              Total DOFs per cell for snow layers
  !> @param[in]     undf_snow             No. of unique DOFs for snow-layer space
  !> @param[in]     map_snow              DOFmap for cells for snow layers
  subroutine snow_layer_density_code(nlayers,                     &
                                     snow_layer_density,          &
                                     snow_layer_ice_mass,         &
                                     snow_layer_liq_mass,         &
                                     snow_layer_thickness,        &
                                     ndf_snow, undf_snow, map_snow)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_snow, undf_snow
    integer(kind=i_def), intent(in) :: map_snow(ndf_snow)

    real(kind=r_def), intent(inout) :: snow_layer_density(undf_snow)
    real(kind=r_def), intent(in)    :: snow_layer_ice_mass(undf_snow)
    real(kind=r_def), intent(in)    :: snow_layer_liq_mass(undf_snow)
    real(kind=r_def), intent(in)    :: snow_layer_thickness(undf_snow)

    integer(kind=i_def) :: df

    do df = 1, ndf_snow
      if (snow_layer_thickness(map_snow(df)) > 0.0_r_def) then
        snow_layer_density(map_snow(df)) =                             &
             (snow_layer_ice_mass(map_snow(df)) +                      &
              snow_layer_liq_mass(map_snow(df))) /                     &
             snow_layer_thickness(map_snow(df))
      else
        snow_layer_density(map_snow(df)) = 0.0_r_def
      end if
    end do

  end subroutine snow_layer_density_code

end module snow_layer_density_kernel_mod

!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Atomic Polarization Tensor calculation by dF/d(E-field) finite differences
!> \author Leo Decking, Hossam Elgabarty
! **************************************************************************************************

MODULE qs_apt_fdiff_types

   USE kinds,                           ONLY: dp
   !USE physcon, ONLY: debye
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gpr_mlff_types'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   PUBLIC :: gpr_mlff_params

   TYPE gpr_mlff_params
      REAL(kind=dp)                                   :: unc = 0.0_dp, &
                                                         unc_threshold = 0.0_dp, &
                                                         radial_cutoff = 3.7_dp, &
                                                         sph_cutoff = 0.0_dp, &
                                                         e_noise = 0.5_dp, & !kj/mol per atom noises!
                                                         f_noise = 0.5_dp    !kj/mol/A
      INTEGER                                         :: n_rad = 12, &
                                                         l_max = 3
      LOGICAL                                         :: do_mlff = .FALSE.
      !REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE     :: gpr_forces, cp_forces
      ! Store internal state of GPR model?
      ! mpi handle? -> driver
   END TYPE gpr_mlff_params

END MODULE qs_apt_fdiff_types

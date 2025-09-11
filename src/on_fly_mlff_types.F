!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief On-the-fly training of machine learning force fields with Gauss. Process Regression (GPR)
!> \author Hossam Elgabarty
! **************************************************************************************************

MODULE on_fly_mlff_types
   USE kinds,                           ONLY: dp, default_path_length
   USE particle_list_types,             ONLY: particle_list_type
   !USE message_passing,                 ONLY: mp_para_env_type

#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'on_fly_mlff_types'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

   !LOGICAL, SAVE                                      :: used_last = .FALSE.


   TYPE prokyon_env_type
      LOGICAL                                            :: do_training, restart
      REAL(kind=dp)                                      :: uncert = 0.0_dp, &
                                                            uncert_threshold = 0.0_dp, &
                                                            update_threshold = 0.0_dp, &
                                                            radial_cutoff = 3.7_dp, &
                                                            e_noise = 5.0_dp, & !kj/mol per atom noise
                                                            f_noise = 10.0_dp   !kj/mol/A per atom noise
      INTEGER                                            :: n_rad = 12, &
                                                            l_max = 3, &
                                                            nskip = 20, &
                                                            natoms
      REAL(kind=dp), DIMENSION(:), POINTER               :: c_forces, coords, uncertainty
      CHARACTER(LEN=default_path_length)                 :: restart_file !, input_file
      !TYPE(force_env_type), POINTER                      :: force_env => NULL()
      TYPE(particle_list_type), POINTER                  :: particles => NULL()
      !TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
      ! Store internal state of GPR model?
      ! mpi handle? -> driver
   END TYPE prokyon_env_type

   PUBLIC :: prokyon_env_type, prokyon_env_release

CONTAINS

!! **************************************************************************************************
!!> \brief Clean up and send exit signal to prokyon
!!> \param ...
!!> \param filename ...
!! **************************************************************************************************
   SUBROUTINE prokyon_env_release(on_fly_mlff_env)
      TYPE(prokyon_env_type), POINTER                      :: on_fly_mlff_env

      CHARACTER(LEN=*), PARAMETER                          :: routineN = 'prokyon_release'


      ! Where to call this function?
      PRINT *, "release MLFF"
      ! send prokyon exit signal
      ! Deallocate force_env%on_fly_mlff_env = NULL()
      IF (ASSOCIATED(on_fly_mlff_env)) THEN
         IF (ASSOCIATED(on_fly_mlff_env%coords)) THEN
            DEALLOCATE (on_fly_mlff_env%coords)
         END IF
         IF (ASSOCIATED(on_fly_mlff_env%c_forces)) THEN
            DEALLOCATE (on_fly_mlff_env%c_forces)
         END IF
         IF (ASSOCIATED(on_fly_mlff_env%uncertainty)) THEN
            DEALLOCATE (on_fly_mlff_env%uncertainty)
         END IF
      END IF

   END SUBROUTINE prokyon_env_release

END MODULE on_fly_mlff_types

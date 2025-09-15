!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to Prokyon, a library for Gaussian Process Regression (GPR)-based on-the-fly
!> \brief Machine Learning Force Fields
!> \brief prokyon github link
!> \author Hossam Elgabarty
! **************************************************************************************************
MODULE on_fly_mlff_methods
   USE on_fly_mlff_types,               ONLY: prokyon_env_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE ISO_C_BINDING,                   ONLY: C_CHAR, C_DOUBLE,&
                                              C_INT, C_NULL_CHAR,&
                                              C_NULL_PTR, C_PTR,&
                                              C_LOC
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE kinds,                           ONLY: dp, int_4
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get, &
                                              section_vals_get_subs_vals
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get, &
                                              qs_subsys_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE particle_types,                  ONLY: particle_type
   USE particle_list_types,             ONLY: particle_list_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get, cp_subsys_type
   USE cell_types,                      ONLY: cell_type

   ! TODO activate once we need to pass the mpi handle to prokyon:
   !USE message_passing,                 ONLY: mp_para_env_type

   ! TODO activate once the restart and input files are implemented:
   !USE kinds,                           ONLY: default_path_length,&
   !                                           default_string_length

   ! We can store atomic GPR uncertainties in the atprop_type, but so far,
   ! we don't need this info inside cp2k as all printing is handled by prokyon
   !USE atprop_types,                    ONLY: atprop_type

   ! TODO prokyon should output traj in pdb format, and write the atomic
   !      uncertainty in the B-factor column


#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'on_fly_mlff_methods'

   PUBLIC :: prokyon_predict, prokyon_train, prokyon_env_init, prokyon_env_release

   INTERFACE

      INTEGER(C_INT) FUNCTION prokyon_train_forces(coords, forces, energy, natoms)  BIND(C, NAME='prokyon_train_forces')
         USE ISO_C_BINDING, ONLY : C_INT, C_DOUBLE
         REAL(C_DOUBLE), DIMENSION(*)                       :: coords
         REAL(C_DOUBLE), DIMENSION(*)                       :: forces
         REAL(C_DOUBLE), VALUE                              :: energy
         INTEGER(C_INT), VALUE                              :: natoms
      END FUNCTION  prokyon_train_forces

      INTEGER(C_INT) FUNCTION prokyon_atoms(labels, natoms) BIND(C, NAME='prokyon_atoms')
         USE ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_PTR
         TYPE(C_PTR), DIMENSION(*)                          :: labels
         INTEGER(C_INT), VALUE                              :: natoms
      END FUNCTION prokyon_atoms

   END INTERFACE



CONTAINS

!! **************************************************************************************************
!!> \brief Initialization
!!> \param prokyon_env
!!> \param qs_env
!!> \param input section
!!> \author Hossam Elgabarty
!! **************************************************************************************************
   SUBROUTINE prokyon_env_init(on_fly_mlff_env, qs_env, force_section)
      TYPE(prokyon_env_type), POINTER                        :: on_fly_mlff_env
      TYPE(section_vals_type), POINTER                       :: force_section, mlff_section
      INTEGER                                                :: natoms, iatom, i, ret
      CHARACTER(LEN=6,KIND=C_CHAR), ALLOCATABLE, TARGET      :: names(:)
      TYPE(C_PTR), DIMENSION(:), ALLOCATABLE                 :: c_names
      TYPE(qs_subsys_type), POINTER                          :: subsys
      TYPE(cell_type), POINTER                               :: cell
      REAL(dp)                                               :: c_cell(10), conv
      CHARACTER(LEN=2)                                       :: element_symbol
      TYPE(particle_list_type), POINTER                      :: particles
      TYPE(qs_environment_type), INTENT(IN), POINTER         :: qs_env
      ! TODO send mpi handle to prokyon
      !TYPE(mp_para_env_type), POINTER                       :: para_env => NULL()


      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prokyon_init'


      PRINT *, "start init MLFF"

      NULLIFY (on_fly_mlff_env, particles)

      ALLOCATE(on_fly_mlff_env)
      !TODO nullify and assign mpi handle
      !NULLIFY (para_env)

      CALL get_qs_env(qs_env=qs_env, subsys=subsys, cell=cell, natom=natoms)
      CALL qs_subsys_get(subsys, particles=particles)
      on_fly_mlff_env%particles => particles
      on_fly_mlff_env%natoms = natoms

      ALLOCATE (on_fly_mlff_env%coords(3 * natoms))
      ALLOCATE (on_fly_mlff_env%c_forces(3 * natoms))
      ALLOCATE (on_fly_mlff_env%uncertainty(3 * natoms))

      mlff_section => section_vals_get_subs_vals(force_section, "OnFly")

      CALL section_vals_val_get(mlff_section, "E_NOISE", r_val=on_fly_mlff_env%e_noise)
      CALL section_vals_val_get(mlff_section, "F_NOISE", r_val=on_fly_mlff_env%f_noise)
      CALL section_vals_val_get(mlff_section, "TRAIN", l_val=on_fly_mlff_env%do_training)
      CALL section_vals_val_get(mlff_section, "CUTOFF_R", r_val=on_fly_mlff_env%radial_cutoff)
      CALL section_vals_val_get(mlff_section, "PREDICTION_THRESHOLD", r_val=on_fly_mlff_env%uncert_threshold)
      CALL section_vals_val_get(mlff_section, "UPDATE_THRESHOLD", r_val=on_fly_mlff_env%update_threshold)
      CALL section_vals_val_get(mlff_section, "N_RAD", i_val=on_fly_mlff_env%n_rad)
      CALL section_vals_val_get(mlff_section, "L_MAX", i_val=on_fly_mlff_env%l_max)
      CALL section_vals_val_get(mlff_section, "N_SKIP", i_val=on_fly_mlff_env%nskip)

      !TODO implement restart files
      !CALL section_vals_val_get(mlff_section, "RESTART_FILE", c_val=on_fly_mlff_env%restart_file)

      !TODO implement input files
      !CALL section_vals_val_get(mlff_section, "INPUT_FILE", c_val=on_fly_mlff_env%input_file)

      !Prepare list of atoms
      !TODO: for QMMM, one should use the particle list from qmmm_env%qs_env
      !NOTE: prokyon should strip any training empty spaces in elements with single-character names
      ALLOCATE (names(natoms))
      ALLOCATE (c_names(natoms))
      conv = cp_unit_from_cp2k(1.0_dp, "angstrom")
      DO iatom = 1, natoms
         CALL get_atomic_kind(atomic_kind=on_fly_mlff_env%particles%els(iatom)%atomic_kind, &
                              element_symbol=element_symbol)
         names(iatom) = element_symbol//C_NULL_CHAR
         c_names(iatom) = C_LOC(names(iatom))
      END DO

      c_cell = 0.0_dp
      DO i = 1, 3
         IF (cell%perd(i) == 1) THEN
             c_cell(10) = c_cell(10) + 1.0_dp
         END IF
         c_cell((i-1)*3+1:(i-1)*3+3)  = cell%hmat(1:3, i)*conv
      END DO

      ! TODO initialize flare
      ! CALL prokyon_init(n_rad, l_max, radial_cutoff, natoms, c_names, c_cell, e_noise, f_noise, &
      !                   uncert_threshold, update_threshold, n_skip, do_training, &
      !                   restart, restart_file, mpi_handle)

      PRINT *, on_fly_mlff_env%natoms
      PRINT *, particles%n_els
      PRINT *, SIZE(particles%els)
      ret = prokyon_atoms(c_names, on_fly_mlff_env%natoms)
      DEALLOCATE (c_names)
      DEALLOCATE (names)
      PRINT *, "end init MLFF"

   END SUBROUTINE prokyon_env_init


!! **************************************************************************************************
!!> \brief Ask prokyon to make a prediction
!!> \brief Prokyon can always refuse by sending use_mlff_force = .FALSE.
!!> \param prokyon_env
!!> \param mlff_energy: predicted energy
!!> \param use_mlff_force: logical flag
!!> \author Hossam Elgabarty
!! **************************************************************************************************
   SUBROUTINE prokyon_predict(on_fly_mlff_env, mlff_energy, use_mlff_force)
      TYPE(prokyon_env_type), POINTER                 :: on_fly_mlff_env
      REAL(KIND=dp)                                       :: mlff_energy, rconv
      LOGICAL, INTENT(OUT)                                :: use_mlff_force
      TYPE(particle_list_type), POINTER                   :: particles
      INTEGER                                             :: iparticle

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prokyon_predict'

      NULLIFY(particles)
      particles => on_fly_mlff_env%particles

      !TODO merge train and predict subroutines into one subroutine?

      rconv = cp_unit_from_cp2k(1.0_dp, "angstrom")
      DO iparticle = 1, particles%n_els
         on_fly_mlff_env%coords((iparticle-1)*3+1:(iparticle-1)*3+3) = &
                                                    particles%els(iparticle)%r(1:3)*rconv
      END DO

      !TODO make sure returned energy has cp2k units
      mlff_energy = 0.0_dp
      use_mlff_force = .FALSE.

      !TODO call prokyon predict function
      !CALL prokyon_predict_forces(on_fly_mlff_env%coords, on_fly_mlff_env%c_forces, &
      !                            on_fly_mlff_env%uncertainty, mlff_energy, use_mlff_force)

      !TODO copy c_forces to on_fly_mlff_env%particles%els(i)%f(1:3)

      PRINT *, "predict MLFF"
   END SUBROUTINE prokyon_predict

!! **************************************************************************************************
!!> \brief Send training set to prokyon
!!> \param prokyon_env
!!> \param energy
!!> \author Hossam Elgabarty
!! **************************************************************************************************
   SUBROUTINE prokyon_train(on_fly_mlff_env, energy)
      TYPE(prokyon_env_type), POINTER                     :: on_fly_mlff_env
      REAL(KIND=dp)                                       :: energy !, fconv
      REAL(KIND=dp)                                       :: rconv !, fconv
      INTEGER                                             :: iparticle, ret
      TYPE(particle_list_type), POINTER                   :: particles

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prokyon_send_training_set'

      PRINT *, "Train MLFF"


      NULLIFY(particles)
      particles => on_fly_mlff_env%particles

      !fconv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_string))
      rconv = cp_unit_from_cp2k(1.0_dp, "angstrom")

      DO iparticle = 1, particles%n_els
         on_fly_mlff_env%coords((iparticle-1)*3+1:(iparticle-1)*3+3) = &
                                           particles%els(iparticle)%r(1:3)*rconv
         on_fly_mlff_env%c_forces((iparticle-1)*3+1:(iparticle-1)*3+3) = &
                                           particles%els(iparticle)%f(1:3)!*fconv
      END DO

      iparticle = particles%n_els
      ret = prokyon_train_forces(on_fly_mlff_env%coords, on_fly_mlff_env%c_forces, energy, iparticle)
      !!IF (ret /= 0) THEN
      !!   CPABORT("libprokyon library returned an error. Aborting.")

      !TODO return forces or particles_list in cp2k units
      particles%els(1)%r(1) = on_fly_mlff_env%coords(1) ! *fconv

      PRINT *, "End train MLFF"

   END SUBROUTINE prokyon_train

!! **************************************************************************************************
!!> \brief Clean up and send exit signal to prokyon
!!> \param prokyon_env
!!> \author Hossam Elgabarty
!! **************************************************************************************************
   SUBROUTINE prokyon_env_release(on_fly_mlff_env)
      TYPE(prokyon_env_type), POINTER                      :: on_fly_mlff_env

      CHARACTER(LEN=*), PARAMETER                          :: routineN = 'prokyon_release'


      PRINT *, "release MLFF"
      ! TODO send prokyon exit signal

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

!! **************************************************************************************************
!!> \brief Flatten a 3*natoms array to 1D
!!> \param array of particles
!!> \param 1D array of doubles
!!> \author Hossam Elgabarty
!! **************************************************************************************************
!~    SUBROUTINE flatten_array(particles, coords)
!~       TYPE(particle_list_type), POINTER                   :: particles
!~       INTEGER                                             :: iparticle
!~       REAL(kind=dp), DIMENSION(:), INTENT(IN), POINTER    :: coords


!~       DO iparticle = 1, on_fly_mlff_env%natoms
!~          on_fly_mlff_env%coords((iparticle-1)*3+1:(iparticle-1)*3+3) = &
!~                                            on_fly_mlff_env%particles%els(iparticle)%r(1:3)*rconv
!~          on_fly_mlff_env%c_forces((iparticle-1)*3+1:(iparticle-1)*3+3) = &
!~                                            on_fly_mlff_env%particles%els(iparticle)%f(1:3)!*fconv
!~       END DO
!~    END SUBROUTINE flatten_array

END MODULE on_fly_mlff_methods

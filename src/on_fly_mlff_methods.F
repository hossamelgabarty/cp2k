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
   !USE force_env_types,                 ONLY: force_env_get, force_env_type
   ! in principle, we can store atomic GPR uncertainties here:
   !USE atprop_types,                    ONLY: atprop_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE kinds,                           ONLY: dp
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
   !USE message_passing,                 ONLY: mp_para_env_type

   !USE kinds,                           ONLY: default_path_length,&
   !                                           default_string_length
   !USE machine,                         ONLY: m_flush

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'on_fly_mlff_methods'

   PUBLIC :: prokyon_predict, prokyon_train, prokyon_env_init, prokyon_env_release

CONTAINS

! **************************************************************************************************
!> \brief Initialization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE prokyon_env_init(qs_env, force_section)
      !TYPE(force_env_type), INTENT(IN), POINTER              :: force_env
      TYPE(prokyon_env_type), POINTER                        :: on_fly_mlff_env
      TYPE(section_vals_type), POINTER                       :: force_section, mlff_section
      INTEGER                                                :: natoms, iatom, i
      CHARACTER(LEN=6,KIND=C_CHAR), ALLOCATABLE, TARGET      :: names(:)
      TYPE(C_PTR), DIMENSION(:), ALLOCATABLE                 :: c_names
      TYPE(qs_subsys_type), POINTER                          :: subsys
      !TYPE(atprop_type), POINTER                            :: atprop_env
      TYPE(cell_type), POINTER                               :: cell
      REAL(dp)                                               :: c_cell(10), conv
      CHARACTER(LEN=2)                                       :: element_symbol
      !TODO I probably don't need the particle_list_type at all
      TYPE(particle_list_type), POINTER                      :: particles
      TYPE(particle_type), DIMENSION(:), POINTER             :: particle_set
      TYPE(qs_environment_type), INTENT(IN), POINTER         :: qs_env
      !TYPE(mp_para_env_type), POINTER                       :: para_env => NULL()

      !REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)              :: atoms_coord



      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prokyon_init'


      PRINT *, "start init MLFF"

      NULLIFY (on_fly_mlff_env, particles)
      !NULLIFY (para_env)

      ALLOCATE(on_fly_mlff_env)
      !TODO assign mpi handle

      CALL get_qs_env(qs_env=qs_env, subsys=subsys, particle_set=particle_set, cell=cell, natom=natoms)
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
      CALL section_vals_val_get(mlff_section, "RESTART_FILE", c_val=on_fly_mlff_env%restart_file)

      !TODO implement input files
      !CALL section_vals_val_get(mlff_section, "INPUT_FILE", c_val=on_fly_mlff_env%input_file)

      ! establish connection with flare

      !Prepare list of atoms
      !TODO: for QMMM, one should use the particle list from qmmm_env%qs_env
      !NOTE: prokyon should strip any training empty spaces elements with single-character names

      ALLOCATE (names(natoms))
      ALLOCATE (c_names(natoms))
      !ALLOCATE (atoms_coord(3 * natoms))
      conv = cp_unit_from_cp2k(1.0_dp, "angstrom")
      !I can also use on_fly_mlff_env%particles%els(iparticle)%atomic_kind
      DO iatom = 1, natoms
         CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                              element_symbol=element_symbol)
         names(iatom) = element_symbol//C_NULL_CHAR
         c_names(iatom) = C_LOC(names(iatom))
      END DO

      ! Alternative, one long string, but not sure about C_PTR => allocatable character(LEN=:)
      !CHARACTER(LEN=:), ALLOCATABLE                      :: atom_list
      !TYPE(C_PTR)                                        :: atoms_string
      !ALLOCATE (CHARACTER(LEN=natoms*2+1) :: atom_list)
      !DO index = 1, natoms, 2
      !   atom_list(index:index+1) = element_symbol(1:2)
      !END DO
      !atom_list(natoms*2+1) = C_NULL_CHAR
      !atoms_string = C_LOC(atom_list) ! I am not sure whether this seamlessly works

      c_cell = 0.0_dp
      DO i = 1, 3
         IF (cell%perd(i) == 1) THEN
             c_cell(10) = c_cell(10) + 1.0_dp
         END IF
         c_cell((i-1)*3+1:(i-1)*3+3)  = cell%hmat(1:3, i)*conv
      END DO

      ! CALL prokyon_init(n_rad, l_max, radial_cutoff, natoms, c_names, c_cell, e_noise, f_noise, &
      !                   uncert_threshold, update_threshold, n_skip, do_training, &
      !                   restart, restart_file, mpi_handle)

      DEALLOCATE (c_names)
      DEALLOCATE (names)

      PRINT *, "end init MLFF"

   END SUBROUTINE prokyon_env_init


!! **************************************************************************************************
!!> \brief Ask prokyon to make a prediction
!!> \brief Prokyon can always refuse by sending use_mlff_force = .FALSE.
!!> \param particles
!!> \param forces
!! **************************************************************************************************
   !SUBROUTINE prokyon_predict(force_env, mlff_energy, use_mlff_force)
   SUBROUTINE prokyon_predict(on_fly_mlff_env, mlff_energy, use_mlff_force)
      TYPE(prokyon_env_type), POINTER                 :: on_fly_mlff_env
      !TYPE(force_env_type), INTENT(IN), POINTER          :: force_env
      !TYPE(cp_subsys_type), POINTER                      :: subsys
      REAL(KIND=dp)                                       :: mlff_energy, rconv
      LOGICAL, INTENT(OUT)                                :: use_mlff_force
      INTEGER                                             :: iparticle
      !REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c_forces, uncertainty
      !TYPE(particle_list_type), INTENT(INOUT), POINTER   :: particles
      !REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coords

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prokyon_predict'


      !CALL force_env_get(force_env, subsys=subsys)
      !CALL cp_subsys_get(subsys, particles=particles)

      !ALLOCATE (coords(3 * particles%n_els))
      !ALLOCATE (c_forces(3 * particles%n_els))

      !ALLOCATE (coords(3 * on_fly_mlff_env%natoms))
      !ALLOCATE (c_forces(3 * on_fly_mlff_env%natoms))

      rconv = cp_unit_from_cp2k(1.0_dp, "angstrom")

      !DO iparticle = 1, particles%n_els
      DO iparticle = 1, on_fly_mlff_env%natoms
         on_fly_mlff_env%coords((iparticle-1)*3+1:(iparticle-1)*3+3) = &
                                       on_fly_mlff_env%particles%els(iparticle)%r(1:3)*rconv
      END DO

      !TODO make sure returned energy has cp2k units
      mlff_energy = 0.0_dp
      use_mlff_force = .FALSE.
      !CALL prokyon_predict_forces(on_fly_mlff_env%coords, on_fly_mlff_env%c_forces, &
      !                            on_fly_mlff_env%uncertainty, mlff_energy, use_mlff_force)

      !TODO copy c_forces to on_fly_mlff_env%particles%els(iparticle)%f(1:3)

      !DEALLOCATE (coords)
      !DEALLOCATE (c_forces)

      PRINT *, "predict MLFF"
   END SUBROUTINE prokyon_predict

!! **************************************************************************************************
!!> \brief Send training set to prokyon
!!> \param particles
!!> \param forces
!!> \param energy
!!> \param filename ...
!! **************************************************************************************************
   !SUBROUTINE prokyon_train(force_env, energy)
   SUBROUTINE prokyon_train(on_fly_mlff_env, energy)
      TYPE(prokyon_env_type), POINTER                 :: on_fly_mlff_env
      !TYPE(force_env_type), POINTER                      :: force_env
      REAL(KIND=dp), INTENT(OUT)                           :: energy !, fconv
      REAL(KIND=dp)                                       :: rconv !, fconv
      INTEGER                                             :: iparticle
      !REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: forces
      !TYPE(particle_list_type), POINTER                  :: particles
      !REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: coords

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'prokyon_send_training_set'


      !CALL force_env_get(force_env, subsys=subsys)
      !CALL cp_subsys_get(subsys, particles=particles)

      !ALLOCATE (coords(3 * particles%n_els))
      !ALLOCATE (forces(3 * particles%n_els))
      !fconv = cp_unit_from_cp2k(1.0_dp, TRIM(unit_string))
      rconv = cp_unit_from_cp2k(1.0_dp, "angstrom")

      !TODO create subroutines update_coords and update_forces
      DO iparticle = 1, on_fly_mlff_env%natoms
         on_fly_mlff_env%coords((iparticle-1)*3+1:(iparticle-1)*3+3) = &
                                           on_fly_mlff_env%particles%els(iparticle)%r(1:3)*rconv
         on_fly_mlff_env%c_forces((iparticle-1)*3+1:(iparticle-1)*3+3) = &
                                           on_fly_mlff_env%particles%els(iparticle)%f(1:3)!*fconv
      END DO

      !CALL prokyon_train_forces(coords, forces, energy)
      energy = 0.0_dp
      !CALL prokyon_train_forces(on_fly_mlff_env%coords, on_fly_mlff_env%c_forces, energy)

      !DEALLOCATE (coords)
      !DEALLOCATE (forces)
      PRINT *, "send forces MLFF"

   END SUBROUTINE prokyon_train

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

END MODULE on_fly_mlff_methods

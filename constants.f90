module constants

    implicit none

    ! Units used:
    ! distance → Angstrom
    ! energy   → eV
    ! time     → fs
    ! mass     → eV·fs²/Angstrom²


    integer, parameter              ::  n_atoms_max = 1000, n_steps = 1000, n_neigh_max = 500
    double precision, parameter     ::  convert_mass = 1.0364269d2 ! 1 u = 1.0364269d2 eV.fs²/Angstrom²
    integer, parameter              ::  n_types = 5
    double precision, parameter     ::  mass(n_types) = (/58.6934d0, 63.546d0, 102.9055d0, 106.42d0, 195.084d0/) * convert_mass

    double precision, parameter     ::  kb = 8.617d-5 ! Boltzmann constant (eV/K)
    double precision, parameter     ::  omega = 0.1d0 ! Nose-Hoover frequency (1/fs)
    double precision, parameter     ::  dt = 1.0d0 ! Time step (fs)
    character(len=2), parameter     ::  iel_to_typ(n_types) = (/'Ni', 'Cu', 'Rh', 'Pd', 'Pt' /)

    double precision                ::  tbsma_a(n_types, n_types), tbsma_xi(n_types, n_types), &
                                        tbsma_p(n_types, n_types), tbsma_q(n_types, n_types), &
                                        tbsma_r0(n_types, n_types), tbsma_rc1(n_types, n_types), &
                                        tbsma_rc2(n_types, n_types), &
                                        tbsma_x5(n_types, n_types), tbsma_x4(n_types, n_types), tbsma_x3(n_types, n_types), &
                                        tbsma_a5(n_types, n_types), tbsma_a4(n_types, n_types), tbsma_a3(n_types, n_types)
    logical, parameter              ::  use_thermostat = .true.
    double precision                ::  skin = 0.3d0
endmodule


module variables
    use constants, only: n_atoms_max, n_neigh_max
    implicit none

    integer                         ::  typ(n_atoms_max), n_atoms
    double precision                ::  box(3), pos(3, n_atoms_max), vel(3, n_atoms_max), epot(n_atoms_max), force(3, n_atoms_max)
    double precision                ::  pos_old(3, n_atoms_max)
    double precision                ::  target_temperature, time, Qth
    integer                         ::  neigh(n_neigh_max, n_atoms_max), n_neigh(n_atoms_max)
endmodule
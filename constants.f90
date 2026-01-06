module constants

    implicit none

    ! Units used:
    ! distance → Angstrom
    ! energy   → eV
    ! time     → fs
    ! mass     → eV·fs²/Angstrom²


    integer, parameter              ::  n_atoms_max = 1000, n_neigh_max = 500
    double precision, parameter     ::  convert_mass = 1.0364269d2 ! 1 u = 1.0364269d2 eV.fs²/Angstrom²
    integer, parameter              ::  n_types = 2
    double precision, parameter     ::  mass(n_types) = (/107.8682d0, 58.933195d0/) * convert_mass

    double precision, parameter     ::  kb = 8.617d-5 ! Boltzmann constant (eV/K)
    double precision, parameter     ::  omega = 0.1d0 ! Nose-Hoover frequency (1/fs)
    double precision, parameter     ::  dt = 1.0d0 ! Time step (fs)
    character(len=2), parameter     ::  iel_to_typ(n_types) = (/'Ag', 'Co'/)

    double precision                ::  tbsma_a(n_types, n_types), tbsma_xi(n_types, n_types), &
                                        tbsma_p(n_types, n_types), tbsma_q(n_types, n_types), &
                                        tbsma_r0(n_types, n_types), tbsma_rc1(n_types, n_types), &
                                        tbsma_rc2(n_types, n_types), tbsma_rc2_sq(n_types, n_types), tbsma_rc2_sq_max, &
                                        tbsma_x5(n_types, n_types), tbsma_x4(n_types, n_types), tbsma_x3(n_types, n_types), &
                                        tbsma_a5(n_types, n_types), tbsma_a4(n_types, n_types), tbsma_a3(n_types, n_types)
    logical, parameter              ::  use_thermostat = .true.
    double precision, parameter     ::  skin = 1.0d0
    double precision, parameter     ::  pi = 4.0d0 * atan(1.0d0)

endmodule


module variables
    use constants, only: n_atoms_max, n_neigh_max
    implicit none

    integer                         ::  typ(n_atoms_max), n_atoms, n_steps
    double precision                ::  box(3), pos(3, n_atoms_max), vel(3, n_atoms_max), epot(n_atoms_max), force(3, n_atoms_max)
    double precision                ::  target_temperature, initial_temperature, final_temperature 
    double precision                ::  time, Qth, ekin, eth
    integer                         ::  neigh(n_neigh_max, n_atoms_max), n_neigh(n_atoms_max)

endmodule
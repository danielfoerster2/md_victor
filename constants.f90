module constants

    implicit none

    ! Units used:
    ! distance → Angstrom
    ! energy   → eV
    ! time     → fs
    ! mass     → eV·fs²/Angstrom²


    integer, parameter              ::  N = 5 ! Stop to the Nth neighbors
    integer, parameter              ::  n_cells = 5 ! number of cells in lattice

    double precision, parameter     ::  convert_mass = 1.0364269d2 ! 1 u = 1.0364269d2 eV.fs²/Angstrom²
    double precision, parameter     ::  kb = 8.617d-5 ! Boltzmann constant (eV/K)
    !double precision, parameter     ::  box_dimension = 10.0d0 * n_cells
    !double precision, parameter     ::  box(3) = (/ 1.0d0, 1.0d0, 1.0d0 /) * box_dimension

    character(len=2), parameter     ::  atoms(16) = (/ "Ni", "Cu", "Rh", "Pd", "Ag", "Ir", "Pt", "Au", &
                                            "Al", "Pb", "Ti", "Zr", "Co", "Cd", "Zn", "Mg" /)

end module
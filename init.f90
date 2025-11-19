module init
    
    use constants
    use file_io, only : read_xyz
    implicit none
    contains


    function random_list_atom(n_atoms, atoms, percent)
        ! ========================
        ! Determine the atoms composing the system
        ! 
        ! Parameters :
        ! ------------
        ! n_atoms : integer
        !           Number of atoms
        ! atoms : 1D character array
        !           List of symbol of all atoms recorded in constants
        ! percent : 1D double precision array
        !           Percent of each atoms in the system
        !
        ! Returns :
        ! ---------
        ! 1D character array
        !           List of atoms symbols composing the system
        ! ========================

        implicit none
        
        integer, intent(in)                 ::  n_atoms
        character(len=2), intent(in)        ::  atoms(:)
        double precision, intent(in)        ::  percent(:)
        character(len=2)                    ::  random_list_atom(n_atoms)
        integer                             ::  i, j, cpt(size(atoms)), k = 1
        double precision                    ::  x = 0.0d0

        do i = 1, size(percent)
            cpt(i) = 0
            x = x + n_atoms * percent(i) / 100.0d0
            do j = k, int(x) 
                random_list_atom(j) = atoms(i)
                cpt(i) = cpt(i) + 1
            end do
            k = int(x) + 1
            print*, atoms(i), ":", cpt(i)
        end do
    
    end function random_list_atom


    function var_array(list_atom)
        ! ========================
        ! Fill an array with the caracteristics of all atoms in the system
        ! 
        ! Parameters :
        ! ------------
        ! list_atom : 1D character array
        !           List of atoms composing the system
        !
        ! Returns :
        ! ---------
        ! 2D double precision array
        !           Array of atoms caracteristics
        !           Dim 1 : Number of atoms
        !           Dim 2 : a0, A, xi, p, q, mass , r0
        ! ======================== 

        implicit none

        character(len=2), intent(in)        ::  list_atom(:)
        double precision                    ::  var_array(size(list_atom), 7)
        integer                             ::  i
        
        do i = 1, size(list_atom)
            select case(list_atom(i))
                case ("Ni")
                    var_array(i,1) = 3.523d0 ! a0 
                    var_array(i,2) = 0.0376d0 ! A
                    var_array(i,3) = 1.070d0 ! xi
                    var_array(i,4) = 16.999d0 ! p
                    var_array(i,5) = 1.189d0 ! q
                    var_array(i,6) = 58.6934d0 * convert_mass ! mass

                case ("Cu")
                    var_array(i,1) = 3.615d0
                    var_array(i,2) = 0.0855d0
                    var_array(i,3) = 1.224d0
                    var_array(i,4) = 10.960d0
                    var_array(i,5) = 2.278d0
                    var_array(i,6) = 63.546d0 * convert_mass
                case ("Rh")
                    var_array(i,1) = 3.803d0
                    var_array(i,2) = 0.0629d0
                    var_array(i,3) = 1.660d0
                    var_array(i,4) = 18.450d0
                    var_array(i,5) = 1.867d0
                    var_array(i,6) = 102.9055d0 * convert_mass
                case ("Pd")
                    var_array(i,1) = 3.936d0
                    var_array(i,2) = 0.1746d0
                    var_array(i,3) = 1.718d0
                    var_array(i,4) = 10.867d0
                    var_array(i,5) = 3.742d0
                    var_array(i,6) = 106.42d0 * convert_mass
                case ("Ag")
                    var_array(i,1) = 4.085d0
                    var_array(i,2) = 0.1028d0
                    var_array(i,3) = 1.178d0
                    var_array(i,4) = 10.928d0
                    var_array(i,5) = 3.139d0
                    var_array(i,6) = 107.8682d0 * convert_mass
                case ("Ir")
                    var_array(i,1) = 3.839d0
                    var_array(i,2) = 0.1156d0
                    var_array(i,3) = 2.289d0
                    var_array(i,4) = 16.980d0
                    var_array(i,5) = 2.691d0
                    var_array(i,6) = 192.217d0 * convert_mass
                case ("Pt")
                    var_array(i,1) = 3.924d0
                    var_array(i,2) = 0.2975d0
                    var_array(i,3) = 2.695d0
                    var_array(i,4) = 10.612d0
                    var_array(i,5) = 4.004d0
                    var_array(i,6) = 195.084d0 * convert_mass
                case ("Au")
                    var_array(i,1) = 4.079d0
                    var_array(i,2) = 0.2061d0
                    var_array(i,3) = 1.790d0
                    var_array(i,4) = 10.229d0
                    var_array(i,5) = 4.036d0
                    var_array(i,6) = 196.96657d0 * convert_mass
                case ("Al")
                    var_array(i,1) = 4.050d0
                    var_array(i,2) = 0.1221d0
                    var_array(i,3) = 1.316d0
                    var_array(i,4) = 8.612d0
                    var_array(i,5) = 2.516d0
                    var_array(i,6) = 26.981539d0 * convert_mass
                case ("Pb")
                    var_array(i,1) = 4.951d0
                    var_array(i,2) = 0.0980d0
                    var_array(i,3) = 0.914d0
                    var_array(i,4) = 9.576d0
                    var_array(i,5) = 3.648d0
                    var_array(i,6) = 207.2d0 * convert_mass
                case ("Ti")
                    var_array(i,1) = 2.492d0
                    var_array(i,2) = 0.1519d0
                    var_array(i,3) = 1.8112d0
                    var_array(i,4) = 8.620d0
                    var_array(i,5) = 2.390d0
                    var_array(i,6) = 47.867d0 * convert_mass
                case ("Zr")
                    var_array(i,1) = 3.232d0
                    var_array(i,2) = 0.1934d0
                    var_array(i,3) = 2.2792d0
                    var_array(i,4) = 8.250d0
                    var_array(i,5) = 2.249d0
                    var_array(i,6) = 91.224d0 * convert_mass
                case ("Co")
                    var_array(i,1) = 2.507d0
                    var_array(i,2) = 0.0950d0
                    var_array(i,3) = 1.4880d0
                    var_array(i,4) = 11.604d0
                    var_array(i,5) = 2.286d0
                    var_array(i,6) = 58.933195d0 * convert_mass
                case ("Cd")
                    var_array(i,1) = 2.959d0
                    var_array(i,2) = 0.1420d0
                    var_array(i,3) = 0.8117d0
                    var_array(i,4) = 10.612d0
                    var_array(i,5) = 5.206d0
                    var_array(i,6) = 112.414d0 * convert_mass
                case ("Zn")
                    var_array(i,1) = 2.653d0
                    var_array(i,2) = 0.1477d0
                    var_array(i,3) = 0.8900d0
                    var_array(i,4) = 9.689d0
                    var_array(i,5) = 4.602d0
                    var_array(i,6) = 65.38d0 * convert_mass
                case ("Mg")
                    var_array(i,1) = 3.176d0
                    var_array(i,2) = 0.0290d0
                    var_array(i,3) = 0.4992d0
                    var_array(i,4) = 12.820d0
                    var_array(i,5) = 2.257d0
                    var_array(i,6) = 24.305d0 * convert_mass
                
                end select
                var_array(i,7) = var_array(i,1) / sqrt(2.0d0) ! r0

        end do

    end function var_array


    function generate_lattice(var, basis, n_atoms)
        !========================================
        ! Generate a FCC lattice of atoms from a basis FCC structure
        ! Only work for a system commposed of one type of particle for the moment
        !
        ! Parameters :
        ! ------------
        ! var : 2D array
        !       Properties of atoms
        !       Dim 1 : Atom
        !       Dim 2 : a0, A, xi, p, q, mass, r0
        ! basis : 2D array
        !       Coordinates of an unit cell
        ! n_atoms : integer
        !       Number of atoms
        !
        ! Returns :
        ! ---------
        ! 2D array
        !        x, y, z coordinates of each atoms in the lattice
        !        Dim 1 : Number of atoms
        !        Dim 2 : x, y, z coordinates
        !========================================
        implicit none

        double precision, intent(in)        ::  basis(:,:), var(:,:)
        integer, intent(in)                 ::  n_atoms
        double precision                    ::  generate_lattice(n_atoms, 3), r(3)
        integer                             ::  nx, ny, nz, i, cpt = 0

        do nx = -n_cells/2, n_cells/2
            do ny = -n_cells/2, n_cells/2
                do nz = -n_cells/2, n_cells/2
                    do i = 1, 4
                        cpt = cpt + 1
                        r = (/ nx*var(cpt,1) + basis(i,1), ny*var(cpt,1) + basis(i,2), nz*var(cpt,1) + basis(i,3) /)
                        generate_lattice(cpt,:) = r(:)
                    end do
                end do
            end do
        end do

    end function generate_lattice


    function generate_random(var, n_atoms, box)
        !========================================
        ! Generate a random lattice of atoms
        ! 
        ! Parameters :
        ! ------------
        ! var : 2D array
        !       Properties of atoms
        !       Dim 1 : Atom
        !       Dim 2 : a0, A, xi, p, q, mass, r0
        ! n_atoms : integer
        !       Number of atoms in the lattice
        !
        ! Returns :
        ! ---------
        ! 2D array
        !       Random x, y, z coordinates of each atoms in the lattice
        !       Dim 1 : Number of atom
        !       Dim 2 : x, y, z coordinates
        !========================================
        implicit none
        double precision, intent(in)            ::  var(:,:), box(3)
        integer, intent(in)                     ::  n_atoms
        double precision                        ::  generate_random(n_atoms,3), x(3)
        integer                                 ::  i, j, cpt = 0
        logical                                 ::  res
        do i = 1, n_atoms
            do
                call random_number(x) ! Generate 3 components position randomly between 0 and 1
                x = x * box(1)
                res = .false. ! .false. if the distance between both atoms i and j is greater than r0
                do j = 1, i-1
                    if (norm2(generate_random(j,:) - x) .lt. (var(i,7) + var(j,7)) / 2.0d0 ) then ! If r < r0
                        res = .true. ! The distance is too short
                        exit
                    end if
                end do
                if (.not. res) exit ! Distance is good
            end do
            generate_random(i,:) = x ! Set the positions x to atom i
            cpt = cpt + 1
            print*, cpt, "atoms created over", n_atoms
        end do
        print*
    end function generate_random


end module
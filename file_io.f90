module file_io

    use constants
    implicit none
    contains


    subroutine read_xyz(path, n_atoms, symbol, pos, box)
        !========================================
        ! Read positions from an extended xyz file (.exyz)
        ! 
        ! Parameters :
        ! ------------
        ! n_atoms : integer
        !       Number of atoms
        ! pos : double precision 2D array
        !       Positions of atoms
        ! epot : double precision
        !       Sum of the potential energy of all atoms in the system
        !
        ! Returns :
        ! ---------
        ! None
        !========================================
        implicit none
        character(len=*), intent(in)                        ::  path
        integer, intent(out)                                ::  n_atoms
        character(len=2), allocatable, intent(out)          ::  symbol(:)
        double precision, allocatable, intent(out)          ::  pos(:,:)
        double precision, intent(out)                       ::  box(3)
        character(len=200)                                  ::  line
        double precision                                    ::  lat(9)
        integer                                             ::  i, p1, p2, ios
        open(10, file=path, status='old', action='read')
        read(10,*) n_atoms
        read(10, '(A)') line
        allocate(symbol(n_atoms))
        allocate(pos(n_atoms,3))
        do i = 1, n_atoms
            read(10,*) symbol(i), pos(i,1), pos(i,2), pos(i,3)
        end do
        close(10)
        p1 = index(line, 'Lattice="')
        p1 = p1 + len('Lattice="')
        p2 = index(line(p1:),'"')
        p2 = p1 + p2 - 2
        read(line(p1:p2), *, iostat=ios) lat
        box(1) = lat(1)
        box(2) = lat(5)
        box(3) = lat(9)
    end subroutine


    subroutine save_xyz(list_atom, n_atoms, pos, epot, box)
        !========================================
        ! Save positions in an extended xyz file (.exyz)
        ! 
        ! Parameters :
        ! ------------
        ! n_atoms : integer
        !       Number of atoms
        ! pos : double precision 2D array
        !       Positions of atoms
        ! epot : double precision
        !       Sum of the potential energy of all atoms in the system
        !
        ! Returns :
        ! ---------
        ! None
        !========================================
        implicit none
        integer, intent(in)                 ::  n_atoms
        double precision, intent(in)        ::  pos(n_atoms,3), epot(n_atoms), box(3)
        character(len=2), intent(in)        ::  list_atom(:)
        character(len=25)                   ::  filename
        integer                             ::  i
        integer, save                       ::  cpt=0
        write(filename, '(A, I0, A)') "XYZ/test", cpt, ".exyz"
        open(10, file=filename, status='replace')
        write(10, '(I0)') n_atoms
        write(10, '(A, 9F12.6, A)') 'Lattice="', box(1), 0.0d0, 0.0d0, 0.0d0, box(2), 0.0d0, 0.0d0, 0.0d0, box(3), &
            '" Properties=species:S:1:pos:R:3:epot:R:1'
        do i = 1, n_atoms
            write(10, '(A,4F12.6)') list_atom(i), pos(i,1), pos(i,2), pos(i,3), epot(i)
        end do
        close(10)
        cpt = cpt + 1
    end subroutine


    subroutine save_data(n_atoms, time, epot, ekin, eth, box)
        !========================================
        ! Save data in a csv file
        ! 
        ! Parameters :
        ! ------------
        ! n_atoms : integer
        !       Number of atoms
        ! time : double precision
        !       Step time
        ! epot : double precision
        !       Sum of the potential energy of all atoms in the system
        ! ekin : double precision
        !       Kinetic energy of the system
        ! eth : double precision
        !       Energy of thermostat
        !
        ! Returns :
        ! ---------
        ! None
        !========================================
        implicit none
        double precision, intent(in)    :: time, epot, ekin, eth, box(3)
        integer, intent(in)             :: n_atoms
        double precision                :: etot, volume
        logical, save                   ::  first_call=.true.
        etot = epot + ekin + eth
        volume = box(1)*box(2)*box(3)
        if (first_call) then
            open(10, file="Data/data.csv", status='replace')
            write(10, '(A)')'n_atoms,time,E_pot,E_kin,E_th,E_tot,volume'
            first_call = .false.
        else
            open(10, file="Data/data.csv", status='old', position='append')
        end if        
        write(10, '(I10, ",", 5(ES20.10, ","),ES20.10)') n_atoms, time, epot, ekin, eth, etot, volume
        close(10)
    end subroutine save_data


end module
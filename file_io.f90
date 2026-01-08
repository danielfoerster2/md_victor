module file_io

    implicit none

    contains

    subroutine read_xyz

        use constants, only: iel_to_typ, n_types
        use variables, only: n_atoms, pos, box, typ

        character(len=255)                                  ::  line, file
        character(len=2)                                    ::  symbol
        integer                                             ::  i_atoms, i_el, p1, p2
        double precision                                    ::  lat(9)
        
        call get_environment_variable("file_input", file)
        open(10, file=file, status='old', action='read')
        read(10, *) n_atoms
        read(10, '(A)') line
        do i_atoms = 1, n_atoms
            read(10, *) symbol, pos(1, i_atoms), pos(2, i_atoms), pos(3, i_atoms)
            do i_el = 1, n_types
                if (symbol .eq. iel_to_typ(i_el)) typ(i_atoms) = i_el
            enddo
        enddo
        close(10)
        p1 = index(line, 'Lattice="')
        p1 = p1 + len('Lattice="')
        p2 = index(line(p1:),'"')
        p2 = p1 + p2 -2
        read(line(p1:p2), *) lat
        box(1) = lat(1)
        box(2) = lat(5)
        box(3) = lat(9)
        write(*,*) "File read : ", file

    endsubroutine


    subroutine save_xyz

        use constants, only: iel_to_typ
        use variables, only: n_atoms, pos, epot, box, typ

        character(len=255)                  ::  file
        integer                             ::  i_atom
        logical, save                       ::  first_call = .true.
        
        call get_environment_variable("file_output", file)
        if (first_call) then
            open(10, file=file, status='replace')
            first_call = .false.
        else
            open(10, file=file, status='old', position='append')
        endif
        write(10, '(I0)') n_atoms
        write(10, '(A, 9ES22.15, A)') 'Lattice="', box(1), 0.0d0, 0.0d0, 0.0d0, box(2), 0.0d0, 0.0d0, 0.0d0, box(3), &
            '" Properties=species:S:1:pos:R:3:epot:R:1'
        do i_atom = 1, n_atoms
            write(10, '(A,4(X,ES22.15))') iel_to_typ(typ(i_atom)), pos(1, i_atom), pos(2, i_atom), pos(3, i_atom), epot(i_atom)
        enddo
        close(10)
        write(*,*) "File saved : ", file

    endsubroutine


    subroutine save_data
        
        use variables, only: n_atoms, epot, time, eth, ekin
        use constants, only: mass
        use potential, only: epot_tbsma

        logical, save                   :: first_call=.true.

        if (first_call) then
            open(10, file="data.dat", status='replace')
            write(10, '(A)')'time, E_pot, E_kin, E_th'
            first_call = .false.
        else
            open(10, file="data.dat", status='old', position='append')
        endif
        write(10, *) time, sum(epot(1:n_atoms)), ekin, eth
        close(10)

    endsubroutine


    subroutine test_velocity_distribution
    
        use variables, only: vel, n_atoms

        integer                 ::  i_atom
        logical, save           ::  first_call=.true.

        if (first_call) then
            open(10, file="vel_init.dat", status='replace')
            write(10, '(A)')'vx, vy, vz'
            do i_atom = 1, n_atoms
                write(10, *) vel(1, i_atom), vel(2, i_atom), vel(3, i_atom)
            enddo
            close(10)
            first_call = .false.
        else
            open(10, file="vel_final.dat", status='replace')
            write(10, '(A)')'vx, vy, vz'
            do i_atom = 1, n_atoms
                write(10, *) vel(1, i_atom), vel(2, i_atom), vel(3, i_atom)
            enddo
            close(10)
        endif
    endsubroutine

end module
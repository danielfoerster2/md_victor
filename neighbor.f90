module neighbor

    implicit none

    contains


    subroutine init_neigh_list

        use constants, only: n_atoms_max, tbsma_rc2, n_neigh_max, skin
        use variables, only: pos, box, n_atoms, neigh, box, n_neigh, typ

        double precision                ::  r_max, rij(3), r, rc2
        integer                         ::  i_atom, j_atom

        neigh = 0
        n_neigh = 0
        do i_atom = 1, n_atoms
            do j_atom = i_atom + 1, n_atoms
                rc2 = tbsma_rc2(typ(i_atom), typ(j_atom))
                r_max = rc2 + skin
                rij = pos(:, j_atom) - pos(:, i_atom)
                rij = rij - box * nint(rij/box)
                !r = norm2(rij)
                r = sum(rij**2)
                if (r .le. r_max**2) then
                    n_neigh(i_atom) = n_neigh(i_atom) + 1
                    n_neigh(j_atom) = n_neigh(j_atom) + 1
                    neigh(n_neigh(i_atom), i_atom) = j_atom
                    neigh(n_neigh(j_atom), j_atom) = i_atom
                endif
                if ( (n_neigh(i_atom) .eq. n_neigh_max) .or. (n_neigh(j_atom) .eq. n_neigh_max )) then
                    write(*,*) 'Neighbor list overflow for atom ', i_atom
                    exit
                endif
            enddo
        enddo

    end subroutine
    

    subroutine update_neigh_list

        use constants, only: skin
        use variables, only: pos, pos_old, n_atoms
        
        double precision                ::  dr
        integer                         ::  i_atoms

        do i_atoms = 1, n_atoms
            dr = norm2(pos(:, i_atoms) - pos_old(:, i_atoms))
            if (dr .gt. 0.5d0*skin) then
                call init_neigh_list
                pos_old = pos
                exit
            endif
        enddo
    endsubroutine

endmodule
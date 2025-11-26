module neighbor

    implicit none

    contains


    subroutine init_neigh_list

        use constants, only: n_atoms_max, tbsma_rc2, n_neigh_max, skin
        use variables, only: pos, box, n_atoms, neigh, n_neigh, typ

        double precision                ::  r_max, rij(3), r, rc2, dr
        double precision, save          ::  pos_old(3, n_atoms_max) = 0.0d0
        integer                         ::  i_atom, j_atom
        logical                         ::  update = .false.

        do i_atom = 1, n_atoms
            dr = norm2(pos(:, i_atom) - pos_old(:, i_atom))
            if (dr .gt. 0.5d0*skin) then
                update = .true.
                exit
            endif
        enddo

        if (update) then
            neigh = 0
            n_neigh = 0
            do i_atom = 1, n_atoms
                do j_atom = i_atom + 1, n_atoms
                    rc2 = tbsma_rc2(typ(i_atom), typ(j_atom))
                    r_max = rc2 + skin
                    rij = pos(:, j_atom) - pos(:, i_atom)
                    rij = rij - box * nint(rij/box)
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
            pos_old = pos
        endif
    end subroutine

endmodule
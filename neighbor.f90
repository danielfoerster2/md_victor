module neighbor

    implicit none

    contains


    subroutine init_neigh_list

        use constants, only: n_atoms_max, tbsma_rc2, n_neigh_max, skin, tbsma_rc2_sq_max
        use variables, only: pos, box, n_atoms, neigh, n_neigh, typ

        double precision                ::  r_max, rij(3), r2, rc2, skin_sq, r_max2
        double precision, save          ::  pos_old(3, n_atoms_max) = 0.0d0
        integer                         ::  i_atom, j_atom
        logical                         ::  update

        update = .false.
        skin_sq = 0.25d0*skin**2
        do i_atom = 1, n_atoms
            if (dot_product(pos(:, i_atom) - pos_old(:, i_atom), pos(:, i_atom) - pos_old(:, i_atom)) .gt. skin_sq) then
                update = .true.
                exit
            endif
        enddo

        if (update) then
            r_max2 = (sqrt(tbsma_rc2_sq_max) + skin)**2
            neigh = 0
            n_neigh = 0
            do i_atom = 1, n_atoms
                do j_atom = i_atom + 1, n_atoms
                    rij = pos(:, j_atom) - pos(:, i_atom)
                    rij = rij - box * nint(rij/box)
                    r2 = dot_product(rij, rij)
                    if (r2 .le. r_max2) then
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
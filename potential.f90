module potential

    use constants

    implicit none

    contains
 
    subroutine init_tbsma

        use constants, only: tbsma_a, tbsma_xi, tbsma_p, tbsma_q, tbsma_r0, tbsma_rc1, tbsma_rc2, tbsma_x5, tbsma_x4, tbsma_x3,&
        & tbsma_a5, tbsma_a4, tbsma_a3
        double precision ::  ar, br, cr, ab, bb, cb

        ! Ni parameters (Cleri-Rosato)
        tbsma_a(1, 1)  = 0.0376d0
        tbsma_xi(1, 1) = 1.070d0
        tbsma_p(1, 1)  = 16.999d0
        tbsma_q(1, 1)  = 1.189d0
        tbsma_r0(1, 1) = 3.523/2**0.5d0

        ! Cu parameters (Cleri-Rosato)
        ! tbsma_a(1, 1)  = 0.0855d0
        ! tbsma_xi(1, 1) = 1.224d0
        ! tbsma_p(1, 1)  = 10.960d0
        ! tbsma_q(1, 1)  = 2.278d0
        ! tbsma_r0(1, 1) = 3.615/2**0.5d0

        ! Rh parameters (Cleri-Rosato table)
        ! tbsma_a(1, 1)  = 0.0629d0
        ! tbsma_xi(1, 1) = 1.660d0
        ! tbsma_p(1, 1)  = 18.450d0
        ! tbsma_q(1, 1)  = 1.867d0
        ! tbsma_r0(1, 1) = 3.803d0/2**0.5d0

        tbsma_rc1(1, 1) = 5.0d0**0.5*tbsma_r0(1, 1) + 0.05d0
        tbsma_rc2(1, 1) = tbsma_rc1(1, 1) + 0.3d0

        ar =-tbsma_a(1, 1)*exp(-tbsma_p(1, 1)*(tbsma_rc1(1, 1)/tbsma_r0(1, 1)-1))/(tbsma_rc2(1, 1)-tbsma_rc1(1, 1))**3
        br =-(tbsma_p(1, 1)/tbsma_r0(1, 1))*tbsma_a(1, 1)*exp(-tbsma_p(1, 1)*(tbsma_rc1(1, 1)/tbsma_r0(1, 1)-1))/(tbsma_rc2(1, 1)&
        &-tbsma_rc1(1, 1))**2
        cr =-((tbsma_p(1, 1)/tbsma_r0(1, 1))**2)*tbsma_a(1, 1)*exp(-tbsma_p(1, 1)*(tbsma_rc1(1, 1)/tbsma_r0(1, 1)-1))&
        &/(tbsma_rc2(1, 1)-tbsma_rc1(1, 1))
        ab =-tbsma_xi(1, 1)*exp(-tbsma_q(1, 1)*(tbsma_rc1(1, 1)/tbsma_r0(1, 1)-1))/(tbsma_rc2(1, 1)-tbsma_rc1(1, 1))**3
        bb =-(tbsma_q(1, 1)/tbsma_r0(1, 1))*tbsma_xi(1, 1)*exp(-tbsma_q(1, 1)*(tbsma_rc1(1, 1)/tbsma_r0(1, 1)-1))/(tbsma_rc2(1, 1)&
        &-tbsma_rc1(1, 1))**2
        cb =-((tbsma_q(1, 1)/tbsma_r0(1, 1))**2)*tbsma_xi(1, 1)*exp(-tbsma_q(1, 1)*(tbsma_rc1(1, 1)/tbsma_r0(1, 1)-1))&
        &/(tbsma_rc2(1, 1)-tbsma_rc1(1, 1))

        tbsma_x5(1, 1) = (12*ab-6*bb+cb)/(2*(tbsma_rc2(1, 1)-tbsma_rc1(1, 1))**2)
        tbsma_x4(1, 1) = (15*ab-7*bb+cb)/(tbsma_rc2(1, 1)-tbsma_rc1(1, 1))
        tbsma_x3(1, 1) = (20*ab-8*bb+cb)/2
        tbsma_a5(1, 1) = (12*ar-6*br+cr)/(2*(tbsma_rc2(1, 1)-tbsma_rc1(1, 1))**2)
        tbsma_a4(1, 1) = (15*ar-7*br+cr)/(tbsma_rc2(1, 1)-tbsma_rc1(1, 1))
        tbsma_a3(1, 1) = (20*ar-8*br+cr)/2
    endsubroutine


    subroutine epot_tbsma
        use constants, only: tbsma_a, tbsma_xi, tbsma_p, tbsma_q, tbsma_r0, tbsma_rc1, tbsma_rc2
        use variables, only: pos, box, n_atoms, epot, typ, neigh, n_neigh

        implicit none
        double precision                    ::  rij(3), r
        double precision                    ::  a, xi, p, q, r0, rc1, rc2, x5, x4, x3, a5, a4, a3
        double precision                    ::  e_rep, rho
        integer                             ::  i_atom, j_atom, j_neigh

        do i_atom = 1, n_atoms
            e_rep = 0.0d0
            rho = 0.0d0
            do j_neigh = 1, n_neigh(i_atom)
                j_atom = neigh(j_neigh, i_atom)
                rij = pos(:, j_atom) - pos(:, i_atom)
                rij = rij - box * nint(rij/box)
                r = norm2(rij)

                a = tbsma_a(typ(i_atom), typ(j_atom))
                xi = tbsma_xi(typ(i_atom), typ(j_atom))
                p = tbsma_p(typ(i_atom), typ(j_atom))
                q = tbsma_q(typ(i_atom), typ(j_atom))
                r0 = tbsma_r0(typ(i_atom), typ(j_atom))
                rc1 = tbsma_rc1(typ(i_atom), typ(j_atom))
                rc2 = tbsma_rc2(typ(i_atom), typ(j_atom))

                x5 = tbsma_x5(typ(i_atom), typ(j_atom))
                x4 = tbsma_x4(typ(i_atom), typ(j_atom))
                x3 = tbsma_x3(typ(i_atom), typ(j_atom))
                a5 = tbsma_a5(typ(i_atom), typ(j_atom))
                a4 = tbsma_a4(typ(i_atom), typ(j_atom))
                a3 = tbsma_a3(typ(i_atom), typ(j_atom))
                
                if (r .lt. rc2) then

                    if (r .gt. rc1) then
                        e_rep = e_rep + a5*(r-rc2)**5 + a4*(r-rc2)**4 + a3*(r-rc2)**3
                        rho = rho     + (x5*(r-rc2)**5 + x4*(r-rc2)**4 + x3*(r-rc2)**3)**2
                    else
                        e_rep = e_rep + a * exp(-p * (r/r0 - 1.0d0))
                        rho = rho + xi**2 * exp(-2.0d0 * q * (r/r0 - 1.0d0))
                    endif
                endif
            enddo
            epot(i_atom) = e_rep - sqrt(rho)
        enddo
    endsubroutine


    subroutine force_tbsma
        use constants, only: n_atoms_max
        use variables, only: pos, box, n_atoms, typ, force, neigh, n_neigh

        double precision                    ::  rho(n_atoms_max), rij(3), r
        double precision                    ::  rep, band, tmp
        double precision                    ::  a, xi, p, q, r0, rc1, rc2, x5, x4, x3, a5, a4, a3
        integer                             ::  i_atom, j_atom, j_neigh

        rho(1:n_atoms) = 0.0d0
        do i_atom = 1, n_atoms
            do j_neigh = 1, n_neigh(i_atom)
                j_atom = neigh(j_neigh, i_atom)
                rij = pos(:, j_atom) - pos(:, i_atom)
                rij = rij - box * nint(rij/box)
                r = norm2(rij)
                xi = tbsma_xi(typ(i_atom), typ(j_atom))
                q = tbsma_q(typ(i_atom), typ(j_atom))
                r0 = tbsma_r0(typ(i_atom), typ(j_atom))
                rc1 = tbsma_rc1(typ(i_atom), typ(j_atom))
                rc2 = tbsma_rc2(typ(i_atom), typ(j_atom))
                x5 = tbsma_x5(typ(i_atom), typ(j_atom))
                x4 = tbsma_x4(typ(i_atom), typ(j_atom))
                x3 = tbsma_x3(typ(i_atom), typ(j_atom))

                if (r .lt. rc2) then
                    if (r .gt. rc1) then
                        tmp = (x5*(r-rc2)**5 + x4*(r-rc2)**4 + x3*(r-rc2)**3)**2
                    else
                        tmp = xi**2 * exp(-2.0d0 * q * (r/r0 - 1.0d0))
                    endif
                    rho(i_atom) = rho(i_atom) + tmp
                endif
            enddo
        enddo

        force(:, 1:n_atoms) = 0.0d0
        do i_atom = 1, n_atoms
            do j_neigh = 1, n_neigh(i_atom)
                j_atom = neigh(j_neigh, i_atom)
                if (j_atom .gt. i_atom) cycle
                rij = pos(:, j_atom) - pos(:, i_atom)
                rij = rij - box * nint(rij/box)
                r = norm2(rij)
                a = tbsma_a(typ(i_atom), typ(j_atom))
                xi = tbsma_xi(typ(i_atom), typ(j_atom))
                p = tbsma_p(typ(i_atom), typ(j_atom))
                q = tbsma_q(typ(i_atom), typ(j_atom))
                r0 = tbsma_r0(typ(i_atom), typ(j_atom))
                rc1 = tbsma_rc1(typ(i_atom), typ(j_atom))
                rc2 = tbsma_rc2(typ(i_atom), typ(j_atom))

                x5 = tbsma_x5(typ(i_atom), typ(j_atom))
                x4 = tbsma_x4(typ(i_atom), typ(j_atom))
                x3 = tbsma_x3(typ(i_atom), typ(j_atom))
                a5 = tbsma_a5(typ(i_atom), typ(j_atom))
                a4 = tbsma_a4(typ(i_atom), typ(j_atom))
                a3 = tbsma_a3(typ(i_atom), typ(j_atom))

                if (r .lt. rc2) then
                    if (r .gt. rc1) then
                        rep = 2.0d0*(5*a5*(r-rc2)**4 + 4*a4*(r-rc2)**3 + 3*a3*(r-rc2)**2)
                        band = -(5*x5*(r-rc2)**4 + 4*x4*(r-rc2)**3 + 3*x3*(r-rc2)**2) &
                                &*(x5*(r-rc2)**5 + x4*(r-rc2)**4 + x3*(r-rc2)**3) &
                                &* (1.0d0/sqrt(rho(j_atom)) + 1.0d0/sqrt(rho(i_atom)))
                    else
                        rep =  -2.0d0 * a * p/r0 * exp(-p * (r/r0 - 1.0d0))
                        band =  q * xi**2 / r0 * exp(-2.0d0 * q * (r/r0 - 1.0d0)) &
                                &* (1.0d0/sqrt(rho(j_atom)) + 1.0d0/sqrt(rho(i_atom)))
                    endif
                    force(:, i_atom) = force(:, i_atom) + (rep + band) * rij/r
                    force(:, j_atom) = force(:, j_atom) - (rep + band) * rij/r
                endif
            enddo
        enddo
    endsubroutine



    subroutine test_force_tbsma

        use variables, only: pos, epot, n_atoms, force

        double precision, parameter         ::  h = 1.0d-5
        double precision                    ::  force_num, e1, e2
        integer                             ::  i_atom, k

        call force_tbsma

        do i_atom = 1, n_atoms
            do k = 1, 3
                pos(k, i_atom) = pos(k, i_atom) + h
                call epot_tbsma
                e1 = sum(epot)
                pos(k, i_atom) = pos(k, i_atom) - 2.0d0 * h
                call epot_tbsma
                e2 = sum(epot)
                pos(k, i_atom) = pos(k, i_atom) + h
                force_num = -(e1 - e2) / (2.0d0 * h)
                if (abs(force(k, i_atom) - force_num)/ abs(force_num + 1.0d-5) .gt. 1.0d-2) then
                    write(*, *) "Atom ", i_atom, "Direction ", k
                    write(*, *) "Analytical force :", force(k, i_atom)
                    write(*, *) "Numerical force  :", force_num
                    stop "Force test failed!"
                endif
            enddo
        enddo
    endsubroutine

endmodule

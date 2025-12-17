module potential

    use constants, only: n_types

    implicit none


    integer, parameter      :: n_tab = 50000
    double precision        :: tab_rho(n_types, n_types, n_tab+1), tab_drho(n_types, n_types, n_tab+1),&
                                & tab_drep(n_types, n_types, n_tab+1), inv_rc2_sq_tab


    contains
 
    subroutine init_tbsma

        use constants, only: tbsma_a, tbsma_xi, tbsma_p, tbsma_q, tbsma_r0, tbsma_rc1, tbsma_rc2, tbsma_x5, tbsma_x4, tbsma_x3,&
        & tbsma_a5, tbsma_a4, tbsma_a3, tbsma_rc2_sq_max, tbsma_rc2_sq
        double precision ::  ar, br, cr, ab, bb, cb, a, p, q, xi, r0, rc1, rc2, r, a3, a4, a5, x3, x4, x5, r2

        integer :: i_tab, i_el, j_el

        ! Ni parameters (Cleri-Rosato)
        tbsma_a(1, 1)  = 0.0376d0
        tbsma_xi(1, 1) = 1.070d0
        tbsma_p(1, 1)  = 16.999d0
        tbsma_q(1, 1)  = 1.189d0
        tbsma_r0(1, 1) = 3.523/2**0.5d0
        
        ! Ag parameters
        tbsma_a(6, 6)  = 0.10433d0
        tbsma_xi(6, 6) = 1.194019d0
        tbsma_p(6, 6)  = 10.79d0
        tbsma_q(6, 6)  = 3.19d0
        tbsma_r0(6, 6) = 2.0d0*1.445d0
        tbsma_rc1(6, 6) = 4.08707719d0
        tbsma_rc2(6, 6) = 5.0056268338740553d0

        ! Co parameters
        tbsma_a(7, 7)  = 0.175700d0
        tbsma_xi(7, 7) = 1.843d0
        tbsma_p(7, 7)  = 9.210d0
        tbsma_q(7, 7)  = 2.975d0
        tbsma_r0(7, 7) = 2.0d0*1.25d0
        tbsma_rc1(7, 7) = 3.53553391d0
        tbsma_rc2(7, 7) = 4.3301270189221932d0

        ! AgCo parameters
        tbsma_a(6, 7)  = 0.15202d0
        tbsma_xi(6, 7) = 1.4319d0
        tbsma_p(6, 7)  = 10.001d0
        tbsma_q(6, 7)  = 3.085d0
        tbsma_r0(6, 7) = 1.445d0+1.25d0
        tbsma_rc1(6, 7) = 4.08707719d0
        tbsma_rc2(6, 7) = 4.3301270189221932d0

        tbsma_a(7, 6)  = tbsma_a(6, 7)
        tbsma_xi(7, 6) = tbsma_xi(6, 7)
        tbsma_p(7, 6)  = tbsma_p(6, 7)
        tbsma_q(7, 6)  = tbsma_q(6, 7)
        tbsma_r0(7, 6) = tbsma_r0(6, 7)
        tbsma_rc1(7, 6) = tbsma_rc1(6, 7)
        tbsma_rc2(7, 6) = tbsma_rc2(6, 7)

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

        !tbsma_rc1(1, 1) = 3.0d0**0.5*tbsma_r0(1, 1) + 0.05d0
        !tbsma_rc2(1, 1) = tbsma_rc1(1, 1) + 0.3d0
        !tbsma_rc2_sq(1, 1) = tbsma_rc2(1, 1)**2
        !tbsma_rc2_sq_max = tbsma_rc2(1, 1)**2 ! maxval(tbsma_rc2**2) XXX

        do i_el = 1, n_types
            do j_el = 1, n_types

                tbsma_rc1(i_el, j_el) = 3.0d0**0.5*tbsma_r0(i_el, j_el) + 0.05d0
                tbsma_rc2(i_el, j_el) = tbsma_rc1(i_el, j_el) + 0.3d0
                tbsma_rc2_sq(i_el, j_el) = tbsma_rc2(i_el, j_el)**2
                tbsma_rc2_sq_max = tbsma_rc2(i_el, j_el)**2 ! maxval(tbsma_rc2**2) XXX

                ar =-tbsma_a(i_el, j_el)*exp(-tbsma_p(i_el, j_el)*(tbsma_rc1(i_el, j_el)/tbsma_r0(i_el, j_el)-1))&
                    &/(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))**3
                br =-(tbsma_p(i_el, j_el)/tbsma_r0(i_el, j_el))*tbsma_a(i_el, j_el)*exp(-tbsma_p(i_el, j_el)&
                    &*(tbsma_rc1(i_el, j_el)/tbsma_r0(i_el, j_el)-1))/(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))**2
                cr =-((tbsma_p(i_el, j_el)/tbsma_r0(i_el, j_el))**2)*tbsma_a(i_el, j_el)*exp(-tbsma_p(i_el, j_el)&
                    &*(tbsma_rc1(i_el, j_el)/tbsma_r0(i_el, j_el)-1))/(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))
                ab =-tbsma_xi(i_el, j_el)*exp(-tbsma_q(i_el, j_el)*(tbsma_rc1(i_el, j_el)/tbsma_r0(i_el, j_el)-1))&
                    &/(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))**3
                bb =-(tbsma_q(i_el, j_el)/tbsma_r0(i_el, j_el))*tbsma_xi(i_el, j_el)*exp(-tbsma_q(i_el, j_el)&
                    &*(tbsma_rc1(i_el, j_el)/tbsma_r0(i_el, j_el)-1))/(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))**2
                cb =-((tbsma_q(i_el, j_el)/tbsma_r0(i_el, j_el))**2)*tbsma_xi(i_el, j_el)*exp(-tbsma_q(i_el, j_el)&
                    &*(tbsma_rc1(i_el, j_el)/tbsma_r0(i_el, j_el)-1))/(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))

                tbsma_x5(i_el, j_el) = (12*ab-6*bb+cb)/(2*(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))**2)
                tbsma_x4(i_el, j_el) = (15*ab-7*bb+cb)/(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))
                tbsma_x3(i_el, j_el) = (20*ab-8*bb+cb)/2
                tbsma_a5(i_el, j_el) = (12*ar-6*br+cr)/(2*(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))**2)
                tbsma_a4(i_el, j_el) = (15*ar-7*br+cr)/(tbsma_rc2(i_el, j_el)-tbsma_rc1(i_el, j_el))
                tbsma_a3(i_el, j_el) = (20*ar-8*br+cr)/2
            enddo
        enddo

        inv_rc2_sq_tab = n_tab/ tbsma_rc2_sq_max

        do i_tab=1, n_tab+1
            r2 = i_tab * (tbsma_rc2_sq_max / n_tab)
            r = sqrt(r2)
            a = tbsma_a(1, 1)
            xi = tbsma_xi(1, 1)
            p = tbsma_p(1, 1)
            q = tbsma_q(1, 1)
            r0 = tbsma_r0(1, 1)
            rc1 = tbsma_rc1(1, 1)
            rc2 = tbsma_rc2(1, 1)

            x5 = tbsma_x5(1, 1)
            x4 = tbsma_x4(1, 1)
            x3 = tbsma_x3(1, 1)
            a5 = tbsma_a5(1, 1)
            a4 = tbsma_a4(1, 1)
            a3 = tbsma_a3(1, 1)

            if (r .gt. rc2) then
                tab_rho(1, 1, i_tab) = 0.0d0
                tab_drep(1, 1, i_tab) = 0.0d0
                tab_drho(1, 1, i_tab) = 0.0d0
            elseif (r .gt. rc1) then
                tab_rho(1, 1, i_tab) = (x5*(r-rc2)**5 + x4*(r-rc2)**4 + x3*(r-rc2)**3)**2
                tab_drep(1, 1, i_tab) = 2.0d0*(5*a5*(r-rc2)**4 + 4*a4*(r-rc2)**3 + 3*a3*(r-rc2)**2)/r
                tab_drho(1, 1, i_tab) = -(5*x5*(r-rc2)**4 + 4*x4*(r-rc2)**3 + 3*x3*(r-rc2)**2) &
                                   &*(x5*(r-rc2)**5 + x4*(r-rc2)**4 + x3*(r-rc2)**3) / r
            else
                tab_rho(1, 1, i_tab) = xi**2 * exp(-2.0d0 * q * (r/r0 - 1.0d0))
                tab_drep(1, 1, i_tab) = -2.0d0 * a * p/r0 * exp(-p * (r/r0 - 1.0d0)) / r
                tab_drho(1, 1, i_tab) = q * xi**2 / r0 * exp(-2.0d0 * q * (r/r0 - 1.0d0)) / r
            endif
        enddo
    endsubroutine


    subroutine epot_tbsma
        use constants, only: tbsma_a, tbsma_xi, tbsma_p, tbsma_q, tbsma_r0, tbsma_rc1, tbsma_rc2, tbsma_a3, tbsma_a4, tbsma_a5,&
                    & tbsma_x3, tbsma_x4, tbsma_x5
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


    subroutine force_tbsma_tab
        use constants, only: n_atoms_max, tbsma_rc2_sq_max
        use variables, only: pos, box, n_atoms, typ, force, neigh, n_neigh

        double precision                    ::  rho(n_atoms_max), rij(3), r2, tmp, tmp2(3)
        integer                             ::  i_atom, j_atom, j_neigh, i_tab
        rho(1:n_atoms) = 0.0d0
        do i_atom = 1, n_atoms
            do j_neigh = 1, n_neigh(i_atom)
                j_atom = neigh(j_neigh, i_atom)
                rij = pos(:, j_atom) - pos(:, i_atom)
                rij = rij - box * nint(rij/box)
                r2 = dot_product(rij, rij)

                if (r2 .ge. tbsma_rc2_sq_max) cycle
                tmp =  r2 * inv_rc2_sq_tab
                i_tab = floor(tmp)

                rho(i_atom) = rho(i_atom) + tab_rho(typ(i_atom), typ(j_atom), i_tab) + &
                    &(tab_rho(typ(i_atom), typ(j_atom), i_tab+1) - tab_rho(typ(i_atom), typ(j_atom), i_tab)) * (tmp - i_tab)
            enddo
        enddo

        rho(1:n_atoms) = rho(1:n_atoms)**(-0.5)
        force(:, 1:n_atoms) = 0.0d0
        do i_atom = 1, n_atoms
            do j_neigh = 1, n_neigh(i_atom)
                j_atom = neigh(j_neigh, i_atom)
                if (j_atom .gt. i_atom) cycle
                rij = pos(:, j_atom) - pos(:, i_atom)
                rij = rij - box * nint(rij/box)
                r2 = dot_product(rij, rij)

                if (r2 .ge. tbsma_rc2_sq_max) cycle
                tmp =  r2 * inv_rc2_sq_tab
                i_tab = floor(tmp)

                tmp2 = (tab_drep(typ(i_atom), typ(j_atom), i_tab) + (tab_drep(typ(i_atom), typ(j_atom), i_tab+1) &
                &- tab_drep(typ(i_atom), typ(j_atom), i_tab)) * (tmp - i_tab) + (tab_drho(typ(i_atom), typ(j_atom), i_tab) &
                &+ (tab_drho(typ(i_atom), typ(j_atom), i_tab+1) &
                &- tab_drho(typ(i_atom), typ(j_atom), i_tab)) * (tmp - i_tab))*(rho(j_atom) + rho(i_atom))) * rij

                force(:, i_atom) = force(:, i_atom) + tmp2
                force(:, j_atom) = force(:, j_atom) - tmp2
            enddo
        enddo
    endsubroutine

    subroutine force_tbsma
        use constants, only: tbsma_a, tbsma_xi, tbsma_p, tbsma_q, tbsma_r0, tbsma_rc1, tbsma_rc2, tbsma_x5, tbsma_x4, tbsma_x3,&
        & tbsma_a5, tbsma_a4, tbsma_a3, n_atoms_max, tbsma_rc2_sq
        use variables, only: pos, box, n_atoms, typ, force, neigh, n_neigh

        double precision                    ::  rho(n_atoms_max), rij(3), r, r2
        double precision                    ::  rep, band, tmp
        double precision                    ::  a, xi, p, q, r0, rc1, rc2, x5, x4, x3, a5, a4, a3, tmp2(3)
        integer                             ::  i_atom, j_atom, j_neigh

        rho(1:n_atoms) = 0.0d0
        do i_atom = 1, n_atoms
            do j_neigh = 1, n_neigh(i_atom)
                j_atom = neigh(j_neigh, i_atom)
                rij = pos(:, j_atom) - pos(:, i_atom)
                rij = rij - box * nint(rij/box)
                r2 = dot_product(rij, rij)

                if (r2 .lt. tbsma_rc2_sq(typ(i_atom), typ(j_atom))) then
                    r = sqrt(r2)
                    if (r .gt. tbsma_rc1(typ(i_atom), typ(j_atom))) then
                        x5 = tbsma_x5(typ(i_atom), typ(j_atom))
                        x4 = tbsma_x4(typ(i_atom), typ(j_atom))
                        x3 = tbsma_x3(typ(i_atom), typ(j_atom))
                        rc2 = tbsma_rc2(typ(i_atom), typ(j_atom))
                        tmp = (x5*(r-rc2)**5 + x4*(r-rc2)**4 + x3*(r-rc2)**3)**2
                    else
                        xi = tbsma_xi(typ(i_atom), typ(j_atom))
                        q = tbsma_q(typ(i_atom), typ(j_atom))
                        r0 = tbsma_r0(typ(i_atom), typ(j_atom))
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
                r2 = dot_product(rij, rij)

                if (r2 .lt. tbsma_rc2_sq(typ(i_atom), typ(j_atom))) then
                    r = sqrt(r2)
                    if (r .gt. tbsma_rc1(typ(i_atom), typ(j_atom))) then
                        x5 = tbsma_x5(typ(i_atom), typ(j_atom))
                        x4 = tbsma_x4(typ(i_atom), typ(j_atom))
                        x3 = tbsma_x3(typ(i_atom), typ(j_atom))
                        a5 = tbsma_a5(typ(i_atom), typ(j_atom))
                        a4 = tbsma_a4(typ(i_atom), typ(j_atom))
                        a3 = tbsma_a3(typ(i_atom), typ(j_atom))
                        rc2 = tbsma_rc2(typ(i_atom), typ(j_atom))
                        rep = 2.0d0*(5*a5*(r-rc2)**4 + 4*a4*(r-rc2)**3 + 3*a3*(r-rc2)**2)
                        band = -(5*x5*(r-rc2)**4 + 4*x4*(r-rc2)**3 + 3*x3*(r-rc2)**2) &
                                &*(x5*(r-rc2)**5 + x4*(r-rc2)**4 + x3*(r-rc2)**3) &
                                &* (1.0d0/sqrt(rho(j_atom)) + 1.0d0/sqrt(rho(i_atom)))
                    else
                        a = tbsma_a(typ(i_atom), typ(j_atom))
                        xi = tbsma_xi(typ(i_atom), typ(j_atom))
                        p = tbsma_p(typ(i_atom), typ(j_atom))
                        q = tbsma_q(typ(i_atom), typ(j_atom))
                        r0 = tbsma_r0(typ(i_atom), typ(j_atom))
                        rep =  -2.0d0 * a * p/r0 * exp(-p * (r/r0 - 1.0d0))
                        band =  q * xi**2 / r0 * exp(-2.0d0 * q * (r/r0 - 1.0d0)) &
                                &* (1.0d0/sqrt(rho(j_atom)) + 1.0d0/sqrt(rho(i_atom)))
                    endif
                    tmp2 = (rep + band) * rij/r
                    force(:, i_atom) = force(:, i_atom) + tmp2
                    force(:, j_atom) = force(:, j_atom) - tmp2
                endif
            enddo
        enddo
    endsubroutine



    subroutine test_force_tbsma

        use variables, only: pos, epot, n_atoms, force

        double precision, parameter         ::  h = 1.0d-5
        double precision                    ::  force_num, e1, e2
        integer                             ::  i_atom, k

        call force_tbsma_tab

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
                if (abs(force(k, i_atom) - force_num)/ abs(force_num + 1.0d-5) .gt. 1.0d-5) then
                    write(*, *) "Atom ", i_atom, "Direction ", k
                    write(*, *) "Analytical force :", force(k, i_atom)
                    write(*, *) "Numerical force  :", force_num
                endif
            enddo
        enddo
        !stop "Force test failed!"
    endsubroutine

endmodule

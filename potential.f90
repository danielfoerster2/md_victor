module potential

    use constants

    implicit none

    contains
 
    subroutine init_potential

        use constants, only: tbsma_a, tbsma_xi, tbsma_p, tbsma_q, tbsma_r0, tbsma_rc1, tbsma_rc2, tbsma_x5, tbsma_x4, tbsma_x3,&
        & tbsma_a5, tbsma_a4, tbsma_a3
        double precision ::  ar, br, cr, ab, bb, cb

        ! tbsma_a(1, 1) = 0.0376d0
        ! tbsma_xi(1, 1) = 1.070d0
        ! tbsma_p(1, 1) = 16.999d0
        ! tbsma_q(1, 1) = 1.189d0
        ! tbsma_r0(1, 1) = 3.523/2**0.5d0

        tbsma_a(1, 1) = 0.0855d0
        tbsma_xi(1, 1) = 1.224d0
        tbsma_p(1, 1) = 10.960d0
        tbsma_q(1, 1) = 2.278d0
        tbsma_r0(1, 1) = 3.615/2**0.5d0
        tbsma_rc1(1, 1) = 5.0d0**0.5*tbsma_r0(1, 1) + 0.1d0
        tbsma_rc2(1, 1) = tbsma_rc1(1, 1) + 0.2d0

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

    subroutine tbsma
        use constants, only: tbsma_a, tbsma_xi, tbsma_p, tbsma_q, tbsma_r0, tbsma_rc1, tbsma_rc2
        use variables, only: pos, box, n_atoms, epot, typ

        implicit none
        double precision                    ::  rij(3), r
        double precision                    ::  a, xi, p, q, r0, rc1, rc2, x5, x4, x3, a5, a4, a3
        double precision                    ::  e_rep, band
        integer                             ::  i_atom, j_atom

        do i_atom = 1, n_atoms
            e_rep = 0.0d0
            band = 0.0d0
            do j_atom = 1, n_atoms
                if (i_atom .ne. j_atom) then
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

                    x5 = tbsma_rc2(typ(i_atom), typ(j_atom))
                    x4 = tbsma_rc2(typ(i_atom), typ(j_atom))
                    x3 = tbsma_rc2(typ(i_atom), typ(j_atom))
                    a5 = tbsma_a5(typ(i_atom), typ(j_atom))
                    a4 = tbsma_a4(typ(i_atom), typ(j_atom))
                    a3 = tbsma_a3(typ(i_atom), typ(j_atom))

                    if (r .lt. tbsma_rc2(typ(i_atom), typ(j_atom))) then

                        if (r .gt. tbsma_rc1(typ(i_atom), typ(j_atom))) then
                            e_rep = e_rep + x5*(r-rc2)**5 + x4*(r-rc2)**4 + x3*(r-rc2)**3
                            band = band + a5*(r-rc2)**5 + a4*(r-rc2)**4 + a3*(r-rc2)**3
                        else
                            e_rep = e_rep + a * exp(-p * (r/r0 - 1.0d0))
                            band = band + xi**2 * exp(-2.0d0 * q * (r/r0 - 1.0d0))
                        endif
                    endif
                endif
            enddo
            epot(i_atom) = e_rep - sqrt(band)
        enddo
    endsubroutine

endmodule
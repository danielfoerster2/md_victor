module potential

    use constants
    implicit none
    contains
 

    function tbsma(var, pos,box, use_cut_function)
        implicit none
        double precision, intent(in)        ::  pos(:,:), var(:,:), box(3)
        logical, intent(in)                 ::  use_cut_function
        double precision                    ::  tbsma(size(pos, dim=1))
        double precision                    ::  rij(3), r, x
        double precision                    ::  a0, A, xi, p, q, r0, rc1, rc2
        double precision                    ::  rep, band, E_rep, dE_rep, E_band, dE_band
        integer                             ::  i, j
        do i = 1, size(pos, dim=1)
            rep = 0.0d0
            band = 0.0d0
            do j = 1, size(pos, dim=1)
                if (i .ne. j) then
                    a0 = (var(i,1) + var(j,1)) / 2.0d0
                    A = (var(i,2) + var(j,2)) / 2.0d0
                    xi = (var(i,3) + var(j,3)) / 2.0d0
                    p = (var(i,4) + var(j,4)) / 2.0d0
                    q = (var(i,5) + var(j,5)) / 2.0d0
                    r0 = (var(i,7) + var(j,7)) / 2.0d0
                    rij = pos(j,:) - pos(i,:)
                    rij = rij - box * nint(rij/box)
                    r = norm2(rij)
                    !rc2 = (sqrt(dble(N))+ sqrt(dble(N+1)))*r0/2
                    rc2 = (sqrt(dble(N))+ sqrt(dble(N+1)))*var(i,7) / 2.0d0
                    if (r .lt. rc2) then
                        if (.not. use_cut_function) then
                            rep = rep + A * exp(-p * (r /r0 - 1.0d0)) ! A*exp(-p*(r/r0 - 1))
                            band = band + xi**2 * exp(-2.0d0 * q * (r/r0 - 1.0d0)) ! xi**2 * exp(-2q(r/r0 - 1))
                        else
                            rc1 = 0.95d0 * rc2
                            if (r .gt. rc1) then
                                E_rep = A * exp(-p*(rc1/r0 - 1.0d0))
                                dE_rep = -p * A / r0 * E_rep
                                E_band = xi**2 * exp(-2.0d0 * q * (rc1/r0 - 1.0d0))
                                dE_band = -2.0d0 * q / r0 * E_band
                                rep = rep + cut_function(E_rep, dE_rep, r, rc1, rc2)
                                band = band + cut_function(E_band, dE_band, r, rc1, rc2)
                            else
                                rep = rep + A * exp(-p * (r/r0 - 1.0d0))
                                band = band + xi**2 * exp(-2.0d0 * q * (r/r0 - 1.0d0))
                            end if
                        end if
                    end if
                end if
            end do
            tbsma(i) = rep - sqrt(band)
        end do
    end function


    function cut_function(E, dE, r, rc1, rc2)
        implicit none
        double precision, intent(in)    ::  E, dE, r, rc1, rc2
        double precision                ::  cut_function
        double precision                ::  b0, b1, b2, b3, b4, b5
        double precision                ::  h
        h = rc2 - rc1
        b0 = E
        b1 = dE
        b2 = 0.0d0
        b3 = (-10.0d0*E - 6.0d0*dE*h) / (h**3)
        b4 = (15.0d0*E + 8.0d0*dE*h) / (h**4)
        b5 = (-6.0d0*E - 3.0d0*dE*h) / (h**5)
        cut_function = b0 + b1*(r-rc1) + b2*(r-rc1)**2 + b3*(r-rc1)**3 + b4*(r-rc1)**4 + b5*(r-rc1)**5
    end function


    subroutine test_tbsma(pos, epot)
        implicit none
        double precision, intent(in)        ::  pos(:,:), epot(:)
        integer                             ::  i
        do i = 1, size(pos, dim=1)
            print*, "Atom", i, ":", epot(i), "eV"
        end do
    end subroutine


end module
module force

    use constants
    use potential
    implicit none
    contains


    function force_tbsma(var, pos, box, use_cut_function)
        implicit none
        double precision, intent(in)        ::  var(:,:), pos(:,:), box(3)
        logical, intent(in)                 ::  use_cut_function
        double precision                    ::  force_tbsma(size(pos, dim=1), 3)
        double precision                    ::  rho(size(pos, dim=1)), rij(3), r, rkj(3)
        double precision                    ::  E_band, E_rep, dE_band, dE_rep, rc1, rc2
        double precision                    ::  rep, band, fcut_rep, fcut_band, dfcut_band, dfcut_rep
        double precision                    ::  a0, A, xi, p, q, r0
        integer                             ::  i, j, k
        do i = 1, size(pos, dim=1)
            rho(i) = 0.0d0
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
                    rc2 = (sqrt(dble(N))+ sqrt(dble(N+1)))*var(i,7) / 2.0d0
                    if (r .lt.rc2) then
                        if (.not. use_cut_function) then
                            rho(i) = rho(i) + xi**2 * exp(-2.0d0 * q * (r/r0 - 1.0d0))
                        else
                            rc1 = 0.95d0 * rc2
                            if (r .gt. rc1) then
                                E_band = xi**2 * exp(-2.0d0 * q * (rc1/r0 - 1.0d0))
                                dE_band = -2.0d0 * q / r0 * E_band
                                fcut_band = cut_function(E_band, dE_band, r, rc1, rc2)
                                rho(i) = rho(i) + fcut_band
                            else
                                rho(i) = rho(i) + xi**2 * exp(-2.0d0 * q * (r/r0 - 1.0d0))
                            end if
                        end if
                    end if
                end if
            end do
        end do
        force_tbsma = 0.0d0
        do j = 1, size(pos, dim=1)-1
            do k = j+1, size(pos, dim=1)
                a0 = (var(j,1) + var(k,1)) / 2.0d0
                A = (var(j,2) + var(k,2)) / 2.0d0
                xi = (var(j,3) + var(k,3)) / 2.0d0
                p = (var(j,4) + var(k,4)) / 2.0d0
                q = (var(j,5) + var(k,5)) / 2.0d0
                r0 = (var(j,7) + var(k,7)) / 2.0d0
                rkj = pos(j,:) - pos(k,:)
                rkj = rkj - box * nint(rkj/box)
                r = norm2(rkj)
                rc2 = (sqrt(dble(N))+ sqrt(dble(N+1)))*r0/2
                if (r .lt. rc2) then
                    if (.not. use_cut_function) then
                        rep =  2.0d0 * A * p/r0 * exp(-p * (r/r0 - 1.0d0))
                        band = q * xi**2 / r0 * exp(-2.0d0 * q * (r/r0 - 1.0d0)) * (1.0d0/sqrt(rho(j)) + 1.0d0/sqrt(rho(k)))
                    else
                        rc1 = 0.95d0 * rc2
                        if (r .gt. rc1) then
                            E_rep = A * exp(-p*(rc1/r0 - 1.0d0))
                            dE_rep = -p * A / r0 * E_rep
                            dfcut_rep = derivative_cut_function(E_rep, dE_rep, r, rc1, rc2)
                            E_band = xi**2 * exp(-2.0d0 * q * (rc1/r0 - 1.0d0))
                            dE_band = -2.0d0 * q / r0 * E_band
                            dfcut_band = derivative_cut_function(E_band, dE_band, r, rc1, rc2)
                            rep = -dfcut_rep
                            !band = dfcut_band
                            !band = -0.5d0 * ( 1.0d0/sqrt(max(rho(j),1.0d-12)) + 1.0d0/sqrt(max(rho(k),1.0d-12)) ) * dfcut_band
                            band = -0.5d0 * ( 1.0d0/sqrt(rho(j)) + 1.0d0/sqrt(rho(k)) ) * dfcut_band
                        else
                            rep =  2.0d0 * A * p/r0 * exp(-p * (r/r0 - 1.0d0))
                            band = q * xi**2 / r0 * exp(-2.0d0 * q * (r/r0 - 1.0d0)) * (1.0d0/sqrt(rho(j)) + 1.0d0/sqrt(rho(k)))
                        end if
                    end if
                    force_tbsma(j,:) = force_tbsma(j,:) + (rep - band) * rkj/r
                    force_tbsma(k,:) = force_tbsma(k,:) - (rep - band) * rkj/r
                end if
            end do
        end do
    end function


    function force_num_tbsma(var, pos, box, use_cut_function)
        implicit none
        double precision, intent(inout)     ::  var(:,:), pos(:,:)
        logical, intent(in)                 ::  use_cut_function
        double precision, intent(in)        ::  box(3)
        double precision                    ::  force_num_tbsma(size(pos,dim=1), 3)
        double precision                    ::  E_p, E_n, h
        integer                             ::  i, k
        h = 1.0d-4
        do i = 1, size(pos, dim=1)
            do k = 1, 3
                pos(i,k) = pos(i,k) + h
                E_p = sum(tbsma(var, pos, box, use_cut_function))
                pos(i,k) = pos(i,k) - 2.0d0 * h
                E_n = sum(tbsma(var, pos, box, use_cut_function))
                pos(i,k) = pos(i,k) + h
                force_num_tbsma(i,k) = -(E_p - E_n) / (2.0d0 * h)
            end do
        end do
    end function


    function derivative_cut_function(E, dE, r, rc1, rc2)
        implicit none
        double precision, intent(in)        ::  E, dE, r, rc1, rc2
        double precision                    ::  derivative_cut_function
        double precision                    ::  b0, b1, b2, b3, b4, b5
        double precision                    ::  h
        h = rc2 - rc1
        b0 = E
        b1 = dE
        b2 = 0.0d0
        b3 = (-10.0d0*E - 6.0d0*dE*h) / (h**3)
        b4 = (15.0d0*E + 8.0d0*dE*h) / (h**4)
        b5 = (-6.0d0*E - 3.0d0*dE*h) / (h**5)
        derivative_cut_function = b1 + 2.0d0*b2*(r-rc1) + 3.0d0*b3*(r-rc1)**2 + 4.0d0*b4*(r-rc1)**3 + 5.0d0*b5*(r-rc1)**4
    end function


    subroutine test_force_tbsma(var, pos, box, use_cut_function)
        implicit none
        double precision, intent(inout)     ::  var(:,:), pos(:,:)
        double precision, intent(in)        ::  box(3)
        logical                             ::  use_cut_function
        double precision                    ::  F(size(pos, dim=1),3), F_num(size(pos, dim=1),3), precision
        integer                             ::  i
        F = force_tbsma(var, pos, box, use_cut_function)
        F_num = force_num_tbsma(var, pos, box, use_cut_function)
        do i = 1, size(pos, dim=1)
            if (abs(norm2(F(i,:))) .gt. abs(norm2(F_num(i,:)))) then ! To have the biggest in the denominator
                precision = norm2(F_num(i,:)) / norm2(F(i,:)) * 100.0d0
                if (precision .lt. 99.999d0) then
                    print*, "Atom", i
                    print*, "Analytical force :", F(i,:)
                    print*, "Numerical force :", F_num(i,:)
                    print*, precision, "%"
                    print*
                end if
            else
                precision = norm2(F(i,:)) / norm2(F_num(i,:)) * 100.0d0
                if (precision .lt. 99.999d0) then
                    print*, "Atom", i
                    print*, "Analytical force :", F(i,:)
                    print*, "Numerical force :", F_num(i,:)
                    print*, precision, "%"
                    print*
                end if
            end if
        end do
    end subroutine


end module
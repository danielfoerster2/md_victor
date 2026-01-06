module dynamic

    implicit none

    contains

    subroutine v_init
        use constants, only: mass, kb, pi
        use variables, only: n_atoms, vel, initial_temperature, typ, pos, box
        
        double precision                ::  v_cm(3), r_cm(3), e_kin, temperature_init, ang_momentum(3), r(3)
        double precision                ::  inertia_tensor(3, 3), inv_tensor(3, 3), omega(3), v_rotational(3), m
        double precision                ::  u1, u2, u3, u4
        integer                         ::  i_atom

        do i_atom = 1, n_atoms
            call random_number(u1)
            call random_number(u2)
            call random_number(u3)
            call random_number(u4)
            vel(1, i_atom) = sqrt(-2.0d0*log(u1)) * cos(2.0d0*pi*u2)
            vel(2, i_atom) = sqrt(-2.0d0*log(u1)) * sin(2.0d0*pi*u2)
            vel(3, i_atom) = sqrt(-2.0d0*log(u3)) * cos(2.0d0*pi*u4)
        enddo
        v_cm = 0.0d0
        r_cm = 0.0d0
        do i_atom = 1, n_atoms
            v_cm = v_cm + mass(typ(i_atom)) * vel(:, i_atom)
            r_cm = r_cm + mass(typ(i_atom)) * pos(:, i_atom)
        enddo
        m = sum(mass(typ(1:n_atoms)))
        v_cm = v_cm / m
        r_cm = r_cm / m

        ang_momentum = 0.0d0
        inertia_tensor = 0.0d0
        do i_atom = 1, n_atoms
            r(1) = pos(1, i_atom) - r_cm(1)
            r(2) = pos(2, i_atom) - r_cm(2)
            r(3) = pos(3, i_atom) - r_cm(3)
            r = r - box * nint(r/box)
            m = mass(typ(i_atom))

            ang_momentum(1) = ang_momentum(1) + (r(2) * m*vel(3, i_atom) - r(3) * m*vel(2, i_atom))
            ang_momentum(2) = ang_momentum(2) + (r(3) * m*vel(1, i_atom) - r(1) * m*vel(3, i_atom))
            ang_momentum(3) = ang_momentum(3) + (r(1) * m*vel(2, i_atom) - r(2) * m*vel(1, i_atom))

            inertia_tensor(1, 1) = inertia_tensor(1, 1) + m * (r(2)**2 + r(3)**2)
            inertia_tensor(2, 2) = inertia_tensor(2, 2) + m * (r(1)**2 + r(3)**2)
            inertia_tensor(3, 3) = inertia_tensor(3, 3) + m * (r(1)**2 + r(2)**2)
            inertia_tensor(1, 2) = inertia_tensor(1, 2) - m * r(1) * r(2)
            inertia_tensor(1, 3) = inertia_tensor(1, 3) - m * r(1) * r(3)
            inertia_tensor(2, 3) = inertia_tensor(2, 3) - m * r(2) * r(3)
        enddo

        inertia_tensor(2, 1) = inertia_tensor(1, 2)
        inertia_tensor(3, 1) = inertia_tensor(1, 3)
        inertia_tensor(3, 2) = inertia_tensor(2, 3)

        inv_tensor(1, :) = (/ inertia_tensor(2,2)*inertia_tensor(3,3) - inertia_tensor(2,3)**2,&
                        & inertia_tensor(2,3)*inertia_tensor(1,3) - inertia_tensor(1,2)*inertia_tensor(3,3), &
                        & inertia_tensor(1,2)*inertia_tensor(2,3) - inertia_tensor(2,2)*inertia_tensor(1,3) /) 
        inv_tensor(2, :) = (/ inertia_tensor(1,3)*inertia_tensor(2,3) - inertia_tensor(1,2)*inertia_tensor(3,3), &
                        & inertia_tensor(1,1)*inertia_tensor(3,3) - inertia_tensor(1,3)**2, &
                        & inertia_tensor(1,2)*inertia_tensor(1,3) - inertia_tensor(1,1)*inertia_tensor(2,3) /) 
        inv_tensor(3, :) = (/ inertia_tensor(1,2)*inertia_tensor(2,3) - inertia_tensor(1,3)*inertia_tensor(2,2), &
                        & inertia_tensor(1,2)*inertia_tensor(1,3) -  inertia_tensor(1,1)*inertia_tensor(2,3), &
                        & inertia_tensor(1,1)*inertia_tensor(2,2) - inertia_tensor(1,2)**2 /) 
        inv_tensor = 1.0d0 / (inertia_tensor(1,1)*inertia_tensor(2,2)*inertia_tensor(3,3) - inertia_tensor(1,1) &
                    & *inertia_tensor(2,3)**2 - inertia_tensor(3,3)*inertia_tensor(1,2)**2 - inertia_tensor(2,2) &
                    & *inertia_tensor(1,3)**2 + 2*inertia_tensor(1,2)*inertia_tensor(1,3)*inertia_tensor(2,3)) * inv_tensor

        omega = matmul(inv_tensor, ang_momentum)

        e_kin = 0.0d0
        do i_atom = 1, n_atoms
            r(1) = pos(1, i_atom) - r_cm(1)
            r(2) = pos(2, i_atom) - r_cm(2)
            r(3) = pos(3, i_atom) - r_cm(3)
            v_rotational(1) = omega(2)*r(3) - omega(3)*r(2)
            v_rotational(2) = omega(3)*r(1) - omega(1)*r(3)
            v_rotational(3) = omega(1)*r(2) - omega(2)*r(1)
            vel(:, i_atom) = vel(:, i_atom) - v_cm - v_rotational
            e_kin = e_kin + 0.5d0 * mass(typ(i_atom)) * sum(vel(:, i_atom)**2)
        enddo
        temperature_init = 2.0d0 * e_kin / ((3*n_atoms-6) * kb)
        do i_atom = 1, n_atoms
            vel(:, i_atom) = vel(:, i_atom) * sqrt(kb*initial_temperature / mass(typ(i_atom)))
        enddo
        
    endsubroutine


    subroutine verlet_velocity
        use constants, only: mass, dt
        use variables, only: n_atoms, typ, force, vel, pos, ekin
        use potential, only: force_tbsma, epot_tbsma

        integer                             ::  i_atom
        double precision                    ::  a(3)
        logical, save                       ::  first_call = .true.

        if (first_call) then
            call force_tbsma
            first_call = .false.
        endif

        do i_atom = 1, n_atoms
            a = force(:, i_atom)/mass(typ(i_atom))
            pos(:, i_atom) = pos(:, i_atom) + vel(:, i_atom)*dt + a*dt**2/2.0d0
            vel(:, i_atom) = vel(:, i_atom) + a*dt/2.0d0
        enddo

        call force_tbsma
        ekin = 0.0d0
        do i_atom = 1, n_atoms
            vel(:, i_atom) = vel(:, i_atom) + force(:, i_atom)/mass(typ(i_atom)) * dt/2.0d0
            ekin = ekin + 0.5d0 * mass(typ(i_atom)) * sum(vel(:, i_atom)**2)
        enddo

    endsubroutine


    subroutine nose_hoover

        use constants, only: mass, kb, dt, n_atoms_max, omega
        use variables, only: n_atoms, typ, force, vel, target_temperature, pos, Qth, ekin, eth
        use potential, only: force_tbsma, epot_tbsma

        implicit none

        double precision                        ::  g_new, v_thermo_new, ekin_new
        double precision                        ::  accel_new(3, n_atoms), vk(3, n_atoms)
        integer                                 ::  k, i_atom
        double precision, save                  ::  accel(3, n_atoms_max), g
        double precision, save                  ::  v_thermo = 0.0d0, x_thermo = 0.0d0, v_thermo_old
        logical, save                           ::  first_call = .true.

        if (first_call) then
            ekin = 0.0d0
            do i_atom = 1, n_atoms
                ekin = ekin + 0.5d0 * mass(typ(i_atom)) * sum(vel(:, i_atom)**2)
            enddo
            Qth = (3*n_atoms-3)*kb*target_temperature / (omega**2)
            g = (2.0d0*ekin - (3*n_atoms-6)*kb*target_temperature) / Qth
            v_thermo_old = -dt * g

            call force_tbsma
            do i_atom = 1, n_atoms
                accel(:, i_atom) = force(:, i_atom)/mass(typ(i_atom))
            enddo
            first_call = .false.
        endif

        pos(:, 1:n_atoms) = pos(:, 1:n_atoms) + vel(:, 1:n_atoms)*dt + (accel(:, 1:n_atoms)-vel(:, 1:n_atoms)*v_thermo)*0.5d0*dt**2

        x_thermo = x_thermo + v_thermo*dt + 0.5d0*g*dt**2

        call force_tbsma
        do i_atom = 1, n_atoms
            accel_new(:, i_atom) = force(:, i_atom)/mass(typ(i_atom))
        enddo

        v_thermo_new = v_thermo_old + 2.0d0*g*dt

        do k = 1, 50
            vk = 1.0d0 / (1.0d0 + v_thermo_new*0.5d0*dt)*(vel(:, 1:n_atoms) + (accel(:, 1:n_atoms) + accel_new &
                    &- vel(:, 1:n_atoms)*v_thermo)*0.5d0*dt)

            ekin_new = 0.0d0
            do i_atom = 1, n_atoms
                ekin_new = ekin_new + 0.5d0 * mass(typ(i_atom)) * sum(vk(:, i_atom)**2)
            enddo

            g_new = (2.0d0* ekin_new - (3*n_atoms-6)*kb*target_temperature) / Qth
            v_thermo_new = v_thermo + (g + g_new)*0.5d0*dt
            if((abs((ekin-ekin_new)/ekin_new) .lt. 1.0d-10).and.(k.ne.1)) exit
            ekin = ekin_new
        enddo

        v_thermo_old = v_thermo
        v_thermo = v_thermo_new
        vel(:, 1:n_atoms) = vk
        g = g_new
        accel(:, 1:n_atoms) = accel_new

        eth = Qth*v_thermo**2 / 2.0d0 + (3*n_atoms-6)*kb*target_temperature*x_thermo
    endsubroutine



endmodule
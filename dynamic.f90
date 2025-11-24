module dynamic

    implicit none

    contains

    subroutine v_init
        use constants, only: mass, kb
        use variables, only: n_atoms, vel, target_temperature, typ
        double precision                ::  v_cm(3), e_kin, temperature_init
        integer                         ::  i_atom

        call random_number(vel(:, 1:n_atoms))
        v_cm = 0.0d0
        do i_atom = 1, n_atoms
            v_cm = v_cm + mass(typ(i_atom)) * vel(:, i_atom)
        enddo
        v_cm = v_cm / sum(mass(typ(1:n_atoms)))
        e_kin = 0.0d0
        do i_atom = 1, n_atoms
            vel(:, i_atom) = vel(:, i_atom) - v_cm
            e_kin = e_kin + 0.5d0 * mass(typ(i_atom)) * sum(vel(:, i_atom)**2)
        enddo
        temperature_init = 2.0d0 * e_kin / ((3*n_atoms-3) * kb)
        vel(:, 1:n_atoms) = vel(:, 1:n_atoms) * sqrt(target_temperature / temperature_init)
    endsubroutine


    subroutine verlet_velocity
        use constants, only: mass, dt
        use variables, only: n_atoms, typ, force, vel, pos
        use potential, only: force_tbsma

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

        do i_atom = 1, n_atoms
            vel(:, i_atom) = vel(:, i_atom) + force(:, i_atom)/mass(typ(i_atom)) * dt/2.0d0
        enddo

    endsubroutine


    subroutine nose_hoover

        use constants, only: mass, kb, dt
        use variables, only: n_atoms, typ, force, vel, target_temperature, pos, epot
        use potential, only: force_tbsma

        implicit none

        double precision                        ::  g_new, g, v_thermo_new, v_thermo_old, ekin, ekin_new
        double precision                        ::  accel(3, n_atoms), accel_new(3, n_atoms), vk(3, n_atoms), eth, Qth
        integer                                 ::  k, i_atom
        double precision, save                  ::  omega=0.1d0
        double precision, save                  ::  v_thermo = 0.0d0, x_thermo = 0.0d0
        logical, save                           ::  first_call = .true.

        if (first_call) then
            call force_tbsma
            first_call = .false.
        endif

        Qth = (3*n_atoms-3)*kb*target_temperature / omega**2
        ekin = 0.0d0
        do i_atom = 1, n_atoms
            ekin = ekin + 0.5d0 * mass(typ(i_atom)) * sum(vel(:, i_atom)**2)
        enddo
        g = (2.0d0*ekin - (3*n_atoms-3)*kb*target_temperature) / Qth

        v_thermo_old = -dt * g

        do i_atom = 1, n_atoms
            accel(:, i_atom) = force(:, i_atom)/mass(typ(i_atom))
        enddo
        pos(:, 1:n_atoms) = pos(:, 1:n_atoms) + vel(:, 1:n_atoms)*dt + (accel(:, 1:n_atoms)-vel(:, 1:n_atoms)*v_thermo)*0.5d0*dt**2

        call force_tbsma

        do i_atom = 1, n_atoms
            accel_new(:, i_atom) = force(:, i_atom)/mass(typ(i_atom))
        enddo

        x_thermo = x_thermo + v_thermo*dt + 0.5d0*g*dt**2
        v_thermo_new = v_thermo_old + 2.0d0*g*dt

        do k = 1, 50
            vk = 1.0d0 / (1.0d0 + v_thermo_new*0.5d0*dt)*(vel(:, 1:n_atoms) + (accel + accel_new - vel(:, 1:n_atoms)*v_thermo)&
                    &*0.5d0*dt)

            ekin_new = 0.0d0
            do i_atom = 1, n_atoms
                ekin_new = ekin_new + 0.5d0 * mass(typ(i_atom)) * sum(vk(:, i_atom)**2)
            enddo

            g_new = (2.0d0* ekin_new - (3*n_atoms-3)*kb*target_temperature) / Qth
            v_thermo_new = v_thermo + (g + g_new)*0.5d0*dt
            if((abs((ekin-ekin_new)/ekin_new) .lt. 1.0d-10).and.(k.ne.1)) exit
            ekin = ekin_new
        enddo
        g = g_new
        v_thermo_old = v_thermo
        v_thermo = v_thermo_new
        vel(:, 1:n_atoms) = vk
        accel = accel_new
        eth = Qth*v_thermo**2 / 2.0d0 + (3*n_atoms-3)*kb*target_temperature*x_thermo
        write(*,*) 'Energy:', eth+ekin + sum(epot(1:n_atoms))
    endsubroutine



end module
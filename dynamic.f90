module dynamic

    implicit none

    contains

    subroutine v_init
        use constants, only: mass, kb
        use variables, only: n_atoms, vel, target_temperature, typ
        double precision                ::  v_cm(3), e_kin, temperature_init
        integer                         ::  i_atom

        call random_number(vel(:, 1:n_atoms))
        vel(:, 1:n_atoms) = vel(:, 1:n_atoms) - 0.5d0
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
        logical, save                       ::  first_call = .true.

        if (first_call) then
            call force_tbsma
            first_call = .false.
        endif

        do i_atom = 1, n_atoms
            pos(:, i_atom) = pos(:, i_atom) + vel(:, i_atom)*dt + force(:, i_atom)/mass(typ(i_atom)) * dt**2/2.0d0
            vel(:, i_atom) = vel(:, i_atom) + force(:, i_atom)/mass(typ(i_atom)) * dt/2.0d0
        enddo

        call force_tbsma

        do i_atom = 1, n_atoms
            vel(:, i_atom) = vel(:, i_atom) + force(:, i_atom)/mass(typ(i_atom)) * dt/2.0d0
        enddo
    endsubroutine verlet_velocity


    subroutine nose_hoover(var, pos, v, accel, x_thermo, v_thermo, target_temp, Qth, dt, Nf, ekin, box)
        use constants, only: mass, kb
        use variables, only: n_atoms, typ, force
        use potential, only: force_tbsma
        !========================================
        ! Implement the Nose-Hoover thermostat
        ! 
        ! Parameters :
        ! ------------
        ! pos : 2D array
        !       Coordinates of atoms
        ! v : 2D array
        !       Velocities of atoms
        ! accel : 2D array
        !       Acceleration of atoms
        ! x_thermo : double precision
        !       Coordinate from the thermostat
        ! v_thermo : double precision
        !       velocity from the thermostat
        ! target_temp : double precision
        !       Target temperature (K)
        ! Qth : double precision
        !       "Mass" of the thermostat
        ! dt : double precision
        !       Step time
        ! Nf : integer
        !       Degree of freedom
        ! ekin : Double precision
        !       Kinetic energy of the system
        ! 
        ! Returns :
        ! ---------
        ! None
        !========================================
        implicit none
        double precision, intent(inout)         ::  pos(:,:), v(:,:), accel(:,:)
        double precision, intent(inout)         ::  v_thermo, x_thermo, ekin
        double precision, intent(in)            ::  var(:,:), target_temp, Qth, dt, box(3)
        integer, intent(in)                     ::  Nf
        double precision                        ::  g_new, g, v_thermo_new, v_thermo_old, ekin_new
        double precision                        ::  v_k(size(pos, dim=1), 3), accel_new(size(accel, dim=1), 3)
        integer                                 ::  k, i_atom
        ekin_new = 0.5d0*sum(var(:,6)*sum(v**2, dim=2))
        g = (sum(var(:,6)*sum(v**2, dim=2)) - Nf*kb*target_temp) / Qth
        v_thermo_old = -dt * g
        pos = pos + v*dt + (accel - v*v_thermo)*0.5d0*dt**2
        call force_tbsma
        do i_atom = 1, n_atoms
            accel_new(:, i_atom) = force(:, i_atom)/mass(typ(i_atom))
        end do
        x_thermo = x_thermo + v_thermo*dt + 0.5d0*g*dt**2
        v_thermo_new = v_thermo_old + 2.0d0*g*dt
        do k = 1, 50
            v_k = 1.0d0 / (1.0d0 + v_thermo_new*0.5d0*dt)*(v + (accel + accel_new - v*v_thermo)*0.5d0*dt)
            ekin_new = 0.5d0*sum(var(:,6)*sum(v_k**2, dim=2))
            g_new = (sum(var(:,6)*sum(v_k**2, dim=2)) - Nf*kb*target_temp) / Qth
            v_thermo_new = v_thermo + (g + g_new)*0.5d0*dt
            if((abs((ekin-ekin_new)/ekin_new) .lt. 1.0d-10).and.(k.ne.1)) exit
            ekin = ekin_new
        end do
        g = g_new
        v_thermo_old = v_thermo
        v_thermo = v_thermo_new
        v = v_k
        accel = accel_new
    end subroutine nose_hoover



end module
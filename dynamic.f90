module dynamic

    use constants
    implicit none
    contains




    function v_init(var, n_atoms, init_temp, Nf)
        !========================================
        ! Initialize the velocity of atoms
        ! 
        ! Parameters :
        ! ------------
        ! var : double precision 2D array
        !       Properties of each atoms
        ! n_atoms : integer
        !       Number of atoms
        ! init_temp : double precision
        !       Initial temperature
        ! Nf : integer
        !       Degree of freedom
        ! 
        ! Returns :
        ! ---------
        ! double precision 2D array
        !           Initial velocity of each atoms
        !           Dim 1 : Number of atoms
        !           Dim 2 : vx, vy, vz 
        !========================================
        implicit none
        integer, intent(in)             ::  n_atoms, Nf
        double precision, intent(in)    ::  init_temp, var(:,:)
        double precision                ::  v_init(n_atoms,3), v_cm(3), ekin, temperature
        integer                         ::  i
        call random_number(v_init)
        v_init = v_init - 0.5d0
        v_cm = sum(v_init, dim=1) / n_atoms
        do i = 1, n_atoms
            v_init(i,:) = v_init(i,:) - v_cm
        end do
        ekin = 0.5d0 * sum(var(:,6)*sum(v_init**2, dim=2))
        temperature = 2.0d0 * ekin / (Nf * kb)
        v_init = v_init * sqrt(init_temp / temperature)
    end function v_init


    subroutine verlet_velocity(var, pos, v, accel, dt, ekin, box)
        use constants, only: mass
        use variables, only: n_atoms, typ, force
        use potential, only: force_tbsma

        !========================================
        ! Implement the Verlet velocity scheme
        ! 
        ! Parameters :
        ! ------------
        ! var : double precision 2D array
        !       Properties of each atoms
        ! pos : 2D array
        !       Coordinates of atoms
        ! v : 2D array
        !       Velocities of atoms
        ! accel : 2D array
        !       Acceleration of atoms
        ! dt : double precision
        !       Step time
        ! ekin : double precision
        !       Kinetic energy
        !
        ! Returns :
        ! ---------
        ! None
        !========================================
        implicit none
        double precision, intent(inout)     ::  var(:,:), pos(:,:), v(:,:), accel(:,:)
        double precision, intent(inout)     ::  ekin
        double precision, intent(in)        ::  dt, box(3)
        integer                             ::  i_atom
        pos = pos + v*dt + accel*dt**2/ 2.0d0
        v = v + accel * dt/2.0d0
        call force_tbsma
        do i_atom = 1, n_atoms
            accel(:, i_atom) = force(:, i_atom)/mass(typ(i_atom))
        end do
        v = v + accel * dt/2.0d0
        ekin = 0.5d0*sum(var(:,6)*sum(v**2, dim=2))
    end subroutine verlet_velocity


    subroutine nose_hoover(var, pos, v, accel, x_thermo, v_thermo, target_temp, Qth, dt, Nf, ekin, box)
        use constants, only: mass
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
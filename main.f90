program main
    use variables, only: target_temperature, time, n_atoms, epot
    use constants, only: dt, n_steps, use_thermostat, mass
    use potential, only: init_tbsma, test_force_tbsma, epot_tbsma
    use file_io, only: read_xyz, save_xyz, save_data
    use dynamic, only: v_init, verlet_velocity, nose_hoover
    use neighbor, only: init_neigh_list

    implicit none
    integer                         ::  i_step, i_atom

    target_temperature = 300.0d0
    call read_xyz
    call init_tbsma
    call init_neigh_list

    ! Problem with Energy
    call epot_tbsma
    do i_atom = 1, n_atoms
        write(*,*) epot(i_atom)
    enddo

    call v_init
    time = 0.0d0

    do i_step = 1, n_steps
        if (.not. use_thermostat) then
            call verlet_velocity
        else
            call nose_hoover
        endif
        call epot_tbsma
!        call save_data
!        call save_xyz
        time = i_step*dt
        call init_neigh_list
        if (mod(i_step, 1000) .eq. 0) then
            write(*,*) 'Step:', i_step
        endif
    enddo
    call save_xyz

endprogram
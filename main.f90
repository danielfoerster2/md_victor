program main
    use variables, only: target_temperature, time, pos_old
    use constants, only: dt, n_steps, use_thermostat
    use potential, only: init_tbsma, test_force_tbsma, epot_tbsma
    use file_io, only: read_xyz, save_xyz, save_data
    use dynamic, only: v_init, verlet_velocity, nose_hoover
    use neighbor, only: init_neigh_list, update_neigh_list

    implicit none
    integer                         ::  i_step

    target_temperature = 300.0d0
    call read_xyz
    call init_tbsma
    pos_old = 0.0d0
    call init_neigh_list
    call v_init
    time = 0.0d0

    do i_step = 1, n_steps
        if (.not. use_thermostat) then
            call verlet_velocity
        else
            call nose_hoover
        endif

        call save_data
        call save_xyz
        time = i_step*dt
        call update_neigh_list
    enddo

endprogram
program main
    use variables, only: target_temperature, time
    use constants, only: dt, n_steps, use_thermostat
    use potential, only: init_tbsma, test_force_tbsma
    use file_io, only: read_xyz, save_xyz, save_data
    use dynamic, only: v_init, verlet_velocity, nose_hoover

    implicit none
    integer                         ::  i_step

    target_temperature = 300.0d0
    call read_xyz
    call init_tbsma
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
    enddo

endprogram
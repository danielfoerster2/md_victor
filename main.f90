program main
    use variables, only: target_temperature, time, initial_temperature, final_temperature, n_steps
    use constants, only: dt, use_thermostat, mass, kb
    use potential, only: init_tbsma, test_force_tbsma, epot_tbsma
    use file_io, only: read_xyz, save_xyz, save_data, test_velocity_distribution
    use dynamic, only: v_init, verlet_velocity, nose_hoover
    use neighbor, only: init_neigh_list

    implicit none

    integer                         ::  i_step
    character(len=255)              ::  env_var

    call get_environment_variable("n_steps", env_var)
    read(env_var, *) n_steps
    call get_environment_variable("initial_temperature", env_var)
    read(env_var, *) initial_temperature
    call get_environment_variable("final_temperature", env_var)
    read(env_var, *) final_temperature
    
    time = 0.0d0
    target_temperature = (final_temperature - initial_temperature) / (n_steps*dt) * time + initial_temperature
    
    call read_xyz
    call init_tbsma
    call init_neigh_list
    ! call epot_tbsma
    call v_init
    !call test_velocity_distribution
    !call save_data
    write(*,*) "Initialization done"
    do i_step = 1, n_steps
        !time = i_step*dt
        if (.not. use_thermostat) then
            call verlet_velocity
        else
            call nose_hoover
        endif
        call init_neigh_list
        !call epot_tbsma
    
        if (mod(i_step, 100000) .eq. 0) then
            write(*,*) 'Step:', i_step, '/', n_steps
            !call save_xyz
        endif
        !call save_data
        target_temperature = (final_temperature - initial_temperature) / n_steps * dble(i_step) + initial_temperature
    enddo
    write(*,*) "MD simulation done"
    !call test_velocity_distribution
    call epot_tbsma
    call save_xyz
endprogram
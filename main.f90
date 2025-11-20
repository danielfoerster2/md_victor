program main
    use variables
    use constants
    use potential
    use file_io
    use dynamic
    implicit none

    !!! ============ Variables ============ !!!
    
    double precision                ::  x_thermo=0.0d0, v_thermo=0.0d0
    double precision                ::  ekin, eth, Qth
    double precision                ::  omega=0.1d0
    double precision                ::  p_Ni, p_Cu, p_Rh, p_Pd, p_Ag, p_Ir, p_Pt, p_Au
    double precision                ::  p_Al, p_Pb, p_Ti, p_Zr, p_Co, p_Cd, p_Zn, p_Mg
    double precision                ::  percent(n_atoms_max), a0, e_kin
    integer                         ::  i, j, Nf, i_atom
    logical                         ::  random=.true., use_thermostat=.false., init_xyz 
    character(len=20)               ::  filename
    character(len=2), allocatable   ::  list(:)



    target_temperature = 300.0d0
    call read_xyz
    call init_tbsma

    pos(1, 2) = 1.0d0
    do i = 1, 800
        pos(1, 2) = pos(1, 2) + 0.01d0
        call epot_tbsma
        write(55,*) pos(1, 2), sum(epot(1:n_atoms))
    enddo
    stop
    call v_init

    ! ! Initialization Nose-Hoover thermostat
    ! Qth = Nf*kb*target_temp / omega**2
    ! eth = Qth*v_thermo**2 / 2.0d0 + Nf*kb*target_temp*x_thermo

    ! call save_data(n_atoms, 0.0d0, sum(epot), ekin, eth, box)

    ! print*, "Time :", 0.0d0, "Ekin =", ekin, "Epot =", sum(epot), "Eth =", eth, "E =", sum(epot) + ekin + eth, "T =", temperature

    !!! ==================================== !!!
    do i = 1, n_steps

        if (.not. use_thermostat) then
            call verlet_velocity
            call test_force_tbsma
        !else
         !   call nose_hoover(var, pos, v, accel, x_thermo, v_thermo, target_temp, Qth, dt, Nf, ekin, box)
        endif
        !epot = tbsma(var, pos, box, .true.)
        !temperature = 2.0d0 * ekin / (Nf * kb)
  !      eth = Qth*v_thermo**2 / 2.0d0 + Nf*kb*target_temp*x_thermo
 !       print*, "Time :", i*dt, "Ekin =", ekin, "Epot =", sum(epot), "Eth =", eth, "E =", sum(epot) + ekin + eth, "T =", temperature
        
!        call save_data(n_atoms, i*dt, sum(epot), ekin, eth, box)

        call epot_tbsma
        e_kin = 0.0d0
        do i_atom = 1, n_atoms
            e_kin = e_kin + 0.5d0 * mass(typ(i_atom)) * sum(vel(:, i_atom)**2)
        enddo
        write(*,*) sum(epot(1:n_atoms)) + e_kin, sum(epot(1:n_atoms)) , e_kin
        call save_xyz
    enddo

endprogram
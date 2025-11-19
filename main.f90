program main
    
    use potential
    use force
    use file_io
    use dynamic
    implicit none

    !!! ============ Variables ============ !!!
    
    double precision                ::  basis(4,3) ! FCC unit cell
    double precision, allocatable   ::  pos(:,:), epot(:), accel(:,:), v(:,:), var(:,:)
    double precision                ::  x_thermo=0.0d0, v_thermo=0.0d0
    double precision                ::  ekin, eth, Qth
    double precision                ::  dt=1.0d0
    double precision                ::  init_temp=100.0d0, temperature, target_temp=500.0d0
    double precision                ::  omega=0.1d0
    double precision                ::  p_Ni, p_Cu, p_Rh, p_Pd, p_Ag, p_Ir, p_Pt, p_Au
    double precision                ::  p_Al, p_Pb, p_Ti, p_Zr, p_Co, p_Cd, p_Zn, p_Mg
    double precision                ::  percent(n_atoms_max), a0
    integer                         ::  n_atoms, i, n_steps=1000, j, Nf
    logical                         ::  random=.true., use_thermostat=.true., init_xyz 
    character(len=20)               ::  filename
    character(len=2), allocatable   ::  list(:)
    double precision                ::  box(3)


    call read_xyz

    call init_potential
    call tbsma
    call save_xyz
    
    !call test_force_tbsma(var, pos, box, .true.)

    stop
    accel = acceleration(var, pos, box)

    !!! ============ Initialization ============ !!!

    allocate(v(n_atoms, 3))
    Nf = 3*n_atoms - 3
    v = v_init(var, n_atoms,init_temp, Nf)
    ekin = 0.5d0 * sum(var(:,6)*sum(v**2, dim=2))
    temperature = 2.0d0 * ekin / (Nf * kb)
    print*, "Target temperature :", init_temp, "K"
    print*, "Temperature of the system :", temperature, "K"

    ! Initialization Nose-Hoover thermostat
    Qth = Nf*kb*target_temp / omega**2
    eth = Qth*v_thermo**2 / 2.0d0 + Nf*kb*target_temp*x_thermo

    call save_data(n_atoms, 0.0d0, sum(epot), ekin, eth, box)

    print*, "Time :", 0.0d0, "Ekin =", ekin, "Epot =", sum(epot), "Eth =", eth, "E =", sum(epot) + ekin + eth, "T =", temperature

    !!! ==================================== !!!
    do i = 1, n_steps

        if (.not. use_thermostat) then
            call verlet_velocity(var, pos, v, accel, dt, ekin, box)
        else
            call nose_hoover(var, pos, v, accel, x_thermo, v_thermo, target_temp, Qth, dt, Nf, ekin, box)
        end if
        !epot = tbsma(var, pos, box, .true.)
        temperature = 2.0d0 * ekin / (Nf * kb)
        eth = Qth*v_thermo**2 / 2.0d0 + Nf*kb*target_temp*x_thermo
        print*, "Time :", i*dt, "Ekin =", ekin, "Epot =", sum(epot), "Eth =", eth, "E =", sum(epot) + ekin + eth, "T =", temperature
        
        call save_data(n_atoms, i*dt, sum(epot), ekin, eth, box)
        call save_xyz
    end do

end program
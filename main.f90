program main
    
    use init
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
    double precision                ::  percent(size(atoms)), a0
    integer                         ::  n_atoms, i, n_steps=1000, j, Nf
    logical                         ::  random=.true., use_thermostat=.true., init_xyz 
    character(len=20)               ::  filename
    character(len=2), allocatable   ::  list(:)
    double precision                ::  box(3)

    !!! ==================================== !!!

    init_xyz = .false.
    if (init_xyz) then
        call read_xyz("XYZ/test0.exyz", n_atoms, list, pos, box)
        allocate(var(n_atoms, 7))
        var = var_array(list)
    else
        p_Ni = 100.0d0 ! Percentage of Nickel
        p_Cu = 0.0d0 ! Percentage of Copper
        p_Rh = 0.0d0 ! Percentage of Rhodium
        p_Pd = 0.0d0 ! Percentage of Palladium
        p_Ag = 0.0d0 ! Percentage of Silver
        p_Ir = 0.0d0 ! Percentage of Iridium
        p_Pt = 0.0d0 ! Percentage of Platinium
        p_Au = 0.0d0 ! Percentage of Gold
        p_Al = 0.0d0 ! Percentage of Aluminium
        p_Pb = 0.0d0 ! Percentage of Plomb
        p_Ti = 0.0d0 ! Percentage of Titanium
        p_Zr = 0.0d0 ! Percentage of 
        p_Co = 0.0d0 ! Percentage of Cobalt
        p_Cd = 0.0d0 ! Percentage of 
        p_Zn = 0.0d0 ! Percentage of Zinc
        p_Mg = 0.0d0 ! Percentage of Magnesium
        percent =  (/ p_Ni, p_Cu, p_Rh, p_Pd, p_Ag, p_Ir, p_Pt, p_Au, &
        p_Al, p_Pb, p_Ti, p_Zr, p_Co, p_Cd, p_Zn, p_Mg /)
    
        if (abs(sum(percent) - 100.0d0) > 1.0d-6) then
            print*, "Warning: Percentages do not sum to 100!"
            stop
        end if
        if (random) then
            n_atoms = 200 ! Good until 1200 for Pb (not test over)
            box = (/ 1.0d0, 1.0d0, 1.0d0 /) * 10.0d0 * n_cells
            list = random_list_atom(n_atoms, atoms, percent)
            allocate(var(n_atoms, 7))
            var = var_array(list)
            pos = generate_random(var, n_atoms, box)
        else ! CFC structure (not read xyz and not random position)
            ! Not implemented for many types of atoms in the system
            n_atoms = 4*n_cells**3
            box = (/ 1.0d0, 1.0d0, 1.0d0 /) * 10.0d0 * n_cells
            list = random_list_atom(n_atoms, atoms, percent)
            allocate(var(n_atoms, 7))
            var = var_array(list)
            a0 = var(1,1)
            basis(1,:) = [0.0d0, 0.0d0, 0.0d0]
            basis(2,:) = [a0/2.0d0, a0/2.0d0, 0.0d0]
            basis(3,:) = [a0/2.0d0, 0.0d0, a0/2.0d0]
            basis(4,:) = [0.0d0, a0/2.0d0, a0/2.0d0]
            pos = generate_lattice(var, basis, n_atoms)
        end if
    end if
    
    print*, "Positions initialized!"
    print*, "Number of atoms :", n_atoms
    
    !!! ==================================== !!!

    allocate(epot(n_atoms))
    epot = tbsma(var, pos, box, .true.)
    !call test_tbsma(pos, epot)
    if (.not. init_xyz) then
        call save_xyz(list, n_atoms, pos, epot, box)
    end if
    
    !call test_force_tbsma(var, pos, box, .true.)

    allocate(accel(n_atoms,3))
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
        epot = tbsma(var, pos, box, .true.)
        temperature = 2.0d0 * ekin / (Nf * kb)
        eth = Qth*v_thermo**2 / 2.0d0 + Nf*kb*target_temp*x_thermo
        print*, "Time :", i*dt, "Ekin =", ekin, "Epot =", sum(epot), "Eth =", eth, "E =", sum(epot) + ekin + eth, "T =", temperature
        
        call save_data(n_atoms, i*dt, sum(epot), ekin, eth, box)
        call save_xyz(list, n_atoms, pos, epot, box)
    end do

end program
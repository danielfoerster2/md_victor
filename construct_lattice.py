#!/usr/bin/python3

n_cells = 5
a = 4.085
f = open('input.xyz', 'w')
print(n_cells**3*4, file=f)
print(f'Lattice="{n_cells*a} 0.0 0.0 0.0 {n_cells*a} 0.0 0.0 0.0 {n_cells*a}" Properties=species:S:1:pos:R:3', file=f)
for ix in range(n_cells):
    for iy in range(n_cells):
        for iz in range(n_cells):
            print("Ag", ix*a, iy*a, iz*a, file=f)
            print("Ag", (ix+0.5)*a, (iy+0.5)*a, iz*a, file=f)
            print("Ag", ix*a, (iy+0.5)*a, (iz+0.5)*a, file=f)
            print("Ag", (ix+0.5)*a, iy*a, (iz+0.5)*a, file=f)
f.close()

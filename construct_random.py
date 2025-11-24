#!/usr/bin/python3

import numpy as np

n_atoms = 50

f = open('input.xyz', 'w')

box_length = 10
xyz = np.zeros((n_atoms, 3))

print(n_atoms, file=f)
print(f'Lattice="{box_length} 0.0 0.0 0.0 {box_length} 0.0 0.0 0.0 {box_length}" Properties=species:S:1:pos:R:3', file=f)
for i in range(n_atoms):
    not_ok = True
    while not_ok:
        coord = np.random.uniform(0, box_length, 3)
        not_ok = False
        for j in range(i):
            d = coord - xyz[j, :]
            d = d - box_length * np.round(d / box_length)
            if np.linalg.norm(d) < 2.0:
                not_ok = True
                break
    xyz[i, :] = coord
    print("Ni", *coord, file=f)

f.close()


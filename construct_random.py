#!/usr/bin/python3

import numpy as np

n_atoms = 50

f = open('input.xyz', 'w')

box_length = 10

print(n_atoms, file=f)
print(f'Lattice="{box_length} 0.0 0.0 0.0 {box_length} 0.0 0.0 0.0 {box_length}" Properties=species:S:1:pos:R:3', file=f)
for i in range(n_atoms):
    print("Ni", *np.random.uniform(0, box_length, 3), file=f)
f.close()


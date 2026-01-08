#!/usr/bin/python3

import os

import numpy as np

n_atoms = int(os.environ.get("n_atoms"))
composition = float(os.environ.get("composition"))
density = float(os.environ.get("density"))
atom_typ1 = os.environ.get("atom_typ1")
atom_typ2 = os.environ.get("atom_typ2")
file = os.environ.get("file_input")

f = open(file, 'w')

box_length = (n_atoms/density)**(1/3)
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
    if i < composition*n_atoms:
        print(atom_typ1, *coord, file=f)
    else:
        print(atom_typ2, *coord, file=f)

f.close()

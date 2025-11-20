#!/usr/bin/python3

f = open('input.xyz', 'w')

box_length = 100
n_atoms = 2
print(n_atoms, file=f)
print(f'Lattice="{box_length} 0.0 0.0 0.0 {box_length} 0.0 0.0 0.0 {box_length}" Properties=species:S:1:pos:R:3', file=f)
print("Ni", 0, 0, 0, file=f)
print("Ni", 2, 0, 0, file=f)
f.close()

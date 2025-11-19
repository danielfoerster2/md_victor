#!/usr/bin/python3

from ase.cluster.icosahedron import Icosahedron
from ase.cluster.decahedron import Decahedron
from ase.cluster.octahedron import Octahedron
import numpy as np
import sys

particle_type = sys.argv[1]
n = int(sys.argv[2])

symbol = "Pd"
a0 = 3.7

atoms = []
if particle_type == "deca" and n>0 :
    p = n
    q = 1
    r = 0
    atoms = Decahedron(symbol, p, q, r, latticeconstant=a0)
if particle_type == "ico" and n>0 :
    noshells = n
    atoms = Icosahedron(symbol, noshells, latticeconstant=a0)
if particle_type == "toct":
    cutoff = n
    length = 3 * n + 1
    atoms = Octahedron(symbol, length, cutoff, latticeconstant=a0, alloy=False)
if particle_type == "cubo":
    cutoff = n
    length = 2 * n + 1
    atoms = Octahedron(symbol, length, cutoff, latticeconstant=a0, alloy=False)
if particle_type == "oct" and n>0 :
    cutoff = 0
    length = n
    atoms = Octahedron(symbol, length, cutoff, latticeconstant=a0, alloy=False)

num_atoms = len(atoms)
pos = np.zeros((num_atoms, 3))
for i, at in enumerate(atoms):
    pos[i, :] = list(at.position)

typ = np.array([[i] * (num_atoms//5) for i in range(1, 6)]).flatten()
remainder = num_atoms - len(typ)
if remainder > 0:
    typ = np.concatenate((typ, np.array([i for i in range(1, remainder+1)])))
typ = np.random.permutation(typ)

pos -= np.mean(pos, axis=0)

with open("particle.xyz", "w") as file:
    print(num_atoms, file=file)
    print("", file=file)
    for t, p in zip(typ, pos):
        print(t, *p, file=file)


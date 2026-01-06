#!/bin/bash

n_atoms_max=500
n_steps_max=50000
temperature_max=2000

for i in 1 2 3
do
    export n_atoms=$(echo "scale=2; $n_atoms_max*$RANDOM/32767" | bc | awk '{print int($1+0.5)}')
    export atom_typ1="Ag"
    export atom_typ2="Co"

    export n_steps=$(echo "scale=2; $n_steps_max*$RANDOM/32767" | bc | awk '{print int($1+0.5)}')
    export initial_temperature=$(echo "scale=2; $temperature_max*$RANDOM/32767" | bc | awk '{print int($1+0.5)}')
    export final_temperature=$(echo "scale=2; $temperature_max*$RANDOM/32767" | bc | awk '{print int($1+0.5)}')


    #path="$HOME/Documents/MD_Data/Data_${atom_typ1}${atom_typ2}/${n_atoms}Atoms_T=${initial_temperature}K"
    path="$HOME/Desktop/MD_Data/Data_${atom_typ1}${atom_typ2}/${n_atoms}Atoms_T=${initial_temperature}K"
    mkdir -p $path
    export file_input="$path/movie0.xyz"
    export file_output="$path/movie1.xyz"
    touch $file_input $file_output

    echo "Simulation $i"
    echo "Number of atoms : $n_atoms"
    echo "Type of atoms : $atom_typ1 and $atom_typ2"
    echo "Number of simulation steps : $n_steps"
    echo "Initial temperature : $initial_temperature"
    echo "Final temperature : $final_temperature"

    ./construct_random.py
    make run

    param_file="$path/param.txt"
    {
        echo "Number of atoms : $n_atoms"
        echo "Atom types : $atom_typ1, $atom_typ2"
        echo "Number of steps : $n_steps"
        echo "Initial temperature : $initial_temperature K"
        echo "Final temperature : $final_temperature K"
    } > "$param_file"
    
    echo "Simulation $i done"
    echo "Saved in $path"
    echo
done

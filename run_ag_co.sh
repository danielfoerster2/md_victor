#!/bin/bash

n_atoms_max=5000
n_steps_min=100000
n_steps_max=10000000
temperature_max=1000
export final_temperature=300
export density=$(echo 4/5.0^3 | bc -l)
export atom_typ1="Ag"
export atom_typ2="Co"

mkdir -p movies
for i in {1..1000} 
do
    export n_atoms=$(echo "scale=2; $n_atoms_max*$RANDOM/32767" | bc | awk '{print int($1+0.5)}')
    export composition=$(echo "scale=5; $RANDOM/32767" | bc)
    export n_steps=$(echo "scale=2; $n_steps_max*$RANDOM/32767+$n_steps_min" | bc | awk '{print int($1+0.5)}')
    export initial_temperature=$(echo "scale=2; $temperature_max*$RANDOM/32767+1500" | bc | awk '{print int($1+0.5)}')

    rm -rf $1
    mkdir $1
    export file_input="$1/movie0.xyz"
    export file_output="$1/movie1.xyz"

    # To remove after (just for test)
    echo "Number of atoms : $n_atoms"
    echo "Number of steps : $n_steps"
    echo "Initial temperature : $initial_temperature"


    ./construct_random.py
    echo "Particules initialized"
    echo "MD simulation ..."
    make run

    i_sim=$(echo "($1-1)*(1000) + $i" | bc)
    mv $file_output movies/$i_sim.xyz
    rm -rf $1
    echo $i_sim $n_atoms $composition $n_steps $initial_temperature >> data.dat

done

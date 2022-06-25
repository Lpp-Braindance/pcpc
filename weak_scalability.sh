#!/bin/bash
nBodies=$1 #prende in input come primo parametro il numero di bodies 
n_proc=$2  #prende in input come secondo parametro il numero di massimo di n processi
proc=1
size=0
while [[ ((proc -le n_proc)) ]]
do  
    ((size = nBodies * proc))
    mpirun -np $proc nbody_def.out $size -t
    wait
((proc++))
done

#!/bin/bash

rm parallel_* # elimina tutti gli eventuali file creati da un test precedente 

progname=$1  #primo parametro path del programma da eseguire
nBodies=$2   #secondo parametro numero di bodies 
n_proc=$3    #terzo parametro numero massimo di n processi
loop=$4		 #quarto parametro indica se il programma deve essere eseguito in loop
n_iters=1	 #quinto parametro, se il quarto viene specificato, indica il n° esecuzioni del programma
proc=1		 #variabile per incrementare il n° di processi
i=1			 #variabile per iterazione i-esima

if [ "$#" -eq 5 ] # viene specificato il n° di esecuzioni del programma
then 
	n_iters=$5
else
	n_iters=1
fi
while [[ ((i -le n_iters)) ]]
do
proc=1 # reset del numero di processi
while [[ ((proc -le n_proc)) ]]
	do  
		mpirun -np $proc --oversubscribe $progname $nBodies -t
		wait # attesa che il programma termini la sua esecuzione 
		# confronto dei risultati con quelli del programma sequenziale tramite il comando diff 
		# diff -q solo se ci stanno differenze le ristituisce altrimenti restituisce stringa vuota
		DIFF=$(diff -q parallel_$proc parallel_1 ) 
		if [ "$DIFF" != "" ] 
		then
			printf "parallel program np %d \!\!\! ERRORE \!\!\! \n\n" $proc
			exit # termina lo script
		else
			printf "parallel program np %d --> OK\n\n%s" $proc $DIFF
		fi
		wait # aspetta diff
		((proc++)) # aumenta n° processi
	done
	# se viene passato solo -t e senza specificare il n° di itearzioni va in loop infinito
	# dunque non incrementiamo la variabile i (arrestare manualmente)
	if [ "$#" -eq 3 ] || [ "$#" -eq 5 ] 
	then
		((i++)) # iterazione successiva
	fi
done

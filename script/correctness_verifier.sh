#!/bin/bash


nBodies=$1 #prende in input come primo parametro il numero di bodies 
n_proc=$2  #prende in input come secondo parametro il numero di massimo di n processi
#controlla se viene specificato il parametro -t per eseguire il programma che salva i risultati su file
#progname=../nbody_split_allgath.out
progname=../nbody_split.out
if [[ "$#" -eq 3 && $3 == "-t" ]] 
then
    mpirun -np 1 $progname $nBodies -t
fi
#per non mostare l'output del programma, ma solo l'output dei test, specificare -dn
if [[ "$#" -eq 4 && $4 == "-dn" ]] 
then
    mpirun -np 1 $progname $nBodies -t 1>/dev/null
fi
# attendiamo che il programma sequenziale termini la sua esecuzione
wait
printf "sequential program ready\n\n"
# settiamo la variabile proc la quale va a specificare il numero di processi per l'i-esima esecuzione 
proc=2
# testiamo la correttezza del programma eseguendo lo stesso con un numero progressivo di processi 
# cioè eseguiamo il programma con 1,2,3,...,max_nproc 
while [[ ((proc -le n_proc)) ]]
do  
	if [ "$#" -eq 4 ] 
	then
	    mpirun -np $proc $progname $nBodies -t 1>/dev/null
    	else
	    mpirun -np $proc $progname $nBodies -t
	fi
	# attendiamo che il programma termini la sua esecuzione 
	wait
	# confrontiamo il risultati che ha generato e memorizzato nel suo corrispondente file con quelli
	# generati dal file sequenziale tramite il comando diff nel quale specifichiamo il parametro -q
  # che indica a diff di fornire l'output solo se ci sono differenze tra i due file
	DIFF=$(diff -q parallel_$proc parallel_1 )
	if [ "$DIFF" != "" ] 
	then
	    printf "parallel program np %d \!\!\! ERRORE \!\!\! \n\n" $proc
		break
	else
	    printf "parallel program np %d --> OK\n\n%s" $proc $DIFF
	fi
((proc++))
done  

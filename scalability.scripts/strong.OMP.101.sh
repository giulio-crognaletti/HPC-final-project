#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=50:00:00
#PBS -q dssc

cd $PBS_O_WORKDIR

module load openmpi/4.0.3/gnu/9.3.0

./compile

echo $PBS_NODEFILE
echo "STRONG SCALABILITY -- OMP K=101" &> strong.OMP.101.data

for procs in 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1; do

        echo "######################################" &>>strong.OMP.101.data
        echo "Running now on " ${procs} " processors" &>>strong.OMP.101.data
        echo "######################################" &>>strong.OMP.101.data

	export OMP_NUM_THREADS=${procs}

	for reps in 1 2 3; do

		/usr/bin/time ./blur.omp 1 101 101 0.2 "earth-large.pgm" &>>strong.OMP.101.data
	done

done 

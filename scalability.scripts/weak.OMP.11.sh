#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=50:00:00
#PBS -q dssc

cd $PBS_O_WORKDIR

module load openmpi/4.0.3/gnu/9.3.0

./compile

echo "WEAK SCALABILITY -- OMP K=11" &> weak.OMP.11.data

for procs in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24; do

	echo "###################################" &>>weak.OMP.11.data
	echo "Running now on " ${procs} " threads" &>>weak.OMP.11.data
	echo "###################################" &>>weak.OMP.11.data

	export OMP_NUM_THREADS=${procs}

	for reps in 1 2 3; do

		/usr/bin/time ./blur.omp 1 11 11 0.2 ${procs}".pgm" &>>weak.OMP.11.data
	done

done 

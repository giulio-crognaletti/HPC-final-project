#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=50:00:00
#PBS -q dssc

cd $PBS_O_WORKDIR

module load openmpi/4.0.3/gnu/9.3.0

./compile

echo "WEAK SCALABILITY -- MPI K=11" &> weak.MPI.11.data

for procs in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24; do

        echo "######################################" &>>weak.MPI.11.data
	echo "Running now on ",${procs}," processors" &>>weak.MPI.11.data
	echo "######################################" &>>weak.MPI.11.data

	for reps in 1 2 3; do

		/usr/bin/time mpirun --mca btl '^openib' -np ${procs} ./blur.mpi 1 11 11 0.2 ${procs}".pgm" &>>weak.MPI.11.data
	done

done 

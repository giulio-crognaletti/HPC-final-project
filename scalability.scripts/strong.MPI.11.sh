#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l walltime=50:00:00
#PBS -q dssc

cd $PBS_O_WORKDIR

module load openmpi/4.0.3/gnu/9.3.0

./compile

echo $PBS_NODEFILE
echo "STRONG SCALABILITY -- MPI K=11" &> strong.MPI.11.data

for procs in 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1; do

	echo "######################################" &>>strong.MPI.11.data
	echo "Running now on " ${procs} " processors" &>>strong.MPI.11.data
	echo "######################################" &>>strong.MPI.11.data

	for rep in 1 2 3; do

		/usr/bin/time mpirun --mca btl '^openib' -np ${procs} ./blur.mpi 1 11 11 0.2 "earth-large.pgm" &>>strong.MPI.11.data
	done

done 

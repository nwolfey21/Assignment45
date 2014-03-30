amos:
	mpicxx -O3 assn3.c -o assn3

kratos:
	mpicc assn3.c -o assn3

local:
	mpicc -O5 -Wall assn3.c -o assn3

runamos:
	srun --nodes=1 --time 15 --partition small --ntasks=1 --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n 1 -t 1 -r 1

runlocal1:
	mpirun -np 1 ./assn3 -n 1 -t 1 -r 1

runlocal2:
	mpirun -np 2 ./assn3 -n 2 -t 2 -r 1

runlocal4:
	mpirun -np 4 ./assn3 -n 4 -t 4 -r 1

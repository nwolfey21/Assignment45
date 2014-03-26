amos:
	mpicxx -O3 assn3.c -o assn3

runamos:
	srun --nodes=1 --time 15 --partition small --ntasks=1 --overcommit --runjob-opts="--mapping TEDCBA" ./assn3

kratos:
	mpicc assn3.c -o assn3

run1:
	mpicc assn3.c -o assn3
	mpirun -np 1 ./assn3

run2:
	mpicc assn3.c -o assn3
	mpirun -np 2 ./assn3

run4:
	mpicc assn3.c -o assn3
	mpirun -np 4 ./assn3

run8:
	mpicc assn3.c -o assn3
	mpirun -np 8 ./assn3

run12:
	mpicc assn3.c -o assn3
	mpirun -np 12 ./assn3

run16:
	mpicxx -O3 assn3.c -o assn3
	srun --partition small --time 15 --ntasks 16 --nodes 16 ./assn3

run32:
	Mpixlc -O3 assn3.c -o assn3
	mpirun -np 32 ./assn3

run64:
	Mpixlc -O3 assn3.c -o assn3
	mpirun -np 64 ./assn3

run128:
	Mpixlc -O3 assn3.c -o assn3
	mpirun -np 128 ./assn3

run256:
	Mpixlc -O3 assn3.c -o assn3
	mpirun -np 256 ./assn3

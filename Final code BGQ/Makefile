amos:
	mpixlc -O3 assn3.c -o assn3

kratos:
	mpicc assn3.c -o assn3

local:
	mpicc assn3.c -O5 -Wall -lm -o assn3

runamos:
	srun --nodes=64 --time 30 --partition small --ntasks=4096 --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n 64 -t 64 -r 0 -a 1 -b 1

runlocal1:
	mpirun -np 1 ./assn3 -n 1 -t 1 -r 0

runlocal2:
#	rm data/cMatrix.dat
	mpirun -np 2 ./assn3 -n 2 -t 2 -r 0 -a 0 -b 0 -p

runlocal4:
#	rm data/cMatrix.dat
	mpirun -np 4 ./assn3 -n 2 -t 2 -r 2 -a 0 -b 0 -p

runlocal8:
#	rm data/cMatrix.dat
	mpirun -np 8 ./assn3 -n 2 -t 2 -r 4 -a 0 -b 0 -p

Assignment3
============

Compile: 
	module load xl
	make amos

Run 1 submission: 
	srun --partition small --time 15 --ntasks 64 --nodes 32 --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n 32 -t 64

Run Jobs:
	./job.sh

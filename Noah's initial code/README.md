Assignment3
============

Compile: 
	module load xl
	make amos

Run 1 submission: 
	srun --partition small --time 15 --ntasks 64 --nodes 32 --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n 32 -t 64 -r 4

Run Jobs:
	./job.sh

Flags:
	-n Total number of nodes
	-t Total number of tasks
	-r Total numnber of ranks per file. Except 1=all ranks share one file.

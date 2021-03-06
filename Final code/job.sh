#!/bin/bash
#--------------------------------------------------------
# Bash script to automate job submissions to BGQ. Jobs are
# are submitted one after another to BGQ with R2 tasks and
# V1 nodes
#--------------------------------------------------------
#Create Directory if doesn't exist
mkdir -p ./data;

R1=16	#number of cores per node
R2=0

for V1 in 256 128		#number of nodes
do
	for V2 in 1 2 4			#number of tasks per core
	do
		let "R2 = $V1 * $R1 * $V2"
		srun --partition medium --time 15 --ntasks $R2 --nodes $V1 --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n $V1 -t $R2
	done
done
for V1 in 64 32		#number of nodes
do
	for V2 in 1 2 4			#number of tasks per core
	do
		let "R2 = $V1 * $R1 * $V2"
		srun --partition small --time 15 --ntasks $R2 --nodes $V1 --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n $V1 -t $R2
	done
done

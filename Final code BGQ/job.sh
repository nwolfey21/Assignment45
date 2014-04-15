#!/bin/bash
#--------------------------------------------------------
# Bash script to automate job submissions to BGQ. Jobs are
# are submitted one after another to BGQ with R2 tasks and
# V1 nodes
#--------------------------------------------------------
#Create Directory if doesn't exist
mkdir -p ./data;

THREADS=0
NTASKS=0

for NODES in 64		#number of nodes
do
	for RANKS in 64 16 4 1			#number of MPI ranks per node
	do
		for RANKSPERFILE in 0 4 8 32	#number of ranks per file 
		do
			for BLOCKED in 0 1			#Blocked = 1 compact =0
			do
				let "THREADS = 64 / $RANKS"
				let "NTASKS = $RANKS * $NODES"
				srun --partition small --time 15 --ntasks $NTASKS --nodes $NODES --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n $NODES -t $RANKS -a $THREADS -b $BLOCKED -r $RANKSPERFILE
			done
		done	
	done
done

#for NODES in 128 256		#number of nodes
#do
#	for RANKS in 64 16 4 1			#number of MPI ranks per node
#	do
#		for RANKSPERFILE in 0 4 8 32	#number of ranks per file 
#		do
#			for BLOCKED in 1 0			#Blocked = 1 compact =0
#			do
#				let "THREADS = 64 / $RANKS"
#				let "NTASKS = $RANKS * $NODES"
#				srun --partition medium --time 30 --ntasks $NTASKS --nodes $NODES --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n $NODES -t $RANKS -a $THREADS -b $BLOCKED -r $RANKSPERFILE
#			done
#		done	
#	done
#done

#for NODES in 128 256		#number of nodes
#do
#	for RANKS in 64 16 4 1			#number of MPI ranks per node
#	do
#		for RANKSPERFILE in 0 4 8 32	#number of ranks per file 
#		do
#			for BLOCKED in 0 1			#Blocked = 1 compact =0
#				let "THREADS = $64 - RANKS"
#				srun --partition medium --time 30 --ntasks $RANKS --nodes $NODES --overcommit --runjob-opts="--mapping TEDCBA" ./assn3 -n $NODES -t $RANKS -a $THREADS -b $BLOCKED -r $RANKSPERFILE
#			done
#		done	
#	done
#done

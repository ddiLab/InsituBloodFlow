#!/bin/sh
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES * 12))
#CONNOR: modify this file path
EMB_PATH=/home/murphyc/BloodFlow/examples/singleCell
#mpirun -f $COBALT_NODEFILE -n $PROCS $EMB_PATH/cellFlow $EMB_PATH/in.lmp4cell
mpirun -f $COBALT_NODEFILE -n 1 $EMB_PATH/cellFlow $EMB_PATH/in.lmp4cell

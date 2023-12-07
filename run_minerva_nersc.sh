#!/bin/bash
for i in {0..5}
do
  shifter --image=fermilab/fnal-wn-sl7:latest --module=cvmfs -- /bin/bash run_flat.sh $i
  shifter --image=fermilab/fnal-wn-sl7:latest --module=cvmfs -- /bin/bash run_minerva.sh $i
done 

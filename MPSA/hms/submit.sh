#!/bin/sh
for (( num=0; num < $1; num++ )); do
bsub -q linux64 -o ./$3/mpsa-out-$num.txt $2 $num
done

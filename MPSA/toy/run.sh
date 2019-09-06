#!/bin/sh
for (( num=0; num < 3; num++ )); do
./MPSA $num 3 5 200.0 100
done

#!/bin/bash
i=0
while test $i -le 39
do
echo DataTPs[$i]=$((($i + 1) * 10000));
i=$(($i + 1))
done

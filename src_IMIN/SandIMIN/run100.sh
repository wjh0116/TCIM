#!/bin/sh
for i in {1..10}
do
    ./IMIN -dataset dataset/Stanford -k 100 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Twitter -k 100 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Youtube -k 100 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
done
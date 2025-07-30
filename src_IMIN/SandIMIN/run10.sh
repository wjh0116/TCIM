#!/bin/sh
for i in {1..10}
do
    ./IMIN -dataset dataset/Facebook -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/EmailCore -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/EmailAll -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/DBLP -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Epinions -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Stanford -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Twitter -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Wiki-Vote -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Youtube -k 10 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
done
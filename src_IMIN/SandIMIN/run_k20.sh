#!/bin/sh
for i in {1..10}
do
    ./IMIN -dataset dataset/Facebook -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/EmailCore -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/EmailAll -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/DBLP -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Epinions -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Stanford -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Twitter -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Wiki-Vote -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Youtube -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset dataset/Pokec -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    ./IMIN -dataset ../../Data/dataset/Orkut -k 20 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
done
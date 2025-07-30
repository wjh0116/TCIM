#!/bin/sh
for i in {1..5}
do
    ./run -dataset ../SandIMIN/dataset/Youtube -k 100 -rumorNum 10 -algo GR -model IC
done

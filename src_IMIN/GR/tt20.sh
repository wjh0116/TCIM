#!/bin/sh
for i in {1..5}
do
    ./run -dataset ../SandIMIN/dataset/Twitter -k 20 -rumorNum 10 -algo GR -model IC
done




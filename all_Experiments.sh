#!/bin/sh

datasets=(
    EmailCore
    Facebook
    Wiki-Vote
    flixster
    EmailAll
    DBLP
    Twitter
    Stanford
    Youtube
    Pokec
    Orkut
)

selected_datasets=(
    flixster
    Pokec
    Orkut
)

algorithms=(
    SandTCIM # Algorithm: CBFM
    SandTCIM-WS # Algorithm: CBFM-WS
    SandTCIM-PR # Algorithm: CBFM-PR
    GB # Algorithm: GB
)

# Exp1
    # Algorithms proposed in the paper
    for algo in "${algorithms[@]}"; do
        for dataset in "${datasets[@]}"; do
            ./TCIM -dataset ../../dataset/${dataset} -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo ${algo} -T 64 -A 1.0
        done
    done
    # algorithm SandIMIN 
    for dataset in "${datasets[@]}"; do
        ./IMIN -dataset ../../dataset/${dataset} -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 1.0
    done

# Exp2
    # Algorithms proposed in the paper
    for algo in "${algorithms[@]}"; do
        for dataset in "${selected_datasets[@]}"; do
            ./TCIM -dataset ../../dataset/${dataset} -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo ${algo} -T 64 -A 0.8
            ./TCIM -dataset ../../dataset/${dataset} -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo ${algo} -T 64 -A 0.9
            ./TCIM -dataset ../../dataset/${dataset} -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo ${algo} -T 64 -A 1.0
            ./TCIM -dataset ../../dataset/${dataset} -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo ${algo} -T 64 -A 1.1
            ./TCIM -dataset ../../dataset/${dataset} -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo ${algo} -T 64 -A 1.2
        done
    done
    # algorithm SandIMIN 
    for dataset in "${selected_datasets[@]}"; do
        ./IMIN -dataset ../../dataset/${dataset} -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 0.8
        ./IMIN -dataset ../../dataset/${dataset} -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 0.9
        ./IMIN -dataset ../../dataset/${dataset} -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 1.0
        ./IMIN -dataset ../../dataset/${dataset} -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 1.1
        ./IMIN -dataset ../../dataset/${dataset} -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 1.2
    done

# Exp3
    # The running command is the same as Exp1, only the input rumorSet_10.txt needs to be modified (according to different IM algorithms, where the code of HIST is available online, and DEG is easy to achieve). 

# Exp4
    # Algorithms proposed in the paper
    for algo in "${algorithms[@]}"; do
        for dataset in "${datasets[@]}"; do
            ./TCIM_node_specific -dataset ../../dataset/${dataset} -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo ${algo} -T 64 -A 1.0
        done
    done
    # algorithm SandIMIN 
    for dataset in "${datasets[@]}"; do
        ./IMIN_node_specific -dataset ../../dataset/${dataset} -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 1.0
    done

# Exp5
    # To output the results for the case study, we need to perform the following steps, using the CBFM algorithm as an example:
    # Uncomment lines 70, 106, and 107 in Sandwich.h.
    # Save the changes and run make to recompile the code.
    # For each dataset, generate an additional file named real_voting_times.txt, which stores the actual voting time of each user.
    # Run the following command:
    ./TCIM_node_specific -dataset ../../dataset/Facebook -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo SandTCIM -T 64 -A 1.0
    ./TCIM_node_specific -dataset ../../dataset/Twitter -k 20 -rumorNum 10 -epsilon 0.1 -gamma 0.1 -algo SandTCIM -T 64 -A 1.0
    ./IMIN_node_specific -dataset ../../dataset/Facebook -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 1.0
    ./IMIN_node_specific -dataset ../../dataset/Twitter -k 20 -rumorNum 10  -epsilon 0.1 -gamma 0.1 -beta 0.1  -T 64 -A 1.0

# Exp6

    datasets_exp6=(
    Stanford
    Twitter
    Youtube
    Pokec
    Orkut
    )

    # Algorithm: GR
    for dataset in "${datasets_exp6[@]}"; do
        ./GR -dataset ../../dataset/${dataset} -k 100 -rumorNum 10 -algo GR -model IC
    done

    # Algorithm: SandIMIN
    for dataset in "${datasets_exp6[@]}"; do
        ./SandIMIN -dataset ../../dataset/${dataset} -k 100 -rumorNum 10 -algo SandIMIN -epsilon 0.1 -gamma 0.1 -beta 0.1
    done

    # Algorithm: CBFM (SandIMIN++)
    for dataset in "${datasets_exp6[@]}"; do
        ./CBFM -dataset ../../dataset/${dataset} -k 100 -rumorNum 10 -algo SandIMIN++ -epsilon 0.1 -gamma 0.1 -beta 0.1
    done



# Exp7
    # algorithm SA
    OMP_NUM_THREADS=1 ./namm -i dataset/emailcore -fakeseeds seeds/ec.seeds -fakeinf seeds/ec.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/ec
    OMP_NUM_THREADS=1 ./namm -i dataset/facebook -fakeseeds seeds/fb.seeds -fakeinf seeds/fb.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/fb
    OMP_NUM_THREADS=1 ./namm -i dataset/wiki -fakeseeds seeds/wiki.seeds -fakeinf seeds/wiki.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/wiki
    OMP_NUM_THREADS=1 ./namm -i dataset/flixster -fakeseeds seeds/fs.seeds -fakeinf seeds/fs.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/fs
    OMP_NUM_THREADS=1 ./namm -i dataset/emailall -fakeseeds seeds/ea.seeds -fakeinf seeds/ea.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/ea
    OMP_NUM_THREADS=1 ./namm -i dataset/twitter -fakeseeds seeds/tt.seeds -fakeinf seeds/tt.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/tt
    OMP_NUM_THREADS=1 ./namm -i dataset/stanford -fakeseeds seeds/sf.seeds -fakeinf seeds/sf.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/sf
    OMP_NUM_THREADS=1 ./namm -i dataset/youtube -fakeseeds seeds/yt.seeds -fakeinf seeds/yt.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/yt
    OMP_NUM_THREADS=1 ./namm -i dataset/pokec -fakeseeds seeds/pk.seeds -fakeinf seeds/pk.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/pk
    OMP_NUM_THREADS=1 ./namm -i dataset/orkut -fakeseeds seeds/ok.seeds -fakeinf seeds/ok.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/ok
    # algorithm SA+
    OMP_NUM_THREADS=1 ./namm2 -i dataset/emailcore -fakeseeds seeds/ec.seeds -fakeinf seeds/ec.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/ec
    OMP_NUM_THREADS=1 ./namm2 -i dataset/facebook -fakeseeds seeds/fb.seeds -fakeinf seeds/fb.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/fb
    OMP_NUM_THREADS=1 ./namm2 -i dataset/wiki -fakeseeds seeds/wiki.seeds -fakeinf seeds/wiki.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/wiki
    OMP_NUM_THREADS=1 ./namm2 -i dataset/flixster -fakeseeds seeds/fs.seeds -fakeinf seeds/fs.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/fs
    OMP_NUM_THREADS=1 ./namm2 -i dataset/emailall -fakeseeds seeds/ea.seeds -fakeinf seeds/ea.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/ea
    OMP_NUM_THREADS=1 ./namm2 -i dataset/twitter -fakeseeds seeds/tt.seeds -fakeinf seeds/tt.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/tt
    OMP_NUM_THREADS=1 ./namm2 -i dataset/stanford -fakeseeds seeds/sf.seeds -fakeinf seeds/sf.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/sf
    OMP_NUM_THREADS=1 ./namm2 -i dataset/youtube -fakeseeds seeds/yt.seeds -fakeinf seeds/yt.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/yt
    OMP_NUM_THREADS=1 ./namm2 -i dataset/pokec -fakeseeds seeds/pk.seeds -fakeinf seeds/pk.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/pk
    OMP_NUM_THREADS=1 ./namm2 -i dataset/orkut -fakeseeds seeds/ok.seeds -fakeinf seeds/ok.inf -k 20 -epsilon 0.1 -delta 0.01 -o results/ok


# Code for the CBFM algorithm

This project implements the CBFM algorithm for the following paper:

- "Time-Critical Influence Minimization via Node Blocking", 2026.

---

## Dataset

All the datasets in this paper can be downloaded in [SNAP](https://snap.stanford.edu/).

---

## Compile

```bash
make
```


## Run

- **Arguments**
  - `dataset`: path to the dataset directory
  - `k`: budget 
  - `T`: deadline 
  - `rumorNum`: number of seed nodes of misinformation
  - `algo`: algorithm 
  - `epsilon`, `gamma`: the parameters


- **Example**
  ```bash
  ./TCIM -dataset dataset/Facebook -k 20 -rumorNum 10 -T 64 -algo SandTCIM -epsilon 0.2 -gamma 0.2

# Time-Critical Influence Minimization via Node Blocking

This repository contains the code implementation for the paper "Time-Critical Influence Minimization via Node Blocking".



## Directory Structure

```
├── all_Experiments.sh          # Script to run all experiments in the paper
├── dataset/  
│   ├── flixster/                   
│   └── Facebook/            
├── src_IMIN/                   # algorithms for IMIN problem
│   ├── README                
│   ├── CBFM/                  
│   ├── GR/                     
│   └── SandIMIN/             
├── src_MM/                     # algorithms for MM problem
│   ├── README                     
│   └── SA_SA+/                
├── src_TCIM/                   # algorithms for TCIM problem
│   ├── README                
│   ├── CBFM/                  
│   ├── CBFM_node_specific/    
│   ├── SandIMIN/        
│   └── SandIMIN_node_specific/ 
├──────────────────────────────   
```


## Datasets

The following datasets can be downloaded from [SNAP](https://snap.stanford.edu/):
- EmailCore, Facebook, Wiki-Vote
- EmailAll, DBLP
- Twitter, Stanford, Youtube
- Pokec, Orkut

**Note:** The flixster dataset is already included in the `dataset/` directory.



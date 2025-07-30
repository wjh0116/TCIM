

# TCIM problem


## Algorithms

- **CBFM**: Unified setting of CBFM
- **CBFM_node_specific**: Node-specific setting of CBFM
- **SandIMIN**: Unified setting of SandIMIN
- **SandIMIN_node_specific**: Node-specific setting of SandIMIN

## Compilation

```bash
make
```


### Parameters

- `dataset`: Path to the dataset directory
- `k`: Budget (number of nodes to block)
- `T`: Time deadline
- `rumorNum`: Number of rumor seed nodes
- `algo`: Algorithm name (SandTCIM, SandTCIM-WS, SandTCIM-PR, GB)
- `epsilon`, `gamma`: Algorithm parameters
- `A`: parameter in the weight function

### Example Commands

For running examples, please refer to the `all_Experiments.sh` script in the root directory.

### Note for the node-specific setting implementation

An additional file containing $x_v$ for each node is required. 


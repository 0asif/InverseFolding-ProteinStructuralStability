#!/bin/bash

python3 protein_mpnn_run.py \
        --pdb_path inputs/spike_af3.pdb \
        --out_folder outputs/spike/ \
        --num_seq_per_target 1000000 \
        --sampling_temp "0.1" \
        --seed 42 \
        --save_probs 1 \
        --save_score 1 \
        

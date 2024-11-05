This package contains the information required to run NMA analysis and compiling data
1) AlphaFold2 (AF2) models/structures
2) NMR ensembles that are mapped to AF2 models (one-to-one mapping), in case of multiple pairs, add one-to-many pairs to list_pairs (see test_example/AF2/NMR_analysis/list_pairs.csv)

Each analysis 1) and 2) requires two main modules to run
- python environment containing WEBnma (installed on Linux/OS): install webnma using "https://github.com/reuter-group/webnma3"

- MDAnalysis 

Important Note: The MDAnalysis module used for testing the code is 'MDAnalysis/2.2.0-foss-2022a'
This module is not compatible with WEBnma python environment. Therefore, 
to run python scripts:
1. "runwebnma_python_single.py" or "runwebnma_python_perstate.py"
It is necessary to do "module purge MDAnalysis/2.2.0-foss-2022a", in case it is loaded, otherwise WEBnma will throw error.

for doing other analysis: do:  "module load MDAnalysis/2.2.0-foss-2022a"
2) compute_fluctuations.py / NMR_compute_flucs.py
3) compute_contacts.py
4) truncate_structures.py
5) NMR_analysis.py


The USAGE of each file is provided at the bottom of each file with decription of the file.

****************************************************************************************

The workflow to use the script is provided in workflow_NMA.pdf

                                WORKFLOW (test_example)
test_example: contains the examples that were tested. 

****************************************************************************************

                                  NMA on AF2 MODELS
                                  

1) IF YOU HAVE MULTIPLE STRUCTURES TO RUN WEBNMA

**[pdb_files]: directory containing .pdb files of AF2 models (non-truncated)

%%%%%%%%%%%%%%%%%%%%%%%% source activate WEBNMA_conda_env  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

step 1) Run webnma (runwebnma_python_single.py) on non-truncated AF2 models **[pdb_files]  > save output (modes.txt) in **[non_trunc_webnmaoutput]

%%%%%%%%%%%%%%%%%%% step 2) module load MDAnalysis/2.2.0-foss-2022a %%%%%%%%%%%%%%%%%%%

step 3) Compute fluctuations (compute_fluctuations.py) on modes.txt of non-truncated AF2 models > save output .csv fluctuations in **[non_trunc_flucs]

step 4) Compute contacts (compute_contacts.py) on .pdb files **[pdb_files] of AF2 models (non-truncated) > save output .csv files in **[contacts_nontrunc]

step 5) Truncate AF2 models based on contacts (truncate_structure.py) > save output in nc_info.csv (dataframe containing truncation info) and the truncated pdb files in **[trunc_pdbfiles]

%%%%%%%%%%%%%%%%%%% step 6) module purge MDAnalysis/2.2.0-foss-2022a %%%%%%%%%%%%%%%%%%%

step 7) Run webnma (runwebnma_python_single.py) on truncated AF2 models **[trunc_pdbfiles]  > save output (modes.txt) in **[trunc_webnmaoutput]

%%%%%%%%%%%%%%%%%%% step 6) module load MDAnalysis/2.2.0-foss-2022a %%%%%%%%%%%%%%%%%%%

step 8) Compute fluctuations (compute_fluctuations.py) on modes.txt of truncated AF2 models **[trunc_pdbfiles] > save output .csv fluctuations in **[trunc_flucs]

________________________________________________________________________________________
1) IF YOU HAVE A SINGLE STRUCTURE TO RUN WEBNMA

**[single]: This is the folder with a single .pdb file which is used to run the code when there is a single structure to run.

Follow the above step 1) to step 8) in the same order as above.
The usage for a single structure is provided in each file for the test example.


*****************************************************************************************
                                        NMA on NMR ENSEMBLES

pre-requisite: 
-mapping_df.csv, the .csv file containing a list of UniprotID, bmrbID, pdbID in exactly the same column names.
-the binary file of STRIDE installed on your machine. (change the path of the binary to your path in NMR_analysis.py)


%%%%%%%%%%%%%%%%%%% step 1) module load MDAnalysis/2.2.0-foss-2022a %%%%%%%%%%%%%%%%%%%

step 1) run NMR_analysis.py [see the usage in the file] > output is saved in the folders below
**[nmr_state_models]: the directory containing pair_1/state_1.pdb, state_2.pdb, and state_N.pdb

**[nmr_ss_states] : the directory containing secondary structures for each state of the pair_1/ss_state_1.mol, ss_state_2.mol, ss_state_N.mol

**[nmr_models] : the .pdb files of downloaded nmr models corresppnding to AF2 models

list_pairs.csv : a list of uniprotID_pdbID (A1Y2K1_2mrj) pairs is created which can be used to run webnma in next step.

%%%%%%%%%%%%%%%%%%% step 2) module purge MDAnalysis/2.2.0-foss-2022a %%%%%%%%%%%%%%%%%%%

step 3) run webnma on NMR ensembles (runwebnma_python_perstate.py) on **[NMR_models] > save output of modes.txt in **[nmr_state_models/uniprotID_pdbID/webnma_states]

%%%%%%%%%%%%%%%%%%% step 3) module load MDAnalysis/2.2.0-foss-2022a %%%%%%%%%%%%%%%%%%%

step 4) run fluctuations (NMR_compute_flucs.py) on **[nmr_state_models/uniprotID_pdbID/webnma_states] > save output in [nmr_state_models/uniprotID_pdbID/all_flucs.csv]







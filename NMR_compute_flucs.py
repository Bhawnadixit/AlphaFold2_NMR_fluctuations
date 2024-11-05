import argparse
import pandas as pd
import NMR_analysis

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run WEBnma on NMR ensembles.')
    
    parser.add_argument('--in_file', help='Path of file containing Uniprot-bmrb-pdb mapping')
    parser.add_argument('--out_path', help='Path to output directory for saving NMA data of NMR ensembles')
    args = parser.parse_args()
    # initialize the function by providing the Dataframe containing uniprotID, bmrbID, and PDB_ID (if there are no bmrbID, the lines in the code for bmrbid can be commented out)
    
    mapping_dataframe = pd.read_csv(args.in_file, index_col=0, header=0)
    print(mapping_dataframe)
    N = NMR_analysis.NMR_analysis(mapping_dataframe, args.out_path)
    
    # create a list of triplets ([uid1, bid1, pid1], [uid2, bid2, pid2],.. )
    N.mapped_ids()
    
    # compute fluctuations on NMR ensembles for each state and compile them into a single dataframe.
    # we will call it in a separate script
    N.compute_flucs_perstate()
    
    
# USAGE EXAMPLE
# python NMR_compute_flucs.py --in_file test_example/AF2/mapping_df.csv --out_path test_example/AF2/NMR_analysis
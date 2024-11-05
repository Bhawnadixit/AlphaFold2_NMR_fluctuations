from Bio import PDB
import MDAnalysis as mda
import os 
import pandas as pd
import numpy as np
import fluctuations
import re 
import argparse

class NMR_analysis:
    def __init__(self, df, nmr_analysispath):
        self.df = df       # df containing 'UniprotID', 'bmrbID', and 'pdbID'
        self.nmr_analysispath = nmr_analysispath # HOME PATH FOR NMR ANALYSIS
        self.mapped_ids_list = None # ID triplet containing UniprotID_BMRBID_PDBID
        
        # CREATE RELEVANT DIRECTORIES
        
        # path where all the NMR models will be saved in a folder named as (UniprotID)
        self.nmr_modelpath = os.path.join(self.nmr_analysispath, 'nmr_models')
        os.makedirs(self.nmr_modelpath, exist_ok=True)
        
        # another directory where NMR models will be split into each state and saved into a folder named as UniprotID_PDBID
        self.nmr_state_models = os.path.join(self.nmr_analysispath, 'nmr_state_models')         
        os.makedirs(self.nmr_state_models, exist_ok=True)
        
        # another directory where ss downloaded from stride for per-state of NMR models will be saved into a folder named as UniprotID_PDBID
        self.nmr_ss_states = os.path.join(self.nmr_analysispath, 'nmr_ss_states')        
        os.makedirs(self.nmr_ss_states, exist_ok=True)
        
        
    def mapped_ids(self):
        # Provides a list containing triplets [(UniprotID, 'bmrbID', 'pdbID')]
        df_nmrsubset = self.df.dropna(subset='pdbID')  # This step is for to remove any Uniprot IDs for which there are no PDB IDs.
        # print(df_nmrsubset)
        uniprot_ids = df_nmrsubset['UniprotID'].unique()
        mapped_ids_list = []
        for n, x in enumerate(uniprot_ids):
            uid = x   # uniprot_ID
            bid = self.df[self.df['UniprotID'] == uid]['bmrbID'].unique()[0]   # bmrb_ID
            pid = self.df[self.df['UniprotID'] == uid]['pdbID'].unique()[0]
            print(uid, bid, pid)
            mapped_ids_list.append([uid, bid, pid]) # f'{uid}_{bid}_{pid}'
        self.mapped_ids_list = mapped_ids_list
        return self.mapped_ids_list
    
    def download_PDBfile(self):
        # Downloads the .pdb files via biopython's pdb module in specific directories
        for x in self.mapped_ids_list:
            uid = x[0]
            bid = x[1]
            pid = x[2]
            #CREATE PATH FOR UNIPROT ID TO SAVE PDB WITHIN THE ID
            download_dir = os.path.join(self.nmr_modelpath, uid)
            if not os.path.exists(download_dir):
                os.makedirs(download_dir)
                
            # List of PDB IDs from your dataframe
            pdb_ids = [pid]
            # print(pdb_ids)
            # DOWNLOAD PDB STRUCTURE
            
            pdb_list = PDB.PDBList()
            pdb_list.download_pdb_files(pdb_ids, file_format='pdb', pdir=download_dir)
            print(f'{pid} is downloaded for {uid}-{bid}.')
        
    def split_NMRensemble(self):
        # function to split the NMR .pdb files containing the ensemble into a single state .pdb file
        # all files in a single NMR ensemble are saved in the directory known as nmr_state_models/UniprotID_PDBID/state_1.pdb, state_n.pdb
        for x in self.mapped_ids_list:
            uid = x[0]
            bid = x[1]
            pid = x[2]

            fname = f'{uid}_{pid}' # dir_name

            in_path = os.path.join(os.path.join(self.nmr_modelpath, uid), f'pdb{pid}.ent') # the input path of downloaded NMR model
            out_path = os.path.join(self.nmr_state_models, fname)    # the output path of per-state .pdb file of each NMR ensemble
            # print(out_path)

            if not os.path.exists(out_path):
                os.makedirs(out_path)

            nmr = mda.Universe(in_path, in_path)
            print(fname, nmr.trajectory)

            nmr_protein = nmr.select_atoms('protein')

            for n, ts in enumerate(nmr.trajectory, 1):
                output_filename = os.path.join(out_path, f"{fname}_{n}.pdb")  # Adjusted to start frame numbering from 1
                with mda.Writer(output_filename, multiframe=True) as pdb:
                    pdb.write(nmr_protein)
                    print(f'The frame no. {n} for {pid} is written.')
                        
    def ss_NMRensemble(self, stride_binary_path):
        # function to compute secondary structure code for each .pdb state file for the NMR ensemble, the output files are .mol files saved in a separate directory.
        # load the binary file from pre-installed stride or provide path after installing stride
        
        for x in self.mapped_ids_list:
            uid = x[0]
            bid = x[1]
            pid = x[2]

            fname = f'{uid}_{pid}' # dir_name

            in_path = os.path.join(self.nmr_state_models, fname)
            out_path = os.path.join(self.nmr_ss_states, fname)
            # print(out_path)

            print(x)
            
            if not os.path.exists(out_path):
                os.makedirs(out_path)

            for n, state in enumerate(os.listdir(in_path), 1):
                if state.endswith('.pdb'):
                    state_name = state
                    state_inpath = os.path.join(in_path, state_name)
                    if os.path.exists(state_inpath):
                        print(state_inpath)
                        state_out = state_name.split('.')[0]
                        state_outpath = os.path.join(out_path, f'{state_out}.mol')
                        # print(state_outpath)
                        os.system(f'{stride_binary_path} -m {state_inpath} > {state_outpath}')
    
    def ss_states_singledf(self):
        # combines .mol files from secondary structure of each state into a single dataframe
        for x in self.mapped_ids_list:
            uid = x[0]
            bid = x[1]
            pid = x[2]

            fname = f'{uid}_{pid}' # dir_name

            in_path = os.path.join(self.nmr_state_models, fname)
            out_path = os.path.join(self.nmr_ss_states, fname)

            ss_states = {}
            for n, state in enumerate(os.listdir(in_path), 1):
                if state.endswith('.pdb'):
                    state_name = state.split('.')[0]
                    state_num = state.split('.')[0].split('_')[-1]
                    # print(state_name)
                    state_inpath = os.path.join(in_path, state)
                    state_outpath = os.path.join(out_path, f'{state_name}.mol')
                    if os.path.exists(state_outpath):
                        # print(state_outpath)
                        with open(state_outpath) as file:
                            lines = file.readlines()
                            # print(lines)
                            start_line = None
                            for jline, line in enumerate(lines):
                                if '|---Residue---|' in line:
                                    start_line = jline
                                    # break
            #                 # Skip rows before the structure details line
                        skiprows = start_line + 1 if start_line is not None else 0

                        # Define the column names from the header line
                        header_line = lines[start_line]
                        column_names = header_line.split()[1:-1]

                        df = {}
                        # Split each line and append as a row to the DataFrame
                        for n, line in enumerate(lines[start_line + 1:]):
                            row_data = line.split()
                            # seq_resname.append(row_data[1])
                            # ss_oneletter.append(row_data[])
                            # print(row_data)
                            df.update({n: row_data})
                        ss_df = pd.DataFrame(df).T.iloc[:,1:-1]
                        # print(ss_df)

                        ss_df.columns = ['resname', 'chainid', 'resid1', 'resid2', 'stride_nmr_ss', 'stride_nmr_ss_full_name', 
                                         'Phi', 'Psi', 'Area']
                        # print(ss_df)

        #                     print(ss_df['stride_af_ss'])
        #                     print(state_name)
                        ss_states.update({f'stride_ss (NMR_state_{state_num})': ss_df['stride_nmr_ss']})
            ss_states_df = pd.DataFrame(ss_states)
            ss_states_df[['resname', 'resid1',  'chainid']] = ss_df[['resname', 'resid1',  'chainid']]
            # print(ss_states_df)
            sorted_labels = sorted(list(ss_states_df.columns), key=self.sort_key)

            # Print the sorted labels
            # print(sorted_labels)

            reordered_ssdf = ss_states_df[[c for c in sorted_labels if not c.startswith('stride_ss')] + [c for c in sorted_labels if c.startswith('stride_ss')]]
            print(reordered_ssdf)
                
            
            reordered_ssdf.to_csv(os.path.join(in_path, f'{fname}_ss.csv'))
            print(f'ss of all states is combined into a single df and saved for {fname}')


    
    def sort_key(self, label):
        # function to reorder the column labels containing state_1, state_5, state_3 > state_1, state_3, state_5
        
        # Match labels with the pattern: prefix followed by (NMR_state_number)
        match = re.match(r"(.+?)\s*\(NMR_state_(\d+)\)", label)
        if match:
            # Return a tuple: (prefix in lowercase, numeric suffix as integer)
            return (match.group(1).strip().lower(), int(match.group(2)))
        else:
            # For labels that don't match the pattern, return a tuple: (label in lowercase, None)
            return (label.lower(), None)

            
    def compute_flucs_perstate(self):
         # To compute fluctuations from the WEBnma modes.txt for each state model from the NMR ensemble
        # MAKE SURE TO RUN WEBnma (runwebnma_python_perstate.py) before calling the function
        
        for x in self.mapped_ids_list:
            uid = x[0]
            bid = x[1]
            pid = x[2]

            fname = f'{uid}_{pid}' # dir_name
            print(fname)

            in_path = os.path.join(self.nmr_state_models, fname)
            out_path = os.path.join(in_path, 'webnma_states')

            try:
                flucs_dict = {}
                for n, state in enumerate(os.listdir(out_path), 1):
                    state_name_pdb = f'{state}'
                    state_num = state.split('.')[0].split('_')[-1]
                    state_name_modes = os.path.join(out_path, state)
                    # print(state_name_modes)

                    path_protein = state_name_modes+'/modes.txt'              # check path correct*****
                    
                    if os.path.exists(path_protein):
                        # print(path_protein)
                        evals, sequence, evecs, nres = fluctuations.read_modes_2(path_protein, doprint=False)
                        # print(evals)
                    #                 # rename variables
                        es = evals              # eigenvalues
                        rs = sequence           # amino-acid sequence 
                        modes = evecs           # eigenvectors
                        # print(es, modes, rs)

                        # double check
                        np.sum(modes**2, axis=1) ## not normed in rows
                        np.sum(modes**2, axis=0) ## normed in columns
                        # so every column is a mode
                        # print(i, np.sum(modes[:,0]**2))

                        nres = len(rs)
                        nmodes = modes.shape[1]
                        # print(nres)
                        # print(nres)
                        # Calculate RMSF (31/1/24) [RMSF in angstrom]
                        fluc_rmsf = fluctuations.calc_fluc_2(es, rs, modes, [], donorm = "rmsf")
                        # print(fluc_rmsf)
                        sequence = [i.replace('1.', 'A.') for i in sequence]
                        # print(sequence)
                        seq_id = [re.findall(r'\d+', i.split('.')[1])[0] for i in sequence]


                        # seq_id = [re.findall(r'\d+', x)[0] for x in sequence]
                        seq_label = [re.findall(r'\D+', i)[0].split('.')[1] for i in sequence]
                        chain_id = [re.findall(r'\D+', i)[0].split('.')[0] for i in sequence]
                        seq_oneletter = [fluctuations.swapped_aa_dict.get(i) for i in seq_label]
                        # print(seq_label, chain_id, seq_oneletter)

                        fluc_df = pd.DataFrame([seq_id, seq_oneletter, seq_label, chain_id, fluc_rmsf]).T
                        fluc_df.columns = ['seqIndex', 'seqLabel', 'residue', 'chain id NMA' , f'RMSF (NMR_state_{state_num})']
                        # print(fluc_df)
                        flucs_dict.update({f'RMSF (NMR_state_{state_num})': fluc_rmsf})
                flucs_dict.update({'seqIndex': seq_id, 'seqLabel': seq_oneletter, 'residue': seq_label, 'chain id NMA' : chain_id})
 
                # print(pd.DataFrame(flucs_dict))
                all_flucs_df = pd.DataFrame(flucs_dict)
                # print(all_flucs_df)
                # print(list(all_flucs_df.columns))
                # Sort the column labels using the custom sort key
                sorted_labels = sorted(list(all_flucs_df.columns), key=self.sort_key)

                # Print the sorted labels
                # print(sorted_labels)
                
                reordered_flucsdf = all_flucs_df[[c for c in sorted_labels if not c.startswith('RMSF')] + [c for c in sorted_labels if c.startswith('RMSF')]]
                print(reordered_flucsdf)
                
                                                 
                reordered_flucsdf.to_csv(os.path.join(in_path, f'{fname}_flucs_allstates.csv'))
#                 print(f'fluctuations of all states is combined into a single df and saved for {fname}')
            except:
                pass

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Run WEBnma on NMR ensembles.')
    
    parser.add_argument('--in_file', help='Path of file containing Uniprot-bmrb-pdb mapping')
    parser.add_argument('--out_path', help='Path to output directory for saving NMA data of NMR ensembles')
    args = parser.parse_args()
    # initialize the function by providing the Dataframe containing uniprotID, bmrbID, and PDB_ID (if there are no bmrbID, the lines in the code for bmrbid can be commented out)
    
    mapping_dataframe = pd.read_csv(args.in_file, index_col=0, header=0)
    print(mapping_dataframe)
    N = NMR_analysis(mapping_dataframe, args.out_path)
    
    # create a list of triplets ([uid1, bid1, pid1], [uid2, bid2, pid2],.. )
    N.mapped_ids()
    
    # download the relevant PDB files corresponding to UniprotID of AlphaFold2 models 
    N.download_PDBfile()
    
    # split the NMR ensemble into single .pdb files (states)
    N.split_NMRensemble()
    
    # compute the secondary structure on each state of NMR ensemble
    N.ss_NMRensemble(stride_binary_path = '$VSC_SCRATCH_VO_USER/NMA/stride_pdb/stride')
    
    # combine the secondary structure of all states into a single dataframe
    N.ss_states_singledf()
    
    # Run WEBnma (WEBnma dependencies are not compatible with MDAnalysis), so make sure MDanalysis module is purged before running WEBnma
    # Create a .csv list of pairs to run webnma
    list_pairs = pd.Series([f'{x}_{y}' for x, y in zip(mapping_dataframe['UniprotID'], mapping_dataframe['pdbID'])])
    list_pairs.to_csv(os.path.join(args.out_path, 'ids_list_pairs.csv'))
    # os.system('python runwebnma_python_perstate.py')
    
    # compute fluctuations on NMR ensembles for each state and compile them into a single dataframe.
    # we will call it in a separate script [NMR_compute_flucs.py]
    # N.compute_flucs_perstate()
    
# USAGE EXAMPLE    
# python NMR_analysis.py --in_file test_example/AF2/mapping_df.csv --out_path test_example/AF2/NMR_analysis


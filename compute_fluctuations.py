import os 
import pandas as pd
import numpy as np
import fluctuations
import re
import argparse
import MDAnalysis as mda

def compute_flucs(in_path, out_path, pdb_path):
    # print(in_path)
    if os.path.isdir(pdb_path):
        for n, a in enumerate(os.listdir(in_path), 1):
            if not a.startswith('.ipynb'):
                a_name_pdb = a
                a_name_modes = f'{os.path.join(in_path, a)}'
                print(a_name_modes)

                path_protein = os.path.abspath(os.path.join(a_name_modes, 'modes.txt'))              # check path correct*****
                if os.path.exists(path_protein):
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
                    # Calculate RMSF (31/1/24) [RMSF in angstrom]
                    fluc_rmsf = fluctuations.calc_fluc_2(es, rs, modes, [], donorm = "rmsf")

                    sequence = [x.replace('1.', 'A.') for x in sequence]
                    # print(sequence)
                    seq_id = [re.findall(r'\d+', x.split('.')[1])[0] for x in sequence]


                    # seq_id = [re.findall(r'\d+', x)[0] for x in sequence]
                    seq_label = [re.findall(r'\D+', x)[0].split('.')[1] for x in sequence]
                    chain_id = [re.findall(r'\D+', x)[0].split('.')[0] for x in sequence]
                    seq_oneletter = [fluctuations.swapped_aa_dict.get(x) for x in seq_label]
                    pdb_id_col = [a_name_pdb for x in sequence]
                    # print(seq_label, chain_id, seq_oneletter)

                    structure = mda.Universe(os.path.join(pdb_path, a_name_pdb+'.pdb'))
                    # print(structure)
                    plddt = [residue.tempfactor for residue in structure.select_atoms('protein and name CA')]

                    fluc_df = pd.DataFrame([seq_id, seq_oneletter, seq_label, chain_id, fluc_rmsf, pdb_id_col, plddt]).T

                    fluc_df.columns = ['seqIndex', 'seqLabel', 'residue', 'chain id NMA' , 'RMSF', 'ID', 'pLDDT']
                    print(fluc_df)

                    if not os.path.exists(out_path):
                        os.makedirs(out_path)
                    # print(out_path+'/'+a+'_flucs_NMA.csv')
                    fluc_df.to_csv(out_path+'/'+a+'_flucs_NMA.csv')

        print('DONE')
    elif pdb_path.endswith('.pdb'):
        file_name = os.path.basename(pdb_path)  # Extracts 'file.pdb'
        a_name_pdb = os.path.splitext(file_name)[0]  # Extracts 'file'
        a_name_modes = f'{os.path.join(in_path, a_name_pdb)}'
        print(a_name_modes)
        path_protein = os.path.abspath((os.path.join(a_name_modes, 'modes.txt')))           # check path 
        
        if os.path.exists(path_protein):
            # print(path_protein)

            evals, sequence, evecs, nres = fluctuations.read_modes_2(path_protein, doprint=False)
            # print(evals)
                        # rename variables
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
            # Calculate RMSF (31/1/24) [RMSF in angstrom]
            fluc_rmsf = fluctuations.calc_fluc_2(es, rs, modes, [], donorm = "rmsf")

            sequence = [x.replace('1.', 'A.') for x in sequence]
            # print(sequence)
            seq_id = [re.findall(r'\d+', x.split('.')[1])[0] for x in sequence]


            # seq_id = [re.findall(r'\d+', x)[0] for x in sequence]
            seq_label = [re.findall(r'\D+', x)[0].split('.')[1] for x in sequence]
            chain_id = [re.findall(r'\D+', x)[0].split('.')[0] for x in sequence]
            seq_oneletter = [fluctuations.swapped_aa_dict.get(x) for x in seq_label]
            pdb_id_col = [a_name_pdb for x in sequence]
            # print(seq_label, chain_id, seq_oneletter)

            structure = mda.Universe(pdb_path)
            # print(structure)
            plddt = [residue.tempfactor for residue in structure.select_atoms('protein and name CA')]

            fluc_df = pd.DataFrame([seq_id, seq_oneletter, seq_label, chain_id, fluc_rmsf, pdb_id_col, plddt]).T

            fluc_df.columns = ['seqIndex', 'seqLabel', 'residue', 'chain id NMA' , 'RMSF', 'ID', 'pLDDT']
            print(fluc_df)

            if not os.path.exists(out_path):
                os.makedirs(out_path)
            # print(out_path+'/'+a+'_flucs_NMA.csv')
            fluc_df.to_csv(out_path+'/'+a_name_pdb+'_flucs_NMA.csv')
            print('DONE')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute fluctuations from webnma modes')
    parser.add_argument('--in_path', help='Path to input directory containing webnma modes file')
    parser.add_argument('--out_path', help='Path to output directory for saving fluctuations computed from webnma')
    parser.add_argument('--pdb_path', help='Path to input directory containing .pdb files / or a single .pdb file')
    args = parser.parse_args()

    compute_flucs(args.in_path, args.out_path, args.pdb_path)
    
# _____________________________________________________________________
# USAGE EXAMPLE: 
#  STEP 1: NON-TRUNCATED ALPHAFOLD2 MODELS

# EXAMPLE 1, WHEN THERE ARE MULTIPLE FILES
# python compute_fluctuations.py -in_path test_example/AF2/non_trunc_webnmaoutput --out_path test_example/AF2/non_trunc_flucs -pdb_path test_example/AF2/pdb_files/

# EXAMPLE 2, WHEM THERE IS A SINGLE .PDB FILE
# python compute_fluctuations.py --in_path test_example/AF2/single/non_trunc_webnmaoutput --out_path test_example/AF2/single/non_trunc_flucs --pdb_path test_example/AF2/single/A7ZTU8.pdb

# _____________________________________________________________________

# STEP 2: TRUNCATED ALPHAFOLD2 MODELS

# EXAMPLE 1, WHEN THERE ARE MULTIPLE FILES
# python compute_fluctuations.py --in_path test_example/AF2/trunc_webnmaoutput --out_path test_example/AF2/trunc_flucs --pdb_path test_example/AF2/trunc_pdbfiles

# EXAMPLE 2, WHEM THERE IS A SINGLE .PDB FILE
# python compute_fluctuations.py --in_path test_example/AF2/single/trunc_webnmaoutput --out_path test_example/AF2/single/trunc_flucs --pdb_path test_example/AF2/single/trunc_pdbfiles
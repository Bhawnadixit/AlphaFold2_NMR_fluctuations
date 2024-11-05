import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
import argparse
import shutil


def compute_trunc(in_path, out_path):
    if os.path.isdir(in_path):
        contacts_dict = {i: pd.read_csv(os.path.join(in_path, i), index_col=0) for i in os.listdir(in_path) if not i.startswith('.ipynb')}
        
    redo = []
    trunc_dict = {}
    for x in list(contacts_dict.keys()):
        if not x.startswith('.ipynb'):
            # print(x)
            x = x.split('.')[0]
            cons = contacts_dict[x+'.csv']
            cons.columns = ['contacts']
            # print(cons)
            nterm = cons
            cterm = cons.iloc[::-1]
            # print(nterm, cterm)
            # print(cterm)

            # Set the threshold to 13 contacts, if first residue or last residue os greater than 13, we dont cut it, otherwise, find the residue to cut for N and C.
            threshold = 13
            if nterm['contacts'].iloc[0] >= threshold and cterm['contacts'].iloc[0] >= threshold:
                # print(x)
                print(x, 'Nothing to cut')
                # continue


            ## NOW IF FIRST RESIDUE OF N-TERM IS >= THRESHOLD, WE ONLY CUT C-TERM RIGHT BEFORE THE THRESHOLD RESIDUE
            elif nterm['contacts'].iloc[0] >= threshold and cterm['contacts'].iloc[0] < threshold: 
                first_occurrence_c = cterm.index[cterm['contacts'] >= threshold][0]
                c_ind = list(cterm.index).index(first_occurrence_c) - 1
                c_residue = cterm.index[c_ind]
                original_c_index = len(cons.index) - 1 - c_ind
                protein_length = cons.shape[0]
                # print(c_ind)
                # print(x)
                # print(f'protein length {protein_length}')
                # print(f'C-term residue to cut: {c_residue}, and its index {original_c_index}')
                # print('NO N terminal to cut @@@@@@@@@@@@@@@@@@@@@@@C')

                trunc_dict.update({x: [None, None, original_c_index, c_residue, protein_length]})


            ## ELSE IF FIRST RESIDUE OF C-TERM IS >= THRESHOLD, WE ONLY CUT N-TERM RIGHT BEFORE THE THRESHOLD RESIDUE
            elif cterm['contacts'].iloc[0] >= threshold and nterm['contacts'].iloc[0] < threshold:
                first_occurrence_n = nterm.index[nterm['contacts'] >= threshold][0]
                n_ind = list(nterm.index).index(first_occurrence_n) - 1
                n_residue = nterm.index[n_ind]
                protein_length = cons.shape[0]

                # print(x)
                # print(f'protein length {protein_length}')
                # print(f'N-term residue to cut: {n_residue}, and its index {n_ind}')
                # print('NO C terminal to cut N@@@@@@@@@@@@@@@@@@@@@@@@@')
                trunc_dict.update({x: [n_ind, n_residue, None, None, protein_length]})


            ## ELSE IF NEITHER OF THE FIRST 2 CONDITIONS ARE MET WE SCAN BOTH N- AND C-TERMINI FOR CUTTING
            elif (nterm['contacts'].iloc[0] < threshold and cterm['contacts'].iloc[0] < threshold) and ((nterm['contacts'] >= threshold).any() and (cterm['contacts'] >= threshold).any()):
                # print(x)
                # print(nterm['contacts'])
                # Now we look at where to cut N- and C-, and calculate how far apart they are, as they could have lower than 13 average contacts.            
                first_occurrence_n = nterm.index[nterm['contacts'] >= threshold][0] 
                first_occurrence_c = cterm.index[cterm['contacts'] >= threshold][0]

                n_ind = list(nterm.index).index(first_occurrence_n)  - 1         # index of item before first occurance in n-term
                c_ind = list(cterm.index).index(first_occurrence_c)  - 1         # index of item before first occurance in c-term
                n_residue = nterm.index[n_ind]                                     # residue at that index
                c_residue = cterm.index[c_ind]


                redo.append(x)
                # print(c_residue)
                original_c_index = len(cons.index) - 1 - c_ind                     # this is original index in sequence, since we reversed the C-term before search
                protein_length = cons.shape[0]
                # print(n_ind, original_c_index, protein_length)
                # print(x)
                # print(f'protein length {protein_length}')
                # print(f'N-term residue to cut: {n_residue}, and its index {n_ind}')
                # print(f'C-term residue to cut: {c_residue}, and its index {original_c_index}')

                trunc_dict.update({x: [n_ind, n_residue, original_c_index, c_residue, protein_length]})

            else:
            #     try:
                print(x, 'Nothing worked here ££££££££££££££££££££££££££££££££££££')
            #     except:
            #         pass
            # print()


    nc_info = pd.DataFrame(trunc_dict).T
    nc_info.columns = ['N-index', 'N-residue', 'C-index', 'C-residue', 'protein length (aa)']
    # print(nc_info)
    nc_info.to_csv(os.path.join(out_path, 'nc_info.csv'))
    return nc_info
    
def truncate_model(nc_info, pdb_path, out_path):
    
    for a in nc_info.index:
        # print(a)
        if os.path.isdir(pdb_path):
            md_af = mda.Universe(os.path.join(pdb_path, a+'.pdb'))
        elif pdb_path.endswith('.pdb'):
            md_af = mda.Universe(os.path.join(pdb_path))
        
        resid_resname = [(x.resid, x.resname) for x in md_af.select_atoms('name CA')]
        # print(resid_resname)
        n_term_seq = resid_resname
        c_term_seq = resid_resname[::-1]
        # print(n_term_seq[0][0], c_term_seq[0][0])

        trunc_n = nc_info.loc[a]['N-residue']
        trunc_c = nc_info.loc[a]['C-residue']
        # print(trunc_n, trunc_c)

        trunc_pdb_path = os.path.join(out_path, 'trunc_pdbfiles')
        if not os.path.exists(trunc_pdb_path):
            os.makedirs(trunc_pdb_path)
        
        
        if trunc_c == None:
            print(a, 'NO C TERMINAL **********')
            trunc_ind_n= resid_resname[nc_info.loc[a]['N-index']][0]
            trunc_protein = md_af.select_atoms(f'(not resid {n_term_seq[0][0]}:{trunc_ind_n}) and name CA')
            output_pdb_path = os.path.join(trunc_pdb_path, a+'_trunc.pdb')
            trunc_protein.write(output_pdb_path)


            # print(trunc_n, trunc_c )
            # print(n_term_seq[nc_info.loc[a]['N-index']])


        elif trunc_n == None:
            print(a, 'NO N terminal """"""""""""""')
            # print(nc_info.loc[a])
            trunc_ind_c = resid_resname[nc_info.loc[a]['C-index']][0]
            # print(trunc_ind_c)
            trunc_protein = md_af.select_atoms(f'(not resid {trunc_ind_c}:{c_term_seq[0][0]}) and name CA')
            output_pdb_path = os.path.join(trunc_pdb_path, a+'_trunc.pdb')
            trunc_protein.write(output_pdb_path)

        else:
            print(a, 'BOTH TERMINAL FOUND')
            trunc_ind_n= resid_resname[nc_info.loc[a]['N-index']][0]
            trunc_ind_c = resid_resname[nc_info.loc[a]['C-index']][0]
            trunc_protein = md_af.select_atoms(f'(not resid {n_term_seq[0][0]}:{trunc_ind_n}) and (not resid {trunc_ind_c}:{c_term_seq[0][0]}) and name CA')
        # print(trunc_protein)
            output_pdb_path = os.path.join(trunc_pdb_path, a+'_trunc.pdb')
            trunc_protein.write(output_pdb_path)


def check_trunc(nc_info, pdb_path, out_path):
    if os.path.isdir(pdb_path):
        non_trunc = [x for x in os.listdir(pdb_path) if not x.startswith('.ipynb')]
    elif pdb_path.endswith('.pdb'):
        file_name = os.path.basename(pdb_path)  # Extracts 'file.pdb'
        non_trunc = [file_name]
    
    trunc_pdb_path = os.path.join(out_path, 'trunc_pdbfiles')
    if not os.path.exists(trunc_pdb_path):
        os.makedirs(trunc_pdb_path)
    
    trunc = [y.split('_')[0]+'.pdb' for y in os.listdir(trunc_pdb_path) if not y.startswith('.ipynb')]

    ids_without_trunc = [g for g in non_trunc if not g in trunc]
    # print(trunc, ids_without_trunc)
    
    if len(ids_without_trunc) != 0:
        print('ID(s) that do not require truncation and are copied to trunc_files for better data compilation', ids_without_trunc)
        for x in ids_without_trunc:
            pdb_id = x.split('.')[0]
            shutil.copy(os.path.join(pdb_path, x), os.path.join(trunc_pdb_path, f'{pdb_id}_trunc.pdb'))
    else:
        print('No structures left to truncate/move')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='computes the truncation amino acid residues for truncating N and/or C-terminal')
    parser.add_argument('--in_path', help='Path to input directory containing per-residue contacts')
    parser.add_argument('--out_path', help='Path to output directory for saving N/C truncation info')
    parser.add_argument('--pdb_path', help='Path to input directory of PDB structure(s) to truncate')
    args = parser.parse_args()
    
    NC = compute_trunc(args.in_path, args.out_path) 
    truncate_model(NC, args.pdb_path, args.out_path)
    check_trunc(NC, args.pdb_path, args.out_path)
    
    
# USAGE EXAMPLE
# EXAMPLE 1, WHEN THERE ARE MULTIPLE FILES
# python truncate_structure.py --in_path test_example/AF2/contacts_nontrunc --out_path test_example/AF2/truncated_pdbfiles --pdb_path test_example/AF2/pdb_files

# EXAMPLE 2, WHEM THERE IS A SINGLE .PDB FILE
# python truncate_structure.py --in_path test_example/AF2/single/contacts_nontrunc --out_path test_example/AF2/single --pdb_path test_example/AF2/single/A7ZTU8.pdb


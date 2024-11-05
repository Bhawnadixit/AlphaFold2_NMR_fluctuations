import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import argparse


def contacts_within_cutoff(u, group_a, group_b, radius):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)


def compute_contacts(in_path, out_path):
    radius = 10
    if not os.path.exists(out_path):
        os.makedirs(out_path)
        
    if os.path.isdir(in_path):
        for i in os.listdir(in_path):
            print(i)
            if not i.startswith('.ipynb'):
                pdb_id = i.split('.')[0]
                protein = mda.Universe(os.path.join(in_path, i))

                ca = protein.select_atoms('name CA')
                n_ca = len(ca)

                seq_con = {}
                for name, ids in zip(ca.resnames, ca.resids):
                    sel1 = ca.select_atoms(f'resid {ids}')
                    cont = contacts_within_cutoff(protein, sel1, ca, radius)
                    seq_con.update({name + '_' + format(ids): cont[0][1]})

                seq_condf = pd.Series(seq_con)
                seq_condf.to_csv(os.path.join(out_path,  pdb_id + '.csv'))
        print('Done')
        
    elif in_path.endswith('.pdb'):
        file_name = os.path.basename(in_path)  # Extracts 'file.pdb'
        pdb_id = os.path.splitext(file_name)[0]  # Extracts 'file'
        protein = mda.Universe(in_path)
        ca = protein.select_atoms('name CA')
        n_ca = len(ca)

        seq_con = {}
        for name, ids in zip(ca.resnames, ca.resids):
            sel1 = ca.select_atoms(f'resid {ids}')
            cont = contacts_within_cutoff(protein, sel1, ca, radius)
            seq_con.update({name + '_' + format(ids): cont[0][1]})

        seq_condf = pd.Series(seq_con)
        seq_condf.to_csv(os.path.join(out_path,  pdb_id + '.csv'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute contacts from .pdb files')
    parser.add_argument('--in_path', help='Path to input directory containing .pdb files')
    parser.add_argument('--out_path', help='Path to  directory for saving contacts computed from .pdb files')
    args = parser.parse_args()

    compute_contacts(args.in_path, args.out_path)
    
# USAGE EXAMPLE
# EXAMPLE 1, WHEN THERE ARE MULTIPLE FILES
# python compute_contacts.py --in_path test_example/AF2/pdb_files/ --out_path test_example/AF2/contacts_nontrunc

# EXAMPLE 2, WHEM THERE IS A SINGLE .PDB FILE
# python compute_contacts.py --in_path test_example/AF2/single/A7ZTU8.pdb --out_path test_example/AF2/single/contacts_nontrunc


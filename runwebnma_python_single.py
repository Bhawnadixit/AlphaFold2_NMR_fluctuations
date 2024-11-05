import os
import argparse
import pandas as pd

# Requires WEBnma installed on the machine within a linux environment

def run_webnma(in_path, out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    if os.path.isdir(in_path):
        invalid_nma = []
        for num, a in enumerate(os.listdir(in_path), 1):
            if not a.startswith('.ipynb'):
                print(num, a)
                pdb_id = a

                if a.endswith('.pdb'):
                    pdb_id = a.split('.')[0]
                    in_pdb_path = os.path.join(in_path, f'{a}')

                    out_path_a = os.path.join(out_path, pdb_id)
                    if not os.path.exists(out_path_a):
                        os.makedirs(out_path_a)
                    try:
                        print(f'webnma mode {in_pdb_path} -p {out_path_a}')
                        os.system(f'webnma mode {in_pdb_path} -p {out_path_a}')
                    except:
                        invalid_nma.append(a)
                        pass

        print(invalid_nma)
        pd.Series(invalid_nma).to_csv(os.getcwd()+'/invalid_NMAerror_ids.csv')
        
    elif in_path.endswith('.pdb'):
        invalid_nma = []
        file_name = os.path.basename(in_path)  # Extracts 'file.pdb'
        pdb_id = os.path.splitext(file_name)[0]  # Extracts 'file'
        in_pdb_path = in_path
        out_path_a = os.path.join(out_path, pdb_id)
        if not os.path.exists(out_path_a):
            os.makedirs(out_path_a)
        try:
            print(f'webnma mode {in_pdb_path} -p {out_path_a}')
            os.system(f'webnma mode {in_pdb_path} -p {out_path_a}')
            print(f'SUCCESS for {pdb_id}')
        except:
            invalid_nma.append(a)
            pass
        print('The proteins with invalid NMA:', invalid_nma)
        pd.Series(invalid_nma).to_csv(os.getcwd()+'/invalid_NMAerror_ids.csv')

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run webnma on protein structures from either a directory of structures or a single structure.')
    parser.add_argument('--in_path', help='Path to input directory containing protein structures (.pdb files) or path of a single .pdb file')
    parser.add_argument('--out_path', help='Path to output directory for webnma modes')
    args = parser.parse_args()
    
    
    run_webnma(args.in_path, args.out_path)
    

# _____________________________________________________________________
# USAGE EXAMPLE: 
#  STEP 1: NON-TRUNCATED ALPHAFOLD2 MODELS

# EXAMPLE 1, WHEN THERE ARE MULTIPLE FILES
# python runwebnma_python_single.py --in_path test_example/AF2/pdb_files/ --out_path test_example/AF2/non_trunc_webnmaoutput

# EXAMPLE 2, WHEM THERE IS A SINGLE .PDB FILE
#  python runwebnma_python_single.py --in_path test_example/AF2/single/A7ZTU8.pdb --out_path test_example/AF2/single/non_trunc_webnmaoutput

# _____________________________________________________________________

# STEP 2: TRUNCATED ALPHAFOLD2 MODELS

# EXAMPLE 1, WHEN THERE ARE MULTIPLE FILES
# python runwebnma_python_single.py --in_path test_example/AF2/trunc_pdbfiles/ --out_path test_example/AF2/trunc_webnmaoutput

# EXAMPLE 2, WHEM THERE IS A SINGLE .PDB FILE
# python runwebnma_python_single.py --in_path test_example/AF2/single/trunc_pdbfiles/ --out_path test_example/AF2/single/trunc_webnmaoutput 1 A7ZTU8_trunc.pdb


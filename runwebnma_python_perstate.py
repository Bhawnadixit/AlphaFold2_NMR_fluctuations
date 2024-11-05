import os
import argparse
import pandas as pd

# Requires WEBnma installed on the machine within a linux environment

def run_webnma(in_path, pairs_list):
    invalid_nmadict = {}
    for num, a in enumerate(pairs_list, 1):
        print(num, a)
        pdb_id = a.split('_')[1]
        # in_pdb_path = os.path.join(in_path, f'{a.split("_")[0]}/pdb{pdb_id}.ent')
        
        out_state_path = os.path.join(os.path.join(in_path, a), 'webnma_states')
        if not os.path.exists(out_state_path):
            os.makedirs(out_state_path)
        invalid_list = []
        for ts in os.listdir(os.path.join(in_path, a)):
            if ts.endswith('.pdb'):
                # print(ts)
                state_name_pdb = ts
                mname = ts.split('.')[0]
                state_name_modes = os.path.join(out_state_path, mname)
                if not os.path.exists(state_name_modes):
                    os.makedirs(state_name_modes)

                try:
                    print(f'webnma mode {os.path.join(os.path.join(in_path, a), state_name_pdb)} -p {state_name_modes}')
                    os.system(f'webnma mode {os.path.join(os.path.join(in_path, a), state_name_pdb)} -p {state_name_modes}')
                except:
                    invalid_list.append(state_name_pdb)
                    pass
        invalid_nmadict.update({a: invalid_list})
    pd.DataFrame(invalid_nma).to_csv(os.path.join(os.path.join(in_path, a), 'invalid_NMAerror_ids.csv'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run webnma on protein structures')
    parser.add_argument('--in_path', help='Path to input directory containing protein structures')
    parser.add_argument('--ids_list', help='List of sequences for webnma', nargs='+')
    parser.add_argument('--ids_csv', help='CSV file containing list of sequences for webnma')
    args = parser.parse_args()

    if args.ids_list:
        pairs_list = args.ids_list
    elif args.ids_csv:
        df = pd.read_csv(args.ids_csv, index_col=0)
        pairs_list = df['0'].tolist()
    else:
        raise ValueError("Either --pairs_list or --pairs_list_csv must be provided.")
    print(pairs_list)
    run_webnma(args.in_path, pairs_list)

    
    
# Usage example: 
# > if input is a list
# python run_webnma.py /path/to/input_directory /path/to/output_directory --ids_list Q52KI8_1mp1 B0FYL5_1xyk P17497_2mp5 Q15121_1n3k P36238_2klm P36655_3pfu P56959_2lcw P35615_2mq9

# > if input is a csv
# python run_webnma.py /path/to/input_directory /path/to/output_directory --ids_csv /path/to/ids.csv
# example
# python runwebnma_python_perstate.py --in_path test_example/AF2/NMR_analysis/nmr_state_models --ids_csv test_example/AF2/NMR_analysis/list_pairs.csv

# NOTE: the IDs used here are uniprot_pdb pairs, but can be any id, given that it has existing folder
# input_dir
# --id1
# --id2
# output_dir = input_dir+webnma_folder
# --id1/modes.txt
# --id2/modes.txt

# The WEBnma command can be edited to retrieve anything other than modes.txt
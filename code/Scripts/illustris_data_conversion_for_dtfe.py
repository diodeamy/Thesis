import h5py
import numpy as np
import os
import illustris_python as il
import sys

sys.path.append('..')

from shared.utilities import get_redshift_dictionary

base_path = '/Users/users/nastase/PROJECT/'
output_dir = '/Users/users/nastase/PROJECT/DATA/files_for_dtfe'
os.makedirs(output_dir, exist_ok=True)

dictionary = get_redshift_dictionary()
snapshot_numbers = dictionary.keys()
# print(f'The keys are: {dictionary.keys()}')

# snapshot_numbers = dictionary['snapshot']

def extract_data(base_path, snapshot):
    positions, masses = il.snapshot.loadSubset(base_path, snapshot, 'dm', ['Coordinates', 'MassTable'])
        
    return positions, masses
    
def save_to_dtfe_format(positions, masses, output_file):
    data = np.column_stack(positions, masses)
    np.savetxt(output_file, data, fmt='%.6f', delimiter=' ')
    
for snapshot in snapshot_numbers:

    positions, masses = extract_data(base_path, snapshot)
    output_file = os.path.join(output_dir, f'snapshot_{snapshot}.dat')
    save_to_dtfe_format(positions, masses, output_file)
    
print('Data extraction and conversion complete')
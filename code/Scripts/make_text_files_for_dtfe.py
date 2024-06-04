import illustris_python as il
import numpy as np
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

sys.path.append('..')

from shared.utilities import get_redshift_dictionary

base_path = "/Users/users/nastase/PROJECT/"
output_path = "/Users/users/nastase/PROJECT/DATA/files_for_dtfe/"

def load_data(snapshot_number):
    data = il.snapshot.loadSubset(base_path, snapshot_number, 'dm', ["Coordinates", "Velocities"])
    data["Masses"] = np.ones(len(data["Coordinates"]))
    
    return data

def process_snapshot(snapshot_number):
    print(f"Start processing {snapshot_number}")
    data = load_data(snapshot_number)
    convert_data_to_tex_even_faster(data, snapshot_number, output_path)
    
    print(f"Finish processing {snapshot_number}")

def convert_data_to_tex_even_faster(data, snapshot_number, path):
    
    particle_number = data["count"]
    coordinates = data["Coordinates"]
    masses = data["Masses"]
    
    box_min = 0
    box_max = 75000
    
    filename = f"{path}/file_{snapshot_number}.txt"
    header = f"{particle_number}\n{box_min:.6f} {box_max:.6f} {box_min:.6f} {box_max:.6f} {box_min:.6f} {box_max:.6f}\n"
    all_data = np.hstack((coordinates, masses.reshape(-1, 1)))

    
    # Convert all data to formatted strings
    lines = "\n".join(f"{x:.6f} {y:.6f} {z:.6f} {m:.6f}" for x, y, z, m in all_data ) 
    
    with open(filename, 'w') as f:
        f.write(header)
        f.write(lines)
        f.write('\n')  # Add a newline at the end
        
def main():
    
    dictionary = get_redshift_dictionary()
    snapshot_numbers = dictionary['snapshots']
#     snapshot_numbers = [135]

#     ##this is for more I/O bound tasks    
#     with ThreadPoolExecutor() as executor:
#         executor.map(process_snapshot, snapshot_numbers)
    
    ##this is for more other kinds of tasks? 
    with ProcessPoolExecutor() as executor:
        executor.map(process_snapshot, snapshot_numbers)

if __name__ == "__main__":
    main()
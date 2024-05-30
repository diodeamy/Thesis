import illustris_python as il
import numpy as np
import sys

sys.path.append('..')

from shared.utilities import get_redshift_dictionary

base_path = "/Users/users/nastase/PROJECT/"
output_path = "/Users/users/nastase/PROJECT/DATA/files_for_dtfe/"

def load_data(snapshot_number):
    data = il.snapshot.loadSubset(base_path, snapshot_number, 'dm', ["Coordinates", "Velocities"])
    data["Masses"] = np.ones(len(data["Coordinates"]))
    
    return data

def convert_data_to_tex_even_faster(data, snapshot_number, path):
    
    particle_number = data["count"]
    coordinates = data["Coordinates"]
    velocities = data["Velocities"]
    masses = data["Masses"]
    
    box_min = 0
    box_max = 75000
    
    filename = f"{path}/file_{snapshot_number}.txt"
    header = f"{particle_number}\n{box_min:.6f} {box_max:.6f} {box_min:.6f} {box_max:.6f} {box_min:.6f} {box_max:.6f}\n"
    
    all_data = np.hstack((coordinates, velocities, masses.reshape(-1, 1)))
    
    # Convert all data to formatted strings
    lines = "\n".join(f"{x:.6f} {y:.6f} {z:.6f} {vx:.6f} {vy:.6f} {vz:.6f} {m:.6f}" for x, y, z, vx, vy, vz, m in all_data)
    
    with open(filename, 'w') as f:
        f.write(header)
        f.write(lines)
        f.write('\n')  # Add a newline at the end
        
def main():
    
    dictionary = get_redshift_dictionary()
    
    for snapshot_number in dictionary["snapshots"]:
        data = load_data(snapshot)
        convert_data_to_tex_even_faster(data, snapshot, output_path)
        
if __name__ == "__main__":
    main()
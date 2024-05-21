import illustris_python as il
import numpy as np
from joblib import Parallel, delayed
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl

base_path = "/Users/users/nastase/PROJECT/"

def get_redshift_dictionary():
    
    snapshot_redshift_pairs = [
        (135, 0.00),
        (133, 0.02),
        (131, 0.05),
        (129, 0.07),
        (127, 0.10),
        (125, 0.13),
        (123, 0.15),
        (122, 0.17),
        (118, 0.23),
        (114, 0.29),
        (110, 0.36),
        (106, 0.44),
        (102, 0.52),
        (98, 0.62),
        (94, 0.73),
        (90, 0.85),
        (86, 0.99), 
        (82, 1.11),
        (78, 1.30),
        (74, 1.53),
        (70, 1.82),
        (66, 2.21),
        (64, 2.44),
        (62, 2.73),
        (60, 3.01),
        (58, 3.28),
        (56, 3.71),
        (54, 4.01),
        (53, 4.18),
        (51, 4.66),
        (49, 5.00),
        
    ]
    
    snapshot_redshifts = {snapshot: redshift for snapshot, redshift in snapshot_redshift_pairs}
    
    return snapshot_redshifts


def get_snapshot_cluster_coordinates(base_snapshot_id: int, subhalo_id: int, test_snapshot_ids: list[int]):
    """
    Params:
    `snapshot_id`: ID of the reference snapshot
    `subhalo_id`: ID of the subhalo. This is dependent on the `snapshot_id`
    `snapshot_ids`: List of IDs

    Returns:
    `test_coordinates`: Dict where the key is the `snapshot_id` and the value is the coordinates for the cluster in that snapshot

    ## Example
    ```python
    # Compute coordinates
    coordinates = get_snapshot_cluster_coordinates(base_snapshot_id=54, subhalo_id=13, test_snapshot_ids=[60])

    # Generate graph
    fig, ax = plt.subplots()

    ax.scatter(coordinates[54][:,0], coordinates[54][:,1], c='blue', alpha=0.03, label='z = 54')
    ax.scatter(coordinates[60][:,0], coordinates[60][:,1], c='red', alpha=0.03, label='z = 60')

    plt.xlabel('x [ckpc/h]')
    plt.ylabel('y [ckpc/h]')
    plt.legend()
    ```
    """    

    subhalo = il.snapshot.loadSubhalo(base_path, base_snapshot_id, subhalo_id, 'DM')
    subhalo_particle_ids = subhalo['ParticleIDs']
    subhalo_coordinates = subhalo['Coordinates']

    test_coordinates = {base_snapshot_id: subhalo_coordinates}
    
    def load_and_match_coordinates(test_snapshot_id):
        coords_ids = il.snapshot.loadSubset(base_path, test_snapshot_id, 'dm', ['Coordinates', 'ParticleIDs'])
        coordinates = get_coordinates_for_particleIDs(coords_ids, subhalo_particle_ids)
        
        return test_snapshot_id, coordinates
    
    results = Parallel(n_jobs=-1)(delayed(load_and_match_coordinates)(test_snapshot_id) for test_snapshot_id in test_snapshot_ids)
    
    for snapshot_id, coordinates in results:
        test_coordinates[snapshot_id] = coordinates

#     for test_snapshot_id in test_snapshot_ids:
#       coords_ids = il.snapshot.loadSubset(base_path, test_snapshot_id, 'dm', ['Coordinates', 'ParticleIDs'])
#       coordinates = get_coordinates_for_particleIDs(coords_ids, subhalo_particle_ids)
#       test_coordinates[test_snapshot_id] = coordinates

    return test_coordinates


def get_coordinates_for_particleIDs(data, particle_ids):
    
    coordinates = data['Coordinates']
    ids = data['ParticleIDs']
    
    sorted_indices = np.argsort(ids)
    sorted_ids = ids[sorted_indices]
    
    pos = np.searchsorted(sorted_ids, particle_ids)
    matched_indices = sorted_indices[pos]
    
    valid_mask = sorted_ids[pos] == particle_ids
    valid_indices = matched_indices[valid_mask]
    
#     id_to_index = {particle_id: index for index, particle_id in enumerate(ids)}
#     particle_indices = [id_to_index.get(particle_id) for particle_id in particle_ids]
#     valid_indices = [index for index in particle_indices if index is not None]
    
    if len(valid_indices) != len(particle_ids):
        print("Warning: Some particle IDs were not found.")
    
    return coordinates[valid_indices]

def load_results(path,base_snapshot_id, subhalo_id):
    filename = f"{path}/coordinates_base_snapshot_{base_snapshot_id}_subhalo_{subhalo_id}.pickle"
    
    print(f"Reading your stuffy stuff from: {filename}")
    with open(filename, "rb") as f:
        python_results = pickle.load(f)
        
    return python_results

def load_all_results(path):
    filename = f"{path}/all_results.pickle"
    
    print(f"Reading ALL stuffy stuff from: {filename}")
    with open(filename, "rb") as f:
        python_all_results = pickle.load(f)
        
    return python_all_results

def plot_2D_histogram(coordinates_dict, coordinate_pair = 'xy', save_path=None):
    """
    plot 2D histograms of particle coordinates of a particular subhalo, for a specific snapshot
    
    Params:
    'coordonates_dict (dict)': Dictionary containing the particle coordinates for a specific cluster
                        Keys are snapshot numbers, values are arrays  of coordinates (x, y, z)
    'coordinate_pair (str)': Pair of coordinates to be plotted. Options: 'xy', 'xz', 'yz'
    'save_path (str)': Path to save the figures. If None, figures will not be saved
    
    
    """
    
    valid_coordinate_pairs = ['xy', 'xz', 'yz']
    if coordinate_pair not in valid_coordinate_pairs:
        raise ValueError("Invalid coordinate pair! Please choose from 'xy', 'xz' or 'yz'.")
  
    
    if coordinate_pair == "xy":
        all_x_coords = [coords[:, 0] for coords in coordinates_dict.values()]
        all_y_coords = [coords[:, 1] for coords in coordinates_dict.values()]
        
    elif coordinate_pair == "xz":
        all_x_coords = [coords[:, 0] for coords in coordinates_dict.values()]
        all_y_coords = [coords[:, 2] for coords in coordinates_dict.values()]        
    else:
        all_x_coords = [coords[:, 1] for coords in coordinates_dict.values()]
        all_y_coords = [coords[:, 2] for coords in coordinates_dict.values()]
        
    all_x_coords = np.concatenate(all_x_coords)
    all_y_coords = np.concatenate(all_y_coords)
    
    x_min, x_max = np.min(all_x_coords), np.max(all_x_coords)
    y_min, y_max = np.min(all_y_coords), np.max(all_y_coords)
    
    for snapshot, coords in coordinates_dict.items():
        plt.figure(figsize=(8,6))
        plt.hist2d(coords[:, 0] if coordinate_pair[0] == 'x' else coords[:,1], 
                   coords[:, 1] if coordinate_pair[1] == 'y' else coords[:,2]
                   , norm=mpl.colors.LogNorm(), bins = 512)
        plt.xlabel("X" if coordinate_pair[0] == 'x' else "Y")
        plt.ylabel("Y" if coordinate_pair[1] == 'y' else "Z")
        plt.title(f"histogram for snapshot {snapshot}")
        plt.colorbar(label="Counts")
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        
        if save_path:
            fig_filename = f"{save_path}/snapshot_{snapshot}.png"
            plt.savefig(fig_filename)
            print(f"Figure saved! File: {fig_filename}")

#     all_x_coords = []
#     all_y_coords = []
    
#     for snapshot, coords in coordinates_dict.items():
#         if coordinate_pair == 'xy':
#             x_coords = coords[:, 0]
#             y_coords = coords[:, 1]
#             xlabel, ylabel = 'X [ckpc/h]', 'Y [ckpc/h]'
#         elif coordinate_pair == 'xz':
#             x_coords = coords[:, 0]
#             y_coords = coords[:, 2]
#             xlabel, ylabel = 'X [ckpc/h]', 'Z [ckpc/h]'
#         else:
#             x_coords = coords[:, 1]
#             y_coords = coords[:, 2]
#             xlabel, ylabel = 'Y [ckpc/h]', 'Z [ckpc/h]'
            
#         all_x_coords.extend(x_coords)
#         all_y_coords.extend(y_coords)
      
#         plt.figure(figsize=(8,6))
#         plt.hist2d(x_coords, y_coords, norm=mpl.colors.LogNorm(), bins = 512)
#         plt.xlabel(xlabel)
#         plt.ylabel(ylabel)
#         plt.title(f"{xlabel} - {ylabel} histogram for snapshot {snapshot}")
#         plt.colorbar(label="Counts")
        
#         plt.xlim(np.min(all_x_coords), np.max(all_x_coords))
#         plt.ylim(np.min(all_y_coords), np.max(all_y_coords))
        
#         if save_path:
#             fig_filename = f"{save_path}/snapshot_{snapshot}.png"
#             plt.savefig(fig_filename)
#             print(f"Figure saved! File: {fig_filename}")
            
        plt.close()
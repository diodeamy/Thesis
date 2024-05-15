import illustris_python as il
import numpy as np

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

    test_coordinates = {}

    test_coordinates[base_snapshot_id] = subhalo_coordinates

    for test_snapshot_id in test_snapshot_ids:
      coords_ids = il.snapshot.loadSubset(base_path, test_snapshot_id, 'dm', ['Coordinates', 'ParticleIDs'])
      coordinates = get_coordinates_for_particleIDs(coords_ids, subhalo_particle_ids)
      test_coordinates[test_snapshot_id] = coordinates

    return test_coordinates


def get_coordinates_for_particleIDs(data, particle_ids):
    
    coordinates = data['Coordinates']
    ids = data['ParticleIDs']
    
    id_to_index = {particle_id: index for index, particle_id in enumerate(ids)}
    particle_indices = [id_to_index.get(particle_id) for particle_id in particle_ids]
    valid_indices = [index for index in particle_indices if index is not None]
    
    if len(valid_indices) != len(particle_indices):
        print("Warning: Some particle IDs were not found.")
    
    return coordinates[valid_indices]
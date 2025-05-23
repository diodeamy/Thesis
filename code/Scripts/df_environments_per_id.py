import numpy as np
import os
import sys
import shutil
nexus_path = os.path.abspath("/Users/users/nastase/Applications/NEXUS_python/NEXUS_1.0/")
sys.path.append(nexus_path)
bin_path = os.path.abspath("/Users/users/nastase/Applications/NEXUS_python/NEXUS_1.0/bin/")
sys.path.append(bin_path)
from density import readDensityHeader, writeDensityData
import re
import illustris_python as il
import pickle
import pandas as pd
import logging

gridsize = 512
sidelen = 106.5
data_directory = "/Users/users/nastase/PROJECT/DATA/nexus_outputs/MMF_outputs/"
base_path = "/Users/users/nastase/PROJECT/"


def giveheader(demoden,DTFEden,outputname,gridsize,sidelen):
    #First we read a file with a good header to edit
#     l = readDensityHeader("demo_density.a_den")
    l = readDensityHeader(demoden)
    #And we read the output from DTFE
    den = np.array(np.fromfile(DTFEden,dtype=np.float32))

    #And we put them in the header
    l.gridSize = np.array([gridsize,gridsize,gridsize])
    l.totalGrid = gridsize**3
    l.box = np.array([0.,sidelen,0.,sidelen,0.,sidelen])
    #Finally we write the new header with the density data to a new NEXUS-ready file
    writeDensityData(outputname,l,den)

def readNEXUS(densityFile,gridsize,vel=False,velFile=None):
    # here we wish to extract the 3 digits corresponding ot the snapshot number, such that we can save the files
    # in the same directory without overwriting
    
    pattern = r'(\d{3})\.a_den$'
    
    # Search for the pattern in the densityFile name
    match = re.search(pattern, densityFile)
    
    if match:
        # Extract the 3 digits
        digits = match.group(1)
#         print("Extracted digits: {digits}")
    else:
        print("No matching digits found before '.a_den'")
    
    #This function will read the files generated by runNEXUS and make them into numpy arrays
    #It can also read DTFE velocity files as an extra option
    #The order is as follows, it returns:
    #Densityfield, NEXUS_nodes, NEXUS_filaments, NEXUS_walls, Velocity (optional)
    #The NEXUS arrays can then be used for contours or converted into boolean arrays
    #Where all cells > 0 are part of the corresponding structure
    #Voids are located where the other three components are all 0 
    shape = (gridsize,gridsize,gridsize)
    denfield = np.fromfile(densityFile,dtype=np.float32)
    denfield = np.reshape(denfield[262:-2],shape)
    MMFn = np.fromfile(f"{data_directory}node_{digits}_clean.MMF",dtype=np.int16)
    MMFn = np.reshape(MMFn[527:-1],shape)
    MMFf = np.fromfile(f"{data_directory}fila_{digits}_clean.MMF",dtype=np.int16)
    MMFf = np.reshape(MMFf[527:-1],shape)
    MMFw = np.fromfile(f"{data_directory}wall_{digits}_clean.MMF",dtype=np.int16)
    MMFw = np.reshape(MMFw[527:-1],shape)
    #Correct for the z-axis missallignment of the density field wrt NEXUS
    def maxmean(densf,MMFfila,axisn,ran):
        k = np.zeros(ran)
        for i in range(ran):
            k[i] = np.mean(np.roll(densf,i,axis=axisn)[MMFfila>0])
        return np.argmax(k)
    
    shift = maxmean(denfield,MMFf,2,gridsize)
    denfield = np.roll(denfield,shift,axis=2)
    if vel==True:
        velo = np.fromfile(velFile,dtype=np.float32)
        k = np.reshape(velo,(gridsize,gridsize,gridsize,3))
        k = np.roll(k,shift,axis=2)
        return denfield, MMFn, MMFf, MMFw, k
    else:
        return denfield, MMFn, MMFf, MMFw
    
    
def generateMask(r, Nex):
    """
    Generates a mask array based on the positions in r and the values in Nex.
    
    Parameters:
    r (array-like): A list of vectors, each containing x, y, and z positions.
    Nex (3D array): A 3D numpy array where each dimension corresponds to x, y, or z. Values are either 1 or 0.
    
    Returns:
    mask (array): An array of 1s and 0s corresponding to whether each vector in r lies on a 1 or 0 in Nex.
    """
    # Convert r to a numpy array for easier manipulation
    r = np.array(r)
    
    # Ensure r is of shape (n, 3)
    if r.shape[1] != 3:
        raise ValueError("Each vector in r must have exactly 3 coordinates (x, y, z)")
    
    # Extract x, y, z coordinates from r
    x_coords = r[:, 0].astype(int)
    y_coords = r[:, 1].astype(int)
    z_coords = r[:, 2].astype(int)
    
    # Ensure the coordinates are within the bounds of Nex
    
    if (x_coords >= Nex.shape[0]).any() or (y_coords >= Nex.shape[1]).any() or (z_coords >= Nex.shape[2]).any():
        
        raise ValueError("Some coordinates in r are out of bounds of the Nex array")
    
    # Generate the mask by indexing into Nex
    mask = Nex[x_coords, y_coords, z_coords]
    
    return mask

def getMorphologicalIDs(MMFf, MMFw, MMFn, positions_grid, ids_array):
    """
    returns:
    
    f_ids: ids of particles in filaments
    w_ids: ids of particles in walls
    n_ids: ids of particles in nodes
    v_ids: ids of particles in voids
    """
    
    
    maskMMFf = np.where(MMFf>0, 1, 0)
    maskMMFw = np.where(MMFw>0, 1, 0)
    maskMMFn = np.where(MMFn>0, 1, 0)
    maskMMFv = np.where((MMFf == 0) & (MMFw == 0) & (MMFn == 0), 1, 0)
    
    maskf = generateMask(positions_grid, maskMMFf)
    maskw = generateMask(positions_grid, maskMMFw)
    maskn = generateMask(positions_grid, maskMMFn)
    maskv = generateMask(positions_grid, maskMMFv)
    
    f_ids = ids_array[maskf == 1]
    w_ids = ids_array[maskw == 1]
    n_ids = ids_array[maskn == 1]
    v_ids = ids_array[maskv == 1]
    
    return f_ids, w_ids, n_ids, v_ids

def determine_particle_label(id_val, filaments, walls, nodes, voids):
    if id_val in filaments:
        return 'f'
    elif id_val in walls:
        return 'w'
    elif id_val in nodes:
        return 'n'
    elif id_val in voids:
        return 'v'
    else:
        return None  # In case the ID is not in any array
    
    
def main():

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    
    if len(sys.argv) < 2:
        print("Usage: python df_environments_per_id.py snapshotnumber1 ...")
        sys.exit(1)
        
    snapnum = sys.argv[1]
    logger.info(f"Creating dataframe with particle IDs and their associated environment at snapshot number: {snapnum}")
    
    logger.info("Reading illustris snapshot")
    snap = il.snapshot.loadSubset(base_path, snapnum, 'dm', ["Coordinates", "ParticleIDs"])
    positions_array = snap["Coordinates"]
    ids_array = snap["ParticleIDs"]
    
    logger.info(f"Normalising illustris particle positions to fit to ({gridsize}, {gridsize}, {gridsize}) shaped obj")
    
    positions_grid = ((positions_array - np.min(positions_array)) / (np.max(positions_array) - np.min(positions_array)) * gridsize) - 1
    
    tempfile = f"/Users/users/nastase/PROJECT/DATA/nexus_outputs/MMF_outputs/output_{snapnum}.a_den"

    logger.info("Using readNEXUS to find environments...")
    denfield, MMFn, MMFf, MMFw = readNEXUS(tempfile,gridsize)

    logger.info("Obtaining IDs in each environment")
    filaments, walls, nodes, voids = getMorphologicalIDs(MMFf, MMFw, MMFn, positions_grid, ids_array)
    
    df = pd.DataFrame(index = ids_array)
    
    logger.info("Creating set of IDs for filaments")
    filament_ids_set = set(filaments)
    logger.info("Creating set of IDs for walls")
    walls_ids_set = set(walls)
    logger.info("Creating set of IDs for nodes")
    nodes_ids_set = set(nodes)
    logger.info("Creating set of IDs for voids")    
    voids_ids_set = set(voids)
    
    logger.info("Creating df with environment for each particle (ID as index)")
    df[f"particle_type_{snapnum}"] = df.index.map(lambda x: determine_particle_label(x, 
                                                filament_ids_set, 
                                                walls_ids_set, 
                                                nodes_ids_set, 
                                                voids_ids_set)
            )
    
    df_dir = f"/Users/users/nastase/PROJECT/DATA/nexus_outputs/particle_selecta/particle_types_dataframe_{snapnum}"
    logger.info(f"Saving df as pickle to {df_dir}")
    
    with open(f'{df_dir}.pickle', 'wb') as handle:
        pickle.dump(df, handle, protocol=4)
     
    logger.info(f"finished writing pickle file for {snapnum} at {df_dir}")
    
if __name__ == "__main__":
    main()
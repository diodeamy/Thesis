import sys
import pickle
import os
import glob
import illustris_python as il
import numpy as np

sys.path.append('..')

from shared.DTFE_utilities_nonparallel import DTFE, GRF, Zeldovich

    
L = 455
sigma = 8
gamma = 1
base_path =  "/Users/users/nastase/PROJECT/"

min_value = 0
max_value = 75000

def save_results(path, data):
    print(f"I'm saving the results and stuff to: {path}")
    filename = f"{path}.pickle"
    
    with open(filename, "wb") as f:
        pickle.dump(data, f)
        
def get_snapshot_number():
    while True:
        snapshot_number = input("Please enter the snapshot number for which you would like to find tesselations: ")
        
        if snapshot_number.isdigit():
            return snapshot_number
        else:
            print("Invalid input: please enter a numeric snapshot number")
            
def check_snapshot_exists(snapshot_number, base_path):
    pattern = f"{base_path}/snapdir_*{snapshot_number}"
    
    matching_dirs = glob.glob(pattern)
    
    if matching_dirs:
        print("Matching snapshot directory found:")
        for directory in matching_dirs:
            print(directory)
        return True
    else:
        print("No matching snapshot found in your directory! :(")
        return False
    
def get_xmin(min_value, max_value):
    while True:
        xmin = input("please state the lower 'X' bound for the area you want tesselated: ")
        
        if xmin.isdigit():
            xmin = int(xmin)
            if min_value <= xmin <= max_value:
                return xmin
            else:
                print(f"the bound is out of bounds! :)) please choose something between {min_value} - {max_value}.")
def get_xmax(min_value, max_value):
    while True:
        xmax = input("please state the upper 'X' bound for the area you want tesselated: ")
        
        if xmax.isdigit():
            xmax = int(xmax)
            if min_value <= xmax <= max_value:
                return xmax
            else:
                print(f"the bound is out of bounds! :)) please choose something between {min_value} - {max_value}.")                

            
def get_ymin(min_value, max_value):
    while True:
        ymin = input("please state the lower 'Y' bound for the area you want tesselated: ")
        
        if ymin.isdigit():
            ymin = int(ymin)
            if min_value <= ymin <= max_value:
                return ymin
            else:
                print(f"the bound is out of bounds! :)) please choose something between {min_value} - {max_value}.")

def get_ymax(min_value, max_value):
    while True:
        ymax = input("please state the upper 'X' bound for the area you want tesselated: ")
        
        if ymax.isdigit():
            ymax = int(ymax)
            if min_value <= ymax <= max_value:
                return ymax
            else:
                print(f"the bound is out of bounds! :)) please choose something between {min_value} - {max_value}.")
                
                
def get_zmin(min_value, max_value):
    while True:
        zmin = input("please state the lower 'Z' bound for the area you want tesselated: ")
        
        if zmin.isdigit():
            zmin = int(zmin)
            if min_value <= zmin <= max_value:
                return zmin
            else:
                print(f"the bound is out of bounds! :)) please choose something between {min_value} - {max_value}.")

def get_zmax(min_value, max_value):
    while True:
        zmax = input("please state the upper 'Z' bound for the area you want tesselated: ")
        
        if zmax.isdigit():
            zmax = int(zmax)
            if min_value <= zmax <= max_value:
                return zmax
            else:
                print(f"the bound is out of bounds! :)) please choose something between {min_value} - {max_value}.")
                
                
def load_datapoints(snapshot_number, xmin, xmax, ymin, ymax, zmin, zmax):
    dm_data = il.snapshot.loadSubset(base_path,snapshot_number,'dm',['Coordinates', 'Velocities'])
    dm_pos_all = dm_data['Coordinates']
    dm_vel_all = dm_data['Velocities']
    
    
    x_filter = (dm_pos_all[:,0] >= xmin) & (dm_pos_all[:,0] <= xmax)
    y_filter = (dm_pos_all[:,1] >= ymin) & (dm_pos_all[:,1] <= ymax)  
    z_filter = (dm_pos_all[:,2] >= zmin) & (dm_pos_all[:,2] <= zmax)    
    
    total_filter = x_filter & y_filter & z_filter
   
    dm_pos = dm_pos_all[total_filter].astype(np.float64)
    dm_vel = dm_vel_all[total_filter].astype(np.float64)
    
    return dm_pos, dm_vel


def main():
    
    snapshot_number = get_snapshot_number()
    check_snapshot_exists(snapshot_number, base_path)
    
    xmin = get_xmin(min_value, max_value)
    xmax = get_xmax(min_value, max_value)
    ymin = get_ymin(min_value, max_value)
    ymax = get_ymax(min_value, max_value)
    zmin = get_zmin(min_value, max_value)
    zmax = get_zmax(min_value, max_value)
    
    dm_pos, dm_vel = load_datapoints(snapshot_number, xmin, xmax, ymin, ymax, zmin, zmax)
        
    m = np.ones(len(dm_pos))
    
    dtfe = DTFE(dm_pos, dm_vel, m)
    
    save_results(f"/Users/users/nastase/PROJECT/DATA/Data/W6_DTFE/persisted_dtfe_{snapshot_number}", dtfe)
    
if __name__ == "__main__":
    main()
    
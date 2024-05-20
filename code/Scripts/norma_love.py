
import sys
import pickle

sys.path.append('..')

from shared.utilities import get_coordinates_for_particleIDs, get_redshift_dictionary, get_snapshot_cluster_coordinates


# def compute stuff(base_snapshot_id:int, subhalo_id:int, test_snapshot_ids:list[int]):
#     print("I'm busy computing... stuff")
        
#     redshift_snapshot_dict = get_redshift_dictionary()
#     snapshot_values = list(redshift_snapshot_dict.keys())
#     coordinates_n_0 = get_snapshot_cluster_coordinates()
    
def save_results(path, data, base_snapshot_id, subhalo_id):
    print(f"I'm saving the results and stuff to: {path}")
    filename = f"{path}/coordinates_base_snapshot_{base_snapshot_id}_subhalo_{subhalo_id}.pickle"
    
    with open(filename, "wb") as f:
        pickle.dump(data, f)
        
        
def main():
    
    #define some stuff
    redshift_snapshot_dict = get_redshift_dictionary()
    snapshot_values = list(redshift_snapshot_dict.keys())
#     coordinates_n_0 = get_snapshot_cluster_coordinates()
    
    #compute your stuff
    print("I'm busy computing... stuff")
    argument_sets = [
    {"base_snapshot_id": 49,
     "subhalo_id": 0,
     "test_snapshot_ids":snapshot_values
    },
    {"base_snapshot_id" :49,
     "subhalo_id": 30,
     "test_snapshot_ids":snapshot_values    
    },
    {"base_snapshot_id" :49,
     "subhalo_id":24 ,
     "test_snapshot_ids":snapshot_values    
    },
    {"base_snapshot_id" :49,
     "subhalo_id": 40,
     "test_snapshot_ids":snapshot_values    
    },
    {"base_snapshot_id" :49,
     "subhalo_id": 6,
     "test_snapshot_ids":snapshot_values    
    },
    {"base_snapshot_id" :135,
     "subhalo_id": 0,
     "test_snapshot_ids":snapshot_values    
    },
    {"base_snapshot_id" :135,
     "subhalo_id" : 574,
     "test_snapshot_ids":snapshot_values    
    },
    {"base_snapshot_id" :135,
     "subhalo_id" : 961,
     "test_snapshot_ids":snapshot_values    
    },
    {"base_snapshot_id" :135,
     "subhalo_id" : 1661,
     "test_snapshot_ids":snapshot_values    
    },
    {"base_snapshot_id" :135,
     "subhalo_id" : 1890,
     "test_snapshot_ids":snapshot_values    
    }
    ]
    
    pickle_results = {}
    
    for args in argument_sets:
        base_snapshot_id = args["base_snapshot_id"]
        subhalo_id = args["subhalo_id"]
        test_snapshot_ids = args["test_snapshot_ids"]
        
        coordinates = get_snapshot_cluster_coordinates(base_snapshot_id, subhalo_id, test_snapshot_ids)
    
        #save the results
        result_key = f"snapshot_{base_snapshot_id}_subhalo{subhalo_id}"
                
        pickle_results[result_key] = coordinates
        
        save_results("/Users/users/nastase/GitBub/Thesis/code/Data/W5_subhalo_movements", coordinates, base_snapshot_id, subhalo_id)
        
    all_results_path = "/Users/users/nastase/GitBub/Thesis/code/Data/W5_subhalo_movements/all_results.pickle"
    with open(all_results_path, "wb") as f:
        pickle.dump(pickle_results, f)
    
    print(f"All your stuff was saved to: {all_results_path}")
    
if __name__ == "__main__":
    main()
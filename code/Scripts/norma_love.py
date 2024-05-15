import json
from shared.utilities import get_coordinates_for_particleIDs, get_redshift_dictionary, get_snapshot_cluster_coordinates


# def compute stuff(base_snapshot_id:int, subhalo_id:int, test_snapshot_ids:list[int]):
#     print("I'm busy computing... stuff")
        
#     redshift_snapshot_dict = get_redshift_dictionary()
#     snapshot_values = list(redshift_snapshot_dict.keys())
#     coordinates_n_0 = get_snapshot_cluster_coordinates()
    
def save_results(path, data):
    print(f"I'm saving the results and stuff to: {path}")
    
    with open(f"{path}/data.json", "w") as f:
        json.dup(data, f)
        
def main():
    
    #define some stuff
    redshift_snapshot_dict = get_redshift_dictionary()
    snapshot_values = list(redshift_snapshot_dict.keys())
    coordinates_n_0 = get_snapshot_cluster_coordinates()
    
    #compute your stuff
    argument_sets = [
    {"base_snapshot_id": 49,
     "subhalo_id": 0,
     "test_snapshot_ids":snapshot_values}
#     },
#     {"base_snapshot_id" :49,
#      "subhalo_id": 30,
#      "test_snapshot_ids":snapshot_values    
#     },
#     {"base_snapshot_id" :49,
#      "subhalo_id":24 ,
#      "test_snapshot_ids":snapshot_values    
#     },
#     {"base_snapshot_id" :49,
#      "subhalo_id": 40,
#      "test_snapshot_ids":snapshot_values    
#     },
#     {"base_snapshot_id" :49,
#      "subhalo_id": 6,
#      "test_snapshot_ids":snapshot_values    
#     },
#     {"base_snapshot_id" :135,
#      "subhalo_id": 0,
#      "test_snapshot_ids":snapshot_values    
#     },
#     {"base_snapshot_id" :135,
#      "subhalo_id" : 574,
#      "test_snapshot_ids":snapshot_values    
#     },
#     {"base_snapshot_id" :135,
#      "subhalo_id" : 961,
#      "test_snapshot_ids":snapshot_values    
#     },
#     {"base_snapshot_id" :135,
#      "subhalo_id" : 1661,
#      "test_snapshot_ids":snapshot_values    
#     },
#     {"base_snapshot_id" :135,
#      "subhalo_id" : 1890,
#      "test_snapshot_ids":snapshot_values    
#     }
    ]
    
    #save the results
    json_results = {}
    
    
    
    
if __name__ = "__main__":
    main()
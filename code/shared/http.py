import requests

def get(path, params = None):
    headers = {"api-key":"a7c770df3fb83d41a97bb145b5c9eda4"}

    r = requests.get(path, params=params, headers=headers)
    
    r.raise_for_status()
    
    if r.headers['content-type'] == 'application/json':
        return r.json()

    if 'content-disposition' in r.headers:
        filename = r.headers['Content-Disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string
    
    return r

def get_coordinates_for_particleIDs(data, particle_ids):
    coordinates = data['Coordinates']
    ids = data['ParticleIDs']
    
    id_to_index = {particle_id: index for index, particle_id in enumerate(ids)}
    particle_indices = [id_to_index.get(particle_id) for particle_id in particle_ids]
    valid_indices = [index for index in particle_indices if index is not None]
    
    if len(valid_indices) != len(particle_indices):
        print("Warning: Some particle IDs were not found.")
    
    return coordinates[valid_indices]
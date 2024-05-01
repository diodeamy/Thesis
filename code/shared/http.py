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

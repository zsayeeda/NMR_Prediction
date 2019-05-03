

import os, sys
import requests
import urllib, json 
 
def AlatisCall(parentFolder):
    for dirName, subdirs, fileList in os.walk(parentFolder):
        print('Scanning %s...' % dirName)
        for filename in fileList:
            if filename.endswith( ".sdf" ):
                print filename
                URL = 'http://alatis.nmrfam.wisc.edu/upload'
                input_file_path = dirName+"/"+filename
                format_ = 'sdf'
                project_2_to_3 = 'on'
                add_hydrogens = 'on'
                files = {'infile': open(input_file_path, 'r')}
                data = {'format': format_, 'response_type': 'json', 'project_2_to_3': project_2_to_3, 'add_hydrogens': add_hydrogens }
                r = requests.post(URL, data=data, files=files)
                json_result = r.json()
                print(json_result['html_url'] )
                molecule_url = json_result['html_url']
                map_text_url = molecule_url+"/"+"map.txt"
                #response = urllib.urlopen(map_text_url)
                #data = response.read()
                f = requests.get(molecule_url)
                f = requests.get(map_text_url)
                #f = requests.post(map_text_url, data = data)
                #f_result = f.json()
                print("map url:{}".format(map_text_url))
                print(f.text)
                map_file = open(dirName+"/"+"map.txt", "w")
                map_file.write(f.text)
                map_file.close()









if __name__ == '__main__':
    if len(sys.argv) > 1:
        folders = sys.argv[1:]
        for i in folders:
            print("Folder name is :{}".format(i))
            if os.path.exists(i):
                AlatisCall(i)
            else:
                print('%s is not a valid path, please verify' % i)
                sys.exit()
    else:
        print('Usage: python FilenameFinder.py folder or python FilenameFinder.py folder1 folder2 folder3')

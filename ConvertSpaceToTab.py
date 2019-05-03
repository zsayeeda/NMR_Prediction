

import os, sys
import requests
import urllib, json 
 
def ToTab(parentFolder, file_name):
    for dirName, subdirs, fileList in os.walk(parentFolder):
        #print('Scanning %s...' % dirName)
        for filename in fileList:
            if filename == file_name:
                print("provided file name:{}".format(file_name))
                new_file_name = dirName+"/"+"tab_"+filename
                print("new writable file name:{}".format(new_file_name))
                new_file = open(new_file_name, "w")
                with open(dirName+"/"+filename, 'r') as fp:
                    text = fp.readlines()
                    new_file.write(text[0])
                    new_file.write(text[1])
                    for i in range(2,len(text)):
                        print("in shift file line:{} is:{}".format(i, text[i]))
                        if text[i].find(" ") > 0:
                            items = text[i].split(" ")
                            new_file.write(items[0]+"\t"+items[2])
                        else:
                            items = text[i].split("\t")
                            new_file.write(items[0]+"\t"+items[1])
            
                new_file.close()






if __name__ == '__main__':
    if len(sys.argv) > 1:
        print sys.argv
        folders = sys.argv[1]
        print("input folders:{}".format(folders))
        file = sys.argv[2]
        print("Folder name is :{}".format(folders))
        if os.path.exists(folders):
                ToTab(folders,file)
        else:
            print('%s is not a valid path, please verify' % folders)
            sys.exit()
    else:
        print('Usage: python ConvertSpaceToTab.py folder or python ConvertSpaceToTab.py  folder1 filename ')

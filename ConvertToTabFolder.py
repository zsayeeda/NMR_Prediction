

import os, sys
import requests
import urllib, json 
 
def ToTab(parentFolder):
    for dirName, subdirs, fileList in os.walk(parentFolder):
        #print('Scanning %s...' % dirName)
        for file in fileList:
            if file.endswith( ".txt" ):
                print ("name of the file is:{}".format(file))
                new_file_name = dirName+"/"+"tab_"+file
                print("new writable file name:{}".format(new_file_name))
                new_file = open(new_file_name, "w")
                with open(dirName+"/"+file, 'r') as fp:
                    text = fp.readlines()
                    new_file.write(text[0])
                    new_file.write(text[1])
                    for i in range(2,len(text)):
                        print("in shift file line:{} is:{}".format(i, text[i]))
                        if text[i].find(" ") > 0:
                            print("reached space")
                            items = text[i].split(" ")
                            len_of_items = len(items)
                            new_file.write(items[0]+"\t"+items[len_of_items-1])
                        else:
                            print("reached tab")
                            items = text[i].split("\t")
                            len_of_items = len(items)
                            new_file.write(items[0]+"\t"+items[len_of_items-1])
            
                new_file.close()






if __name__ == '__main__':
    if len(sys.argv) > 1:
        print sys.argv
        folders = sys.argv[1]
        print("input folders:{}".format(folders))
        if os.path.exists(folders):
                ToTab(folders)
        else:
            print('%s is not a valid path, please verify' % folders)
            sys.exit()
    else:
        print('Usage: python ConvertSpaceToTabFolder.py folder ')

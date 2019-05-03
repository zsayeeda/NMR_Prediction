

import os, sys
import requests
import urllib, json 
 
def AssignmentToShiftFile(parentFolder, file_name):
    for dirName, subdirs, fileList in os.walk(parentFolder):
        #print('Scanning %s...' % dirName)
        for filename in fileList:
            if filename == file_name:
                print("provided file name:{}".format(file_name))
                new_file_name = dirName+"/"+"new_"+filename
                print("new writable file name:{}".format(new_file_name))
                new_file = open(new_file_name, "w")
                with open(dirName+"/"+filename, 'r') as fp:
                    text = fp.readlines()
                    new_file.write(text[0])
                    new_file.write(text[1])
                    for i in range(2,len(text)):
                        #print("in shift file line:{} is:{}".format(i, text[i]))
                        items = text[i].split(" ")
                        atom_temp = int(items[0])
                        atom = atom_temp + 1
                        #print("atom after increasing to 1:{}".format(atom))
                        atom_in_assignment, shift_in_assignment = AssignmentParse(dirName+"/assignment.txt")
                        new_map, old_map = MapWithBMRB(dirName+"/map.txt")
                        index_of_assignment_atom = atom_in_assignment.index(atom)
                        index_of_new_map = old_map.index(atom)
                        chemical_shift = shift_in_assignment[index_of_assignment_atom]
                        mapped_atom = new_map[index_of_new_map]
                        new_file.write(str(atom_temp)+"\t"+chemical_shift+"\t"+str(mapped_atom)+"\n")
                new_file.close()




def AssignmentParse( file_name):
    atom = []
    shift =[]
    with open(file_name, 'r') as fp:
        text = fp.readlines()
        for i in range(1,len(text)):
            #print("in assignment file:{}".format(text[i]))
            items = text[i].split("\t")
            atom.append(int(items[0]))
            shift_temp = items[1].split("\n")
            shift.append(shift_temp[0])
    #print type(atom[0])
    #print shift
    return atom, shift

def MapWithBMRB(file_name):
    new_map = []
    old_map = []
    with open(file_name, 'r') as fp:
        text = fp.readlines()
        for i in range(1,len(text)):
            items = text[i].split("\t")
            new_map.append(int(items[0]))
            old_map_temp = items[2].split("\n")
            old_map.append(int(old_map_temp[0]))
    print("new:{}".format(new_map))
    print("old:{}".format(old_map))
    return new_map, old_map








if __name__ == '__main__':
    if len(sys.argv) > 1:
        #print sys.argv
        folders = sys.argv[1]
        print("input folders:{}".format(folders))
        file = sys.argv[2]
        print("Folder name is :{}".format(folders))
        if os.path.exists(folders):
                AssignmentToShiftFile(folders,file)
        else:
            print('%s is not a valid path, please verify' % folders)
            sys.exit()
    else:
        print('Usage: python FilenameFinder.py folder or python FilenameFinder.py folder1 filename ')

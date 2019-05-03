
#  To run this script use python FilenameFinder.py /folder1 ./folder2.

# FilenameFinder.py
import os, sys
 
def findDup(parentFolder):
    for dirName, subdirs, fileList in os.walk(parentFolder):
        print('Scanning %s...' % dirName)
        for filename in fileList:
            print filename


 
if __name__ == '__main__':
    if len(sys.argv) > 1:
        folders = sys.argv[1:]
        for i in folders:
            print("Folder name is :{}".format(i))
            if os.path.exists(i):
                findDup(i)
            else:
                print('%s is not a valid path, please verify' % i)
                sys.exit()
    else:
        print('Usage: python FilenameFinder.py folder or python FilenameFinder.py folder1 folder2 folder3')

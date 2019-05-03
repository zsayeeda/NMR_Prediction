
#  To run this script use python particular_line_finder.py /folder1 ./folder2.

# dupFinder.py
import os, sys

 
 
if __name__ == '__main__':
    print(" To run this script use python particular_line_finder.py /folder1 ./folder2.")
    if len(sys.argv) > 1:
        non_proton_nmr_files = {}
        folder = sys.argv[1]
        if os.path.exists(folder):
            for dirName, subdirs, fileList in os.walk(folder):
                print('Scanning %s...' % dirName)
                word = "ADDRESS"
                for filename in fileList:
                    #print("filename is :{}".format(filename))
                    path = os.path.join(dirName, filename)
                    print("filename with total path is:{}".format(path))
                    # Calculate hash. Check for particular line
                    f = open(path)
                    f.seek(0)
                    for line in f:
                        print(line)
                        if word in line:
                            non_proton_nmr_files["non_proton_nmr_file"] = path
                            os.remove(path)
                            print("file: {} has been deleted".format(path))
                            break
               

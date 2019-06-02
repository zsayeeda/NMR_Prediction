import csv
import os, sys
import operator





def GroupDataset(parentFolder):
  dataset_array = []
  for dirName, subdirs, fileList in os.walk(parentFolder):
     for filename in fileList:
      #if filename.endswith( ".txt" ):
      if filename.endswith( ".sdf" ):
        splitted_filename = filename.split(".sdf")
        hmdb_id = splitted_filename[0]
        print("in dataset hmdb_id:{}".format(hmdb_id))
        #with open("/Users/zinatsayeeda/anaconda3/envs/rdkit/class.csv", "r") as fp: #class.csv is the file I parsed using ruby. Parse is done in HMDB metabocad
        with open("/Users/zinatsayeeda/anaconda3/envs/rdkit/class_holdout.csv", "r") as fp:
          lines = fp.readlines()
          for line in lines:
            splitted_line = line.split("\t")
            hmdb_id_in_class = splitted_line[0]
            hmdb_id_class = splitted_line[1].split("\n")[0]
            if hmdb_id == hmdb_id_in_class:
              dataset_array.append({"id":hmdb_id, "klass": hmdb_id_class})
              print("inside dataset array:{}".format(dataset_array))
  #with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_2nd_priority/molecule_class.csv', 'w') as csv_file: #change path and file name here too
  with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/test_1/molecule_class_holdout.csv', 'w') as csv_file: #change path and file name here too
    writer = csv.writer(csv_file)
    for i in range(len(dataset_array)):
      data_row = []
      data_row.append(dataset_array[i]["id"])
      data_row.append(dataset_array[i]["klass"])
      writer.writerow(data_row)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        folders = sys.argv[1]
        print("input folders:{}".format(folders))
        if os.path.exists(folders):
                GroupDataset(folders)
        else:
            print('%s is not a valid path, please verify' % folders)
            sys.exit()
    else:
        print('Usage: python GroupDataset.py folder name with path')




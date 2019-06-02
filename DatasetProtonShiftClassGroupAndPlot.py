
# this script take the hmdb id  and chemical shift from my original training dataset. Then it takes corresponding taxonomy from my created csv file that is with
# molecule and class (already processed from the file ruby parse created). Altogether another csv file is created with hmdb id, chemical shift and taxonomy
#after that the graph is created to see how the classes are scattered




import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os,sys
import seaborn as sns
import csv



def DatasetProtonShiftClassGroupAndPlot():
  dataset_file = pd.read_csv("/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_csv/molecule_shift_class_merge.csv" , header = None)
  data_for_graph = pd.read_csv("/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_csv/molecule_shift_class_merge.csv" , sep=',')
  #print(" data is read and dataset file is:{}".format(dataset_file))
  y_data =  dataset_file.iloc[:,-2].values # this is the chemical shift values
  y_data = np.array(y_data)
  hmdb_ids =  dataset_file.iloc[:,-3].values # this is the chemical shift values
  hmdb_ids = np.array(hmdb_ids)
  group =  dataset_file.iloc[:,-1].values # this is the chemical group values
  group = np.array(group)
  print("group:{}".format(group))
  x_label= np.array(np.unique(y_data, return_counts=True)).T # x_label returns uniq values in the dataset depending on what u are looking for as uniq set
  print("uniq classes:{} and length:{}".format(x_label, len(x_label)))
  sns.set(font_scale=0.5)
  countplt=sns.countplot(x='CHEMICAL_SHIFT_CLASS', data=data_for_graph, palette ='hls')
  plt.title('Visualization of Chemical Shift Classes', fontsize=18)
  plt.xlabel('Chemical Shift Class', fontsize=16)
  plt.ylabel('Frequency', fontsize=16)
  #sns.scatterplot(x='GROUP', y='HMDB_ID',hue='CHEMICAL_SHIFT_CLASS',data=data_for_graph)
  plt.show()

  ##### this part is to take the HMDBID for uniq chemical shift values####
  dataset_array =[]
  for i in range(len(x_label)):
    hm_id =[]
    chemical_group =[]
    for j in range(len(y_data)):
      if x_label[i][0] == y_data[j]:
        hm_id.append(hmdb_ids[j])
        chemical_group.append(group[j])
    dataset_array.append({"chemical_shift_class":x_label[i][0], "HMDB_IDS": hm_id,"klass": chemical_group})
  with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_csv/shiftclass_id_group.csv', 'w') as csv_file: #change path and file name here too
    writer = csv.writer(csv_file)
    for i in range(len(dataset_array)):
      data_row = []
      data_row.append(dataset_array[i]["chemical_shift_class"])
      data_row.append(dataset_array[i]["HMDB_IDS"])
      data_row.append(dataset_array[i]["klass"])
      writer.writerow(data_row)









def GroupDataset(parent_file):
  dataset_array = []
  training_dataset = pd.read_csv("/Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_2nd_priority_megred.csv", sep=',')  #change path here
  training_dataset_hmdb_id = training_dataset.iloc[1:, -2]
  training_dataset_chemical_shift = training_dataset.iloc[1:, -3]
  training_dataset_hmdb_id = np.array(training_dataset_hmdb_id)
  training_dataset_chemical_shift = np.array(training_dataset_chemical_shift)
  training_dataset_chemical_shift = training_dataset_chemical_shift.astype(float)
  print("training dataset hmdb id:{}".format(training_dataset_hmdb_id))
  print("training dataset chemical shift:{}".format(training_dataset_chemical_shift))
  for i in range(len(training_dataset_hmdb_id)):
    with open(parent_file, "r") as fp:
      lines = fp.readlines()
      for j in range(1,len(lines)):
        splitted_line = lines[j].split(",")
        hmdb_id_in_class = splitted_line[0]
        hmdb_id_class = splitted_line[1].split("\r")[0]
        if training_dataset_hmdb_id[i] == hmdb_id_in_class:
          dataset_array.append({"id":training_dataset_hmdb_id[i], "chemical_shift_class": round(training_dataset_chemical_shift[i],2)*100,"klass": hmdb_id_class})
          print("inside dataset array:{}".format(dataset_array))
  with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_csv/molecule_shift_class_merge.csv', 'w') as csv_file: #change path and file name here too
    writer = csv.writer(csv_file)
    for i in range(len(dataset_array)):
      data_row = []
      data_row.append(dataset_array[i]["id"])
      data_row.append(int(dataset_array[i]["chemical_shift_class"]))
      data_row.append(dataset_array[i]["klass"])
      writer.writerow(data_row)





if __name__ == '__main__':
  if len(sys.argv) > 1:
    folders = sys.argv[1]
    print("input folders:{}".format(folders))
    if os.path.exists(folders):
      #GroupDataset(folders) # i called it first time and then using the created file. if you want to create new file, omit this comment
      DatasetProtonShiftClassGroupAndPlot()
    else:
      print('%s is not a valid path, please verify' % folders)
      sys.exit()
  else:
    print('Usage: python DatasetProtonShiftClassGroupAndPlot.py file name with path ') #give the file name I created with HMDB id and taxonomy already from the ruby parse(molecule_class_merge.csv)
    # like /Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_csv/molecule_class_merge.csv






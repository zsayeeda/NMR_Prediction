import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os,sys
import seaborn as sns
import csv

  
def ClassVsNoOfGroup(file_name):

  dataset_file = pd.read_csv(file_name , header = None)
  data_for_graph = pd.read_csv(file_name , sep=',')
  sns.set(font_scale=0.5)
  countplt=sns.catplot(x='chemical_shift_class', y= 'NO_GROUP', data=data_for_graph, kind ='bar')
  plt.title('Visualization of Chemical Shift Classes Vs No. of Chem Group', fontsize=10)
  plt.xlabel('Chemical Shift Class', fontsize=10)
  plt.ylabel('No. of the Chem Group', fontsize=10)
  #sns.scatterplot(x='GROUP', y='HMDB_ID',hue='CHEMICAL_SHIFT_CLASS',data=data_for_graph)
  plt.show()



if __name__ == '__main__':
  if len(sys.argv) > 1:
    folders = sys.argv[1]
    print("input folders:{}".format(folders))
    if os.path.exists(folders):
      #GroupDataset(folders) # i called it first time and then using the created file. if you want to create new file, omit this comment
      ClassVsNoOfGroup(folders)
    else:
      print('%s is not a valid path, please verify' % folders)
      sys.exit()
  else:
    print('Usage: python ClassVsNoOfGroup.py file name with path ') #give the file name I created with uniq chemical shift class with hmdb_id/class, no of com/class,
    # no of protorn/class. like /Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_csv/shiftclass_id_group.csv
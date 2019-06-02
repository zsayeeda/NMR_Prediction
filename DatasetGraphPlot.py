import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os,sys
import seaborn as sns



def DatasetGraphPlot(file_name):
  dataset_file = pd.read_csv("/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_csv/molecule_class_merge.csv" , header = None)
  #dataset_file = pd.read_csv("/Users/zinatsayeeda/anaconda3/envs/rdkit/test_1/molecule_class_holdout.csv" , header = None)
  data_for_graph = pd.read_csv("/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_csv/molecule_class_merge.csv" , sep=',')
  #data_for_graph = pd.read_csv("/Users/zinatsayeeda/anaconda3/envs/rdkit/test_1/molecule_class_holdout.csv" , sep=',')
  print(" data is read and dataset file is:{}".format(dataset_file))
  x_data = dataset_file.iloc[: , 0:1].values
  y_data =  dataset_file.iloc[:,-1].values
  x_data = np.array(x_data)
  y_data = np.array(y_data)
  x_label= np.array(np.unique(y_data, return_counts=True)).T
  print("uniq group:{} and length:{}".format(x_label, len(x_label)))
  sns.set(font_scale=0.5)
  #print(data_for_graph.head(n=3))
  countplt=sns.countplot(x='GROUP', data=data_for_graph, palette ='hls')
  plt.title('Type of Molecules in Dataset', fontsize=18)
  plt.xlabel('Chemical Group', fontsize=16)
  plt.ylabel('Frequency', fontsize=16)
  plt.show()



if __name__ == '__main__':
  if len(sys.argv) > 1:
    folders = sys.argv[1]
    print("input folders:{}".format(folders))
    if os.path.exists(folders):
      DatasetGraphPlot(folders)
    else:
      print('%s is not a valid path, please verify' % folders)
      sys.exit()
  else:
    print('Usage: python DatasetGraphPlot.py file name with path')
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import scale
import os, csv, sys



def pca(trainning_dataset):
  dataset = pd.read_csv(trainning_dataset , header = None)
  X = dataset.iloc[1: , 0:-3].values
  X = np.array(X)
  X = X.astype(float)
  threshold = 10**-3
  sd = np.std(X, axis =0) 
  X_new = X[:,sd >= threshold]
  print("shape of X_new:{}".format(X_new.shape))
  # print("first twenty samples:{}".format(X[1:21, [1,2,20]]))
  print("standard daviation of X:{}".format(sd))
  #print("shape of  sd:{}".format(sd.shape))
  X_new = scale(X_new)
  pca = PCA(n_components= X_new.shape[1])
  pca.fit(X_new)
  var= pca.explained_variance_ratio_
  var1=np.cumsum(np.round(pca.explained_variance_ratio_, decimals=4)*100)
  print var1
  plt.plot(var1)
  plt.show()

  ### after plotting variance I got n = 21 (95% variance i am going to calculate)




  pca = PCA(n_components=21)
  pca.fit(X_new)
  X1=pca.fit_transform(X)
  print X1







if __name__ == '__main__':
  if len(sys.argv) > 1:
    folders = sys.argv[1]
    print("input folders:{}".format(folders))
    if os.path.exists(folders):
      pca(folders)
    else:
      print('%s is not a valid path, please verify' % folders)
      sys.exit()
  else:
    print('Usage: python PCA.py file name with path') #python PCA.py /Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_2nd_priority_without_2D.csv


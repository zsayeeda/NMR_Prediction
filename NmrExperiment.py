
import os
import NmrStructure as NmrStructure
from rdkit import Chem
import  GetCDKDescriptors as GetCDKDescriptors
from sklearn.cross_validation import cross_val_score
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, auc
import pickle
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import cross_val_predict, cross_val_score, GridSearchCV
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.preprocessing import scale
import csv, sys


class NmrExperiment:
    #  This constructor function runs experiment on a training set doing 10-fold cross validation
    #    * measuring CV accuracy and validating the model
    saved_model = []
    best_param = None
    def __init__(self, dataset): # if dataset = True run classifier with CDK will be executed
        folder = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/"
        test_folder = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test/"
        #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/whole_training_nmr_1063_instance.csv" # used only those descriptors I have generated using CDK
        #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_priority.csv" # used only those descriptors I have generated using CDK
        #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_2nd_priority_megred.csv" # used only those descriptors I have generated using CDK
        #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_2nd_priority_without_2D.csv" # used only those descriptors I have generated using CDK
        trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_2nd_priority_all_atom.csv" # used only those descriptors I have generated using CDK
        #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_2nd_priority_megred_2.csv"
        #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/whole_training_nmr_3000_plus_instance.csv" # used only those descriptors I have generated using CDK
        # try:
        #     Instances isTrainingSet = (Instances) weka.core.SerializationHelper.read("models/train_classification_6") # fix it
        #     runClassifier(isTrainingSet, True)
        # except FileNotFoundError as fnf_error:
        #     print(fnf_error)
        #NmrStructure = NmrStruc.NmrStructures() 
        print("dataset value is:{}".format(dataset))
        if dataset == "False":
            print("executing False")
            trainingSet = self.buildTrainingClassification(folder)
            self.runClassifier(trainingSet, False)
            #self.runPrediction(test_folder)
        else:
            self.runClassifierWithCDK(trainingSetFile)
            print("executing True")

    def buildTrainingClassification(self, folder):
        feature_factor = 4
        nmr_structures = []  
        nmr_structures = self.getChemicalShifts(folder) # retruns array of nmr structures. assign carbon(H) position and chemical shifts in nmr structures. assign checmical shift classes also which is chemical_shift/100 to get 1000 classes
                                                        # this is the first time NmrStructure.py file is callled 
                                                        # also hmdbid is assigned for the structure           
        print("The NMR STRUCTURES after getting chemical shifts:{}".format(nmr_structures))
        print("\n\nAfter carbon position and chemical shift, now started taking H positions:\n")
        self.getStructures(nmr_structures, folder) # for each structure,H atom those are not connected with O and N are assigned and for all of those H atom position the HMDBIDs are returned
                                                   # do we really need it? cause no chemical shift is calculated here
        print("Length of nmrstructures{}".format(len(nmr_structures)))
        for nmr_str in nmr_structures:
            print("Afte getting structure, atomic descriptor and nearest atoms to H calculation starts for:{}".format(nmr_str.hmdb_id))
            nmr_str.atomic_descriptors = GetCDKDescriptors.getAtomicDescriptor(nmr_str.structure_sdf)
            nmr_str.findNearestAtomToHydrogens(GetCDKDescriptors.getNearestAtoms(nmr_str.structure_sdf))
        print("calculation for all the descriptors for all the files are complete")
        print("#######     dataset creation started #################\n\n")
        values = nmr_structures[0].atomic_descriptors[0]
        print("length of values is:{}".format(len(values)))
        #### The next block is For testing purpose. Delet/comment out the block once test is done ######
        for nmr_str in nmr_structures:
            print("print descriptors for all atoms in  file:{} is:{}".format(nmr_str.hmdb_id,nmr_str.atomic_descriptors))
        #######   Test Print Done ##########
        isTrainingSet = [[]]
        for nmr_str in nmr_structures:
            print("taining data calculation for the file:{}".format(nmr_str.hmdb_id))
            for i in range(len(nmr_str.hydrogen_positions)):
                h_position_in_nmr_str = nmr_str.hydrogen_positions[i]
                print("hydrogen position:{}".format(h_position_in_nmr_str))
                iExample = []
                print("length of atomic descriptos for file:{} is:{}".format(nmr_str.hmdb_id,len(nmr_str.atomic_descriptors[0])))
                print("Printing the arrays for nearest atoms for:{} {}".format(nmr_str.hmdb_id,nmr_str.nearest_atoms))
                #for j in range(len(nmr_str.atomic_descriptors)):
                for j in range(len(values)):
                    #h_position_in_nmr_str = nmr_str.hydrogen_positions[i]
                    #print("hydrogen position:{}".format(h_position_in_nmr_str))
                    iExample.append(nmr_str.atomic_descriptors[int(h_position_in_nmr_str)][j])
                for k in range(feature_factor - 1):
                    nearest_atom_for_all_h = nmr_str.nearest_atoms[1]
                    nearest_atoms_for_i_H_and_k_position = nearest_atom_for_all_h[i][k]
                    print("for {}th H {}th nearest atom is:{}".format(i,k,nearest_atoms_for_i_H_and_k_position))
                    #for j in range(len(nmr_str.atomic_descriptors)):
                    for j in range(len(values)):
                        print("length of atomic descriptor:{}".format(len(nmr_structures[0].atomic_descriptors[0])))
                        #nearest_ato = nmr_str.nearest_atoms[i][k]
                        #print("Printing the arrays for nearest atoms:{}".format(nmr_str.nearest_atoms))
                        #print("The length for nearest atoms:{}".format(len(nmr_str.nearest_atoms[0])))
                        #nearest_ato = nmr_str.nearest_atoms[int(h_position_in_nmr_str)][k]
                        # nearest_atom_for_all_h = nmr_str.nearest_atoms[0][i][k]
                        # nearest_ato = nmr_str.nearest_atoms[0][i][k]
                        #iExample.append(nmr_str.atomic_descriptors[int(nearest_ato)][j])
                        if j != 0:
                            iExample.append(nmr_str.atomic_descriptors[int(nearest_atoms_for_i_H_and_k_position)][j])
                #iExample.append(nmr_str.chemical_shifts[i]) # no need to use chemical shift as features
                #iExample.append(nmr_str.assignShiftClasses[i])
                iExample.append(nmr_str.c_shift_classes[i])
                #print("format of iexample is:{}\nand the length is:{}".format(iExample,len(iExample)))
                isTrainingSet.append(iExample)
                #print("format of training set:{}\n and length is:{}".format(isTrainingSet,len(isTrainingSet)))
        isTrainingSet.pop(0)
        print("The trainigset is:{}".format(isTrainingSet))
        return isTrainingSet


    def runClassifier(self, isTrainingSet, read):
        print("#########    Random forest model train started ##############\n")
        folds = 10
        true_values = [[]]
        predicted_values = [[]]
        #train = None
        #test = None
        hmdb_ids = self.mapInstancesToMolecules()
        print("hmdb ids:{}".format(hmdb_ids))
        dataset = np.array(isTrainingSet)
        print("This is how the dataset looks like after converting to np array:{}".format(dataset))
        X = dataset[:,1:-1]
        X = X.astype(float)
        print("X portion of dataset:{}".format(X))
        y =  dataset[:,-1]
        y = y.reshape(y.shape[0],1)
        y = y.astype(float)
        y = y/100 # creating 1000 classes
        print("y portion of dataset:{} and shape is".format(y,y.shape))
        y_temp = y[:]
        y_temp = y_temp.ravel()
        print("after ravel, y_temp shape is :{}".format(y_temp.shape))
        factor = pd.factorize(y_temp)
        Y = factor[0]
        definitions = factor[1]
        reversefactor = dict(zip(range(len(factor[0])),definitions))
        print("reverse factor looks like:{}".format(reversefactor))
        Y = Y.reshape(Y.shape[0],1)
        print("after factorization  Y data is :{} and shape is {}".format(Y,Y.shape))

        '''
        # implementation of PCA #####
        sd = np.std(X, axis =0) 
        X_new = X[:,sd >= threshold]
        pca = PCA(n_components= X_new.shape[1])
        pca.fit(X_new)
        pca = PCA(n_components=21)
        pca.fit(X_new)
        X1=pca.fit_transform(X) # X1 is the new dataset now
        #### implementation of PCA is done
        '''



        cv = StratifiedKFold(n_splits=folds, random_state=123, shuffle=True)
        fprs, tprs, scores = [], [], []
        model = RandomForestClassifier()
        average_error = 0
        outlier_index = 0
        for (train_index, test_index), i in zip(cv.split(X, Y), range(folds)):
        #for (train_index, test_index), i in zip(cv.split(X1, Y), range(folds)): # before it was X which was without PCA
            #print("checking value of i in fold:{}".format(i))
            X_train, X_test = X[train_index], X[test_index]
            #print(" X_train is:{} and shape is:{}\n X_test is:{} and shape is:{}".format(X_train,X_train.shape,X_test,X_test.shape))
            y_train, y_test = Y[train_index], Y[test_index] # y_test is the true value of prediction
            #print(" y_train is:{} and shape is:{}\n y_test is:{} and shape is:{}".format(y_train,y_train.shape,y_test,y_test.shape))
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            #print(" y_pred is:{} and shape:{}".format(y_pred,y_pred.shape))
            y_pred = np.vectorize(reversefactor.get)(y_pred)
            y_pred = 100 * y_pred 
            y_pred = y_pred.reshape(y_pred.shape[0],1)
            #print(" y_pred is:{} and shape:{}".format(y_pred,y_pred.shape))
            y_test = np.vectorize(reversefactor.get)(y_test)
            y_test = 100 * y_test
            true_values.append(y_test)
            predicted_values.append(y_pred)
            error = 0
            outlier_num = 0
            test_arbitary = 0
            for j in range(len(y_pred)):
                #if not((abs(y_pred[j] - y_test[j]) < 2.0)):
                if ((abs(y_pred[j] - y_test[j]) < 8.00)):
                    #print("no outlier in {}th prediction".format(j))
                    test_arbitary = test_arbitary +1
                else:
                    outlier_num = outlier_num + 1
                    print("true values of the outlier:{}".format(y_test[j]))
                    print("HMDB ID of the outlier is:{}".format(hmdb_ids[outlier_index]))
                error = abs(y_pred[j] - y_test[j]) + error
                outlier_index = outlier_index + 1
                #print("In fold i:{} y_pred is:{}    and   y_test:{}".format(i,y_pred[j],y_test[j]))

            error = error / len(y_pred)
            print("avg error in {}th fold:{}".format(i,error))  
            print("total number of outliers in {} fold:{}".format(i,outlier_num))
            self.saved_model.append({"error": error, "model": model})
            average_error = error + average_error

        average_error = average_error/folds
        print("after all fold the average error:{}".format(average_error))
        #print(self.saved_model)
        self.saved_model.sort()
        #print("\nafter sort:{}".format(self.saved_model))
    
        #### plot all the values

   
    def mapInstancesToMolecules(self):
        folder = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/"
        nmr_structures = []
        hmdb_ids = []
        nmr_structures = self.getChemicalShifts(folder)
        self.getStructures(nmr_structures, folder)
        for nmr_str in nmr_structures:
            for f in nmr_str.chemical_shifts:
                hmdb_ids.append(nmr_str.hmdb_id)
        return hmdb_ids

    
  
   
    def createNmrStructures(self, folder):
        structures = []
        print("reached create structure")
        for fileEntry in os.listdir(folder):
            print("file entry in create structure:{}".format(fileEntry))
            if fileEntry.endswith( ".sdf" ):
                print("sending to nmrstructure:{}".format(fileEntry.split(".")[0]))
                structure = NmrStructure(None,None,None,fileEntry.split(".")[0])
                structures.append(structure)
        print("structure creation done")
        return structures

    
    def getChemicalShifts(self, folder):
        structures = []
        file_names = []
        for fileEntry in os.listdir(folder):
            if fileEntry.endswith( ".txt" ):
                print("File entry in getchemicalshift method:{}".format(fileEntry))
                carbon_position = []
                chemical_shift = []
                file_names.append(fileEntry)
                self.readChemicalShifts(folder+fileEntry, carbon_position, chemical_shift) #return carbon position and chemical shifts. Here the carbon positions are actually H atom position
                structure = NmrStructure.NmrStructure(carbon_position, chemical_shift, "", fileEntry) # no need fileEntry.split(".")[0]. Here the carbon positions are actually H atom position
                structure.assignShiftClasses(structure.chemical_shifts)
                structures.append(structure)
        return structures

    
    def readChemicalShifts(self, file, c_pos, c_shift): 
        #reader = None
        try:
            print("start with file: {} to calculate carbon position and chemical shift from readchemicalshift method".format(file))
            with open(file, 'r') as fp:
                text = fp.readlines()
            parse_shifts = False
            for i in range(len(text)):
                items = text[i].split("\t")
                if parse_shifts:
                    if text[i] == "":
                        continue 
                    if "Table" in text[i]:
                        parse_shifts = False
                        continue 
                    #print("The file with file path: {}".format(file))
                    if items[0].strip() == "No." or ("M" in items[2]):
                        continue 
                    c_pos.append(items[1].strip())
                    if items[2] != None and items[2] != "":
                        c_shift.append(items[2].strip())
                    else:
                        c_shift.append(items[3].strip())
                    #print("The file with file path: {}".format(file))
                if len(items) == 0 or len(items) == 1:
                    continue 
                if "Atom" in items[1].strip():
                    parse_shifts = True
            print("carbon position:{}".format(c_pos))
            print("chemical shift:{}".format(c_shift))
        except:
            print("Some error occured")

    
    def getStructures(self, structures, folder):
        hmdb_ids = [[]]
        for nmr_s in structures:
            file = folder+nmr_s.hmdb_id+".sdf"
            print("\nStarted to get H atoms from the getstructures method, from the file:{}".format(file))
            try:
                text = Chem.SDMolSupplier(file)
                print("The mols:{}".format(text))
                nmr_s.structure_sdf = text  # the mols are sent as input
                print("mols inside nmrstructure:{}".format(nmr_s.structure_sdf))
                print("Is chemical shift true/false?:{}".format(nmr_s.has_chemical_shift))
                #if nmr_s.has_chemical_shift == False: # should it be false??? i don't think so. need to check
                if nmr_s.has_chemical_shift == False:
                    Print("Now taking H position from getHydrogenAtoms method")
                    nmr_s.hydrogen_positions = GetCDKDescriptors.getHydrogenAtoms(nmr_s.structure_sdf) #The mols are sent as input.H atom's position number's are returned where H atoms are connected only with C atoms
                    nmr_s.has_chemical_shift = True
                print("The hydrogen positions in nmrstructure:{}".format(nmr_s.hydrogen_positions))
                for s in nmr_s.hydrogen_positions:
                    hmdb_id = []
                    print("hmdb  id is:{}".format(nmr_s.hmdb_id))
                    hmdb_id_for_atom = nmr_s.hmdb_id
                    #hmdb_ids.append(nmr_s.hmdb_id)
                    hmdb_id.append(hmdb_id_for_atom)
                    Print("HMDBID addition done")
                hmdb_ids.append(hmdb_id)
            except:
                print("Error Occured")
        hmdb_ids.pop(0)
        return hmdb_ids

# This function is for testing purpose to check what result the random forest provides with the training set I have created seperately using CDK
    def runClassifierWithCDK(self, trainingSetFile):
        print("#########    Random forest model train started ##############\n")
        folds = 5
        true_values = [[]]
        predicted_values = [[]]
        #train = None
        #test = None
    #    hmdb_ids = self.mapInstancesToMolecules()
    #    print("hmdb ids:{}".format(hmdb_ids))
    #    dataset = np.array(isTrainingSet)
        dataset = pd.read_csv(trainingSetFile , header = None)
        print("This is how the dataset looks like after converting to np array:{}".format(dataset))
    #    X = dataset[1:,1:-1]
        hmdb_ids = dataset.iloc[1:,-2].values # it was -1 before for xuan's dataset
        hmdb_ids = hmdb_ids.reshape(hmdb_ids.shape[0],1)
        hydrogen_positions = dataset.iloc[1:,-1].values
        hydrogen_positions = hydrogen_positions.reshape(hydrogen_positions.shape[0],1)
        X = dataset.iloc[1: , 0:-3].values #it was 0:-1 before]
        X = X.astype(float)
        print("X portion of dataset:{}".format(X))
    #    y =  dataset[:,-1]
        y =  dataset.iloc[1:,-3].values # it was -1 before
        y = y.reshape(y.shape[0],1)
        y = y.astype(float)
    #    y = y/100
        print("y portion of dataset:{} and shape is".format(y,y.shape))
    ##    y_temp = y[:]
    ##    y_temp = y_temp.ravel()
    ##    print("after ravel, y_temp shape is :{}".format(y_temp.shape))
    ##    factor = pd.factorize(y_temp)
    ##    Y = factor[0]
    ##   definitions = factor[1]
    ##    reversefactor = dict(zip(range(len(factor[0])),definitions))
    ##    Y = Y.reshape(Y.shape[0],1)
        Y = self.classLabel(y)
        print("Y after classLabel:{}".format(Y))
        Y = Y.reshape(Y.shape[0],1)
        Y = Y.astype(int)
    ##    print("after factorization  Y data is :{} and shape is {}".format(Y,Y.shape))
        print("after converting to class label  Y data is :{} and shape is {}".format(Y,Y.shape))

        '''
        #### grid searcg start #### ignore it for now
        print("####grid search started ######")
        Y_temp = Y
        gsc = GridSearchCV(estimator=RandomForestClassifier(),param_grid={'max_depth': range(3,10),'n_estimators': (10, 50, 100, 500, 1000),},cv=5, scoring='neg_mean_absolute_error', verbose=0,n_jobs=-1)
        grid_result = gsc.fit(X, Y_temp.ravel())
        print("gsc done")
        best_params = grid_result.best_params_
        self.best_param = best_params
        print("best_params={}".format(self.best_param))
        #p_result_from_gsc = grid_result.predict(X)
        #print("gsc predicted result:{} and type:{}".format(p_result_from_gsc,type(p_result_from_gsc)))
        cv_results = grid_result.cv_results_.keys()
        print("gsc cv results keys:{}".format(cv_results))
        print("gsc cv results :{}".format(grid_result.cv_results_))
        print("split0_train_score:{}         split0_test_score:{}".format(cv_results['split0_train_score'],cv_results['split0_test_score']))
        print("split1_train_score:{}         split1_test_score:{}".format(cv_results['split1_train_score'],cv_results['split1_test_score']))
        print("split2_train_score:{}         split2_test_score:{}".format(cv_results['split2_train_score'],cv_results['split2_test_score']))
        print("split3_train_score:{}         split3_test_score:{}".format(cv_results['split3_train_score'],cv_results['split3_test_score']))
        print("split4_train_score:{}         split4_test_score:{}".format(cv_results['split4_train_score'],cv_results['split4_test_score']))
        print("std_train_score:{}         std_test_score:{}".format(cv_results['std_train_score'],cv_results['std_test_score']))
        print("mean_train_score:{}         mean_test_score:{}".format(cv_results['mean_train_score'],cv_results['mean_test_score']))
        print("param_max_depth:{}         param_n_estimators:{}".format(cv_results['param_max_depth'],cv_results['param_n_estimators']))
        print("params:{}         rank_test_score:{}".format(cv_results['params'],cv_results['rank_test_score']))
        print("####grid search ends ######")

        #### grid search end ###
        '''


        '''
        # implementation of PCA #####
        threshold = 10**-3
        sd = np.std(X, axis =0) 
        X_new = X[:,sd >= threshold]
        pca = PCA(n_components= X_new.shape[1])
        pca.fit(X_new)
        pca = PCA(n_components=21)
        pca.fit(X_new)
        X1=pca.fit_transform(X_new) # X1 is the new dataset now
        #### implementation of PCA is done
        '''
        
        



        #cv = StratifiedKFold(n_splits=folds, random_state=123, shuffle=True)
        cv = StratifiedKFold(n_splits=folds, random_state=None, shuffle=False)
        fprs, tprs, scores = [], [], []
        #model = RandomForestClassifier(max_depth=self.best_param["max_depth"], n_estimators=self.best_param["n_estimators"])
        model = RandomForestClassifier()
        average_error = 0
        outlier_index = 0
        for (train_index, test_index), i in zip(cv.split(X, Y), range(folds)): #change here. uncommit it when PCA is not used
        #for (train_index, test_index), i in zip(cv.split(X1, Y), range(folds)): # it was X before implementing the PCA analysis
            #print("checking value of i in fold:{}".format(i))
            X_train, X_test = X[train_index], X[test_index]
            #print(" X_train is:{} and shape is:{}\n X_test is:{} and shape is:{}".format(X_train,X_train.shape,X_test,X_test.shape))
            y_train, y_test = Y[train_index], Y[test_index] # y_test is the true value of prediction
            #print(" y_train is:{} and shape is:{}\n y_test is:{} and shape is:{}".format(y_train,y_train.shape,y_test,y_test.shape))
            hmdb_ids_in_each_fold_train = hmdb_ids[train_index]
            hmdb_ids_in_each_fold_test = hmdb_ids[test_index]
            hydrogen_positions_in_each_fold_train = hydrogen_positions[train_index]
            hydrogen_positions_in_each_fold_test = hydrogen_positions[test_index]
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            #print(" y_pred is:{} and shape:{}".format(y_pred,y_pred.shape))
    ##        y_pred = np.vectorize(reversefactor.get)(y_pred)
            y_pred = y_pred.astype(float)
            #y_pred = y_pred/float(10) # for 100 class
            y_pred = y_pred/float(100)
        #    y_pred = 10 * y_pred 
    ##        y_pred = y_pred.reshape(y_pred.shape[0],1)
            #print(" y_pred is:{} and shape:{}".format(y_pred,y_pred.shape))
    ##        y_test = np.vectorize(reversefactor.get)(y_test)
        #    y_test = 10 * y_test
            y_test = y_test.astype(float)
            #y_test = y_test/float(10) # for 100 class
            y_test = y_test/float(100)
            true_values.append(y_test)
            predicted_values.append(y_pred)
            error = 0
            outlier_num = 0
            test_arbitary = 0
            meanAbsError = mean_absolute_error(y_test, y_pred) # from sklearn
            print("Mean Absolute Error from fold:{} is:{}".format(i,meanAbsError))
            for j in range(len(y_pred)):
                #print("in fold {} y_true = {} and y_pred ={}".format(i,y_test[j],y_pred[j]))
                if ((abs(round(y_pred[j],2) - y_test[j]) < 8.0)): # it was 20.0. don't know whether that was correct
    ##            if ( (abs(round(y_pred[j],3) - round(y_test[j],3)) < 20)):
                    #print("no outlier in {}th prediction".format(j))
                    test_arbitary = test_arbitary +1
                else:
                    outlier_num = outlier_num + 1
                    print("true values of the outlier:{} and predicted value of the outlier is:{}, HMDB ID is{} and H position is:{}".format(y_test[j],round(y_pred[j],2),hmdb_ids_in_each_fold_test[j],hydrogen_positions_in_each_fold_test[j]))
                    #print("HMDB ID of the outlier is:{}".format(hmdb_ids[outlier_index]))
                error = abs(round(y_pred[j],2) - y_test[j]) + error
                outlier_index = outlier_index + 1
                #print("In fold i:{} y_pred is:{}    and   y_test:{}".format(i,y_pred[j],y_test[j]))

            error = error / len(y_pred)
            print("avg error in {}th fold:{}".format(i,error))  
            print("total number of outliers in {} fold:{}".format(i,outlier_num))
            self.saved_model.append({"error": error, "model": model})
            average_error = error + average_error

        average_error = average_error/folds
        print("after all fold the average error:{}".format(average_error))
        #print(self.saved_model)
        self.saved_model.sort()
        #print("\nafter sort:{}".format(self.saved_model))


    def classLabel(self, y_value):
        print("received y_value is:{}".format(y_value))
        cls = []
        #for i in range(1,101): # for 100 class
        for i in range(1,1001): #for 1000 class
            #cls.append({"value": float(i)/float(10), "class": i}) # for 100 class
            cls.append({"value": float(i)/float(100), "class": i})
        print("array of has for cls is :{}".format(cls))
        y_value_temp = y_value[:]
        print("y_value_temp = {}".format(y_value_temp))
        y_value_classlabel = np.zeros(len(y_value_temp))
        print("initially y_value_classlabel is :{}".format(y_value_classlabel))
        for i in range(len(y_value_temp)):
            y_value_temp[i] = round(y_value_temp[i],2) # changed it from 1 to 2
        cls_value = [cls[i]["value"] for i in range(len(cls))]
        print("cls_value:{}".format(cls_value))
        cls_class = [cls[i]["class"] for i in range(len(cls))]
        print("cls_class:{}".format(cls_class))
        cls_value_tuple = tuple(cls_value)
        print("cls_value_tuple:{}".format(cls_value_tuple))
        cls_class_tuple = tuple(cls_class)
        print("cls_class_tuple:{}".format(cls_class_tuple))
        cls_class_value_dict = dict(zip(cls_class_tuple,cls_value_tuple))
        print("cls_class_value_dict:{}".format(cls_class_value_dict))
        #y_value_classlabel = [k for k, v in cls_class_value_dict.items() if v == y_value_temp[i]]
        for i in range(len(y_value_temp)):
            for k, v in cls_class_value_dict.items():
                #print("k = {} and v={}".format(k,v))
                if v == y_value_temp[i]:
                    #print("y_value_temp_{} = {}".format(i,y_value_temp[i]))
                    y_value_classlabel[i] = k
                    #print("y_value_classlabel_{} = {}".format(i,y_value_classlabel[i]))
                    break
                else:
                    #y_value_classlabel[i] = 101.0 # for 100 class
                    y_value_classlabel[i] = 1001.0 
        print("y value has been transferred into class:{}".format(y_value_classlabel)) 
        return y_value_classlabel



    # def buildTestClassification(self, folder, nmr_structures):
    #     """ generated source for method buildTestClassification """
    #     feature_factor = 4
    #     try:
    #         for nmr_str in nmr_structures:
    #             print("test file id:{}".format(nmr_str.hmdb_id))
    #             nmr_str.atomic_descriptors = GetCDKDescriptors.getAtomicDescriptor(nmr_str.structure_sdf)
    #             nmr_str.findNearestAtomToHydrogens(GetCDKDescriptors.getNearestAtoms(nmr_str.structure_sdf))
    #     except:
    #         print("error in calculation of H position, nearest atom position in test")
    #     values = nmr_structures[0].atomic_descriptors[0]
    #     print("length of values is:{}".format(len(values)))
    #     #### The next block is For testing purpose. Delet/comment out the block once test is done ######
    #     for nmr_str in nmr_structures:
    #         print("print descriptors for all atoms in  file:{} is:{}".format(nmr_str.hmdb_id,nmr_str.atomic_descriptors))
    #     #######   Test Print Done ##########
    #     isTestSet = [[]]
    #     for nmr_str in nmr_structures:
    #         print("taining data calculation for the file:{}".format(nmr_str.hmdb_id))
    #         for i in range(len(nmr_str.hydrogen_positions)):
    #             h_position_in_nmr_str = nmr_str.hydrogen_positions[i]
    #             print("hydrogen position:{}".format(h_position_in_nmr_str))
    #             iExample = []
    #             print("length of atomic descriptos for file:{} is:{}".format(nmr_str.hmdb_id,len(nmr_str.atomic_descriptors[0])))
    #             print("Printing the arrays for nearest atoms for:{} {}".format(nmr_str.hmdb_id,nmr_str.nearest_atoms))
    #             #for j in range(len(nmr_str.atomic_descriptors)):
    #             for j in range(len(values)):
    #                 #h_position_in_nmr_str = nmr_str.hydrogen_positions[i]
    #                 #print("hydrogen position:{}".format(h_position_in_nmr_str))
    #                 iExample.append(nmr_str.atomic_descriptors[int(h_position_in_nmr_str)][j])
    #             for k in range(feature_factor - 1):
    #                 nearest_atom_for_all_h = nmr_str.nearest_atoms[1]
    #                 nearest_atoms_for_i_H_and_k_position = nearest_atom_for_all_h[i][k]
    #                 print("for {}th H {}th nearest atom is:{}".format(i,k,nearest_atoms_for_i_H_and_k_position))
    #                 #for j in range(len(nmr_str.atomic_descriptors)):
    #                 for j in range(len(values)):
    #                     print("length of atomic descriptor:{}".format(len(nmr_structures[0].atomic_descriptors[0])))
    #                     #nearest_ato = nmr_str.nearest_atoms[i][k]
    #                     #print("Printing the arrays for nearest atoms:{}".format(nmr_str.nearest_atoms))
    #                     #print("The length for nearest atoms:{}".format(len(nmr_str.nearest_atoms[0])))
    #                     #nearest_ato = nmr_str.nearest_atoms[int(h_position_in_nmr_str)][k]
    #                     # nearest_atom_for_all_h = nmr_str.nearest_atoms[0][i][k]
    #                     # nearest_ato = nmr_str.nearest_atoms[0][i][k]
    #                     #iExample.append(nmr_str.atomic_descriptors[int(nearest_ato)][j])
    #                     if j != 0:
    #                         iExample.append(nmr_str.atomic_descriptors[int(nearest_atoms_for_i_H_and_k_position)][j])
    #             #iExample.append(nmr_str.chemical_shifts[i])
    #             #iExample.append(nmr_str.assignShiftClasses[i])
    #             iExample.append(nmr_str.c_shift_classes[i])
    #             #print("format of iexample is:{}\nand the length is:{}".format(iExample,len(iExample)))
    #             isTestSet.append(iExample)
    #             #print("format of training set:{}\n and length is:{}".format(isTrainingSet,len(isTrainingSet)))
    #     isTestSet.pop(0)
    #     print("The trainigset is:{}".format(isTrainingSet))
    #     return isTestSet


    # def runPrediction(self,folder):
    #     print("##################Prediction for holdout test starts here###############\n")
    #     try:
    #         # LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
    #         # RandomForest modesl = (RandomForest) weka.core.SerializationHelper.read("models/classification.model_1");
    #         structures = self.createNmrStructures(folder)
    #         print("structures created in holdout")
    #         hmdb_ids =  self.getStructures(structures, folder)
    #         structures = self.getChemicalShifts(folder)
    #         print("chemical shift calculated in holdout")
    #         #model_saved_in_train_phase = exp.saved_model[6]
    #         model_saved_in_train_phase = saved_model[6]
    #         print("holdout section took error:{}  and model:{}".format(model_saved_in_train_phase["error"],model_saved_in_train_phase["model"]))
    #         model = model_saved_in_train_phase["model"]
    #         test = self.buildTestClassification(folder, structures)
    #         dataset = np.array(test)
    #         X = dataset[:,1:-1]
    #         X = X.astype(float)
    #         print("X portion of dataset:{}".format(X))
    #         y =  dataset[:,-1]
    #         y = y.reshape(y.shape[0],1)
    #         y = y.astype(float)
    #         print("y portion of dataset:{} and shape is".format(y,y.shape))
    #         y_temp = y[:]
    #         y_temp = y_temp.ravel()
    #         factor = pd.factorize(y_temp)
    #         Y = factor[0]
    #         Y = Y.reshape(Y.shape[0],1)
    #         definitions = factor[1]
    #         reversefactor = dict(zip(range(len(factor[0])),definitions))
    #         Y = Y.reshape(Y.shape[0],1)
    #         print("after factorization  Y data is :{} and shape is {}".format(Y,Y.shape))
    #         average_error = 0
    #         outlier_index = 0
    #         X_test = X
    #         y_test = Y
    #         y_pred = model.predict(X_test)
    #         y_pred = np.vectorize(reversefactor.get)(y_pred)
    #         y_pred = 100 * y_pred 
    #         y_pred = y_pred.reshape(y_pred.shape[0],1)
    #         y_test = np.vectorize(reversefactor.get)(y_test)
    #         y_test = 100 * y_test
    #         error = 0
    #         outlier_num = 0
    #         test_arbitary = 0
    #         for j in range(len(y_pred)):
    #             if ((abs(y_pred[j] - y_test[j]) < 20)):
    #                 #print("no outlier in {}th prediction".format(j))
    #                 test_arbitary = test_arbitary +1
    #             else:
    #                 outlier_num = outlier_num + 1
    #                 print("true values of the outlier:{}".format(y_test[j]))
    #                 print("HMDB ID of the outlier is:{}".format(hmdb_ids[outlier_index]))
    #             error = abs(y_pred[j] - y_test[j]) + error
    #             outlier_index = outlier_index + 1
    #         error = error / len(y_pred)
    #         print("avg error in holdout test:{}".format(error))  
    #         print("total number of outliers:{}".format(outlier_num))
    #     except:
    #         print("error occur in hold out test")

    



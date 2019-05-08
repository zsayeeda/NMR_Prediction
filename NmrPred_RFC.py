
import os
import NmrExperiment as NmrExp
#import matlabcontrol
import sys
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import roc_curve, auc
import pickle
from sklearn.metrics import mean_absolute_error


# if __name__ == '__main__':
    
#     folder = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/"
#     testData = "/Users/zinatsayeeda/anaconda3/envs/rdkit/whole_test_nmr.csv"
#     print("Start")
#     print("Start with:{}".format(sys.argv[1]))
#     exp = NmrExp.NmrExperiment(sys.argv[1])
#     #runPrediction(folder);
#     runPredictionWithCDK(testData)
#     print("End")


    #  RunPrediction function builds a classifier from a saved Training Set. It then builds a testing instance
    #    * from a test folder. It then runs predict function to calculate ppm values for all Hydrogen atoms in 
    #    * for all the test molecules
    #    
    # def runPrediction(self, folder):
    #     """ generated source for method runPrediction """
    #     try:
    #         # LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
    #         # RandomForest modesl = (RandomForest) weka.core.SerializationHelper.read("models/classification.model_1");
    #         structures = NmrExperiment.createNmrStructures(folder)
    #         hmdb_ids = NmrExperiment.getStructures(structures, folder)
    #         model.buildClassifier(isTrainingSet)
    #         for nmr_str in structures:
    #             for h_pos in nmr_str.hydrogen_positions:
    #                 #  Class is a double from 1 to 100 classes. So divide by 10 to get shift
    #                 nmr_str.chemical_shifts.add(Float.valueOf(label))
    #                 print hmdb_ids.get(i) + "\t" + label
    #                 i += 1
    #         callMatlab(structures)
    #     except Exception as e:
    #         e.printStackTrace()

    # @classmethod
    # def callMatlab(cls, structures):
    #     """ generated source for method callMatlab """
    #     options = MatlabProxyFactoryOptions.Builder().setUsePreviouslyControlledSession(True).build()
    #     factory = MatlabProxyFactory(options)
    #     proxy = factory.getProxy()
    #     proxy.eval("cd ..")
    #     proxy.eval("cd matlab")
    #     for nmr_str in structures:
    #         proxy.eval("inter.coupling.scalar = []")
    #         while i < len(nmr_str.hydrogen_positions):
    #             if i == len(nmr_str.hydrogen_positions) - 1:
    #                 isotopes = isotopes + "'1H'}"
    #                 ppms = ppms + String.valueOf(nmr_str.chemical_shifts.get(i)) + "}"
    #                 proxy.eval(coupling)
    #             else:
    #                 isotopes = isotopes + "'1H', "
    #                 ppms = ppms + String.valueOf(nmr_str.chemical_shifts.get(i)) + ", "
    #             while j < len(nmr_str.hydrogen_positions):
    #                 proxy.eval(coupling)
    #                 j += 1
    #             i += 1
    #         print ppms
    #         proxy.eval(isotopes)
    #         proxy.eval(ppms)
    #         proxy.feval("create_nmr1H_plot")
    #     proxy.disconnect()

def runPrediction(folder):
       # """ generated source for method runPrediction """
    try:
            # LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
            # RandomForest modesl = (RandomForest) weka.core.SerializationHelper.read("models/classification.model_1");
        structures = exp.createNmrStructures(folder)
        hmdb_ids = exp.getStructures(structures, folder)
        structures = exp.getChemicalShifts(folder)
        model_saved_in_train_phase = exp.saved_model[6]
        model = model_saved_in_train_phase["model"]
        test = exp.buildTestClassification(folder, structures)
        dataset = np.array(isTrainingSet)
        X = dataset[:,1:-1]
        X = X.astype(float)
        print("X portion of dataset:{}".format(X))
        y =  dataset[:,-1]
        y = y.reshape(y.shape[0],1)
        y = y.astype(float)
        print("y portion of dataset:{} and shape is".format(y,y.shape))
        y_temp = y[:]
        y_temp = y_temp.ravel()
        factor = pd.factorize(y_temp)
        Y = factor[0]
        Y = Y.reshape(Y.shape[0],1)
        definitions = factor[1]
        reversefactor = dict(zip(range(len(factor[0])),definitions))
        Y = Y.reshape(Y.shape[0],1)
        print("after factorization  Y data is :{} and shape is {}".format(Y,Y.shape))
        average_error = 0
        outlier_index = 0
        X_test = X
        y_test = Y
        y_pred = model.predict(X_test)
        y_pred = np.vectorize(reversefactor.get)(y_pred)
        y_pred = 100 * y_pred 
        y_pred = y_pred.reshape(y_pred.shape[0],1)
        y_test = np.vectorize(reversefactor.get)(y_test)
        y_test = 100 * y_test
        error = 0
        outlier_num = 0
        test_arbitary = 0
        for j in range(len(y_pred)):
            if ((abs(y_pred[j] - y_test[j]) < 20)):
                    #print("no outlier in {}th prediction".format(j))
                test_arbitary = test_arbitary +1
            else:
                outlier_num = outlier_num + 1
                print("true values of the outlier:{}".format(y_test[j]))
                print("HMDB ID of the outlier is:{}".format(hmdb_ids[outlier_index]))
            error = abs(y_pred[j] - y_test[j]) + error
            outlier_index = outlier_index + 1
        error = error / len(y_pred)
        print("avg error in holdout test:{}".format(error))  
        print("total number of outliers:{}".format(outlier_num))
    except:
        print("error occur in holdout prediction")



def runPredictionWithCDK(trainingSetFile, testData):
       # """ generated source for method runPrediction """
    try:
            # LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
            # RandomForest modesl = (RandomForest) weka.core.SerializationHelper.read("models/classification.model_1");
            # structures = exp.createNmrStructures(folder)
            # hmdb_ids = exp.getStructures(structures, folder)
            # structures = exp.getChemicalShifts(folder)
        print("\n\n ######### \n\n Started Prediction\n\n#########")
    #    model_saved_in_train_phase = exp.saved_model[0]
    #    model = model_saved_in_train_phase["model"]
    ##### fit the model with training set. no need for cross validation #########
        model = RandomForestClassifier()
        dataset_training = pd.read_csv(trainingSetFile , header = None)
        print("training data read and dataset_training:{}".format(dataset_training))
        X_training = dataset_training.iloc[1: , 0:-3].values # it was -1 before
        X_training = X_training.astype(float)
        y_training =  dataset_training.iloc[1:,-3].values # it was -1 before
        y_training = y_training.reshape(y_training.shape[0],1)
        y_training = y_training.astype(float)
    ##    y_training = y_training/100
    ##    print("before factorization  y_training data is :{} and shape is {}".format(y_training,y_training.shape))
    ##    y_temp_training = y_training[:]
    ##    y_temp_training = y_temp_training.ravel()
    ##    factor_training = pd.factorize(y_temp_training)
    ##    Y_training = factor_training[0]
    ##    print("after factorization  Y_training data is :{} and shape is {}".format(Y_training,Y_training.shape))
    ##    definitions_training = factor_training[1]
    ##    reversefactor_training = dict(zip(range(len(factor_training[0])),definitions_training))
    ##    print("Printing reverse factor training:{}".format(reversefactor_training))
    ##    Y_training = Y_training.reshape(Y_training.shape[0],1)
        Y_training = exp.classLabel(y_training)
        print("Y_training after classLabel:{}".format(Y_training))
        Y_training = Y_training.reshape(Y_training.shape[0],1)
        Y_training = Y_training.astype(int)
        model.fit(X_training, Y_training)
        print("training done")
    #    print("Saved model's average error was:{}".format(model_saved_in_train_phase["error"]))
    #    print("Saved best model is:{}".format(model))
            #test = exp.buildTestClassification(folder, structures)
    #### reading the data for prediction ###############
    ###################################################
        dataset = pd.read_csv(testData , header = None)
        print("This is how the dataset looks like after converting to np array:{}".format(dataset))
        hmdb_ids = dataset.iloc[1:,-2].values # it was -1 before for xuan's dataset
        hmdb_ids = hmdb_ids.reshape(hmdb_ids.shape[0],1)
        hydrogen_positions = dataset.iloc[1:,-1].values
        hydrogen_positions = hydrogen_positions.reshape(hydrogen_positions.shape[0],1)
        X = dataset.iloc[1: , 0:-3].values # it was -2 before for xuan's dataset
        X = X.astype(float)
        print("X portion of dataset:{}".format(X))
        y =  dataset.iloc[1:,-3].values # it was -1 before for xuan's dataset
        y = y.reshape(y.shape[0],1)
        y = y.astype(float)
    ##    y = y/100
    ##    print("y portion of dataset:{} and shape is".format(y,y.shape))
    ##    y_temp = y[:]
    ##    y_temp = y_temp.ravel()
    ##    print("after ravel, y_temp shape is :{}".format(y_temp.shape))
    ##    factor = pd.factorize(y_temp)
    ##    Y = factor[0]
    ##    definitions = factor[1]
    ##    reversefactor = dict(zip(range(len(factor[0])),definitions))
    ##    Y = Y.reshape(Y.shape[0],1)
    ##    print("after factorization  Y data is :{} and shape is {}".format(Y,Y.shape))
        Y = exp.classLabel(y)
        print("Y after classLabel:{}".format(Y))
        Y = Y.reshape(Y.shape[0],1)
        Y = Y.astype(int)


            # dataset = np.array(isTrainingSet)
            # X = dataset[:,1:-1]
            # X = X.astype(float)
            # print("X portion of dataset:{}".format(X))
            # y =  dataset[:,-1]
            # y = y.reshape(y.shape[0],1)
            # y = y.astype(float)
            # print("y portion of dataset:{} and shape is".format(y,y.shape))
            # y_temp = y[:]
            # y_temp = y_temp.ravel()
            # factor = pd.factorize(y_temp)
            # Y = factor[0]
            # Y = Y.reshape(Y.shape[0],1)
            # definitions = factor[1]
            # reversefactor = dict(zip(range(len(factor[0])),definitions))
            # Y = Y.reshape(Y.shape[0],1)
            # print("after factorization  Y data is :{} and shape is {}".format(Y,Y.shape))
        average_error = 0
        outlier_index = 0
        X_test = X
        y_test = Y
        y_pred = model.predict(X_test)
        print("Prediction phase passed successfully")
        print("Prediction result is:{}".format(y_pred))
        print("True result is:{}".format(y_test))
        print("Prediction result size is:{}".format(y_pred.shape))
    #    y_pred = np.vectorize(reversefactor.get)(y_pred)
    ##    y_pred = np.vectorize(reversefactor_training.get)(y_pred)
    ##    print("y_pred reverse factorization is done:{}".format(y_pred))
    #    y_pred = 10 * y_pred 
        y_pred = y_pred.astype(float)
        y_pred = y_pred/float(100)
        y_pred = y_pred.reshape(y_pred.shape[0],1)
    ##    print("Prediction result after reverse factor is:{}".format(y_pred))
    ##    y_test = np.vectorize(reversefactor.get)(y_test)
    ##    print("y_test after factorzation is:{}".format(y_test))
    #    y_test = 10 * y_test
        y_test = y_test.astype(float)
        y_test = y_test/float(100)
        error = 0
        outlier_num = 0
        test_arbitary = 0
        for j in range(len(y_pred)):
            print("Printing prediction:{}------true:{}, HMDB_ID:{} and H position is:{}".format(y_pred[j], y_test[j],hmdb_ids[j],hydrogen_positions[j]))
            if abs(y_pred[j] - y_test[j]) <= 8.00: # it was 20.0
                print("no outlier in {}th prediction".format(j))
                test_arbitary = test_arbitary +1
            else:
                outlier_num = outlier_num + 1
                print("true values of the outlier:{}".format(y_test[j]))
                    #print("HMDB ID of the outlier is:{}".format(hmdb_ids[outlier_index]))
            error = abs(y_pred[j] - y_test[j]) + error
            outlier_index = outlier_index + 1
        error = error / len(y_pred)
        print("avg error in holdout test:{}".format(error))  
        print("total number of outliers:{}".format(outlier_num))
    except:
        print("error occur in holdout prediction")


if __name__ == '__main__':
    
    folder = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/"
    #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/whole_training_nmr_1063_instance.csv"
    #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/whole_training_nmr_3000_plus_instance.csv"
    #trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_priority.csv"
    trainingSetFile = "/Users/zinatsayeeda/anaconda3/envs/rdkit/training_nmr_1st_2nd_priority_megred.csv"
    #testData = "/Users/zinatsayeeda/anaconda3/envs/rdkit/whole_test_nmr.csv"
    testData = "/Users/zinatsayeeda/anaconda3/envs/rdkit/holdout_nmr.csv"
    print("Start")
    print("Start with:{}".format(sys.argv[1]))
    exp = NmrExp.NmrExperiment(sys.argv[1])
    #runPrediction(folder);
    runPredictionWithCDK(trainingSetFile,testData)
    print("End")




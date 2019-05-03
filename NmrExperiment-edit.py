
import os
from NmrStructure import NmrStructure 
from  GetCDKDescriptors import GetCDKDescriptors
from sklearn.ensemble import RandomForest
from sklearn.cross_validation import cross_val_score
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, auc


class NmrExperiment:
    #  This constructor function runs experiment on a training set doing 10-fold cross validation
    #    * measuring CV accuracy and validating the model
    def __init__(self):
        folder = "dataset/"
        # try:
        #     Instances isTrainingSet = (Instances) weka.core.SerializationHelper.read("models/train_classification_6") # fix it
        #     runClassifier(isTrainingSet, True)
        # except FileNotFoundError as fnf_error:
        #     print(fnf_error)
        trainingSet = buildTrainingClassification(folder)
        runClassifier(trainingSet, False)

    def buildTrainingClassification(self, folder):
        feature_factor = 4
        nmr_structures = []  
        nmr_structures = getChemicalShifts(folder) # retruns array of nmr structures
        getStructures(nmr_structures, folder) # for each structure, there are several H atom position and for all of those H atom position the HMDBIDs are returned
        for nmr_str in nmr_structures:
            print(nmr_str.hmdb_id)
            nmr_str.atomic_descriptors = GetCDKDescriptors.getAtomicDescriptor(nmr_str.structure_sdf, "")
            nmr_str.findNearestAtomToHydrogens(GetCDKDescriptors.getNearestAtoms(nmr_str.structure_sdf))
        values = nmr_structures[0].atomic_descriptors
        isTrainingSet = [[]]
        for nmr_str in nmr_structures:
            for i in range(len(nmr_str.hydrogen_positions)):
                iExample = []
                for j in range(len(nmr_str.atomic_descriptors)):
                    h_position_in_nmr_str = nmr_str.hydrogen_positions[i]
                    iExample.append(nmr_str.atomic_descriptors[h_position_in_nmr_str][j])
                for k in range(feature_factor - 1):
                    for j in range(len(nmr_str.atomic_descriptors)):
                        nearest_ato = nmr_str.nearest_atoms[i][k]
                        iExample.append(nmr_str.atomic_descriptors[nearest_ato][j])
                iExample.append(nmr_str.chemical_shifts[i])
                iExample.append(nmr_str.assignShiftClasses[i])
                isTrainingSet.append(iExample)
        return isTrainingSet.pop(0)


    def runClassifier(self, isTrainingSet, read):
        folds = 10
        true_values = [[]]
        predicted_values = [[]]
        #train = None
        #test = None
        hmdb_ids = mapInstancesToMolecules()
        dataset = np.array(isTrainingSet)
        X = dataset[:,0:-1]
        y =  dataset[:,-1]
        y = y.reshape(Y.shape[0],1)
        cv = StratifiedKFold(n_splits=folds, random_state=123, shuffle=True)
        fprs, tprs, scores = [], [], []
        model = RandomForestClassifier()
        average_error = 0
        outlier_index = 0
        for (train_indes, test_index), i in zip(cv.split(X, y), range(folds)):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index] # y_test is the true value of prediction
            model.fit(X_train, y_train)
            y_pred = model.predict(X_test)
            true_values.append(y_test)
            predicted_values.append(y_pred)
            error = 0
            outlier_num = 0

            for j in range(y_pred.size):
                if !(abs(y_pred[j] - y_test[j]) < 2.0):
                    outlier_num = outlier_num + 1
                    print("true values of the outlier:{}".format(y_test[j])
                    print("HMDB ID of the outlier is:{}".format(hmdb_ids[outlier_index])
                error = abs(y_pred[j] - y_test[j]) + error
                outlier_index = outlier_index + 1

            error = error / y_pred.size()
            print("error in {} fold:{}".format(i,error))  
            println("total number of outliers in {} fold:{}".format(i,outlier_num))
            average_error = error + average_error

        average_error = average_error/folds
        print("after all fold the average error:{}".format(average_error))
    
        #### plot all the values

   
    def mapInstancesToMolecules(self):
        folder = "dataset/"
        nmr_structures = []
        hmdb_ids = []
        nmr_structures = getChemicalShifts(folder)
        getStructures(nmr_structures, folder)
        for nmr_str in nmr_structures:
            for f in nmr_str.chemical_shifts:
                hmdb_ids.append(nmr_str.hmdb_id)
        return hmdb_ids

    
    

   
    # def buildTestClassification(cls, folder, nmr_structures):
    #     """ generated source for method buildTestClassification """
    #     feature_factor = 4
    #     try:
    #         for nmr_str in nmr_structures:
    #             print nmr_str.hmdb_id
    #             nmr_str.atomic_descriptors = GetCDKDescriptors.getAtomicDescriptor(nmr_str.structure_sdf, "")
    #             nmr_str.findNearestAtomToHydrogens(GetCDKDescriptors.getNearestAtoms(nmr_str.structure_sdf))
    #     except Exception as e:
    #         e.printStackTrace()
    #     values = nmr_structures.get(0).atomic_descriptors
    #     attributes = ArrayList()
    #     i = 0
    #     while i < feature_factor * len(values):
    #         attributes.add(Attribute(String.valueOf(i)))
    #         i += 1
    #     fv = FastVector(100)
    #     df = DecimalFormat("#.#")
    #     i = 1
    #     while i <= 100:
    #         fv.addElement(df.format(i / 10))
    #         i += 1
    #     attributes.add(Attribute("Class", fv))
    #     wekaAttributes = FastVector(feature_factor * len(values))
    #     for a in attributes:
    #         wekaAttributes.addElement(a)
    #     isTestSet = Instances("Rel", wekaAttributes, 3000)
    #     isTestSet.setClassIndex(feature_factor * len(values))
    #     for nmr_str in nmr_structures:
    #         while i < len(nmr_str.hydrogen_positions):
    #             print nmr_str.hydrogen_positions.get(i)
    #             while j < len(nmr_str.atomic_descriptors):
    #                 iExample.setValue(wekaAttributes.elementAt(j), nmr_str.atomic_descriptors.get(j)[Integer.valueOf(nmr_str.hydrogen_positions.get(i))])
    #                 j += 1
    #             while k < feature_factor - 1:
    #                 while j < len(nmr_str.atomic_descriptors):
    #                     iExample.setValue(wekaAttributes.elementAt(j + (k + 1) * len(values)), nmr_str.atomic_descriptors.get(j)[Integer.valueOf(nmr_str.nearest_atoms.get(i).get(k))])
    #                     j += 1
    #                 k += 1
    #             # iExample.setValue((Attribute)wekaAttributes.elementAt(feature_factor*len(values)), );
    #             isTestSet.add(iExample)
    #             i += 1
    #     return isTestSet

   #""" generated source for method createNmrStructures """
    def createNmrStructures(self, folder):
        structures = []
        for fileEntry in os.listdir(folder):
            if fileEntry.endswith( ".sdf" ):
                structure = NmrStructure(fileEntry.split(".")[0])
                structures.append(structure)
        return structures

    
    def getChemicalShifts(self, folder):
        structures = []
        file_names = []
        for fileEntry in os.listdir(folder):
            if fileEntry.endswith( ".txt" ):
                carbon_position = []
                chemical_shift = []
                file_names.append(fileEntry)
                readChemicalShifts(folder+fileEntry, carbon_position, chemical_shift)
                structure = NmrStructure(carbon_position, chemical_shift, "", fileEntry) # no need fileEntry.split(".")[0]
                structure.assignShiftClasses(structure.chemical_shifts)
                structures.append(structure)
        return structures

    
    def readChemicalShifts(self, file, c_pos, c_shift): #why carbon positions are showing three times at the begining
        #reader = None
        try:
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
                    print("The file with file path: {}".format(file))
                    if items[0].strip() == "No." or ("M" in items[2]):
                        continue 
                    c_pos.append(items[1].strip())
                    if items[2] != None and items[2] != "":
                        c_shift.append(items[2].strip())
                    else:
                        c_shift.append(items[3].strip())
                    print("The file with file path: {}".format(file))
                if len(items) == 0 or len(items) == 1:
                    continue 
                if "Atom" in items[1].strip():
                    parse_shifts = True
        except FileNotFoundError as e:
            print(e)
        except IOError as e:
            print(e)
        finally:
            try:
                if fp != None:
                    fp.close()
            except IOError as e:
                print(e)

    
    def getStructures(self, structures, folder):
        hmdb_ids = []
        for nmr_s in structures:
            file = nmr_s.hmdb_id + ".sdf"
            try:
                text = Chem.SDMolSupplier(file)
                nmr_s.structure_sdf = text
                if nmr_s.has_chemical_shift == False:
                    nmr_s.hydrogen_positions = GetCDKDescriptors.getHydrogenAtoms(nmr_s.structure_sdf) #H atom's position number's are returned
                    nmr_s.has_chemical_shift = True
                for s in nmr_s.hydrogen_positions:
                    hmdb_ids.append(nmr_s.hmdb_id)
            except FileNotFoundError as e:
                print(e)
            except IOError as e:
                print(e)
        return hmdb_ids



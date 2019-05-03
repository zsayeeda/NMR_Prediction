
""" generated source for module NmrPred """
import os
from NmrExperiment import NmrExperiment
# import java.lang.String

# import java.io

# import java.util

# import java.util.regex

# import java.nio.file_

# import weka.core

# import weka.classifiers.functions

# import weka.classifiers.functions.supportVector

# import weka.classifiers.Evaluation

# import weka.classifiers.evaluation.Prediction

# import weka.classifiers.trees

# import weka.classifiers

import matlabcontrol

if __name__ == '__main__':
    
    folder = "test/"
    #runPrediction(folder);
    exp = NmrExperiment();


    #  RunPrediction function builds a classifier from a saved Training Set. It then builds a testing instance
    #    * from a test folder. It then runs predict function to calculate ppm values for all Hydrogen atoms in 
    #    * for all the test molecules
    #    
    def runPrediction(self, folder):
        """ generated source for method runPrediction """
        try:
            # LinearRegression model = (LinearRegression) weka.core.SerializationHelper.read("models/regression.model_3d");
            # RandomForest modesl = (RandomForest) weka.core.SerializationHelper.read("models/classification.model_1");
            structures = NmrExperiment.createNmrStructures(folder)
            hmdb_ids = NmrExperiment.getStructures(structures, folder)
            model.buildClassifier(isTrainingSet)
            for nmr_str in structures:
                for h_pos in nmr_str.hydrogen_positions:
                    #  Class is a double from 1 to 100 classes. So divide by 10 to get shift
                    nmr_str.chemical_shifts.add(Float.valueOf(label))
                    print hmdb_ids.get(i) + "\t" + label
                    i += 1
            callMatlab(structures)
        except Exception as e:
            e.printStackTrace()

    @classmethod
    def callMatlab(cls, structures):
        """ generated source for method callMatlab """
        options = MatlabProxyFactoryOptions.Builder().setUsePreviouslyControlledSession(True).build()
        factory = MatlabProxyFactory(options)
        proxy = factory.getProxy()
        proxy.eval("cd ..")
        proxy.eval("cd matlab")
        for nmr_str in structures:
            proxy.eval("inter.coupling.scalar = []")
            while i < len(nmr_str.hydrogen_positions):
                if i == len(nmr_str.hydrogen_positions) - 1:
                    isotopes = isotopes + "'1H'}"
                    ppms = ppms + String.valueOf(nmr_str.chemical_shifts.get(i)) + "}"
                    proxy.eval(coupling)
                else:
                    isotopes = isotopes + "'1H', "
                    ppms = ppms + String.valueOf(nmr_str.chemical_shifts.get(i)) + ", "
                while j < len(nmr_str.hydrogen_positions):
                    proxy.eval(coupling)
                    j += 1
                i += 1
            print ppms
            proxy.eval(isotopes)
            proxy.eval(ppms)
            proxy.feval("create_nmr1H_plot")
        proxy.disconnect()



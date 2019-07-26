#!/jython2.7.0/bin/jython.exe

import sys
import glob
import os,sys
#jarfiles = glob.glob("C:\\mystuff\\bin\\java\\*.jar")
#for jar in jarfiles:
sys.path.append("/Users/zinatsayeeda/anaconda3/envs/rdkit/cdk-2.2.jar")

import java.io
from org.openscience import cdk

def molToCDK(filename):
    reader = cdk.io.MDLV2000Reader(java.io.FileReader(filename))
    molecule = cdk.AtomContainer()
    # print("type of molecule:{}".format(molecule))
    #print("length of molecule is :{}".format(len(molecule)))
    return reader.read(molecule)

count_NONCARBONPROTON = 0
count_CONJUGATED = 0
count_NONROTATABLE = 0
count_ALIPHATIC = 0
count_AROMATIC = 0
count_file = 0
excluded_compound =[]
with open("/Users/zinatsayeeda/anaconda3/envs/rdkit/2D_Holdout_data.tsv", "r") as fp:
    lines = fp.readlines()
    for line in lines:
        excluded_compound.append(line.split("\n")[0])

#parentFolder = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority"
#parentFolder = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority/test"
parentFolder = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_with_all_atoms/"
for dirName, subdirs, fileList in os.walk(parentFolder):
    for file_name in fileList:
        if file_name.endswith( ".sdf" ):
            file_name_temp = file_name.split(".sdf")[0]
            if file_name_temp not in excluded_compound:
                print("file name is:{}".format(file_name))
                count_file = count_file +1
#filename = 'Conformer3D_CID_16078.sdf'
#filename = "/Users/zinatsayeeda/anaconda3/envs/rdkit/HMDB0000001.sdf"
                filename = parentFolder+"/"+file_name
                mol = molToCDK(filename)
                mol.setProperty("FILENAME",filename)

#set varisou properties
                cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol)

#sets both ring and arom flags
                arom_model = cdk.aromaticity.Aromaticity(cdk.aromaticity.ElectronDonation.daylight(),
                            cdk.graph.Cycles.or(cdk.graph.Cycles.all(), cdk.graph.Cycles.all(6))
                            )
                arom_model.apply(mol)

#also set non aromatic rings
                cdk.graph.Cycles.markRingAtomsAndBonds(mol)

#if no hydrogens then maybe?
                cdk.tools.manipulator.AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol)

                testarom = cdk.qsar.descriptors.atomic.IsProtonInAromaticSystemDescriptor()
                testconj = cdk.qsar.descriptors.atomic.IsProtonInConjugatedPiSystemDescriptor()

#print(mol.getProperty("FILENAME"))

                for a in mol.atoms():

                    printstr = []
    
                    a_index = mol.indexOf(a)+1
                    a_symbol = a.getSymbol()
                    print("index:{}".format(a_index))
                    print("atom symbol:{}".format(a_symbol.upper()))
                    a_isarom = a.isAromatic()
                    a_inring = a.isInRing()
                    a_descr_arom = testarom.calculate(a,mol).getValue().intValue()
                    a_descr_conj = testconj.calculate(a,mol).getValue().booleanValue()
                    print(str(a_index)+" "+a_symbol+" "+"aromatic:"+str(a_descr_arom))
                    print(str(a_index)+" "+a_symbol+" "+"conj:"+str(a_descr_conj))
                    
                    
                    printstr.append(str(a_index))
                    printstr.append(a_symbol)
    
                    if a.getSymbol() == "H":
    
                        a_bonded = mol.getConnectedAtomsList(a)
                        # print("connected atoms:{}".format(a_bonded))
                        for a2 in a_bonded: #only should be one atom bonded to H
                            a2_index = mol.indexOf(a2)+1
                            a2_symbol = a2.getSymbol()
                            print("connect atom symbol:{}".format(a2_symbol))
                            a2_isarom = a2.isAromatic()
                            a2_inring = a2.isInRing()

                        if a2_symbol != "C":
                            printstr.append("NONCARBONPROTON")
                            count_NONCARBONPROTON = count_NONCARBONPROTON + 1
                        elif a_descr_arom == 1: #or if a2_isAromatic...
                            #print("the value of is proton in aromatic system:{}".format(a_descr_arom))
                            printstr.append("AROMATIC")
                            count_AROMATIC = count_AROMATIC +1
                        elif a_descr_conj:
                            #print("the value of is proton in conjugated pi system:{}".format(a_descr_arom))
                            printstr.append("CONJUGATED")
                            count_CONJUGATED = count_CONJUGATED + 1
                        elif a2_inring:
                            printstr.append("NONROTATABLE")
                            count_NONROTATABLE = count_NONROTATABLE +1
                        else:
                            printstr.append("ALIPHATIC")
                            count_ALIPHATIC = count_ALIPHATIC +1
    
                #else:
                    #printstr.append("NONPROTON")
    
                    print(" ".join(printstr))
print("NONCARBONPROTON:{}".format(count_NONCARBONPROTON))
print("AROMATIC:{}".format(count_AROMATIC))
print("CONJUGATED:{}".format(count_CONJUGATED))
print("NONROTATABLE:{}".format(count_NONROTATABLE))
print("ALIPHATIC:{}".format(count_ALIPHATIC))
print("total molecule:{}".format(count_file))


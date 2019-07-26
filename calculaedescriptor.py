#!/jython2.7.0/bin/jython.exe

import sys
import glob
import os,sys
import operator
import csv

#jarfiles = glob.glob("C:\\mystuff\\bin\\java\\*.jar")
#for jar in jarfiles:
sys.path.append("/Users/zinatsayeeda/anaconda3/envs/rdkit/cdk-2.2.jar")

import java.io
from org.openscience import cdk
import javax.vecmath.Point3d;

def molToCDK(filename):
    reader = cdk.io.MDLV2000Reader(java.io.FileReader(filename))
    molecule = cdk.AtomContainer()
    #print("type of molecule:{}".format(molecule))
    #print("length of molecule is :{}".format(len(molecule)))
    return reader.read(molecule)


def getHydrogenAtom(raw_mol_file): #didn't want to keep the same parameter name as getAtomicDescriptor parameter
    print("Started taking H atom from mol file and excluding H atoms those are connected with O and N:")
    hydrogen_positions = []
    mol = raw_mol_file
    # atoms = mol.atoms()
    # print("this is how atom looks:{}".format(atoms))
    for atom in mol.atoms():
        atom_symbol = atom.getSymbol()
        # print("atom symbol:{}".format(atom_symbol))
            
        if atom_symbol.upper() == "H":
            connected_atoms = mol.getConnectedAtomsList(atom)
            for single_connected_atom in connected_atoms:
                single_connected_atom_symbol = single_connected_atom.getSymbol()
                #print("connected atom symbol:{}".format(single_connected_atom_symbol))
                if (single_connected_atom_symbol.upper() != "O" and single_connected_atom_symbol.upper() != "N"):
                    hydrogen_positions.append(mol.indexOf(atom))
                
    print("So the hydrogen positions connected with C only:{}".format(hydrogen_positions))
    return hydrogen_positions


def getNearestAtoms(transferred_mol_file): #didn't want to keep the same parameter name as getAtomicDescriptor parameter
    #print("Started calculating nearest atoms for all atoms")
    mol = transferred_mol_file
    atom_distances = [[]]

    atom_count = mol.getAtomCount()
    #print("total number of atoms:{}".format(atom_count))
        #i = 0
        #for atom1 in atoms: #### why this iteration is not stopping after 1st iteration
    for i in range(atom_count): #cdk.interfaces.IAtom.getPoint3d()
            #print("started with atom:{}".format(atom1))
        # print("started with atom:{}".format(mol.getAtom(i)))
        # print("i increases to:{}".format(i))
        distances = []  # distances represent the distances for one atom only. that means one--> other atoms
        j = 0
            #for atom2 in atoms:
        for j in range(atom_count): 
            if i == j:
                # print("i = j reached")
                distances.append(99999.12)
            else:
                atom1_point = mol.getAtom(i).getPoint3d()
                atom2_point = mol.getAtom(j).getPoint3d()
                dst = atom1_point.distance(atom2_point)
                # print("distance from i:{} to j:{} is:{}".format(i,j,dst))
                distances.append(dst)
                #j = j+1
        # print("distance from i:{} to all the distances:{}".format(i,distances))
            # indices =[]
            # d = distances[:]
            # d.sort()
            # print("sorted distance:{}".format(d))
            # print(" before sorting the distances are:{}".format(distances))
            # d_list = distances
            # for x in range(len(distances)):
            #     value_of_d_x = d[x]
            #     index = d_list.index(value_of_d_x)
            #     indices.append(index)
            #     if x !=0
            # print("indices of sorted distance:{}".format(indices))
            # atom_distances.append(indices)
            # #i = i+1
        index_of_distances = [ x for x in range(len(distances))]
        index_of_distances_tuple = tuple(index_of_distances)
        values_of_distances_tuple = tuple(distances)
        dic_with_index_value_of_distances = dict(zip(index_of_distances_tuple,values_of_distances_tuple))
        sorted_dic_with_index_value_of_distances = sorted(dic_with_index_value_of_distances.items(), key=operator.itemgetter(1))
        indices = [v[0] for v in sorted_dic_with_index_value_of_distances]
        # print("indices of sorted distance:{}".format(indices))
        atom_distances.append(indices)
        # print("i increases to:{}".format(i))
    atom_distances.pop(0)
    # print("calculated nearest atoms for all atoms:{}".format(atom_distances))
    # print("the size of calculated nearest atoms for all atoms.{}".format(len(atom_distances)))
    return atom_distances


def ReadShiftFile1(file_name, h_position):
  #print("h_position is:{} and length:{}".format(h_position,len(h_position)))
  H_position_and_chemicalshift_in_shift_file = [] # this returns after checking "O" and "N" as neighbor atoms
  hydrogen_positions_in_mol_file = []
  hydrogen = []
  shift = []
  with open(file_name, "r") as fp:
    lines = fp.readlines()
    for i in range(2,len(lines)):
      splited_line = lines[i].split("\t")
      len_of_splitted_line = len(splited_line)
      # print("H:{}".format(splited_line[0].strip()))
      # print("S:{}".format(splited_line[len_of_splitted_line-1].strip()))
      hydrogen.append(int(splited_line[0].strip()))
      shift.append(float(splited_line[len_of_splitted_line-1].strip()))
    # print("hydrogen:{}, shif:{}".format(hydrogen, shift))
  for i in range(len(h_position)):
    if h_position[i] in hydrogen:
      index_of_hydrogen_array = hydrogen.index(h_position[i])
      H_position_and_chemicalshift_in_shift_file.append({"H_position": hydrogen[index_of_hydrogen_array], "chemical_shift": shift[index_of_hydrogen_array]})
  # print("returning H_position_and_chemicalshift_in_shift_file:{}".format(H_position_and_chemicalshift_in_shift_file))
  return H_position_and_chemicalshift_in_shift_file





def getAtomicDescriptor(mol_file):
    Headers = ["AtomdDegree","AtomHybridization","AtomHybridizationVSEPR","AtomValence","BondsToAtom","CovalentRadius","DistanceToAtom","EffectiveAtomPolarizability","InductiveAtomicHardness","InductiveAtomicSoftness","IpAtomicHOSE","IpAtomicLearning","IsProtonInAromaticSystem","IsProtonInConjugatedPiSystem","PartialPiCharge","PartialSigmaCharge","PartialTChargeMMFF94","PartialTChargePEOE","PeriodicTablePosition","PiElectronegativity","ProtonAffinityHOSE","ProtonTotalPartial","RDFProtonDescriptorG3R","RDFProtonDescriptorGDR","RDFProtonDescriptorGHR","RDFProtonDescriptorGHRTopol","RDFProtonDescriptorGSR","SigmaElectronegativity","StabilizationPlusCharge","VdWRadius","ChemicalShift","HMDBId","HydrogenPosition"]
    descriptors = [[]]
    #descriptors.append(Header)
    for a in mol.atoms():
        a_index = mol.indexOf(a)
        a_symbol = a.getSymbol()
        print("index:{} and symbol:{}".format(a_index,a_symbol.upper()))
        
        values = []
        atom_degree_descriptor = cdk.qsar.descriptors.atomic.AtomDegreeDescriptor().calculate(a,mol).getValue()
        # atom_hybridization_descriptor = cdk.qsar.descriptors.atomic.AtomHybridizationDescriptor().calculate(a,mol).getValue()
        # atom_hybridization_VSEPR_descriptor = cdk.qsar.descriptors.atomic.AtomHybridizationVSEPRDescriptor().calculate(a,mol).getValue()
        # atom_valence_descriptor = cdk.qsar.descriptors.atomic.AtomValenceDescriptor().calculate(a,mol).getValue()
        # bonds_to_atom_descriptor = cdk.qsar.descriptors.atomic.BondsToAtomDescriptor().calculate(a,mol).getValue()
        # covalent_radius_descriptor = cdk.qsar.descriptors.atomic.CovalentRadiusDescriptor().calculate(a,mol).getValue()
        # distance_to_atom_descriptor = cdk.qsar.descriptors.atomic.DistanceToAtomDescriptor().calculate(a,mol).getValue()
        # effective_atom_polarizability_descriptor = cdk.qsar.descriptors.atomic.EffectiveAtomPolarizabilityDescriptor().calculate(a,mol).getValue()
        # inductive_atomic_hardness_descriptor = cdk.qsar.descriptors.atomic.InductiveAtomicHardnessDescriptor().calculate(a,mol).getValue()
        # inductive_atomic_softness_descriptor = cdk.qsar.descriptors.atomic.InductiveAtomicSoftnessDescriptor().calculate(a,mol).getValue()
        # ip_atomic_HOSE_descriptor = cdk.qsar.descriptors.atomic.IPAtomicHOSEDescriptor().calculate(a,mol).getValue()
        # ip_atomic_learning_descriptor = cdk.qsar.descriptors.atomic.IPAtomicLearningDescriptor().calculate(a,mol).getValue()
        # is_proton_in_aromatic_system_descriptor = cdk.qsar.descriptors.atomic.IsProtonInAromaticSystemDescriptor().calculate(a,mol).getValue()
        is_proton_in_conjugated_pi_system_descriptor = cdk.qsar.descriptors.atomic.IsProtonInConjugatedPiSystemDescriptor().calculate(a,mol).getValue().booleanValue()
        if str(is_proton_in_conjugated_pi_system_descriptor) == "False":
            is_proton_in_conjugated_pi_system_descriptor = 0
        else:
            is_proton_in_conjugated_pi_system_descriptor = 1
        # partial_pi_charge_descriptor = cdk.qsar.descriptors.atomic.PartialPiChargeDescriptor().calculate(a,mol).getValue()
        # partial_sigma_charge_descriptor = cdk.qsar.descriptors.atomic.PartialSigmaChargeDescriptor().calculate(a,mol).getValue()
        # partial_t_charge_MMFF94_descriptor = cdk.qsar.descriptors.atomic.PartialTChargeMMFF94Descriptor().calculate(a,mol).getValue()
        # partial_t_charge_PEOE_descriptor = cdk.qsar.descriptors.atomic.PartialTChargePEOEDescriptor().calculate(a,mol).getValue()
        # periodic_table_position_descriptor = cdk.qsar.descriptors.atomic.PeriodicTablePositionDescriptor().calculate(a,mol).getValue()

        # hAdder = cdk.tools.HydrogenAdder()
        # hAdder.addExplicitHydrogensToSatisfyValency(mol)
                    
        # pi_electronegativity_descriptor = cdk.qsar.descriptors.atomic.PiElectronegativityDescriptor().calculate(a,mol).getValue()
        # proton_affinity_HOSE_descriptor = cdk.qsar.descriptors.atomic.ProtonAffinityHOSEDescriptor().calculate(a,mol).getValue()
        # proton_total_partial_charge_descriptor = cdk.qsar.descriptors.atomic.ProtonTotalPartialChargeDescriptor().calculate(a,mol).getValue()
        # RDF_proton_descriptor_G3R = cdk.qsar.descriptors.atomic.RDFProtonDescriptor_G3R().calculate(a,mol).getValue()
        # RDF_proton_descriptor_GDR = cdk.qsar.descriptors.atomic.RDFProtonDescriptor_GDR().calculate(a,mol).getValue()
        # RDF_proton_descriptor_GHR = cdk.qsar.descriptors.atomic.RDFProtonDescriptor_GHR().calculate(a,mol).getValue()
        # RDF_proton_descriptor_GHR_topol = cdk.qsar.descriptors.atomic.RDFProtonDescriptor_GHR_topol().calculate(a,mol).getValue()
        # RDF_proton_descriptor_GSR = cdk.qsar.descriptors.atomic.RDFProtonDescriptor_GSR().calculate(a,mol).getValue()
        # sigma_electronegativity_descriptor = cdk.qsar.descriptors.atomic.SigmaElectronegativityDescriptor().calculate(a,mol).getValue()
        stabilization_plus_charge_descriptor = cdk.qsar.descriptors.atomic.StabilizationPlusChargeDescriptor().calculate(a,mol).getValue()
        VdWRadius_descriptor = cdk.qsar.descriptors.atomic.VdWRadiusDescriptor().calculate(a,mol).getValue()
        
        # descriptors.append([atom_degree_descriptor,atom_hybridization_descriptor,atom_hybridization_VSEPR_descriptor,atom_valence_descriptor,bonds_to_atom_descriptor,covalent_radius_descriptor,distance_to_atom_descriptor,effective_atom_polarizability_descriptor,inductive_atomic_hardness_descriptor,inductive_atomic_softness_descriptor,ip_atomic_HOSE_descriptor,ip_atomic_learning_descriptor,is_proton_in_aromatic_system_descriptor,is_proton_in_conjugated_pi_system_descriptor,partial_pi_charge_descriptor,partial_sigma_charge_descriptor,partial_t_charge_MMFF94_descriptor,partial_t_charge_PEOE_descriptor,periodic_table_position_descriptor,pi_electronegativity_descriptor,proton_affinity_HOSE_descriptor,proton_total_partial_charge_descriptor,RDF_proton_descriptor_G3R,RDF_proton_descriptor_GDR,RDF_proton_descriptor_GHR,RDF_proton_descriptor_GHR_topol,RDF_proton_descriptor_GSR,sigma_electronegativity_descriptor,stabilization_plus_charge_descriptor,VdWRadius_descriptor])
        descriptors.append([atom_degree_descriptor,is_proton_in_conjugated_pi_system_descriptor,stabilization_plus_charge_descriptor,VdWRadius_descriptor])
    #temp_heardes = ["atom_degree_descriptor","is_proton_in_conjugated_pi_system_descriptor","stabilization_plus_charge_descriptor","VdWRadius_descriptor"]
    descriptors[0] = Headers
    # print("descriptors beforr pop:{}".format(descriptors))
    #descriptors[0] = temp_heardes
        #print(str(a_index)+" "+a_symbol+"  atom_degree="+str(atom_degree_descriptor))
            #print(str(a_index)+" "+a_symbol+"  atom_hybridization="+str( atom_hybridization_descriptor))
            #print(str(a_index)+" "+a_symbol+"  atom_hybridization_VSEPR="+str(atom_hybridization_VSEPR_descriptor))
            #print(str(a_index)+" "+a_symbol+"  atom_valence="+str(atom_valence_descriptor))
            #print(str(a_index)+" "+a_symbol+"  bonds_to_atom_descriptor="+str(bonds_to_atom_descriptor))
            #print(str(a_index)+" "+a_symbol+"  covalent_radius_descriptor="+str(covalent_radius_descriptor))
            #print(str(a_index)+" "+a_symbol+"  distance_to_atom_descriptor="+str(distance_to_atom_descriptor)) # check
            #print(str(a_index)+" "+a_symbol+"  effective_atom_polarizability_descriptor="+str(effective_atom_polarizability_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  inductive_atomic_hardness_descriptor="+str(inductive_atomic_hardness_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  inductive_atomic_softness_descriptor="+str(inductive_atomic_softness_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  ip_atomic_HOSE_descriptor="+str(ip_atomic_HOSE_descriptor)) # all are zero
            #print(str(a_index)+" "+a_symbol+"  ip_atomic_learning_descriptor="+str(ip_atomic_learning_descriptor)) #check the values
            #print(str(a_index)+" "+a_symbol+"  is_proton_in_aromatic_system_descriptor="+str(is_proton_in_aromatic_system_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  is_proton_in_conjugated_pi_system_descriptor="+str(is_proton_in_conjugated_pi_system_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  partial_pi_charge_descriptor="+str(partial_pi_charge_descriptor)) #check the calculation. all the values for H is zero
            #print(str(a_index)+" "+a_symbol+"  partial_sigma_charge_descriptor="+str(partial_sigma_charge_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  partial_t_charge_MMFF94_descriptor="+str(partial_t_charge_MMFF94_descriptor))
            #print(str(a_index)+" "+a_symbol+"  partial_t_charge_PEOE_descriptor="+str(partial_t_charge_PEOE_descriptor))
            #print(str(a_index)+" "+a_symbol+"  periodic_table_position_descriptor="+str(periodic_table_position_descriptor))
            # print(str(a_index)+" "+a_symbol+"  pi_electronegativity_descriptor="+str(pi_electronegativity_descriptor)) # check the calcultation
            # print(str(a_index)+" "+a_symbol+"  proton_affinity_HOSE_descriptor="+str(proton_affinity_HOSE_descriptor))
            # print(str(a_index)+" "+a_symbol+"  proton_total_partial_charge_descriptor="+str(proton_total_partial_charge_descriptor))# check the calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_G3R="+str(RDF_proton_descriptor_G3R)) # check calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_GDR="+str(RDF_proton_descriptor_GDR)) # check calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_GHR="+str(RDF_proton_descriptor_GHR)) #check the calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_GHR_topol="+str(RDF_proton_descriptor_GHR_topol)) #check the calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_GSR="+str(RDF_proton_descriptor_GSR)) #check the calsultation
            # print(str(a_index)+" "+a_symbol+"  sigma_electronegativity_descriptor="+str(sigma_electronegativity_descriptor)) #check the calculation
            # print(str(a_index)+" "+a_symbol+"  stabilization_plus_charge_descriptor="+str(stabilization_plus_charge_descriptor))#check the calculation
            # print(str(a_index)+" "+a_symbol+"  VdWRadius_descriptor="+str(VdWRadius_descriptor))
    return descriptors


            #print(str(a_index)+" "+a_symbol+"  atom_degree="+str(atom_degree_descriptor))
            #print(str(a_index)+" "+a_symbol+"  atom_hybridization="+str( atom_hybridization_descriptor))
            #print(str(a_index)+" "+a_symbol+"  atom_hybridization_VSEPR="+str(atom_hybridization_VSEPR_descriptor))
            #print(str(a_index)+" "+a_symbol+"  atom_valence="+str(atom_valence_descriptor))
            #print(str(a_index)+" "+a_symbol+"  bonds_to_atom_descriptor="+str(bonds_to_atom_descriptor))
            #print(str(a_index)+" "+a_symbol+"  covalent_radius_descriptor="+str(covalent_radius_descriptor))
            #print(str(a_index)+" "+a_symbol+"  distance_to_atom_descriptor="+str(distance_to_atom_descriptor)) # check
            #print(str(a_index)+" "+a_symbol+"  effective_atom_polarizability_descriptor="+str(effective_atom_polarizability_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  inductive_atomic_hardness_descriptor="+str(inductive_atomic_hardness_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  inductive_atomic_softness_descriptor="+str(inductive_atomic_softness_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  ip_atomic_HOSE_descriptor="+str(ip_atomic_HOSE_descriptor)) # all are zero
            #print(str(a_index)+" "+a_symbol+"  ip_atomic_learning_descriptor="+str(ip_atomic_learning_descriptor)) #check the values
            #print(str(a_index)+" "+a_symbol+"  is_proton_in_aromatic_system_descriptor="+str(is_proton_in_aromatic_system_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  is_proton_in_conjugated_pi_system_descriptor="+str(is_proton_in_conjugated_pi_system_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  partial_pi_charge_descriptor="+str(partial_pi_charge_descriptor)) #check the calculation. all the values for H is zero
            #print(str(a_index)+" "+a_symbol+"  partial_sigma_charge_descriptor="+str(partial_sigma_charge_descriptor)) #check the calculation
            #print(str(a_index)+" "+a_symbol+"  partial_t_charge_MMFF94_descriptor="+str(partial_t_charge_MMFF94_descriptor))
            #print(str(a_index)+" "+a_symbol+"  partial_t_charge_PEOE_descriptor="+str(partial_t_charge_PEOE_descriptor))
            #print(str(a_index)+" "+a_symbol+"  periodic_table_position_descriptor="+str(periodic_table_position_descriptor))
            # print(str(a_index)+" "+a_symbol+"  pi_electronegativity_descriptor="+str(pi_electronegativity_descriptor)) # check the calcultation
            # print(str(a_index)+" "+a_symbol+"  proton_affinity_HOSE_descriptor="+str(proton_affinity_HOSE_descriptor))
            # print(str(a_index)+" "+a_symbol+"  proton_total_partial_charge_descriptor="+str(proton_total_partial_charge_descriptor))# check the calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_G3R="+str(RDF_proton_descriptor_G3R)) # check calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_GDR="+str(RDF_proton_descriptor_GDR)) # check calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_GHR="+str(RDF_proton_descriptor_GHR)) #check the calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_GHR_topol="+str(RDF_proton_descriptor_GHR_topol)) #check the calculation
            # print(str(a_index)+" "+a_symbol+"  RDF_proton_descriptor_GSR="+str(RDF_proton_descriptor_GSR)) #check the calsultation
            # print(str(a_index)+" "+a_symbol+"  sigma_electronegativity_descriptor="+str(sigma_electronegativity_descriptor)) #check the calculation
            # print(str(a_index)+" "+a_symbol+"  stabilization_plus_charge_descriptor="+str(stabilization_plus_charge_descriptor))#check the calculation
            # print(str(a_index)+" "+a_symbol+"  VdWRadius_descriptor="+str(VdWRadius_descriptor))








if __name__ == '__main__':
    if len(sys.argv) > 1:
        folders = sys.argv[1]
        number_of_nearest_atom = int(sys.argv[2])
        print("input folders:{}".format(folders))
        if os.path.exists(folders):
                print("folder exists")
        else:
            print('%s is not a valid path, please verify' % folders)
            sys.exit()
    else:
        print('Usage: jython calculatedescriptor.py folder name with path: jython calculatedescriptor.py /Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority/test 4(which is 3 nearest neighbor')

    #temp_heardes = ["atom_degree_descriptor","is_proton_in_conjugated_pi_system_descriptor","stabilization_plus_charge_descriptor","VdWRadius_descriptor","Chemical_shift","HMDB_ID","Hydrogen_position"]
    count_file = 0
    isTrainingSet = [[]]
    excluded_compound =[]
    with open("/Users/zinatsayeeda/anaconda3/envs/rdkit/2D_Holdout_data.tsv", "r") as fp:
        lines = fp.readlines()
        for line in lines:
            excluded_compound.append(line.split("\n")[0])

    parentFolder = folders
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
                    cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol)
                    arom_model = cdk.aromaticity.Aromaticity(cdk.aromaticity.ElectronDonation.daylight(),
                            cdk.graph.Cycles.or(cdk.graph.Cycles.all(), cdk.graph.Cycles.all(6))
                            )
                    arom_model.apply(mol)
                    cdk.graph.Cycles.markRingAtomsAndBonds(mol)
                    cdk.tools.manipulator.AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol)
                    hmdb_id = file_name_temp
                    features = getAtomicDescriptor(mol)
                    print("features:{} and length:{}".format(features,len(features)))
                    hydrogen_positions_in_mol_file = getHydrogenAtom(mol)
                    print("after parsing mol file H positions are:{}".format(hydrogen_positions_in_mol_file))
                    nearest_H = getNearestAtoms(mol)
                    print("nearest atoms:{}".format(nearest_H))
                    # print("nearest atoms length:{}".format(len(nearest_H)))
                    # print("taining data calculation for the file:{}".format(hmdb_id))
                    chem_file_path = parentFolder+"/"+hmdb_id+".txt"
                    H_position_and_chemicalshift_in_shift_file = ReadShiftFile1(chem_file_path,hydrogen_positions_in_mol_file)
                    values = len(features[1])
                    feature_factor = number_of_nearest_atom
                    # print("features1{}".format(features[1]))
                    for i in range(len(hydrogen_positions_in_mol_file)):
                        iExample = []
                        h_position_in_nmr_str = hydrogen_positions_in_mol_file[i]
                        print("hydrogen position in nmr str:{}".format(h_position_in_nmr_str))
                        for j in range(values):
                            iExample.append(features[int(h_position_in_nmr_str)+1][j])
                        print("iExample after{} th H:{}".format(h_position_in_nmr_str, iExample))
                        for k in range(feature_factor - 1):
                            nearest_atom_for_i_H_and_k_position = nearest_H[int(h_position_in_nmr_str)][k]
                            print("for {}th H {}th nearest atom is:{}".format(h_position_in_nmr_str,k,nearest_atom_for_i_H_and_k_position))
                            for j in range(values):
                                iExample.append(features[int(nearest_atom_for_i_H_and_k_position)+1][j])
                            print(" i example for nearest atom:{}, is:{}".format(nearest_atom_for_i_H_and_k_position, iExample))
                            
                        h_p = h_position_in_nmr_str
                        for d in H_position_and_chemicalshift_in_shift_file:
                            if  d["H_position"] == h_p:
                                iExample.append(d["chemical_shift"])
                                iExample.append(hmdb_id)
                                iExample.append(h_p)
                        isTrainingSet.append(iExample)
    isTrainingSet[0] = features[0]
    print("final training set:{}".format(isTrainingSet))
    with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority/test/training_nmr_1st_2nd_priority_without_2D.csv', 'w') as csv_file: #change path and file name here too
        writer = csv.writer(csv_file)
        for i in range(len(isTrainingSet)):
            data_row = isTrainingSet[i]
            writer.writerow(data_row)
    print("total number of file:{}".format(count_file))








               
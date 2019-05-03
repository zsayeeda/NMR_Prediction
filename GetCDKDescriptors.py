## from rdkit it calculates atomic features ##
from rdkit import Chem
import os
import numpy as np
from scipy.spatial import distance
import operator


def getAtomicDescriptor(sdf):
    mols = sdf
    print("read {} mols from the given file".format(len(mols)))
    descriptors = ['Atoms','GetAtomicNum','GetBonds','GetDegree','GetExplicitValence','GetFormalCharge','GetHybridization','GetImplicitValence','GetIsAromatic','GetIsotope','GetNoImplicit','GetNumExplicitHs','GetNumImplicitHs','GetNumRadicalElectrons','GetTotalDegree','GetTotalNumHs','IsInRing']
    values = [[]]
    i = 0
    for mol in mols:
        mol = Chem.AddHs(mol)
        atoms = mol.GetAtoms()
        for atom in atoms:
            atomic_num = atom.GetAtomicNum()
            bonds = len(atom.GetBonds())
            get_degree = atom.GetDegree()
            exp_valence = atom.GetExplicitValence()
            formal_charge = atom.GetFormalCharge()
            hyberdization = 0
            hyber = atom.GetHybridization().name
            if hyber == "SP":
                hyberdization = 2
            elif hyber == "SP2":
                hyberdization = 3
            elif hyber == "SP3":
                hyberdization = 4
            impli_valence = atom.GetImplicitValence()
            bool_aromatic = atom.GetIsAromatic()
            if bool_aromatic == "True":
                is_aromatic = 1
            else: 
                is_aromatic = 0
            get_iso = atom.GetIsotope()
            impli = atom.GetNoImplicit()
            num_of_implicit = 0
            if impli == "True":
                num_of_implicit = 1
            else:
                num_of_implicit = 0
            num_of_exp_h = atom.GetNumExplicitHs()
            num_of_impli_h = atom.GetNumImplicitHs()
            num_of_radical_elec = atom.GetNumRadicalElectrons()
            total_degree = atom.GetTotalDegree()
            total_num_h = atom.GetTotalNumHs()
            ring = atom.IsInRing()
            is_in_ring = 0
            if ring == "True":
                is_in_ring = 1
            else:
                is_in_ring = 0
            symbol = atom.GetSymbol().upper()
            values.append([symbol,atomic_num,bonds,get_degree,exp_valence,formal_charge,hyberdization,impli_valence,is_aromatic,get_iso,num_of_implicit,num_of_exp_h,num_of_impli_h,num_of_radical_elec,total_degree,total_num_h,is_in_ring])
    values[0] = descriptors
    print("printing descriptors:\n\n")
    print(values)
    return values

    #  Find all relevant hydrogen atoms in molecule
   
def getHydrogenAtoms(sdf): # it returns the 
    print("Started taking H atom from getcdkdescriptor.py file and excluding H atoms those are connected with O and N:")
    hydrogen_positions = []
    mols = sdf
    if len(mols) != 0:
        for mol in mols:
            mol = Chem.AddHs(mol) # we need to add H atoms in the mol
            for atom in mol.GetAtoms():
                #print("before testing 'O' and 'N' atoms, the atom is:{} and index is:{}".format(atom,atom.GetIdx()))
                #if ((atom.GetSmarts()).upper() != "O") and ((atom.GetSmarts()).upper() != "N"):
                if (atom.GetSymbol().upper()) == "H":
                    connected_atoms = atoms.GetNeighbors()
                    if ((connected_atoms[0].GetSymbol()).upper() != "O" and (connected_atoms[0].GetSymbol()).upper() != "N"):
                        hydrogen_positions.append(atom.GetIdx())
                # if ((atom.GetSymbol()).upper() != "O") and ((atom.GetSmarts()).upper() != "N"):
                #     print("after testing 'O' and 'N' atoms, the atom is:{}".format(atom))
                #     hydrogen_positions.append(atom.GetIdx())
    print("So the hydrogen positions connected with C only:{}".format(hydrogen_positions))
    return hydrogen_positions

def getNearestAtoms(sdf):
    print("Started calculating nearest atoms for all atoms")
    mols = sdf
    atom_distances = [[]]
    for mol in mols:
        mol = Chem.AddHs(mol)
        atoms = mol.GetAtoms()
        #i = 0
        #for atom1 in atoms: #### why this iteration is not stopping after 1st iteration
        for i in range(len(atoms)):
            #print("started with atom:{}".format(atom1))
            print("started with atom:{}".format(atoms[i]))
            print("i increases to:{}".format(i))
            distances = []  # distances represent the distances for one atom only. that means one--> other atoms
            j = 0
            #for atom2 in atoms:
            for j in range(len(atoms)): 
                if i == j:
                    print("i = j reached")
                    distances.append(99999.12)
                else:
                    conf = mol.GetConformer()
                    atom1_point = conf.GetAtomPosition(i)
                    atom2_point = conf.GetAtomPosition(j)
                    firstPoint = (atom1_point.x,atom1_point.y,atom1_point.z)
                    secondPoint = (atom2_point.x,atom2_point.y,atom2_point.z)
                    dst = distance.euclidean(firstPoint, secondPoint)
                    print("distance from i:{} to j:{} is:{}".format(i,j,dst))
                    distances.append(dst)
                #j = j+1
            print("distance from i:{} to all the distances:{}".format(i,distances))
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
            print("indices of sorted distance:{}".format(indices))
            atom_distances.append(indices)
            print("i increases to:{}".format(i))
    atom_distances.pop(0)
    print("calculated nearest atoms for all atoms:{}".format(atom_distances))
    print("the size of calculated nearest atoms for all atoms.{}".format(len(atom_distances)))
    return atom_distances







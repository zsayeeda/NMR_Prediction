from rdkit import Chem
import os
import numpy as np
from scipy.spatial import distance
import operator

def main():
  #file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_1/log_out_2.txt"
  file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_3/log_out.txt"
  file_text = open(file_path).read()
  file_blocks = file_text.split("file entry in create structure:")
  dataset_array = []
  for i in range(1, len(file_blocks)): # iterate over each block
    array_of_lines_in_each_block = file_blocks[i].split("\n") # each block is splitted. fetching the required values for each block starts from here.
    hmdb_id = []
    hydrogen_position_string = []
    hydrogen_position = []
    hmdb_id.append(array_of_lines_in_each_block[0].split(".sdf")[0])
    print("block:{} has hmbd_id:{}".format(i,array_of_lines_in_each_block[0].split(".sdf")[0]))
    index_of_starts = array_of_lines_in_each_block.index('*********************************')
    for j in range(2,index_of_starts):
      hydrogen_position_string.append(array_of_lines_in_each_block[j])
      print("block:{} has H positions:{}".format(i,array_of_lines_in_each_block[j]))
    hydrogen_position = [int(m) for m in hydrogen_position_string]
    starting_of_index_of_descriptors_in_block_lines_array = index_of_starts+1
    descriptors_final_values = [[]] #this is the final descriptors array
    for d in range(starting_of_index_of_descriptors_in_block_lines_array, len(array_of_lines_in_each_block)-1): # for each H
      raw_descriptors = array_of_lines_in_each_block[d].split(",")
      raw_descriptors.pop()
      descriptors_string_values = []
      descriptors = []
      for k in range(len(raw_descriptors)):
        if k == 0:
          descriptors_string_values.append(raw_descriptors[k].strip().split("[")[1])
        else:
          descriptors_string_values.append(raw_descriptors[k].strip())
      descriptors = [float(l) for l in descriptors_string_values]
      print("block:{} has descriptors:{}".format(i,descriptors))
      descriptors_final_values.append(descriptors) 
    descriptors_final_values.pop(0)
    print("For block:{}, HMDB IDs:{}, H Positions:{} and Descriptor Values:{}".format(i, hmdb_id,hydrogen_position,descriptors_final_values))
    chem_file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_1/"+array_of_lines_in_each_block[0].split(".sdf")[0]
    print("chem file path:{}".format(chem_file_path))
    H_position_and_chemicalshift_in_shift_file = ReadShiftFile(chem_file_path)
    print("in chemical shift file{}".format(H_position_and_chemicalshift_in_shift_file))
    for i in range(len(hydrogen_position)):
      h_p = hydrogen_position[i]
      for d in H_position_and_chemicalshift_in_shift_file:
        if  d["H_position"] == h_p:
          dataset_array.append({"descriptor": descriptors_final_values[i], "chemical_shift": d["chemical_shift"], "HMDB_ID": hmdb_id[0], "hydrogen_position": d["H_position"]})
  print("Final Dataset:{}".format(dataset_array))


def ReadShiftFile(file_name):
  H_position_and_chemicalshift_in_shift_file = [] # this returns after checking "O" and "N" as neighbor atoms
  hydrogen_positions_in_mol_file = []
  with open(file_name+".txt", "r") as fp:
    lines = fp.readlines()
    sdf_file_name = file_name+".sdf"
    print("sdf file name:{}".format(sdf_file_name))
    sdf = Chem.SDMolSupplier(sdf_file_name)
    mols = sdf
    if len(mols) != 0:
      for mol in mols:
        mol = Chem.AddHs(mol)
        atoms = mol.GetAtoms()
        print("Atoms:{}".format(atoms))
        for atom in mol.GetAtoms():
          print("Atom symbol:{}".format(atom.GetSymbol().upper()))
          if atom.GetSymbol().upper() == "H":
            connected_atoms = atom.Getneighbors()
            if ((connected_atoms[0].GetSymbol()).upper() != "O" and (connected_atoms[0].GetSymbol()).upper() != "N"):
              hydrogen_positions_in_mol_file.append(atom.GetIdx())
    for i in range(2,len(lines)):
      splitted_line = lines[i].split("\t") 
      len_of_splitted_line = len(splited_line)
      print("H:{}".format(splited_line[0].strip()))
      print("S:{}".format(splited_line[len_of_splitted_line-1].strip()))
      hydrogen = int(splited_line[0].strip())
      shift = float(splited_line[len_of_splitted_line-1].strip())
      if hydrogen in  hydrogen_positions_in_mol_file:
        H_position_and_chemicalshift_in_shift_file.append({"H_position": hydrogen, "chemical_shift": shift})
  print("returning H_position_and_chemicalshift_in_shift_file:{}".format(H_position_and_chemicalshift_in_shift_file))
  return H_position_and_chemicalshift_in_shift_file



  if __name__ == '__main__':
    print("Start")
    main()
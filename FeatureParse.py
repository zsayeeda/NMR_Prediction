# do not delete new line at the end of log_out file"
#always use this file for parsing (until get-descriptor is fixed)
from rdkit import Chem
import os
import numpy as np
from scipy.spatial import distance
import operator
import csv

def main():
  #file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_1/log_out_2.txt"
  #file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_3/log_out.txt"
  #file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_4/log_out.txt"
  #file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_5/log_out.txt"
  #file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_priority/log_out.txt" # change here the file name
  #file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_2nd_priority/log_out.txt"
  #file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority/log_out.txt"
  file_path = "//Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_with_all_atoms/log_out.txt"
  file_text = open(file_path).read()
  file_blocks = file_text.split("file entry in create structure:")
  print("file blocks are:{}".format(file_blocks))
  dataset_array = []
  header = ["5YQ8","D6M9","SAR9","ZQ50","1OGV","UVGB","TWWV","FZUL","LDBJ","3O99","C0TG","U5NM","K03U","2FZG","68AK","O54T","XXYJ","ZT3G","KCS9","JV9V","L8O7","08W0","3O2D","QAID","OVCT","YAQ5","C77C","8K3B","HTZ9","G4N0","67JR","N8KG","9ZJ9","4CX3","NAD7","X3KZ","8YQA","TV2E","1B1U","XYFG","AN91","YPYO","OURC","I3T6","GASV","A9XG","35N1","YKSQ","3BMV","YKAQ","WOLC","TALS","Z7E2","VQP0","HPNW","EM8L","G3VD","R8H6","SBQX","GGQ4","ZE85","CPLP","DI39","7T6E","R7OH","6HUB","M1UL","75PM","2J0O","UX3R","CSIF","FXYO","LHG6","IX4T","MCGP","GCJI","1X82","RRH0","Y07Q","WY2M","HNZD","LBZF","LREX","IHWP","FZ2K","ZWGF","A0GN","2S89","RRNR","IUFK","OM9G","ZZ6G","N4PZ","1HI2","BOJ3","2HZR","LBDE","D9QV","VN0S","TQY0","6OP0","TR6T","CLHS","LY8T","3EVY","MK8O","1FO8","1CMK","CDQJ","LKU8","0MRV","UILJ","X9XK","XEQX","KOBD","3GTM","ChemicalShift","HMDBID","HPosition"]
  for i in range(1, len(file_blocks)): # iterate over each block
    array_of_lines_in_each_block = file_blocks[i].split("\n") # each block is splitted. fetching the required values for each block starts from here.
    print("array of lines in each block:{}".format(array_of_lines_in_each_block))
    hmdb_id = []
    hydrogen_position_string = []
    hydrogen_position = []
    hmdb_id.append(array_of_lines_in_each_block[0].split(".sdf")[0])
    print("block:{} has hmbd_id:{}".format(i,array_of_lines_in_each_block[0].split(".sdf")[0]))
    index_of_starts = array_of_lines_in_each_block.index('*********************************')
    print("index of starts:{}".format(index_of_starts))
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
    #chem_file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_1/"+array_of_lines_in_each_block[0].split(".sdf")[0] #change the path here
    #chem_file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_4/"+array_of_lines_in_each_block[0].split(".sdf")[0] #change the path here
    #chem_file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/test_5/"+array_of_lines_in_each_block[0].split(".sdf")[0] #change the path here
    chem_file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_with_all_atoms/"+array_of_lines_in_each_block[0].split(".sdf")[0] #change the path here
    #chem_file_path = "/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority/"+array_of_lines_in_each_block[0].split(".sdf")[0] #change the path here
    print("chem file path:{}".format(chem_file_path))
    H_position_and_chemicalshift_in_shift_file = ReadShiftFile1(chem_file_path,hydrogen_position)
    print("in chemical shift file{}".format(H_position_and_chemicalshift_in_shift_file))
    for i in range(len(hydrogen_position)):
      h_p = hydrogen_position[i]
      for d in H_position_and_chemicalshift_in_shift_file:
        if  d["H_position"] == h_p:
          dataset_array.append({"descriptor": descriptors_final_values[i], "chemical_shift": d["chemical_shift"], "HMDB_ID": hmdb_id[0], "hydrogen_position": d["H_position"]})
  print("Final Dataset:{}".format(dataset_array))
  #with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/test_1/holdout_nmr_2.csv', 'w') as csv_file:
  #with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_priority/training_nmr_1st_priority.csv', 'w') as csv_file: #change path and file name here too
  #with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/test_4/holdout_nmr_test4.csv', 'w') as csv_file: #change path and file name here too
  #with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/test_5/holdout_nmr_test5.csv', 'w') as csv_file: #change path and file name here too
  with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority_with_all_atoms/training_nmr_1st_2nd_priority_all_atom.csv', 'w') as csv_file: #change path and file name here too
  #with open('/Users/zinatsayeeda/anaconda3/envs/rdkit/dataset/dataset_1st_2nd_priority/training_nmr_1st_2nd_priority_without_2D.csv', 'w') as csv_file: #change path and file name here too
    writer = csv.writer(csv_file)
    writer.writerow(header)
    for i in range(len(dataset_array)):
      data_row = [d_value for d_value in dataset_array[i]["descriptor"]]
      data_row.append(dataset_array[i]["chemical_shift"])
      data_row.append(dataset_array[i]["HMDB_ID"])
      data_row.append(dataset_array[i]["hydrogen_position"])
      writer.writerow(data_row)



#don't use it for now    
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
          print("Atom Symbol:{}".format(atom.GetSymbol().upper()))
          if (atom.GetSymbol().upper()) == "H":
            connected_atoms = atom.GetNeighbors()
            if ((connected_atoms[0].GetSymbol()).upper() != "O" and (connected_atoms[0].GetSymbol()).upper() != "N"):
              hydrogen_positions_in_mol_file.append(atom.GetIdx())
    for i in range(2,len(lines)):
      splited_line = lines[i].split("\t")
      len_of_splitted_line = len(splited_line)
      print("H:{}".format(splited_line[0].strip()))
      print("S:{}".format(splited_line[len_of_splitted_line-1].strip()))
      hydrogen = int(splited_line[0].strip())
      shift = float(splited_line[len_of_splitted_line-1].strip())
      if hydrogen in  hydrogen_positions_in_mol_file:
        H_position_and_chemicalshift_in_shift_file.append({"H_position": hydrogen, "chemical_shift": shift})
  print("returning H_position_and_chemicalshift_in_shift_file:{}".format(H_position_and_chemicalshift_in_shift_file))
  return H_position_and_chemicalshift_in_shift_file

#use it for now     
def ReadShiftFile1(file_name, h_position):
  print("h_position is:{} and length:{}".format(h_position,len(h_position)))
  H_position_and_chemicalshift_in_shift_file = [] # this returns after checking "O" and "N" as neighbor atoms
  hydrogen_positions_in_mol_file = []
  hydrogen = []
  shift = []
  with open(file_name+".txt", "r") as fp:
    lines = fp.readlines()
    for i in range(2,len(lines)):
      splited_line = lines[i].split("\t")
      len_of_splitted_line = len(splited_line)
      print("H:{}".format(splited_line[0].strip()))
      print("S:{}".format(splited_line[len_of_splitted_line-1].strip()))
      hydrogen.append(int(splited_line[0].strip()))
      shift.append(float(splited_line[len_of_splitted_line-1].strip()))
    print("hydrogen:{}, shif:{}".format(hydrogen, shift))
  for i in range(len(h_position)):
    if h_position[i] in hydrogen:
      index_of_hydrogen_array = hydrogen.index(h_position[i])
      H_position_and_chemicalshift_in_shift_file.append({"H_position": hydrogen[index_of_hydrogen_array], "chemical_shift": shift[index_of_hydrogen_array]})
  print("returning H_position_and_chemicalshift_in_shift_file:{}".format(H_position_and_chemicalshift_in_shift_file))
  return H_position_and_chemicalshift_in_shift_file








if __name__ == '__main__':
  main()
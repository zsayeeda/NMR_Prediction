 
#generated source for module NmrStructure 

class NmrStructure:
    hydrogen_positions = []
    chemical_shifts = []
    c_shift_classes = []
    nearest_atoms = [[]]
    structure_sdf = ""
    hmdb_id = ""
    has_chemical_shift = False
    atomic_descriptors = [[]]


    def __init__(self, h_pos, c_shifts, sdf, id):
        self.hydrogen_positions = h_pos
        self.chemical_shifts = c_shifts
        self.structure_sdf = sdf
        print("Inside the NMRSTRUCTURE.py file, name of the file entry passed:{}".format(id))
        print(" And length of  hydrogen postions:{} and length of chemical shifts:{}".format(len(self.hydrogen_positions),len(self.chemical_shifts)))
        print("And hydrogen positions:{},  carbon positions:{}".format(self.hydrogen_positions,self.chemical_shifts))
        temp_id = id.split('.txt')
        self.hmdb_id = temp_id[0].strip()
        self.c_shift_classes = []
        if (len(self.hydrogen_positions) != 0 and len(self.chemical_shifts) != 0):
            self.has_chemical_shift = True
            print("From nmrstructure.py file the has_chemical_shif is assigned for the file:{} is {}".format(id,self.has_chemical_shift))
        print("From NMRSTRUCTURE.py got H atom position and chemical shift for the file:{}".format(temp_id[0].strip()))

    def __init___(self, hmid):
        print("received id from create structure in nmrexp file{}".format(hmid))
        self.hmdb_id = hmid
        self.chemical_shifts = []
        print("assigned id:{}".format(self.hmdb_id))

    def assignShiftClasses(self, c_shifts):
        for c_shift in c_shifts:
            #c_shift_number = round(float(c_shift), 1)
            #c_shift_number_class = c_shift_number/10
            c_shift_number = float(c_shift)/100
            #self.c_shift_classes.append(round(c_shift_number,2))  # for now no round figue 
            self.c_shift_classes.append(c_shift_number)

    #  [ [ 0, 4, 2, 1], [ 5, 3, 2, 10] ] => indices of nearest_atoms
    #  Only contains hydrogen atoms but they can contain indices of any atoms, sorted by distance
    def findNearestAtomToHydrogens(self, n_atom):
        print("initially the nearest_atom variable is:{}".format(self.nearest_atoms))
        self.nearest_atoms = [[]]
        print("nearest atoms received:{}".format(n_atom))
        print("the size of n_atom is:{} and print n_atom[0]:{}".format(len(n_atom),n_atom[0]))
        nearest_pos = []
        for h_position in self.hydrogen_positions:
            #nearest_pos = []
            int_value_of_H_position = int(h_position)
            print("The int value of H position is:{}".format(int_value_of_H_position))
            print("nearest atoms of hydrogen position{} is:{}".format(h_position,n_atom[int_value_of_H_position]))
            nearest_pos.append(n_atom[int_value_of_H_position]) 
        print("before adding to self nearest atoms, the nearest_pos is : {}".format(nearest_pos))
        self.nearest_atoms.append(nearest_pos)
            #self.nearest_atoms.append(n_atom[int(h_position)])
        #self.nearest_atoms.pop(0)
        print("final nearest atom for {} file/nmr structure:{}".format(self.hmdb_id,self.nearest_atoms))
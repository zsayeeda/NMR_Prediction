 
""" generated source for module NmrStructure """

class NmrStructure():
    """ generated source for class NmrStructure """
    self.hydrogen_positions = []
    self.chemical_shifts = []
    self.c_shift_classes = []
    self.nearest_atoms = [[]]
    self.structure_sdf = ""
    self.hmdb_id = ""
    self.has_chemical_shift = False
    self.atomic_descriptors = [[]]


    def __init__(self, h_pos, c_shifts, sdf, id):
        self.hydrogen_positions = h_pos
        self.chemical_shifts = c_shifts
        self.structure_sdf = sdf
        temp_id = id.split('.txt')
        self.hmdb_id = temp_id[0].strip()
        self.c_shift_classes = []
        self.has_chemical_shift = True

    def __init___(self, id):
        self.hmdb_id = id
        self.chemical_shifts = []

    def assignShiftClasses(self, c_shifts):
        for c_shift in c_shifts:
            c_shift_number = round(float(c_shift), 1)
            c_shift_number_class = c_shift_number/10
            self.c_shift_classes.appened(c_shift_number)

    #  [ [ 0, 4, 2, 1], [ 5, 3, 2, 10] ] => indices of nearest_atoms
    #  Only contains hydrogen atoms but they can contain indices of any atoms, sorted by distance
    def findNearestAtomToHydrogens(self, n_atom):
        self.nearest_atoms = [[]]
        nearest_pos = []
        for h_position in self.hydrogen_positions:
            nearest_pos = n_atom[h_position] 
            self.nearest_atoms.append(nearest_pos)
        self.nearest_atoms.pop(0)

    # def printDescriptors(self):   
    #     for h_pos in self.nearest_atoms:
    #         for index in h_pos:
    #             for values in self.atomic_descriptors:
    #                 print String.valueOf(values[Integer.valueOf(index)])
    #             print "For atom number ---- " + index